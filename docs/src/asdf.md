# ASDF and metadata

```@meta
CurrentModule = GWCS
```

GWCS.jl reads and writes the WCS metadata Python GWCS stores in ASDF files. [`parse_wcs`](@ref) builds a WCS from a
tagged metadata mapping, such as a node of a file loaded with ASDF.jl, or from tagged YAML text:

```@example asdf
using GWCS

roman = parse_wcs(read(joinpath(pkgdir(GWCS), "test", "fixtures", "roman_wcs.yaml"), String))
available_frames(roman), pixel_to_world(roman, 11.0, 22.0)
```

With ASDF.jl loaded, a WCS in an ASDF file is read from its node:

```julia
using ASDF
using GWCS

file = ASDF.load("roman_cal.asdf"; extensions=true, validate_checksum=false)
wcs = parse_wcs(file["roman"]["meta"]["wcs"])
```

## Supported metadata

The core parser supports:

- GWCS WCS, step, simple and compound bounding boxes, and generic, 2D, celestial, spectral, temporal, Stokes,
  composite, and name-only frames;
- ICRS and other astropy reference frames, kept as stored when they carry more than the frame name (such as a fiducial
  position), structured temporal reference metadata, units, axis classifications, and frame axis order;
- compose, concatenate, elementwise arithmetic, identity, shift, scale, affine, axis remapping, and polynomial
  transforms;
- legacy and current constant transforms, fixed-input compounds, and transform names and input and output names;
- array-, dictionary-, range-, and transform-backed label mappers and region selectors;
- the 27 FITS projections that have ASDF tags (all but ZPN; a conic `sigma` is required, as in FITS), 2D rotations,
  native/celestial rotations, Cartesian and spherical rotation sequences, and spherical/Cartesian conversions;
- explicit inverses;
- astropy FITS WCSs (`tag:astropy.org:astropy/wcs/wcs`), which asdf-astropy writes for an `astropy.wcs.WCS` as its FITS
  header. The header is parsed with [`parse_fits_wcs`](@ref), the stored `pixel_shape` is kept, and `pixel_bounds`
  becomes the bounding box.

Unknown tags, unsupported schema versions, malformed dimensions, and missing required properties raise errors that
name their path in the tree. Built-in parsers register the exact tag versions in the ASDF schema packages instead of
assuming that every future `1.x` version is compatible. Rotation and projection angles given as quantities in deg, rad,
arcmin, arcsec, mas, or hourangle are converted to degrees, and dimensionless projection parameters must have an empty
unit; the parsed transform records each unit, and writing restores it. Other unit-bearing parameters, such as a shift
in pixels, raise an error.

## Writing

[`to_metadata`](@ref) returns the tagged tree of a WCS, transform, frame, bounding box, or label mapper. It can be
inspected, printed with YAML.jl, or parsed back:

```@example asdf
using GWCS: Scale, Shift, parse_transform, to_metadata
using YAML

tree = to_metadata(Shift(1.0) & Scale(2.0))
print(YAML.write(tree))
```

The shift and scale were fused into one affine when combined with `&`; see [What a round trip keeps](@ref).

```@example asdf
parse_transform(tree)
```

With ASDF.jl loaded, [`write_asdf`](@ref) writes a WCS to a file, and `ASDF.write_file` and `save` write GWCS objects
anywhere in a document. ASDF.jl handles the file, its blocks, and every other value in the document, so a matrix there
must be an `ASDF.NDArrayWrapper`:

```@example asdf
using ASDF

path = joinpath(mktempdir(), "roman.asdf")
write_asdf(path, roman)
pixel_to_world(parse_wcs(ASDF.load_file(path; extensions=true)["wcs"]), 11.0, 22.0)
```

Arrays are written as Python GWCS writes them: affine matrices and translations, polynomial coefficients, label maps,
and unmodeled arrays are ndarrays in uncompressed binary blocks, the default of Python asdf. ASDF.jl 2.1.1 cannot write
an array of `String`s as an ndarray, so writing a label map of strings fails. ASDF.jl's writer keeps global state, so
files must not be written from several threads at once.

## What a round trip keeps

Serialization represents the compiled Julia object, not the exact source tree. A round trip keeps frames, annotations
(names, inputs, and outputs), explicit inverses, bounding boxes, and the mathematical transform, with values changing at
most by rounding; the tree itself may be canonicalized.

- [`compose`](@ref) and [`parallel`](@ref) fuse unannotated transforms, so a sequence or concatenation of shifts,
  scales, and affines may be written as one equivalent affine, and an affine larger than 2D is written back as a
  concatenation of its 1D and 2D diagonal blocks. Annotated nodes are never fused, so pipeline files that annotate
  their transforms keep their stored tree.
- Output uses the latest exact schema versions the writer supports, so original tag versions are not kept.
- A two-axis [`CoordinateFrame`](@ref) whose axes are both spatial and that has no reference frame is written as
  `frame2d`, which Python GWCS reads as its `Frame2D`; other `CoordinateFrame`s are written as `frame`. Reading a
  `frame2d` ignores `axes_type`, as Python GWCS does, so its axes are spatial.
- FITS SIP, TPV, and Paper IV TPD distortions have no ASDF tags, so they are written as the standard polynomial
  transforms they equal, which Python GWCS evaluates; their reverse polynomials (SIP `AP`/`BP`, TPD `REV`) only start
  the numerical inverse and are not written.
- Callable closures, callable label mappers, affines of more than two dimensions that are not block diagonal,
  unregistered external transforms, TPV and TPD radial terms, D2IM and CPDIS lookup tables, `-TAB` coordinate arrays,
  and nonlinear spectral algorithms have no ASDF form, and writing them raises an error.

### Preserved properties

Properties GWCS.jl does not model are kept and written back unchanged, on the object they were read on:

- Python GWCS's parameter constraints (`fixed` and `bounds`) and `input_units_equivalencies`;
- a `bounding_box` on any transform but the first step's, whose box is the WCS's;
- unknown properties of transforms, label mappers, frames, steps, and the WCS.

An unmodeled array stored in an ASDF block is read while parsing, so a parsed WCS does not keep the file open, and is
written back as an ndarray with its element type. A composite frame with other properties is rejected, as in Python
GWCS.

### Effects of editing

An edit that keeps a step's frame keeps the step's and the frame's properties, and a replaced transform takes its own
with it. `fix_inputs`, which reduces the input frame, drops those of the first step. A transform written with another
tag family than it was read with, such as an Euler `rotate3d` written as a composition, is written without them,
since they describe the original tag's parameters.

## Registering transforms

External packages can add metadata support for their own [`AbstractTransform`](@ref) subtypes without modifying
GWCS.jl. Dimensions are in the subtype, and inverses, Jacobians, and optimized forms are ordinary method extensions:

```julia
module ExampleTransforms

using GWCS
using GWCS: AbstractTransform, Shift, metadata_value, register_transform!, register_transform_serializer!, required_metadata
using StaticArrays
import GWCS: compile_transform, inverse_transform, transform_dependencies, transform_jacobian

struct OffsetTransform{T<:Real} <: AbstractTransform{1,1}
    offset::T
end

@inline (t::OffsetTransform)(coords::StaticVector{1}) = SVector(coords[1] + t.offset)
inverse_transform(t::OffsetTransform) = OffsetTransform(-t.offset)
transform_jacobian(t::OffsetTransform, coords::StaticVector{1}) = SMatrix{1,1}(one(eltype(coords)))
transform_dependencies(::OffsetTransform) = SMatrix{1,1,Bool}(true)
compile_transform(t::OffsetTransform) = Shift(t.offset)

function parse_offset(node, context)
    offset = metadata_value(required_metadata(node, "offset", context))
    return OffsetTransform(offset)
end

serialize_offset(t, context) = Dict("offset" => t.offset)

function __init__()
    register_transform!("tag:example.org:gwcs/offset", parse_offset; versions=(v"1.0.0", v"1.1.0"), properties=("offset",))
    register_transform_serializer!(OffsetTransform, "tag:example.org:gwcs/offset", serialize_offset; version=v"1.1.0")
end

end
```

Tags use their full namespace and omit the trailing schema version. Prefer `versions` when support is known for exact
schema tags; half-open `min_version`/`max_version` ranges are for external schemas that guarantee compatibility, and
overlapping registrations are rejected. [`register_transform!`](@ref) and [`register_transform_serializer!`](@ref)
describe the other keywords. Register from an extension package's `__init__`, so the process-wide
[`default_transform_registry`](@ref) is populated at run time. For deterministic or restricted conversion, use
`registry = copy(default_transform_registry())`, register into that isolated registry, and pass it to
[`parse_transform`](@ref), `parse_wcs`, or `to_metadata`.

Registration and metadata callbacks run only while constructing a transform. A parsed `WCS` stores concrete compiled
transforms, so evaluation performs no registry lookup or dynamic metadata dispatch. Registered constructors must return
an `AbstractTransform`, and custom inverse and compilation methods are checked for compatible dimensions before the
pipeline is accepted.
