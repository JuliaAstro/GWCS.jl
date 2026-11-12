# Transforms

```@meta
CurrentModule = GWCS
```

Every transform is an [`AbstractTransform`](@ref) with its input and output dimensions in its type, and is called on a
static vector. Transforms combine three ways:

- `first | second` applies `first`, then `second` to its outputs ([`compose`](@ref));
- `first & second` applies them side by side, `first` to the first inputs and `second` to the rest
  ([`parallel`](@ref));
- `first + second`, and likewise `-`, `*`, `/`, and `^`, evaluates both on the same inputs and combines their outputs
  ([`elementwise`](@ref)).

```@example transforms
using GWCS
using GWCS: Scale, Shift
using StaticArrays

offset = Shift(1.0) & Scale(2.0)  # x + 1 and 2y
offset(SVector(3.0, 4.0))
```

```@example transforms
(offset | (Scale(10.0) & Shift(-8.0)))(SVector(3.0, 4.0))
```

```@example transforms
(Shift(1.0) * Scale(2.0))(SVector(3.0))  # (x + 1) * 2x
```

The library includes identity, shift, scale, affine, constant-output, fixed-input, and axis-remapping transforms, one-
and two-dimensional polynomials, spherical and Cartesian conversions, 2D and 3D rotations, and all 28 FITS spherical
projections. The [Transform API](@ref) lists them.

## Projections

Projection algorithms are separate from transform direction. [`IntermediateToNative`](@ref) converts projection-plane
coordinates to native spherical longitude and latitude, and [`NativeToIntermediate`](@ref) is the inverse direction.
Both sides use degrees:

```@example transforms
using GWCS: IntermediateToNative, ZenithalPerspective, inverse_transform

pix2sky = IntermediateToNative(ZenithalPerspective(1.5, 15.0))
sky2pix = inverse_transform(pix2sky)
sky2pix(pix2sky(SVector(1.0, 2.0)))
```

Parameterized projections use the ASDF parameter names of the FITS projection codes. Invalid parameters raise an
`ArgumentError`; coordinates on singular, overlapping, or non-visible projection branches, and non-finite coordinates,
give `NaN` in both outputs. FITS parsing uses the same projection types, so an algorithm has one numeric
implementation whatever its metadata source.

[`celestial_transform`](@ref) combines a projection with the FITS Paper II native-to-celestial rotation for a reference
longitude and latitude. It computes the default `LONPOLE` when omitted and uses `LATPOLE` to choose between general pole
solutions:

```@example transforms
using GWCS: GnomonicProjection, celestial_transform

tangent_plane_to_sky = celestial_transform(GnomonicProjection(), (83.8221, -5.3911))
tangent_plane_to_sky(SVector(0.0, 0.0))  # the reference longitude and latitude
```

It returns a transform rather than a complete WCS or frame: it is the low-level construction used by the FITS parser
and programmatic pipelines.

## Simplification

[`compose`](@ref) and [`parallel`](@ref) simplify as they combine: nested compositions flatten, identities drop out,
and neighboring shifts, scales, and affines fuse into one affine. An [`AnnotatedTransform`](@ref) is
never fused, so a pipeline read from metadata keeps its structure; only the compiled WCS kernel fuses it. `compose`
also recognizes Python GWCS's distortion idiom `Mapping((0, 1, 0, 1)) | Polynomial2D & Polynomial2D`, written
`RemapAxes{2}((1, 2, 1, 2)) | (px & py)`, and builds one two-output polynomial that is written back as the same idiom;
when both polynomials are linear it becomes an affine with an exact inverse.

Arithmetic transforms have no analytical inverse, because combining invertible operands does not in general give a
uniquely invertible function. A square WCS inverts them numerically, or [`WithInverse`](@ref) supplies an inverse.

## Inverses

[`inverse_transform`](@ref) returns a transform's analytical or explicit inverse, or [`NoInverse`](@ref). A WCS whose
steps all have one evaluates their composition in [`world_to_pixel`](@ref). Otherwise, for a WCS with as many world as
pixel axes, `world_to_pixel` solves numerically:

```julia
world_to_pixel(wcs, world; tolerance=1e-8, maxiter=50, initial_guess=(10.0, 20.0), quiet=true)
```

A solve that does not converge throws [`NoConvergence`](@ref), or with `quiet=true` returns its best pixel. The solver
first inverts the pipeline's longest analytically invertible tail, such as the projection and sky rotation of a
TAN-SIP WCS, and iterates only on the rest, as astropy's `all_world2pix` does. A world coordinate no pixel can reach
then gives `fill_value`, or a `DomainError` with `strict=true`, just as with an analytical inverse.
[`numerical_inverse`](@ref) runs the solver directly, even for a WCS with an analytical inverse.

## Callable transforms

[`CallableTransform`](@ref) turns functions of positional coordinates into a transform:

```@example transforms
using GWCS: CallableTransform

forward(x, y) = (x + y, x - y)
backward(u, v) = ((u + v) / 2, (u - v) / 2)
jacobian(x, y) = SMatrix{2,2}(1, 1, 1, -1)

callable = CallableTransform{2,2}(forward; inverse=backward, jacobian)
inverse_transform(callable)(callable(SVector(3.0, 1.0)))
```

The callable, inverse, Jacobian, and dependency declarations are stored in concrete type parameters, so named functions
and type-stable closures remain inferred and allocation-free. Supplied Jacobians propagate through composed and
parallel transforms and are used by numerical inversion. A `DomainError` thrown by the user function (for example Base
`sqrt` of a negative number during a numerical-inverse step) becomes `NaN`, while other errors propagate.

Use `CallableTransform` for local functions, prototypes, and closures over calibration values. Define an
`AbstractTransform` subtype when a transform needs metadata tags, serialization, specialized compilation, or a
reusable package API; see [Registering transforms](@ref).

## Region-selected transforms

[`RegionSelector`](@ref) implements discontinuous transforms such as IFU slices, detector regions, and spectral orders.
A label mapper chooses a label from the full input coordinate, and the selector evaluates the transform for that label:

```@example transforms
using GWCS: ArrayLabelMapper, RegionSelector

labels = [0 2; 1 0]
regions = RegionSelector(ArrayLabelMapper(labels; no_label=0), (
    1 => (Shift(1.0) & Scale(1.0)),
    2 => (Shift(2.0) & Scale(2.0)),
))
regions(SVector(0.0, 1.0)), regions(SVector(1.0, 0.0))
```

Array lookup follows the GWCS pixel convention: `(x, y)` coordinates are zero-based pixel centers, rounded with
`floor(value + 0.5)`, and index the label array as `[y + 1, x + 1]`. [`DictionaryLabelMapper`](@ref) matches numeric
keys within a tolerance, [`RangeLabelMapper`](@ref) uses non-overlapping open intervals, [`TransformLabelMapper`](@ref)
adapts a one-output transform, and [`CallableLabelMapper`](@ref) accepts arbitrary label logic.

A coordinate without a label, or with a label that has no branch, evaluates to `NaN` in every output, so
`pixel_to_world` and `world_to_pixel` treat it as having no result. [`WithInverseLabelMapper`](@ref) pairs forward and
world-side label mappers; a `RegionSelector` has an analytical inverse only when the inverse mapper and every branch
are invertible. Without an inverse mapper, `world_to_pixel` starts Newton from the first branch inverse whose result
lies in that branch's own region, then checks it against the forward selector.

Region selection and compound bounding boxes are separate: the compound box decides whether an input is in the WCS
domain, and the label mapper decides which transform applies there.

## Fixed inputs and annotations

[`fix_inputs`](@ref) holds inputs constant and keeps the result visible to composition, Jacobian propagation, metadata,
and inversion:

```@example transforms
using GWCS: AffineTransform, AnnotatedTransform, fix_inputs

annotated = AnnotatedTransform(AffineTransform([1.0 1.0; 1.0 -1.0]); inputs=("x", "y"), outputs=("u", "v"))
slice = fix_inputs(annotated, "y" => 2.0)
slice(SVector(3.0)), inverse_transform(slice)(SVector(5.0, 1.0))
```

Integer input indices are one-based, and string and symbol keys use the transform's [`input_names`](@ref).
[`annotate_transform`](@ref) attaches a name and input and output names without changing the numerical transform.
`fix_inputs(wcs, ...)` also reduces the input frame, bounding box, and pixel shape; fixing every input gives a
zero-dimensional WCS evaluated with `pixel_to_world(wcs)`. [`ConstantTransform`](@ref) adds a constant output axis
without consuming an input.

## Axis dependencies

[`axis_correlation_matrix`](@ref) returns which world axes depend on which pixel axes, and
[`transform_dependencies`](@ref) the same for a transform's outputs and inputs:

```@example transforms
using GWCS: transform_dependencies

transform_dependencies(Shift(1.0) & Scale(2.0))
```

Composition propagates these boolean matrices, parallel transforms produce block-diagonal matrices, and fixed inputs
remove the corresponding columns. FITS polynomial, lookup, and `-TAB` distortions depend only on the axes their models
read, so unrelated spectral or time axes stay separable. A transform without its own method is treated as dense.
