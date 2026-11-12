# GWCS.jl

```@meta
CurrentModule = GWCS
```

GWCS.jl implements composable generalized world coordinate systems in Julia. A WCS is a pipeline of coordinate frames
joined by transforms, from detector pixels to world coordinates such as sky positions, wavelengths, and times.
Pipelines are normalized and compiled during construction, so evaluation uses concrete, statically sized transforms
without runtime graph traversal. GWCS.jl reads WCSs from FITS headers and from the ASDF metadata Python GWCS writes,
and writes them back to ASDF.

!!! warning "Pixel coordinates are zero-based"
    For every WCS, whatever format it was read from, the center of the first pixel is pixel `0`, as in Python GWCS and
    astropy's `pixel_to_world`. Julia array indices, FITS, and WCSLIB number it `1`. `origin=1` takes and
    returns one-based pixels instead; see [Pixel origin](@ref) for the conversions.

## Installation

```julia-repl
pkg> add https://github.com/cgarling/GWCS.jl
```

FITS support for FITSIO.jl and FITSFiles.jl headers and HDUs, and ASDF reading and writing, load automatically with
those packages.

## Quick start

GWCS.jl exports what loading, evaluating, and inspecting a WCS needs: [`WCS`](@ref), [`parse_wcs`](@ref),
[`parse_fits_wcs`](@ref), [`pixel_to_world`](@ref), [`world_to_pixel`](@ref), the APE 14 axis accessors, and a few
others. Transforms, frames, bounding boxes, label mappers, pipeline editors, and the extension API are public but not
exported, so the examples import them by name.

Transforms compose left to right with `|`, and independent transforms combine with `&`:

```@example quickstart
using GWCS
using GWCS: AffineTransform, Frame2D, Shift, Step

detector = Frame2D(name="detector", axes_names=("x", "y"), unit=("pixel", "pixel"))
world = Frame2D(name="world", axes_names=("u", "v"))

shift = Shift(-10.0) & Shift(-20.0)
linear = AffineTransform([0.1 0.01; -0.02 0.2], [1.0, 2.0])
wcs = WCS(Step(detector, shift | linear), Step(world))

pixel_to_world(wcs, 11.0, 22.0)
```

```@example quickstart
world_to_pixel(wcs, 1.12, 2.38)
```

A WCS read from a file is used the same way:

```@example quickstart
header = Dict(
    "NAXIS" => 2, "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
    "CRPIX1" => 512.0, "CRPIX2" => 512.0,
    "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
    "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
)
sky = parse_fits_wcs(header)
pixel_to_world(sky, 511.0, 511.0)  # FITS pixel CRPIX = (512, 512) is GWCS pixel (511, 511)
```

## Contents

- [Building and evaluating a WCS](pipelines.md): pipelines, evaluation, editing, domains, and frames.
- [Transforms](transforms.md): the transform library, projections, region-selected and callable transforms, and inverses.
- [FITS WCS](fits.md): reading FITS image-header WCSs.
- [ASDF and metadata](asdf.md): reading and writing Python GWCS's ASDF metadata, and registering new transforms.
- [Conventions and differences](conventions.md): coordinate order and the deliberate differences from Python GWCS.
- [Python equivalents](python_equivalents.md) and the [API reference](api/wcs.md).
