# GWCS.jl

[![CI](https://github.com/cgarling/GWCS.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/cgarling/GWCS.jl/actions/workflows/Test.yml?query=branch%3Amain)

GWCS.jl implements generalized world coordinate systems in Julia. A WCS is a pipeline of coordinate frames joined by
transforms, from detector pixels to world coordinates such as sky positions, wavelengths, and times. Pipelines are
compiled when they are built, so evaluating a coordinate is fast and allocation-free. GWCS.jl reads the WCSs in FITS
headers and in the ASDF metadata Python GWCS writes, and writes the pipelines ASDF can represent back to ASDF.

## Installation

```julia-repl
pkg> add https://github.com/cgarling/GWCS.jl
```

FITS and ASDF support load automatically with FITSIO.jl, FITSFiles.jl, or ASDF.jl.

## Conventions

- **Pixels are zero-based**: the center of the first pixel is `0`, as in Python GWCS and astropy. `origin=1` takes and
  returns one-based pixels, as Julia arrays and FITS number them.
- **Coordinates are in transform order**: a frame's `axes_order` labels its axes and never reorders them.
- **A coordinate without a result** gives `fill_value` (default `NaN`), or a `DomainError` with `strict=true`. A
  numerical inverse that does not converge throws `NoConvergence`.
- **Angles are in degrees**: angle parameters and the coordinates that projections and celestial rotations take and
  return are in degrees. Coordinates are in the units their frame records, which `world_axis_units` returns.
  `parse_fits_wcs` converts celestial axes to degrees, or keeps the header's `CUNIT` with `preserve_units=true`, and
  an ASDF frame keeps its stored units, such as arcseconds for a V2/V3 frame.
- **Exports**: GWCS.jl exports what loading and evaluating a WCS needs. Transforms, frames, and the other building
  blocks are public but not exported; import them by name, as in `using GWCS: Shift`.

The documentation describes these conventions and the differences from Python GWCS in detail.

## Usage

From a FITS header as a dictionary:

```julia
using GWCS

header = Dict("CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CRPIX1" => 512.0, "CRPIX2" => 512.0,
    "CRVAL1" => 83.8221, "CRVAL2" => -5.3911, "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4)
wcs = parse_fits_wcs(header)
world = pixel_to_world(wcs, 511.0, 511.0)
pixel = world_to_pixel(wcs, world)
```

From an ASDF file, and back to one:

```julia
using ASDF, GWCS

wcs = parse_wcs(ASDF.load("roman_cal.asdf"; extensions=true)["roman"]["meta"]["wcs"])
write_asdf("wcs.asdf", wcs)
```

From a FITS file with FITSIO.jl:

```julia
using FITSIO, GWCS

wcs = FITS("image.fits") do file
    parse_fits_wcs(file[1])
end
```

From a FITS file with FITSFiles.jl:

```julia
using FITSFiles, GWCS

wcs = parse_fits_wcs(FITSFiles.fits("image.fits"), 1)
```

## Features

- Composable pipelines (`|`, `&`, and elementwise arithmetic) compiled for inferred, allocation-free evaluation, with
  matrix forms for many coordinates at once
- Analytical inverses, and a damped Newton numerical inverse that first inverts the analytically invertible part of
  the pipeline, as astropy's `all_world2pix` does
- Jacobians in both directions, and pixel solid angles for single pixels or whole images
- Shift, scale, affine, axis remapping, constant, fixed-input, polynomial, 2D and 3D rotation, and spherical-Cartesian
  transforms, and all 28 FITS spherical projections
- Region-selected transforms for IFU slices, detector regions, and spectral orders, with array, dictionary, range,
  transform, and callable label mappers
- Bounding boxes and compound bounding boxes, with `in_domain`, `in_image`, and grid sampling
- Celestial, spectral, temporal, Stokes, composite, and generic frames, with the APE 14 axis accessors
- Non-mutating pipeline editing: insert, remove, rename, and replace frames and transforms
- FITS WCS reading: PC, CD, and CROTA linear forms; SIP, TPV, TPD, TNX, and ZPX distortions; D2IM and CPDIS lookup
  tables; `-TAB` coordinate arrays; spectral algorithms and grisms; time axes; and repairs for legacy conventions
- Reading and writing Python GWCS's ASDF metadata, keeping properties GWCS.jl does not model
- An extension API for new transforms, including registering their metadata
- Callable transforms for prototypes and closures
- Validated against Python GWCS, astropy, and WCSLIB

## Acknowledgements

The design of GWCS.jl was inspired by [Python GWCS](https://github.com/spacetelescope/gwcs). The development of the FITS WCS features was aided by validation (see validation/) against [astropy](https://www.astropy.org) and [WCSLIB](https://www.atnf.csiro.au/computing/software/wcs/). The sample headers in `test/fixtures/wcslib` are from WCSLIB and are
distributed under the LGPL v3; see `test/fixtures/wcslib/COPYING.LESSER`.
