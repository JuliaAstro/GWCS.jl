# Conventions and differences

```@meta
CurrentModule = GWCS
```

## Pixel origin

For every WCS, whatever format it was read from, the center of the first pixel is pixel `0`. This matches Python GWCS
and astropy's `pixel_to_world`/`world_to_pixel` (APE 14). It does **not** match Julia array indices, the FITS standard,
or WCSLIB, which number the first pixel center `1`.

- Julia array index `(i, j)`, with `i` along the first pixel axis, is GWCS pixel `(i - 1, j - 1)`. A GWCS pixel `p`
  falls in array index `round.(Int, p) .+ 1`.
- A FITS pixel `p` (WCSLIB, DS9, SourceExtractor `X_IMAGE`) is GWCS pixel `p .- 1`.
  [`parse_fits_wcs`](@ref) applies this shift itself, so `CRPIX` needs no adjustment.
- Pixels from Python GWCS, astropy's APE 14 interface, and photutils are already zero-based.

To pass and receive one-based pixels instead, give `origin=1` to [`pixel_to_world`](@ref), [`world_to_pixel`](@ref),
[`numerical_inverse`](@ref), [`in_domain`](@ref), [`in_image`](@ref) (for `initial_guess`), or
[`grid_from_bounding_box`](@ref) (for `fixed` values and the grid). As with the `origin` argument of astropy.wcs's
`all_pix2world`, it shifts **every** input axis by one and never shifts world coordinates or `fill_value`. Do not use it
for a WCS with an input that is not an array position, such as a grism spectral order or a selector label.

## Coordinates

- **Coordinates are in transform order**, the order of a transform's inputs and outputs. A frame stores its axes in its
  own order, and `axes_order[i]` is the one-based transform position of its axis `i`; it labels coordinates and never
  reorders them. The axis accessors, bounding boxes, `pixel_shape`, and [`fix_inputs`](@ref) indices use transform
  order, while frame fields and constructor keywords use frame order. This is Python GWCS's convention.
- **Matrices hold one coordinate per column** in the evaluation functions.
- **A coordinate without a result gives `NaN`** in every component of a transform's output, and `fill_value` from
  [`pixel_to_world`](@ref) and [`world_to_pixel`](@ref), which throw a `DomainError` instead with `strict=true`.
  Transforms never throw for such coordinates, except where a zero-dimensional result cannot carry `NaN`.
- **Labels are not coordinates.** A label mapper returns a label or `nothing` for no match, so `0` and `""` remain
  ordinary labels.
- **Angles are in degrees.** Angle parameters, such as rotation angles, projection parameters and reference points, and
  the coordinates the projections and celestial rotations take and return, are in degrees. Frames record each axis's
  unit but evaluation does not convert units; Python GWCS can carry astropy units through evaluation instead.

## Deliberate differences from Python GWCS

- **Region misses.** A [`RegionSelector`](@ref) evaluates to `NaN` for a coordinate in no region, so the evaluation
  functions report it as having no result. Python GWCS returns its `undefined_transform_value`; GWCS.jl stores that
  value only for ASDF round trips. Pass the same number as `fill_value` to reproduce Python's output. As in Python GWCS,
  an inverse selector does not inherit a custom undefined value.
- **Region-selector numerical inverse.** Without an inverse label mapper, `world_to_pixel` starts Newton from the first
  branch inverse whose result lies in that branch's own region. Python GWCS starts from a fixed point that can fall in no
  region.
- **Linear polynomial distortions.** `RemapAxes{2}((1, 2, 1, 2)) | (px & py)` with polynomials of degree at most one
  becomes an exact affine with an analytical inverse; Python GWCS inverts it numerically.
- **Frame-to-frame evaluation.** [`GWCS.transform`](@ref transform) applies no bounding box or `fill_value`, while
  Python's `WCS.transform` applies the bounding box by default.
- **Values, not objects.** [`pixel_to_world`](@ref) and [`world_to_pixel`](@ref) take and return numbers, like Python's
  `pixel_to_world_values` and `world_to_pixel_values`; Python's `pixel_to_world` returns astropy coordinate objects.
- **Immutable WCSs.** Editing functions such as [`set_transform`](@ref) return a new WCS, while Python GWCS changes its
  WCS in place.
- **FITS time axes** are in seconds from the reference epoch, whatever `TIMEUNIT` the header uses.

The [Python equivalents](@ref) page lists the Python counterpart of every name that has one.
