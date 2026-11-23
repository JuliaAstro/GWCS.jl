# Building and evaluating a WCS

```@meta
CurrentModule = GWCS
```

A [`WCS`](@ref) is a sequence of [`Step`](@ref)s, each a coordinate frame and the transform from it to the next
frame. The first frame holds pixel coordinates and the last world coordinates:

```@example pipeline
using GWCS
using GWCS: AffineTransform, Frame2D, Shift, Step

detector = Frame2D(name="detector", axes_names=("x", "y"), unit=("pixel", "pixel"))
world = Frame2D(name="world", axes_names=("u", "v"))

shift = Shift(-10.0) & Shift(-20.0)
linear = AffineTransform([0.1 0.01; -0.02 0.2], [1.0, 2.0])
wcs = WCS(Step(detector, shift | linear), Step(world))
```

As in Python GWCS, coordinates are in transform order, and a frame's `axes_order` only labels them: a frame stores its
axes in its own order (a `CelestialFrame`'s is longitude, latitude), and `axes_order[i]` is the one-based transform
position of the frame's axis `i`. The axis accessors ([`axis_names`](@ref), [`pixel_axis_names`](@ref),
[`world_axis_physical_types`](@ref), and so on) report the axes in transform order, as do bounding boxes,
`pixel_shape`, and [`fix_inputs`](@ref) indices. Pixel coordinates are zero-based: parsers encode each format's pixel
origin explicitly at the boundary, and a pipeline built by hand is evaluated exactly as written, with no implicit
offset.

## Evaluation

[`pixel_to_world`](@ref) and [`world_to_pixel`](@ref) accept separate numbers, tuples, vectors, and matrices whose
columns are coordinates:

```@example pipeline
pixel_matrix = [11.0 12.0 13.0; 22.0 22.0 22.0]
pixel_to_world(wcs, pixel_matrix)
```

A coordinate without a result returns `fill_value` (default `NaN`) in every output component. That covers coordinates
outside the bounding box, which is enforced by default when present (pass `with_bounding_box=false` to evaluate
outside it), and coordinates the transforms cannot map, such as the corners of an all-sky projection, the far side of
a TAN projection, or a detector pixel between IFU slices. As with APE 14 WCS objects in Python, whole image grids can
therefore be evaluated directly, and invalid points cost no more than valid ones.

With `strict=true`, any coordinate without a result throws a `DomainError`. Here a point off an all-sky Hammer-Aitoff
projection has no world coordinate:

```@example pipeline
using GWCS: HammerAitoffProjection, IntermediateToNative

allsky = WCS(Step(Frame2D(name="intermediate"), IntermediateToNative(HammerAitoffProjection())), Step(Frame2D(name="native")))
pixel_to_world(allsky, 200.0, 0.0)
```

```@repl pipeline
pixel_to_world(allsky, 200.0, 0.0; strict=true)
```

For an undefined coordinate, the message names the coordinate, the pipeline step, and the transform that is undefined
there. [`NoConvergence`](@ref) from the numerical inverse (controlled by `quiet`) always propagates. Internally every transform returns `NaN` where a coordinate has no image
rather than throwing, so the default and strict modes share one compiled pipeline; the step-by-step diagnosis runs
only when a strict error is reported.

[`pixel_to_world_jacobian`](@ref) returns the Jacobian of `pixel_to_world`, the `world_n_dim × pixel_n_dim` `SMatrix`
of world-coordinate derivatives with respect to the pixel, and [`world_to_pixel_jacobian`](@ref) that of
`world_to_pixel`. Both evaluate the compiled pipeline in one allocation-free pass. A zenithal projection and the sky
rotation after it are differentiated together, so the Jacobian is finite and accurate at and near the reference point,
where the native longitude is undefined.

[`pixel_area`](@ref) returns the solid angle in steradians that a pixel covers on the sky, computed in one pass with
the world coordinate. Given ranges of integer pixel coordinates, `pixel_area(wcs, xs, ys)` returns the areas of the
grid `xs × ys` as a matrix indexed like `data[x, y]`, each from the pixel's four corners, which neighboring pixels
share, so a full Roman detector map takes about 0.2 s.

## Editing a pipeline

`WCS` values are immutable. Editing functions return a newly validated and compiled WCS, leaving the source unchanged.
This splits the detector-to-world step at a new frame: the new edge gets `shift`, and the next edge keeps
`shift | linear`, so [`set_transform`](@ref) restores the original mapping.

```@example pipeline
using GWCS: Scale, insert_frame, insert_transform, remove_frame, rename_frame, set_transform

split = insert_frame(wcs, "detector", shift, Frame2D(name="centered"))
split = set_transform(split, "centered", "world", linear)
available_frames(split)
```

```@example pipeline
pixel_to_world(split, 11.0, 22.0)  # as for wcs
```

```@example pipeline
recalibrated = insert_transform(split, "world", Scale(1.01) & Scale(1.01))
renamed = rename_frame(recalibrated, "centered", "focal")
merged = remove_frame(renamed, "focal")
available_frames(merged), pixel_to_world(merged, 11.0, 22.0)
```

[`set_transform`](@ref) requires adjacent frames. [`insert_transform`](@ref) composes before the named frame by
default or after it with `after=true`; it changes coordinate values without changing frame metadata.
[`insert_frame`](@ref) requires exactly one endpoint to already exist and injects a new edge, so it does not preserve
the overall mapping unless the added transform is an identity. [`remove_frame`](@ref) removes an intermediate frame
while composing its adjacent transforms. [`prepend_transform`](@ref), [`append_transform`](@ref),
[`prepend_frame`](@ref), and [`append_frame`](@ref) are conveniences for the ends of the pipeline. [`get_frame`](@ref)
retrieves a frame stored in the pipeline, and [`get_transform`](@ref) the compiled transform between any two frames,
which [`GWCS.transform`](@ref transform) evaluates. Unlike `pixel_to_world`, `GWCS.transform` applies no bounding box
and no `fill_value`.

Each edit validates and compiles a complete WCS before returning it, so a failed edit throws rather than returning a
partial pipeline. Several edits chain by passing each result to the next, as above. Edits that keep the input frame
keep `pixel_shape` and `bounding_box`. Prepending a new input frame requires explicit replacements, or explicit
`nothing` values to clear existing input-domain metadata, because GWCS.jl cannot infer that metadata in the new
coordinate system.

## Domains and bounding boxes

A WCS's input domain is the pixels it is defined on. A [`BoundingBox`](@ref) limits it to an inclusive interval on
each pixel axis, and [`with_bounding_box`](@ref) returns a copy of a WCS with one:

```@example pipeline
using GWCS: BoundingBox

boxed = with_bounding_box(wcs, BoundingBox(((0.0, 99.0), (0.0, 49.0))))
pixel_to_world(boxed, 200.0, 0.0)  # outside the box
```

[`in_domain`](@ref) answers whether a pixel is in the domain, and [`in_image`](@ref) whether a world coordinate
inverts to a pixel in it:

```@example pipeline
in_domain(boxed, 11.0, 22.0), in_domain(boxed, 200.0, 0.0), in_image(boxed, pixel_to_world(wcs, 200.0, 0.0))
```

`with_bounding_box=false` evaluates or tests outside the box. `in_image` takes the numerical-inverse controls of
`world_to_pixel`, and a solve that does not converge gives `false`. [`with_pixel_shape`](@ref) sets the image size in
the same way.

[`grid_from_bounding_box`](@ref) returns the pixel centers in a box, one per column, ready for `pixel_to_world`:

```@example pipeline
using GWCS: grid_from_bounding_box

grid = with_bounding_box(wcs, BoundingBox(((0.0, 3.0), (0.0, 1.0))))
grid_from_bounding_box(grid)
```

### Unbounded and selected axes

The dimensioned constructor bounds only some axes and leaves the others unconstrained, though a point in the domain
must still be finite:

```@example pipeline
using GWCS: bounded_axes, ignored_axes

box = BoundingBox{3}(((0.0, 4095.0), (0.0, 4095.0)), (1, 2))
bounded_axes(box), ignored_axes(box)
```

A [`CompoundBoundingBox`](@ref) selects a box by the exact value of one or more inputs. It represents detector regions,
slices, and spectral orders whose valid spatial bounds differ:

```@example pipeline
using GWCS: CompoundBoundingBox, select_bounding_box

boxes = Dict(
    (1,) => BoundingBox{3}(((0.0, 2047.0), (0.0, 2047.0)), (1, 2)),
    (2,) => BoundingBox{3}(((2048.0, 4095.0), (0.0, 2047.0)), (1, 2)),
)
domain = CompoundBoundingBox(boxes, (3,))
select_bounding_box(domain, 2)
```

A coordinate whose selector value has no box is outside the domain. A grid of a compound box needs the `selector` of
the box to sample, and `fixed` values for ignored axes. `fix_inputs(wcs, ...)` reduces simple and compound domains,
selecting and collapsing a compound box when its selector inputs are fixed. Without another starting point, numerical
inversion starts from the center of a simple box or of a compound box with one entry; a compound box with several
needs an explicit `initial_guess`, because a world coordinate does not identify a selector branch before solving.

These tests use the input domain, not an approximation of the output footprint.

## Coordinate frames

Frames classify axes with the [`AxisType`](@ref) values `SPATIAL`, `SPECTRAL`, `TIME`, `STOKES`, `PIXEL`, and
`UNDEFINED`; generic frames may also use custom strings. [`SpectralFrame`](@ref), [`TemporalFrame`](@ref), and
[`StokesFrame`](@ref) are one-dimensional, and [`CompositeFrame`](@ref) combines frames without changing the compiled
numeric kernel:

```@example pipeline
using GWCS: CelestialFrame, CompositeFrame, SpectralFrame, TemporalFrame, axis_names

sky = CelestialFrame(name="sky", axes_order=(2, 4))
wave = SpectralFrame(name="wave", axes_order=1, axes_names="lambda", unit="nm")
time = TemporalFrame(name="time", axes_order=3, unit="s", reference_epoch=2460000.5, time_format="jd")
cube = CompositeFrame((sky, wave, time); name="world")
axis_names(cube)
```

Constituent `axes_order` values are positions in the whole composite and must form a complete permutation across it,
so the coordinates above are (lambda, lon, time, lat). [`axis_names`](@ref), [`axis_physical_types`](@ref),
[`axis_units`](@ref), and [`axis_types`](@ref) report frame metadata in that order, as do the corresponding
`pixel_axis_*` and `world_axis_*` WCS accessors. ASDF frames list their per-axis metadata in frame order, as Python
GWCS reads and writes them.

[`EmptyFrame`](@ref) is a frame known only by its name, with no axis metadata. Python GWCS writes such a step frame as a
bare string, and GWCS.jl reads and writes it the same way.

`SpectralFrame` infers common wavelength, frequency, energy, and velocity physical types from unit strings.
`TemporalFrame` keeps reference epoch, scale, format, and location metadata. [`offset_to_time`](@ref) and
[`time_to_offset`](@ref) handle numeric JD and MJD origins for common time units and delegate ordinary arithmetic to
compatible external time values; other numeric time formats require an optional time integration. Rich units,
celestial coordinates, and astronomy time objects remain concerns of optional extensions.

### Celestial systems

[`celestial_frame`](@ref) finds a WCS's celestial frame, and [`celestial_system`](@ref) reports its system, such as
ICRS, FK5 with its equinox, or Galactic. [`with_celestial_system`](@ref) replaces the output frame with one in another
system, rotating the celestial coordinates appropriately, so the WCS gives world coordinates in that system
and `world_to_pixel` takes them. The rotation is folded into the WCS's final sky rotation where it has one,
so it adds nothing to evaluation. `keep_frame=true` keeps the original
frame as an intermediate one instead:

```@example pipeline
using GWCS: celestial_system, celestial_frame, with_celestial_system

galactic = parse_fits_wcs(Dict("CTYPE1" => "GLON-TAN", "CTYPE2" => "GLAT-TAN", "CRPIX1" => 50.0, "CRPIX2" => 50.0,
    "CRVAL1" => 120.0, "CRVAL2" => 20.0, "CDELT1" => -1e-3, "CDELT2" => 1e-3))
icrs = with_celestial_system(galactic, :icrs)
celestial_system(celestial_frame(icrs)), pixel_to_world(icrs, 49.0, 49.0)
```

ICRS, FK5, FK4, and FK4-NO-E at any equinox, Galactic, Supergalactic, and the FITS ecliptic systems can be converted.

#### Conversion details

The conversions use astropy's definitions of each system, validated against astropy, except where the systems are
formally defined otherwise:

| Between | Model |
|---|---|
| ICRS and FK5 J2000 | The frame bias of USNO Circular 179 |
| FK5 equinoxes | IAU 1976 precession (Lieske et al. 1977; ERFA's `pmat76`), the FK5 catalog's |
| FK4 equinoxes | Newcomb's precession |
| FK4-NO-E B1950 and FK5 J2000 | Murray's (1989) rotation, with its correction for FK4's rotation at the obstime; every conversion from or to FK4 passes through it |
| FK4 and FK4-NO-E | The E-terms of aberration at the FK4 equinox, removed and added as SLALIB's `sla_SUBET` and `sla_ADDET` do |
| Galactic and FK5 J2000 | astropy's J2000 north Galactic pole |
| Supergalactic and Galactic | astropy's north supergalactic pole |
| FK5 ecliptic and FK5 at its equinox | The IAU 1980 mean obliquity (ERFA's `obl80`, SLALIB's `sla_ECMAT`) |
| FK4 ecliptic and FK4-NO-E at its equinox | Newcomb's mean obliquity, as SLALIB's `sla_ETRMS` gives it |
| ICRS ecliptic and ICRS | The IAU 2006 mean ecliptic and equinox of J2000 (ERFA's `ecm06`), astropy's `BarycentricMeanEcliptic` |

An FK4 obstime is the time of observation, which Murray's correction needs since FK4's axes rotate slowly against
FK5's. [`parse_fits_wcs`](@ref) takes it from `MJD-OBS`, else `DATE-OBS`, as FITS Paper II (Calabretta & Greisen
2002, sect. 3.1) specifies, and otherwise as the equinox. A FITS ecliptic (`ELON`/`ELAT`) is the mean ecliptic and
equinox of `EQUINOX` in the `RADESYS` system, as Paper II defines it; the FK5 conversion reproduces its example 2 to about 1 mas. An FK4
ecliptic has no E-terms, and Paper II does not apply `RADESYS = 'FK4-NO-E'` to ecliptics, so one is read as FK4's.
ICRS defines no ecliptic, so converting from or to `:ecliptic_icrs` warns that it takes astropy's J2000 mean ecliptic.
Apparent (`RADESYS = 'GAPPT'`) and helioecliptic (`HLON`/`HLAT`) coordinates depend on the time of observation, and
are not converted. astropy has no frame for the FK5 and FK4 ecliptics.

The results differ from astropy's in six ways:

- **FK5 equinoxes other than J2000** are precessed with IAU 1976, the model that defines FK5. astropy precesses them
  with IAU 2006, which differs by about 75 mas at J1975 and J2025 and 150 mas at J2050. A WCS written to ASDF in
  such a frame, and a FITS header with such an equinox, are read the same way by Python GWCS, whose numbers come from
  the pipeline, but astropy's `SkyCoord` conversions from them use IAU 2006.
- **FK4 obstimes from FITS headers.** astropy's `wcs_to_celestial_frame` ignores `MJD-OBS` and `DATE-OBS` for FK4 and
  takes the obstime as the equinox. For a typical archival FK4 image, observed decades after 1950, Murray's correction
  then differs by about 0.1″.
- **FK4 and Galactic coordinates.** Every conversion from or to FK4 passes through FK5 J2000 with Murray's correction
  at the obstime, so its result does not depend on the route. astropy converts FK4 to Galactic coordinates through the
  B1950 north Galactic pole without the obstime, and between FK4 systems by precession alone, so its results through
  different intermediate systems differ. The two agree at an obstime of B1950, and differ by up to about 0.25″ at
  obstimes near B2000.
- **FK4 to FK5.** The conversion removes the E-terms, precesses to B1950, applies Murray's rotation, and precesses
  within FK5, the path astropy's direct FK4-NO-E-to-FK5 transform takes. astropy's `transform_to` from FK4 passes
  instead through an FK4-NO-E frame at the FK5 equinox, because its frame graph gives intermediate frames the target's
  attributes; Newcomb's precession there and back is not exactly the identity, which shifts its result by up to about
  60 µas for equinoxes near J2025.
- **Murray's rotation** is made exactly orthogonal. Its linear obstime correction leaves it up to about 1e-12 rad
  (0.3 µas) from orthogonal at obstimes near B2000, so a conversion keeps an exact inverse.
- **FK5 to FK4** is the exact inverse of FK4 to FK5. astropy evaluates Newcomb's precession in the reverse direction
  instead, which differs by tens of µas for FK4 equinoxes far from B1950.
