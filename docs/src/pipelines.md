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
