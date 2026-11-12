# Frames and domains

```@meta
CurrentModule = GWCS
```

## Frames

```@docs
AbstractCoordinateFrame
CoordinateFrame
Frame2D
CelestialFrame
SpectralFrame
TemporalFrame
StokesFrame
CompositeFrame
EmptyFrame
AxisType
SPATIAL
SPECTRAL
TIME
STOKES
PIXEL
UNDEFINED
frame_dimension
axis_names
axis_physical_types
axis_units
axis_types
offset_to_time
time_to_offset
```

## Bounding boxes

```@docs
BoundingBox
CompoundBoundingBox
bounded_axes
ignored_axes
bounds
selector_axes
selector_ignored
select_bounding_box
grid_from_bounding_box
```

## Region selection

```@docs
RegionSelector
selected_transform
map_label
ArrayLabelMapper
DictionaryLabelMapper
RangeLabelMapper
TransformLabelMapper
CallableLabelMapper
WithInverseLabelMapper
```
