# WCS

```@meta
CurrentModule = GWCS
```

GWCS.jl exports the names in the first three sections of this page. The pipeline functions after them, and the names
on the other API pages, are public but not exported: use them as `GWCS.Step`, or import them with `using GWCS: Step`.

```@docs
GWCS
```

## Reading and writing

```@docs
WCS
parse_wcs
parse_fits_wcs
find_all_fits_wcs
write_asdf
```

## Evaluation

```@docs
pixel_to_world
world_to_pixel
pixel_to_world_jacobian
world_to_pixel_jacobian
in_domain
in_image
pixel_area
```

## WCS properties

```@docs
pixel_n_dim
world_n_dim
pixel_axis_names
world_axis_names
pixel_axis_physical_types
world_axis_physical_types
pixel_axis_units
world_axis_units
pixel_axis_types
world_axis_types
axis_correlation_matrix
input_frame
output_frame
available_frames
with_bounding_box
with_pixel_shape
```

## Pipelines

```@docs
Step
get_frame
get_transform
transform
set_transform
insert_transform
prepend_transform
append_transform
insert_frame
prepend_frame
append_frame
remove_frame
rename_frame
fix_inputs
numerical_inverse
NoConvergence
```
