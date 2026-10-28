"""
    GWCS

Generalized world coordinate systems with composable, compiled transformations.
"""
module GWCS

using ForwardDiff
using LinearAlgebra: Diagonal, I, det, inv
using PrecompileTools: @compile_workload, @setup_workload
using StaticArrays: MMatrix, MVector, SMatrix, SVector, StaticMatrix, StaticVector
using YAML

include("transforms/core.jl")
include("transforms/annotations.jl")
include("transforms/fixed.jl")
include("transforms/callable.jl")
include("transforms/elementary.jl")
include("selectors.jl")
include("transforms/polynomial.jl")
include("transforms/interpolation.jl")
include("transforms/spectral.jl")
include("transforms/geometry.jl")
include("transforms/projections.jl")
include("transforms/composition.jl")
include("transforms/selectors.jl")
include("transforms/dependencies.jl")
include("metadata_registry.jl")
include("metadata_schemas.jl")
include("frames.jl")
include("domains.jl")
include("pipeline.jl")
include("display.jl")
include("fits.jl")
include("api.jl")
include("numerical_inverse.jl")
include("metadata.jl")
include("serialization.jl")
include("precompile.jl")

export AbstractTransform, AffineTransform, AnnotatedTransform, ArithmeticTransform, ArrayLabelMapper, CallableLabelMapper, CallableTransform, CartesianToSpherical,
       AbstractConicProjection, AbstractSphericalProjection, AbstractZenithalProjection, AiryProjection, BonneProjection, COBESphericalCube, ConicEqualArea, ConicEquidistant, ConicOrthomorphic, ConicPerspective, CylindricalEqualArea, CylindricalPerspective,
       AbstractBoundingBox, AbstractCoordinateFrame, AxisType, BoundingBox, CelestialFrame, CompoundBoundingBox, CompositeFrame, EmptyFrame,
       CelestialToNative,
       ComposedTransform, ConstantTransform, CoordinateFrame, DictionaryLabelMapper, FixedInputsTransform,
       BilinearLookupTable, FITSGrismTransform, FITSLogTransform, FITSLookupDistortion, FITSSpectralTransform, FITSTabularTable, FITSTabularTransform, Frame2D, GnomonicProjection, HammerAitoffProjection, IdentityTransform, IntermediateToNative,
       NativeToCelestial, NativeToIntermediate, NoInverse,
       HEALPixPolarProjection, HEALPixProjection, MercatorProjection, MollweideProjection, ParabolicProjection, ParallelTransform, PlateCarreeProjection, PolyconicProjection, Polynomial1D, Polynomial2D, QuadrilateralizedSphericalCube, RangeLabelMapper, RegionSelector, RemapAxes, Rotation2D, Rotation3D,
       PIXEL, SPATIAL, SPECTRAL, STOKES, TIME, UNDEFINED,
       SansonFlamsteedProjection, Scale, Shift, SIPTransform, SlantOrthographic, SlantZenithalPerspective, SpectralFrame, SphericalToCartesian, StereographicProjection, Step, StokesFrame, TangentialSphericalCube,
       TemporalFrame, TPVTransform, WCS, WCSBuilder, ZenithalEqualArea, ZenithalEquidistant, ZenithalPerspective, ZenithalPolynomial,
       TransformLabelMapper,        WithInverse, WithInverseLabelMapper, TransformParseContext, TransformRegistry, TransformWriteContext, available_frames,
       annotate_transform, axis_names, axis_physical_types, axis_types, axis_units, bounded_axes, bounds,
       celestial_transform, compile_label_mapper, compile_transform, compose, default_lonpole, default_transform_registry, degree, elementwise,
       append_frame, append_frame!, append_transform, append_transform!, axis_correlation_matrix, fix_inputs, frame_dimension, get_frame, get_transform,
       grid_from_bounding_box, has_inverse, ignored_axes, input_dimension, input_frame, insert_frame, insert_frame!, insert_transform, insert_transform!,
       in_domain, in_image, input_names, inverse_label_mapper, inverse_transform, label_dependencies, label_mapper_name, label_type, map_label, materialize_metadata, metadata_path, metadata_tag,
       metadata_value, metadata_version, output_dimension, output_frame,
       native_fiducial, offset_to_time, output_names, parallel, pixel_axis_names, pixel_axis_physical_types, pixel_axis_types,
       pixel_area, pixel_axis_units, pixel_n_dim, pixel_to_world, pixel_to_world_jacobian, numerical_inverse, NoConvergence,
       parse_fits_wcs, find_all_fits_wcs, parse_transform, parse_wcs, register_transform!, required_metadata, select_bounding_box, selected_transform, selector_axes, selector_ignored,
       register_transform_serializer!,
       prepend_frame, prepend_frame!, prepend_transform, prepend_transform!, remove_frame, remove_frame!, rename_frame, rename_frame!,
       set_transform, set_transform!, time_to_offset, transform, transform_dependencies, transform_jacobian, transform_name,
       to_metadata, with_bounding_box, with_pixel_shape, world_axis_names, world_axis_physical_types, world_axis_types, world_axis_units, world_n_dim, world_to_pixel, world_to_pixel_jacobian, write_asdf

end
