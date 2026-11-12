"""
    GWCS

Generalized world coordinate systems with composable, compiled transformations.
"""
module GWCS

using ForwardDiff
using LinearAlgebra: Diagonal, I, cross, det, dot, inv, norm
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

export WCS, parse_wcs, parse_fits_wcs, find_all_fits_wcs, write_asdf, pixel_to_world, world_to_pixel, pixel_to_world_jacobian,
       world_to_pixel_jacobian, in_domain, in_image, pixel_area, axis_correlation_matrix, pixel_n_dim, world_n_dim, pixel_axis_names,
       pixel_axis_types, pixel_axis_units, pixel_axis_physical_types, world_axis_names, world_axis_types, world_axis_units,
       world_axis_physical_types, input_frame, output_frame, available_frames, with_bounding_box, with_pixel_shape

# Public API that is not exported; `public` needs Julia 1.11, so it is parsed only there.
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("""
        public AbstractCoordinateFrame, AbstractTransform, AffineTransform, AiryProjection, AnnotatedTransform, ArrayLabelMapper, AxisType,
               BonneProjection, BoundingBox, COBESphericalCube, CallableLabelMapper, CallableTransform, CartesianToSpherical, CelestialFrame,
               CelestialToNative, CompositeFrame, CompoundBoundingBox, ConicEqualArea, ConicEquidistant, ConicOrthomorphic, ConicPerspective,
               ConstantTransform, CoordinateFrame, CylindricalEqualArea, CylindricalPerspective, DictionaryLabelMapper, EmptyFrame, Frame2D,
               GnomonicProjection, HEALPixPolarProjection, HEALPixProjection, HammerAitoffProjection, IdentityTransform, IntermediateToNative,
               MercatorProjection, MollweideProjection, NativeToCelestial, NativeToIntermediate, NoConvergence, NoInverse, PIXEL,
               ParabolicProjection, PlateCarreeProjection, PolyconicProjection, Polynomial1D, Polynomial2D, QuadrilateralizedSphericalCube,
               RangeLabelMapper, RegionSelector, RemapAxes, Rotation2D, Rotation3D, SPATIAL, SPECTRAL, STOKES, SansonFlamsteedProjection, Scale,
               Shift, SlantOrthographic, SlantZenithalPerspective, SpectralFrame, SphericalToCartesian, Step, StereographicProjection, StokesFrame,
               TIME, TangentialSphericalCube, TemporalFrame, TransformLabelMapper, TransformParseContext, TransformRegistry, TransformWriteContext,
               UNDEFINED, WithInverse, WithInverseLabelMapper, ZenithalEqualArea, ZenithalEquidistant, ZenithalPerspective, ZenithalPolynomial,
               annotate_transform, append_frame, append_transform, axis_names, axis_physical_types, axis_types, axis_units, bounded_axes, bounds,
               celestial_transform, compile_transform, compose, default_transform_registry, degree, elementwise, fix_inputs, frame_dimension,
               get_frame, get_transform, grid_from_bounding_box, has_inverse, ignored_axes, input_dimension, input_names, insert_frame,
               insert_transform, inverse_transform, map_label, materialize_metadata, metadata_path, metadata_tag, metadata_value, metadata_version,
               numerical_inverse, offset_to_time, output_dimension, output_names, parallel, parse_transform, prepend_frame, prepend_transform,
               register_transform!, register_transform_serializer!, remove_frame, rename_frame, required_metadata, select_bounding_box,
               selected_transform, selector_axes, selector_ignored, set_transform, time_to_offset, to_metadata, transform, transform_dependencies,
               transform_jacobian, transform_name
        """))
end

end
