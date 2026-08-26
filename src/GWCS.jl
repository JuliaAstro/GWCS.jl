"""
    GWCS

Generalized world coordinate systems with composable, compiled transformations.
"""
module GWCS

using ForwardDiff
using LinearAlgebra: I, det, inv
using StaticArrays: MMatrix, MVector, SMatrix, SVector, StaticMatrix, StaticVector
using YAML

include("transforms/core.jl")
include("transforms/annotations.jl")
include("transforms/fixed.jl")
include("transforms/callable.jl")
include("transforms/elementary.jl")
include("transforms/polynomial.jl")
include("transforms/geometry.jl")
include("transforms/composition.jl")
include("transforms/dependencies.jl")
include("metadata_registry.jl")
include("frames.jl")
include("pipeline.jl")
include("api.jl")
include("numerical_inverse.jl")
include("metadata.jl")

export AbstractTransform, AffineTransform, AnnotatedTransform, ArithmeticTransform, CallableTransform, CartesianToSpherical,
       AbstractCoordinateFrame, AxisType, BoundingBox, CelestialFrame, CompositeFrame,
       CelestialToNative,
       ComposedTransform, ConstantTransform, CoordinateFrame, FixedInputsTransform,
       Frame2D, IdentityTransform,
       NativeToCelestial, NativeToTangentPlane, NoInverse,
       ParallelTransform, Polynomial1D, Polynomial2D, RemapAxes, Rotation3D,
       PIXEL, SPATIAL, SPECTRAL, STOKES, TIME, UNDEFINED,
       Scale, Shift, SpectralFrame, SphericalToCartesian, Step, StokesFrame, TangentPlaneToNative, TemporalFrame, WCS,
       WithInverse, TransformParseContext, TransformRegistry, available_frames,
       annotate_transform, axis_names, axis_physical_types, axis_types, axis_units,
       compile_transform, compose, default_transform_registry, degree, elementwise,
       axis_correlation_matrix, fix_inputs, frame_dimension, get_transform, has_inverse, input_dimension, input_frame,
       input_names, inverse_transform, materialize_metadata, metadata_path, metadata_tag,
       metadata_value, metadata_version, output_dimension, output_frame,
       offset_to_time, output_names, parallel, pixel_axis_names, pixel_axis_physical_types, pixel_axis_types,
       pixel_axis_units, pixel_n_dim, pixel_to_world, numerical_inverse, NoConvergence,
       parse_transform, parse_wcs, register_transform!, required_metadata,
       time_to_offset, transform, transform_dependencies, transform_jacobian, transform_name,
       world_axis_names, world_axis_physical_types, world_axis_types, world_axis_units, world_n_dim, world_to_pixel

end
