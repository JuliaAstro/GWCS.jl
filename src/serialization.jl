const _UNIT_TAG = "tag:stsci.edu:asdf/unit/unit"
const _TIME_TAG = "tag:stsci.edu:asdf/time/time"
const _ICRS_TAG = "tag:astropy.org:astropy/coordinates/frames/icrs"

_tagged_metadata(tag, version, properties) = TaggedMapping("$tag-$version", Dict{Any,Any}(properties))
_metadata_rows(matrix) = [collect(row) for row in eachrow(matrix)]
_metadata_vector(values) = collect(values)

function _metadata_mapping(node, path)
    node isa TaggedMapping || _metadata_error(path, "serializer produced $(typeof(node)), expected tagged mapping")
    return node
end

function _decorate_metadata(node::TaggedMapping, properties)
    values = Dict{Any,Any}(node)
    for (key, value) in properties
        value === nothing ? delete!(values, key) : (values[key] = value)
    end
    return TaggedMapping(node.tag, values)
end

function _metadata_shorthand_tag(tag::AbstractString)
    prefix = "tag:stsci.edu:asdf/"
    return startswith(tag, prefix) ? "!" * chopprefix(tag, prefix) : "!<" * tag * ">"
end

function YAML._print(io::IO, value::Union{TaggedMapping,TaggedSequence}, level::Int=0, ignore_level::Bool=false)
    if isempty(value.value)
        print(io, _metadata_shorthand_tag(value.tag), " ")
        return YAML._print(io, value.value, level, ignore_level)
    end
    print(io, YAML._indent(_metadata_shorthand_tag(value.tag) * "\n", level, ignore_level))
    return YAML._print(io, value.value, level)
end

function YAML._print(io::IO, value::TaggedScalar, level::Int=0, ignore_level::Bool=false)
    print(io, _metadata_shorthand_tag(value.tag), " ")
    return YAML._print(io, value.value, level, ignore_level)
end

"""
    to_metadata(value; registry=default_transform_registry())

Return a canonical tagged metadata tree for a supported transform, frame,
domain, pipeline, or nested document value.
"""
to_metadata(value; registry=default_transform_registry()) = _to_metadata(value, "metadata", registry)

function to_metadata(value, context::TransformWriteContext; path=metadata_path(context))
    return _to_metadata(value, path, context.registry)
end

_to_metadata(value::Nothing, path, registry) = nothing
_to_metadata(value::Union{Bool,Number,AbstractString}, path, registry) = value
_to_metadata(value::Symbol, path, registry) = String(value)
_to_metadata(value::AxisType, path, registry) = string(value)

function _to_metadata(value::TaggedMapping, path, registry)
    properties = Dict{Any,Any}()
    for (key, item) in value
        properties[_to_metadata(key, "$path/key", registry)] = _to_metadata(item, "$path/$key", registry)
    end
    return TaggedMapping(value.tag, properties)
end

function _to_metadata(value::TaggedSequence, path, registry)
    return TaggedSequence(value.tag, [_to_metadata(item, "$path/$index", registry) for (index, item) in enumerate(value)])
end

function _to_metadata(value::TaggedScalar, path, registry)
    return TaggedScalar(value.tag, _to_metadata(value.value, path, registry))
end

function _to_metadata(value::NamedTuple, path, registry)
    return Dict{Any,Any}(String(key) => _to_metadata(item, "$path/$key", registry) for (key, item) in pairs(value))
end

function _to_metadata(value::AbstractDict, path, registry)
    return Dict{Any,Any}(_to_metadata(key, "$path/key", registry) => _to_metadata(item, "$path/$key", registry) for (key, item) in value)
end

function _to_metadata(value::Tuple, path, registry)
    return [_to_metadata(item, "$path/$index", registry) for (index, item) in enumerate(value)]
end

function _to_metadata(value::AbstractArray, path, registry)
    return map(CartesianIndices(value)) do index
        indices = join(Tuple(index), "/")
        _to_metadata(value[index], "$path/$indices", registry)
    end
end

function _to_metadata(transform::AnnotatedTransform, path, registry)
    node = _metadata_mapping(_to_metadata(transform.transform, path, registry), path)
    properties = Pair{String,Any}[]
    transform.name === nothing || push!(properties, "name" => transform.name)
    transform.inputs === nothing || push!(properties, "inputs" => _metadata_vector(transform.inputs))
    transform.outputs === nothing || push!(properties, "outputs" => _metadata_vector(transform.outputs))
    return _decorate_metadata(node, properties)
end

function _to_metadata(transform::WithInverse, path, registry)
    node = _metadata_mapping(_to_metadata(transform.forward, path, registry), path)
    inverse = _to_metadata(transform.backward, "$path/inverse", registry)
    return _decorate_metadata(node, ("inverse" => inverse,))
end

function _to_metadata(transform::AffineTransform{1,1}, path, registry)
    factor = transform.matrix[1, 1]
    offset = transform.offset[1]
    isone(factor) && !iszero(offset) && return _to_metadata(Shift(offset), path, registry)
    transforms = Any[_to_metadata(Scale(factor), "$path/forward/1", registry)]
    iszero(offset) || push!(transforms, _to_metadata(Shift(offset), "$path/forward/2", registry))
    length(transforms) == 1 && return only(transforms)
    return _tagged_metadata("$_ASDF_TRANSFORM_TAG/compose", v"1.4.0", Dict("forward" => transforms))
end

# parallel() fuses Shift, Scale, and 2D affine branches into one block-diagonal affine; astropy has no N-D affine, so write the blocks back as a concatenate.
function _to_metadata(transform::AffineTransform{N,N}, path, registry) where {N}
    N == 2 && return invoke(_to_metadata, Tuple{AbstractTransform,Any,Any}, transform, path, registry)
    blocks = _affine_blocks(transform.matrix)
    blocks === nothing && _metadata_error(path, "a $N-dimensional affine transform can be serialized only when it is block diagonal with 1- and 2-dimensional blocks")
    return _to_metadata(ParallelTransform(Tuple(AffineTransform(transform.matrix[block, block], transform.offset[block]) for block in blocks)), path, registry)
end

function _affine_blocks(matrix::StaticMatrix{N,N}) where {N}
    blocks = UnitRange{Int}[]
    start = 1
    while start <= N
        size = findfirst(k -> start + k - 1 <= N && all(iszero, matrix[start:(start + k - 1), (start + k):N]) && all(iszero, matrix[(start + k):N, start:(start + k - 1)]), 1:2)
        size === nothing && return nothing
        push!(blocks, start:(start + size - 1))
        start += size
    end
    return blocks
end

# compose() fuses this idiom into one leaf; write it back in the form Python GWCS reads.
function _to_metadata(transform::_Polynomial2DPair, path, registry)
    return _to_metadata(ComposedTransform((RemapAxes{2}((1, 2, 1, 2)), ParallelTransform((transform.x, transform.y)))), path, registry)
end

function _to_metadata(transform::ArithmeticTransform, path, registry)
    name = transform.operator isa _AddOperator ? "add" :
        transform.operator isa _SubtractOperator ? "subtract" :
        transform.operator isa _MultiplyOperator ? "multiply" :
        transform.operator isa _DivideOperator ? "divide" :
        transform.operator isa _PowerOperator ? "power" :
        _metadata_error(path, "unsupported arithmetic operator $(typeof(transform.operator))")
    forward = [_to_metadata(transform.left, "$path/forward/1", registry), _to_metadata(transform.right, "$path/forward/2", registry)]
    return _tagged_metadata("$_ASDF_TRANSFORM_TAG/$name", v"1.4.0", Dict("forward" => forward))
end

function _to_metadata(transform::AbstractTransform, path, registry)
    registration = _lookup_transform_serializer(registry, transform, path)
    context = TransformWriteContext(registration.tag, registration.version, String(path), registry)
    properties = registration.serializer(transform, context)
    properties isa AbstractDict || _metadata_error(path, "serializer for $(typeof(transform)) returned $(typeof(properties)), expected an AbstractDict")
    canonical = _to_metadata(properties, path, registry)
    return _tagged_metadata(registration.tag, registration.version, canonical)
end

function _serialize_identity(transform, context)
    return Dict("n_dims" => input_dimension(transform))
end

_serialize_shift(transform, context) = Dict("offset" => transform.offset)
_serialize_scale(transform, context) = Dict("factor" => transform.factor)

function _serialize_affine(transform, context)
    return Dict("matrix" => _metadata_rows(transform.matrix), "translation" => _metadata_vector(transform.offset))
end

function _serialize_constant(transform, context)
    properties = Dict{String,Any}("value" => transform.value)
    metadata_version(context) >= v"1.4.0" && (properties["dimensions"] = input_dimension(transform))
    return properties
end

function _serialize_remap(transform, context)
    return Dict("mapping" => [index - 1 for index in transform.mapping], "n_inputs" => input_dimension(transform))
end

function _serialize_polynomial(transform::Polynomial1D, context)
    properties = Dict{String,Any}("coefficients" => _metadata_vector(transform.coefficients))
    if transform.domain !== nothing
        properties["domain"] = _metadata_vector(transform.domain)
        properties["window"] = _metadata_vector(transform.window)
    end
    return properties
end

function _serialize_polynomial(transform::Polynomial2D, context)
    properties = Dict{String,Any}("coefficients" => _metadata_rows(transform.coefficients))
    if transform.domain !== nothing
        properties["domain"] = [_metadata_vector(interval) for interval in transform.domain]
        properties["window"] = [_metadata_vector(interval) for interval in transform.window]
    end
    return properties
end

_projection_direction(::IntermediateToNative) = "pix2sky"
_projection_direction(::NativeToIntermediate) = "sky2pix"

_projection_properties(::GnomonicProjection) = Dict{String,Any}()
_projection_properties(::StereographicProjection) = Dict{String,Any}()
_projection_properties(::ZenithalEquidistant) = Dict{String,Any}()
_projection_properties(::ZenithalEqualArea) = Dict{String,Any}()
_projection_properties(projection::ZenithalPerspective) = Dict{String,Any}("mu" => projection.mu, "gamma" => projection.gamma)
_projection_properties(projection::SlantZenithalPerspective) = Dict{String,Any}("mu" => projection.mu, "phi0" => projection.phi0, "theta0" => projection.theta0)
_projection_properties(projection::SlantOrthographic) = Dict{String,Any}("xi" => projection.xi, "eta" => projection.eta)
_projection_properties(projection::AiryProjection) = Dict{String,Any}("theta_b" => projection.theta_b)

function _serialize_projection(transform::Union{IntermediateToNative,NativeToIntermediate}, context)
    properties = _projection_properties(transform.projection)
    properties["direction"] = _projection_direction(transform)
    return properties
end

function _serialize_native_rotation(transform::NativeToCelestial, context)
    return Dict("phi" => transform.alpha_p, "theta" => transform.delta_p, "psi" => transform.phi_p, "direction" => "native2celestial")
end

function _serialize_native_rotation(transform::CelestialToNative, context)
    return Dict("phi" => transform.alpha_p, "theta" => transform.delta_p, "psi" => transform.phi_p, "direction" => "celestial2native")
end

function _rotation_xyz_angles(transform::Rotation3D, path)
    matrix = Matrix(transform.matrix)
    T = float(eltype(matrix))
    values = Matrix{T}(matrix)
    all(isfinite, values) || _metadata_error(path, "Rotation3D matrix must be finite to serialize")
    tolerance = T(100) * eps(T)
    isapprox(transpose(values) * values, Matrix{T}(I, 3, 3); atol=tolerance, rtol=tolerance) || _metadata_error(path, "Rotation3D matrix must be orthogonal to serialize")
    isapprox(det(values), one(T); atol=tolerance, rtol=tolerance) || _metadata_error(path, "Rotation3D matrix must have determinant one to serialize")

    horizontal = hypot(values[1, 1], values[2, 1])
    beta = atan(-values[3, 1], horizontal)
    if horizontal > sqrt(eps(T))
        alpha = atan(values[3, 2], values[3, 3])
        gamma = atan(values[2, 1], values[1, 1])
    elseif beta > 0
        alpha = atan(values[1, 2], values[1, 3])
        gamma = zero(T)
    else
        alpha = atan(-values[1, 2], -values[1, 3])
        gamma = zero(T)
    end
    angles = map(angle -> -rad2deg(angle), (alpha, beta, gamma))
    reconstructed = Rotation3D(collect(angles), "xyz")
    isapprox(reconstructed.matrix, transform.matrix; atol=T(10) * tolerance, rtol=T(10) * tolerance) || _metadata_error(path, "Rotation3D matrix could not be represented as a stable xyz rotation sequence")
    return angles
end

function _serialize_rotation(transform::Rotation3D, context)
    angles = _rotation_xyz_angles(transform, metadata_path(context))
    return Dict("angles" => _metadata_vector(angles), "axes_order" => "xyz", "rotation_type" => "cartesian")
end

function _serialize_spherical_cartesian(transform::SphericalToCartesian, context)
    return Dict("transform_type" => "spherical_to_cartesian", "wrap_lon_at" => transform.wrap_lon_at)
end

function _serialize_spherical_cartesian(transform::CartesianToSpherical, context)
    return Dict("transform_type" => "cartesian_to_spherical", "wrap_lon_at" => transform.wrap_lon_at)
end

function _serialize_compound_sequence(transforms::Tuple, tag, version, path, registry)
    left = _to_metadata(first(transforms), "$path/forward/1", registry)
    if length(transforms) == 2
        right = _to_metadata(last(transforms), "$path/forward/2", registry)
    else
        properties = _serialize_compound_sequence(Base.tail(transforms), tag, version, "$path/forward/2", registry)
        right = _tagged_metadata(tag, version, properties)
    end
    return Dict("forward" => Any[left, right])
end

function _serialize_composed(transform::ComposedTransform, context)
    return _serialize_compound_sequence(transform.transforms, metadata_tag(context), metadata_version(context), metadata_path(context), context.registry)
end

function _serialize_parallel(transform::ParallelTransform, context)
    return _serialize_compound_sequence(transform.transforms, metadata_tag(context), metadata_version(context), metadata_path(context), context.registry)
end

function _serialize_fixed_inputs(transform::FixedInputsTransform, context)
    transform.atol === nothing || _metadata_error(metadata_path(context), "fixed-input atol cannot be represented by the supported metadata schema")
    transform.rtol === nothing || _metadata_error(metadata_path(context), "fixed-input rtol cannot be represented by the supported metadata schema")
    keys = Int[]
    values = Any[]
    for (index, source) in enumerate(transform.sources)
        source < 0 || continue
        push!(keys, index - 1)
        push!(values, transform.fixed_values[-source])
    end
    path = metadata_path(context)
    forward = Any[_to_metadata(transform.transform, "$path/forward/1", context.registry), Dict("keys" => keys, "values" => values)]
    return Dict("forward" => forward)
end

function _register_builtin_transform_serializers!()
    registry = default_transform_registry()
    register_transform_serializer!(registry, IdentityTransform, "$_ASDF_TRANSFORM_TAG/identity", _serialize_identity; version=v"1.4.0")
    register_transform_serializer!(registry, Shift, "$_ASDF_TRANSFORM_TAG/shift", _serialize_shift; version=v"1.4.0")
    register_transform_serializer!(registry, Scale, "$_ASDF_TRANSFORM_TAG/scale", _serialize_scale; version=v"1.4.0")
    register_transform_serializer!(registry, AffineTransform{2,2}, "$_ASDF_TRANSFORM_TAG/affine", _serialize_affine; version=v"1.5.0")
    register_transform_serializer!(registry, ConstantTransform{0}, "$_ASDF_TRANSFORM_TAG/constant", _serialize_constant; version=v"1.3.0")
    register_transform_serializer!(registry, ConstantTransform{1}, "$_ASDF_TRANSFORM_TAG/constant", _serialize_constant; version=v"1.6.0")
    register_transform_serializer!(registry, ConstantTransform{2}, "$_ASDF_TRANSFORM_TAG/constant", _serialize_constant; version=v"1.6.0")
    register_transform_serializer!(registry, RemapAxes, "$_ASDF_TRANSFORM_TAG/remap_axes", _serialize_remap; version=v"1.5.0")
    register_transform_serializer!(registry, Polynomial1D, "$_ASDF_TRANSFORM_TAG/polynomial", _serialize_polynomial; version=v"1.3.0")
    register_transform_serializer!(registry, Polynomial2D, "$_ASDF_TRANSFORM_TAG/polynomial", _serialize_polynomial; version=v"1.3.0")
    for (projection_type, tag, version) in (
        (GnomonicProjection, "gnomonic", v"1.4.0"),
        (ZenithalPerspective, "zenithal_perspective", v"1.5.0"),
        (SlantZenithalPerspective, "slant_zenithal_perspective", v"1.4.0"),
        (StereographicProjection, "stereographic", v"1.4.0"),
        (SlantOrthographic, "slant_orthographic", v"1.4.0"),
        (ZenithalEquidistant, "zenithal_equidistant", v"1.4.0"),
        (ZenithalEqualArea, "zenithal_equal_area", v"1.4.0"),
        (AiryProjection, "airy", v"1.4.0")
    )
        register_transform_serializer!(registry, IntermediateToNative{<:projection_type}, "$_ASDF_TRANSFORM_TAG/$tag", _serialize_projection; version)
        register_transform_serializer!(registry, NativeToIntermediate{<:projection_type}, "$_ASDF_TRANSFORM_TAG/$tag", _serialize_projection; version)
    end
    register_transform_serializer!(registry, NativeToCelestial, "$_ASDF_TRANSFORM_TAG/rotate3d", _serialize_native_rotation; version=v"1.5.0")
    register_transform_serializer!(registry, CelestialToNative, "$_ASDF_TRANSFORM_TAG/rotate3d", _serialize_native_rotation; version=v"1.5.0")
    register_transform_serializer!(registry, Rotation3D, "$_ASDF_TRANSFORM_TAG/rotate_sequence_3d", _serialize_rotation; version=v"1.2.0")
    register_transform_serializer!(registry, SphericalToCartesian, "$_GWCS_TAG/spherical_cartesian", _serialize_spherical_cartesian; version=v"1.3.0")
    register_transform_serializer!(registry, CartesianToSpherical, "$_GWCS_TAG/spherical_cartesian", _serialize_spherical_cartesian; version=v"1.3.0")
    register_transform_serializer!(registry, ComposedTransform, "$_ASDF_TRANSFORM_TAG/compose", _serialize_composed; version=v"1.4.0")
    register_transform_serializer!(registry, ParallelTransform, "$_ASDF_TRANSFORM_TAG/concatenate", _serialize_parallel; version=v"1.4.0")
    register_transform_serializer!(registry, FixedInputsTransform, "$_ASDF_TRANSFORM_TAG/fix_inputs", _serialize_fixed_inputs; version=v"1.4.0")
    return nothing
end

function _label_mapper_properties(mapper, path, registry)
    properties = Dict{String,Any}()
    names = input_names(mapper)
    names === nothing || (properties["inputs"] = _metadata_vector(names))
    name = label_mapper_name(mapper)
    name === nothing || (properties["name"] = name)
    return properties
end

function _mapper_input_mapping(mapper, path, registry)
    return _to_metadata(RemapAxes{input_dimension(mapper)}(mapper.mapping), "$path/inputs_mapping", registry)
end

function _to_metadata(mapper::ArrayLabelMapper, path, registry)
    properties = _label_mapper_properties(mapper, path, registry)
    properties["mapper"] = _metadata_rows(mapper.labels)
    properties["inputs_mapping"] = _mapper_input_mapping(mapper, path, registry)
    mapper.no_label === nothing || (properties["no_label"] = mapper.no_label)
    return _tagged_metadata("$_GWCS_TAG/label_mapper", v"1.3.0", properties)
end

function _to_metadata(mapper::TransformLabelMapper, path, registry)
    properties = _label_mapper_properties(mapper, path, registry)
    properties["mapper"] = _to_metadata(mapper.transform, "$path/mapper", registry)
    mapper.no_label === nothing || (properties["no_label"] = mapper.no_label)
    return _tagged_metadata("$_GWCS_TAG/label_mapper", v"1.3.0", properties)
end

function _to_metadata(mapper::DictionaryLabelMapper, path, registry)
    mapper.rtol == 1e-5 || _metadata_error(path, "dictionary label-mapper rtol $(mapper.rtol) cannot be represented by the supported schema")
    properties = _label_mapper_properties(mapper, path, registry)
    properties["inputs_mapping"] = _mapper_input_mapping(mapper, path, registry)
    properties["atol"] = mapper.atol
    properties["mapper"] = Dict("labels" => _metadata_vector(mapper.keys),
        "models" => [_to_metadata(transform, "$path/mapper/models/$index", registry) for (index, transform) in enumerate(mapper.transforms)])
    return _tagged_metadata("$_GWCS_TAG/label_mapper", v"1.3.0", properties)
end

function _to_metadata(mapper::RangeLabelMapper, path, registry)
    properties = _label_mapper_properties(mapper, path, registry)
    properties["inputs_mapping"] = _mapper_input_mapping(mapper, path, registry)
    properties["mapper"] = Dict("labels" => [_metadata_vector(interval) for interval in mapper.intervals],
        "models" => [_to_metadata(transform, "$path/mapper/models/$index", registry) for (index, transform) in enumerate(mapper.transforms)])
    return _tagged_metadata("$_GWCS_TAG/label_mapper", v"1.3.0", properties)
end

function _to_metadata(mapper::CallableLabelMapper, path, registry)
    _metadata_error(path, "CallableLabelMapper cannot be serialized canonically")
end

function _to_metadata(mapper::WithInverseLabelMapper, path, registry)
    _metadata_error(path, "WithInverseLabelMapper cannot be represented by the supported metadata schema")
end

function _to_metadata(mapper::AbstractLabelMapper, path, registry)
    _metadata_error(path, "no metadata serializer is registered for label mapper $(typeof(mapper))")
end

function _serialize_region_selector(transform::RegionSelector, context)
    path = metadata_path(context)
    inputs = something(input_names(transform), _asdf_input_names(transform))
    outputs = something(output_names(transform), _asdf_output_names(transform))
    selector = Dict("labels" => _metadata_vector(transform.labels),
        "transforms" => [_to_metadata(value, "$path/selector/transforms/$index", context.registry) for (index, value) in enumerate(transform.transforms)])
    return Dict("inputs" => _metadata_vector(inputs), "outputs" => _metadata_vector(outputs),
        "label_mapper" => _to_metadata(transform.label_mapper, "$path/label_mapper", context.registry),
        "selector" => selector, "undefined_transform_value" => transform.undefined_value)
end

function _axis_type_metadata(value)
    return string(value)
end

_unit_metadata(value) = isempty(value) ? "" : TaggedScalar("$_UNIT_TAG-1.0.0", String(value))

function _reference_metadata(value, path, registry)
    value === nothing && return nothing
    value === :ICRS && return _tagged_metadata(_ICRS_TAG, v"1.1.0", Dict("frame_attributes" => Dict()))
    return _to_metadata(value, path, registry)
end

function _frame_properties(frame, path, registry)
    properties = Dict{String,Any}("name" => frame.name, "axes_names" => _metadata_vector(frame.axes_names),
        "axes_order" => [value - 1 for value in frame.axes_order],
        "axis_physical_types" => _metadata_vector(frame.axis_physical_types),
        "unit" => [_unit_metadata(value) for value in frame.unit])
    reference = _reference_metadata(frame.reference_frame, "$path/reference_frame", registry)
    reference === nothing || (properties["reference_frame"] = reference)
    return properties
end

function _to_metadata(frame::CoordinateFrame, path, registry)
    properties = _frame_properties(frame, path, registry)
    properties["naxes"] = frame_dimension(frame)
    properties["axes_type"] = [_axis_type_metadata(value) for value in frame.axis_types]
    return _tagged_metadata("$_GWCS_TAG/frame", v"1.2.0", properties)
end

_to_metadata(frame::EmptyFrame, path, registry) = frame.name

function _to_metadata(frame::CelestialFrame, path, registry)
    return _tagged_metadata("$_GWCS_TAG/celestial_frame", v"1.2.0", _frame_properties(frame, path, registry))
end

function _to_metadata(frame::SpectralFrame, path, registry)
    properties = _frame_properties(frame, path, registry)
    properties["reference_position"] = frame.reference_position
    return _tagged_metadata("$_GWCS_TAG/spectral_frame", v"1.2.0", properties)
end

function _temporal_reference_metadata(frame::TemporalFrame, path, registry)
    frame.reference_frame === :ICRS && return _reference_metadata(:ICRS, path, registry)
    frame.reference_frame isa Union{TaggedMapping,TaggedSequence,TaggedScalar} && return _to_metadata(frame.reference_frame, path, registry)
    frame.reference_epoch === nothing && return _reference_metadata(frame.reference_frame, path, registry)
    properties = Dict{String,Any}("value" => _to_metadata(frame.reference_epoch, "$path/value", registry))
    isempty(frame.time_scale) || (properties["scale"] = frame.time_scale)
    isempty(frame.time_format) || (properties["format"] = frame.time_format)
    frame.location === nothing || (properties["location"] = _to_metadata(frame.location, "$path/location", registry))
    return _tagged_metadata(_TIME_TAG, v"1.4.0", properties)
end

function _to_metadata(frame::TemporalFrame, path, registry)
    properties = _frame_properties(frame, path, registry)
    reference = _temporal_reference_metadata(frame, "$path/reference_frame", registry)
    reference === nothing ? delete!(properties, "reference_frame") : (properties["reference_frame"] = reference)
    return _tagged_metadata("$_GWCS_TAG/temporal_frame", v"1.2.0", properties)
end

function _to_metadata(frame::StokesFrame, path, registry)
    properties = Dict{String,Any}("name" => frame.name, "axes_names" => _metadata_vector(frame.axes_names),
        "axes_order" => [value - 1 for value in frame.axes_order],
        "axis_physical_types" => _metadata_vector(frame.axis_physical_types))
    return _tagged_metadata("$_GWCS_TAG/stokes_frame", v"1.0.0", properties)
end

function _to_metadata(frame::CompositeFrame, path, registry)
    frames = [_to_metadata(value, "$path/frames/$index", registry) for (index, value) in enumerate(frame.frames)]
    return _tagged_metadata("$_GWCS_TAG/composite_frame", v"1.0.0", Dict("name" => frame.name, "frames" => frames))
end

function _to_metadata(frame::AbstractCoordinateFrame, path, registry)
    _metadata_error(path, "no metadata serializer is registered for coordinate frame $(typeof(frame))")
end

function _transform_axis_names(transform)
    return something(input_names(transform), _asdf_input_names(transform))
end

function _metadata_axis_name(names, axes_order, native_axis)
    return names[axes_order[native_axis]]
end

function _serialize_bounding_box(box::BoundingBox, names, axes_order, path, registry; include_ignored=true)
    intervals = Dict{Any,Any}()
    for (index, native_axis) in enumerate(box.axes)
        intervals[_metadata_axis_name(names, axes_order, native_axis)] = _metadata_vector(box.intervals[index])
    end
    properties = Dict{String,Any}("intervals" => intervals, "order" => "F")
    if include_ignored
        ignored = [_metadata_axis_name(names, axes_order, axis) for axis in ignored_axes(box)]
        isempty(ignored) || (properties["ignore"] = ignored)
    end
    return _tagged_metadata("$_ASDF_TRANSFORM_TAG/property/bounding_box", v"1.2.0", properties)
end

function _serialize_bounding_box(box::CompoundBoundingBox, names, axes_order, path, registry)
    selector_args = [Dict("argument" => _metadata_axis_name(names, axes_order, box.selector_axes[index]), "ignore" => box.selector_ignored[index]) for index in eachindex(box.selector_axes)]
    entries = Any[]
    for (index, entry) in enumerate(box.boxes)
        nested = _serialize_bounding_box(entry.second, names, axes_order, "$path/cbbox/$index/bbox", registry; include_ignored=false)
        push!(entries, Dict("key" => _metadata_vector(entry.first), "bbox" => nested))
    end
    properties = Dict{String,Any}("selector_args" => selector_args, "cbbox" => entries, "order" => "F")
    ignored = [_metadata_axis_name(names, axes_order, axis) for axis in box.globally_ignored]
    isempty(ignored) || (properties["ignore"] = ignored)
    return _tagged_metadata("$_ASDF_TRANSFORM_TAG/property/compound_bounding_box", v"1.2.0", properties)
end

function _to_metadata(box::AbstractBoundingBox, path, registry)
    names = ntuple(index -> "x$(index - 1)", input_dimension(box))
    return _serialize_bounding_box(box, names, ntuple(identity, input_dimension(box)), path, registry)
end

function _serialize_step(step::Step, path, registry; bounding_box=nothing, input_order=nothing)
    properties = Dict{String,Any}("frame" => _to_metadata(step.frame, "$path/frame", registry))
    if step.transform !== nothing
        transform = _metadata_mapping(_to_metadata(step.transform, "$path/transform", registry), "$path/transform")
        if bounding_box !== nothing
            names = _transform_axis_names(step.transform)
            box = _serialize_bounding_box(bounding_box, names, input_order, "$path/transform/bounding_box", registry)
            transform = _decorate_metadata(transform, ("bounding_box" => box,))
        end
        properties["transform"] = transform
    end
    return _tagged_metadata("$_GWCS_TAG/step", v"1.3.0", properties)
end

function _to_metadata(step::Step, path, registry)
    return _serialize_step(step, path, registry)
end

function _metadata_pixel_shape(wcs::WCS)
    order = input_frame(wcs).axes_order
    return [wcs.pixel_shape[findfirst(==(position), order)::Int] for position in 1:length(order)]
end

function _to_metadata(wcs::WCS, path, registry)
    steps = Any[]
    order = input_frame(wcs).axes_order
    for (index, step) in enumerate(wcs.steps)
        box = index == 1 ? wcs.bounding_box : nothing
        push!(steps, _serialize_step(step, "$path/steps/$index", registry; bounding_box=box, input_order=order))
    end
    properties = Dict{String,Any}("name" => wcs.name, "steps" => steps)
    wcs.pixel_shape === nothing || any(isnothing, wcs.pixel_shape) || (properties["pixel_shape"] = _metadata_pixel_shape(wcs))
    return _tagged_metadata("$_GWCS_TAG/wcs", v"1.4.0", properties)
end

function _to_metadata(value, path, registry)
    _metadata_error(path, "value of type $(typeof(value)) cannot be serialized canonically")
end

function write_asdf(filename, value; kwargs...)
    throw(ArgumentError("load ASDF.jl before calling write_asdf"))
end

_register_builtin_transform_serializers!()
register_transform_serializer!(default_transform_registry(), RegionSelector, "$_GWCS_TAG/regions_selector", _serialize_region_selector; version=v"1.3.0")
