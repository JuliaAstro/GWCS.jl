"""Mapping node that retains an external metadata tag."""
struct TaggedMapping{D<:AbstractDict} <: AbstractDict{Any,Any}
    tag::String
    value::D
end

Base.getindex(node::TaggedMapping, key) = node.value[key]
Base.get(node::TaggedMapping, key, default) = get(node.value, key, default)
Base.haskey(node::TaggedMapping, key) = haskey(node.value, key)
Base.length(node::TaggedMapping) = length(node.value)
Base.iterate(node::TaggedMapping, state...) = iterate(node.value, state...)

"""Sequence node that retains an external metadata tag."""
struct TaggedSequence{V<:AbstractVector} <: AbstractVector{Any}
    tag::String
    value::V
end

Base.size(node::TaggedSequence) = size(node.value)
Base.getindex(node::TaggedSequence, index::Int) = node.value[index]
Base.IndexStyle(::Type{<:TaggedSequence}) = IndexLinear()

"""Scalar node that retains an external metadata tag."""
struct TaggedScalar{T}
    tag::String
    value::T
end

_tag(node::TaggedMapping) = node.tag
_tag(node::TaggedSequence) = node.tag
_tag(node::TaggedScalar) = node.tag
_tag(node) = throw(ArgumentError("metadata node $(typeof(node)) does not retain a tag"))

_unwrap(node::TaggedScalar) = node.value
_unwrap(node) = node
_materialize_array(value) = value

"""Unwrap a scalar value retained by a tagged metadata adapter."""
metadata_value(value) = _unwrap(value)

"""Materialize an inline or lazy metadata array."""
materialize_metadata(value) = _materialize_array(value)

function _construct_tagged_node(constructor, node)
    if node isa YAML.MappingNode
        return TaggedMapping(node.tag, YAML.construct_mapping(Dict{Any,Any}, constructor, node))
    elseif node isa YAML.SequenceNode
        return TaggedSequence(node.tag, YAML.construct_sequence(constructor, node))
    end
    return TaggedScalar(node.tag, YAML.construct_scalar(constructor, node))
end

function _load_tagged_yaml(metadata::AbstractString)
    constructors = copy(YAML.default_yaml_constructors)
    constructors[nothing] = _construct_tagged_node
    return YAML.load(IOBuffer(metadata), constructors; dicttype=Dict{Any,Any})
end

function _tag_name(node)
    tag = _tag(node)
    matched = match(r"/([^/]+)-(\d+\.\d+\.\d+)$", tag)
    matched === nothing && throw(ArgumentError("unsupported metadata tag $tag"))
    VersionNumber(matched.captures[2]).major == 1 || throw(ArgumentError("unsupported metadata tag version $tag"))
    return matched.captures[1]
end

function _metadata_error(path, message)
    throw(ArgumentError("$message at $path"))
end

struct _NoMetadataDefault end
const _NO_METADATA_DEFAULT = _NoMetadataDefault()

function _required(node, key, path)
    haskey(node, key) || _metadata_error(path, "missing required property $key")
    return node[key]
end

"""Return a required mapping property or raise a path-aware metadata error."""
function required_metadata(node, key, context::TransformParseContext)
    return _required(node, key, metadata_path(context))
end

_property_path(context::TransformParseContext, key) = "$(metadata_path(context))/$key"

function _metadata_property(node, key, context; default=_NO_METADATA_DEFAULT)
    haskey(node, key) && return node[key]
    default isa _NoMetadataDefault && _metadata_error(metadata_path(context), "missing required property $key")
    return default
end

function _metadata_string(node, key, context; default=_NO_METADATA_DEFAULT)
    value = metadata_value(_metadata_property(node, key, context; default))
    value isa AbstractString || _metadata_error(_property_path(context, key), "expected a string")
    return String(value)
end

function _metadata_integer(node, key, context; default=_NO_METADATA_DEFAULT)
    value = metadata_value(_metadata_property(node, key, context; default))
    value isa Integer && !(value isa Bool) || _metadata_error(_property_path(context, key), "expected an integer")
    return Int(value)
end

function _metadata_enum(node, key, allowed, context; default=_NO_METADATA_DEFAULT)
    value = _metadata_string(node, key, context; default)
    value in allowed || _metadata_error(_property_path(context, key), "expected one of $(join(allowed, ", ")), got $value")
    return value
end

function _numeric_metadata(value, path)
    unwrapped = metadata_value(value)
    unwrapped isa Real && !(unwrapped isa Bool) || _metadata_error(path, "expected a unitless numeric value")
    return unwrapped
end

function _numeric_property(node, key, context)
    return _numeric_metadata(_metadata_property(node, key, context), _property_path(context, key))
end

function _array_metadata(value, path)
    materialized = materialize_metadata(value)
    if materialized isa AbstractDict
        name = try
            _tag_name(materialized)
        catch error
            error isa ArgumentError || rethrow()
            nothing
        end
        if name == "ndarray"
            data = _required(materialized, "data", path)
            shape = haskey(materialized, "shape") ? Tuple(Int(_unwrap(axis)) for axis in materialized["shape"]) : nothing
            return data, shape
        end
    end
    return materialized, nothing
end

function _numeric_vector(value, path)
    materialized, declared_shape = _array_metadata(value, path)
    materialized isa AbstractVector || _metadata_error(path, "expected a one-dimensional numeric array")
    values = map(enumerate(materialized)) do (index, item)
        item isa AbstractVector && _metadata_error("$path/$index", "expected a numeric value")
        _numeric_metadata(item, "$path/$index")
    end
    declared_shape === nothing || declared_shape == (length(values),) || _metadata_error("$path/shape", "declared shape $declared_shape does not match vector length $(length(values))")
    isempty(values) && return Float64[]
    return collect(promote(values...))
end

function _numeric_matrix(value, path)
    materialized, declared_shape = _array_metadata(value, path)
    matrix = if materialized isa AbstractMatrix
        map(CartesianIndices(materialized)) do index
            _numeric_metadata(materialized[index], "$path/$(index[1])/$(index[2])")
        end
    elseif materialized isa AbstractVector
        isempty(materialized) && _metadata_error(path, "expected a nonempty numeric matrix")
        all(row -> row isa AbstractVector, materialized) || _metadata_error(path, "expected a two-dimensional numeric array")
        width = length(first(materialized))
        width > 0 || _metadata_error(path, "expected a nonempty numeric matrix")
        all(row -> length(row) == width, materialized) || _metadata_error(path, "matrix rows must have equal lengths")
        rows = map(enumerate(materialized)) do (row_index, row)
            map(enumerate(row)) do (column_index, item)
                _numeric_metadata(item, "$path/$row_index/$column_index")
            end
        end
        permutedims(reduce(hcat, rows))
    else
        _metadata_error(path, "expected a two-dimensional numeric array")
    end
    declared_shape === nothing || declared_shape == size(matrix) || _metadata_error("$path/shape", "declared shape $declared_shape does not match matrix size $(size(matrix))")
    return matrix
end

function _metadata_axis_names(node, key, count, boundary, context)
    path = _property_path(context, key)
    values = node[key]
    values isa AbstractVector || _metadata_error(path, "expected an array of strings")
    names = map(enumerate(values)) do (index, value)
        unwrapped = metadata_value(value)
        unwrapped isa AbstractString || _metadata_error("$path/$index", "expected a string")
        String(unwrapped)
    end
    try
        return _axis_names(names, count, boundary)
    catch error
        error isa ArgumentError || error isa DimensionMismatch || rethrow()
        _metadata_error(path, sprint(showerror, error))
    end
end

function _tuple_strings(values)
    return Tuple(String(_unwrap(value)) for value in values)
end

function _tuple_ints(values; zero_based=false)
    offset = zero_based ? 1 : 0
    return Tuple(Int(_unwrap(value)) + offset for value in values)
end

function _coefficient_array(value)
    materialized = _materialize_array(value)
    materialized isa AbstractMatrix && return Array(materialized)
    materialized isa AbstractVector || throw(ArgumentError("polynomial coefficients must be an array"))
    isempty(materialized) && return Float64[]
    first(materialized) isa AbstractVector || return collect(materialized)
    rows = map(collect, materialized)
    all(row -> length(row) == length(first(rows)), rows) || throw(ArgumentError("polynomial coefficient rows must have equal lengths"))
    return permutedims(reduce(hcat, rows))
end

function _domain_window(node)
    domain = get(node, "domain", nothing)
    window = get(node, "window", nothing)
    if domain === nothing || all(isnothing, domain)
        return nothing, nothing
    end
    window === nothing && throw(ArgumentError("polynomial metadata provides a domain without a window"))
    if first(domain) isa AbstractVector
        return Tuple(Tuple(pair) for pair in domain), Tuple(Tuple(pair) for pair in window)
    end
    return Tuple(domain), Tuple(window)
end

function _parse_compose(node, context)
    values = required_metadata(node, "forward", context)
    path = metadata_path(context)
    transforms = map(enumerate(values)) do (index, value)
        parse_transform(value; path="$path/forward/$index", registry=context.registry)
    end
    return compose(transforms...)
end

function _parse_concatenate(node, context)
    values = required_metadata(node, "forward", context)
    path = metadata_path(context)
    transforms = map(enumerate(values)) do (index, value)
        parse_transform(value; path="$path/forward/$index", registry=context.registry)
    end
    inputs = _asdf_parallel_names(transforms, _asdf_input_names)
    outputs = _asdf_parallel_names(transforms, _asdf_output_names)
    return AnnotatedTransform(parallel(transforms...); inputs, outputs)
end

function _parse_arithmetic(node, context, operator)
    path = metadata_path(context)
    values = required_metadata(node, "forward", context)
    values isa AbstractVector || _metadata_error("$path/forward", "expected an array containing two transforms")
    length(values) == 2 || _metadata_error("$path/forward", "arithmetic metadata must contain exactly two transforms")
    left = parse_transform(values[1]; path="$path/forward/1", registry=context.registry)
    right = parse_transform(values[2]; path="$path/forward/2", registry=context.registry)
    try
        return elementwise(operator, left, right)
    catch error
        error isa DimensionMismatch || rethrow()
        _metadata_error(path, sprint(showerror, error))
    end
end

_parse_add(node, context) = _parse_arithmetic(node, context, +)
_parse_subtract(node, context) = _parse_arithmetic(node, context, -)
_parse_multiply(node, context) = _parse_arithmetic(node, context, *)
_parse_divide(node, context) = _parse_arithmetic(node, context, /)
_parse_power(node, context) = _parse_arithmetic(node, context, ^)

_parse_shift(node, context) = Shift(metadata_value(required_metadata(node, "offset", context)))
_parse_scale(node, context) = Scale(metadata_value(required_metadata(node, "factor", context)))

function _parse_identity(node, context)
    dimensions = _metadata_integer(node, "n_dims", context; default=1)
    dimensions > 0 || _metadata_error(_property_path(context, "n_dims"), "identity dimension must be positive")
    return IdentityTransform{dimensions}()
end

function _parse_affine(node, context)
    matrix = _numeric_matrix(_metadata_property(node, "matrix", context), _property_path(context, "matrix"))
    size(matrix) == (2, 2) || _metadata_error(_property_path(context, "matrix"), "affine matrix must be 2 x 2")
    translation = haskey(node, "translation") ? _numeric_vector(node["translation"], _property_path(context, "translation")) : zeros(eltype(matrix), 2)
    length(translation) == 2 || _metadata_error(_property_path(context, "translation"), "affine translation must have length 2")
    return AffineTransform(matrix, translation)
end

function _parse_gnomonic(node, context)
    direction = _metadata_enum(node, "direction", ("pix2sky", "sky2pix"), context; default="pix2sky")
    return direction == "pix2sky" ? TangentPlaneToNative() : NativeToTangentPlane()
end

function _parse_constant(node, context)
    value = metadata_value(required_metadata(node, "value", context))
    metadata_version(context) < v"1.4.0" && return ConstantTransform(value)
    dimensions = Int(metadata_value(required_metadata(node, "dimensions", context)))
    dimensions in (1, 2) || _metadata_error(metadata_path(context), "constant dimensions must be 1 or 2")
    return ConstantTransform(value; n_inputs=dimensions)
end

_asdf_input_names(transform::AnnotatedTransform) = something(input_names(transform), _asdf_input_names(transform.transform))
_asdf_input_names(transform::WithInverse) = _asdf_input_names(transform.forward)
_asdf_input_names(transform::ComposedTransform) = _asdf_input_names(first(transform.transforms))
_asdf_input_names(transform::ArithmeticTransform) = something(input_names(transform), _asdf_input_names(transform.left))
_asdf_output_names(transform::AnnotatedTransform) = something(output_names(transform), _asdf_output_names(transform.transform))
_asdf_output_names(transform::WithInverse) = _asdf_output_names(transform.forward)
_asdf_output_names(transform::ComposedTransform) = _asdf_output_names(last(transform.transforms))
_asdf_output_names(transform::ArithmeticTransform) = something(output_names(transform), _asdf_output_names(transform.left))

function _asdf_parallel_names(transforms, accessor)
    branch_names = map(accessor, transforms)
    flattened = Tuple(Iterators.flatten(branch_names))
    length(unique(flattened)) == length(flattened) && return flattened
    return Tuple(name * string(branch - 1) for (branch, names) in enumerate(branch_names) for name in names)
end

_asdf_input_names(transform::ParallelTransform) = _asdf_parallel_names(transform.transforms, _asdf_input_names)
_asdf_output_names(transform::ParallelTransform) = _asdf_parallel_names(transform.transforms, _asdf_output_names)

function _asdf_input_names(transform::FixedInputsTransform{Nfree}) where {Nfree}
    names = _asdf_input_names(transform.transform)
    return ntuple(index -> names[transform.free_indices[index]], Nfree)
end

_asdf_output_names(transform::FixedInputsTransform) = _asdf_output_names(transform.transform)

function _asdf_input_names(transform::AbstractTransform{N}) where {N}
    N == 1 && return ("x",)
    N == 2 && return ("x", "y")
    return ntuple(index -> "x$(index - 1)", N)
end

function _asdf_output_names(transform::AbstractTransform{Nin,Nout}) where {Nin,Nout}
    Nout == 1 && return ("y",)
    Nout == 2 && return ("x", "y")
    return ntuple(index -> "y$(index - 1)", Nout)
end

function _parse_fixed_inputs(node, context)
    path = metadata_path(context)
    forward = required_metadata(node, "forward", context)
    length(forward) == 2 || _metadata_error(path, "fixed-input metadata must contain a transform and fixed values")
    transform = parse_transform(forward[1]; path="$path/forward/1", registry=context.registry)
    specification = forward[2]
    specification isa AbstractDict || _metadata_error("$path/forward/2", "fixed-input specification must be a mapping")
    keys = _required(specification, "keys", "$path/forward/2")
    values = _required(specification, "values", "$path/forward/2")
    length(keys) == length(values) || _metadata_error("$path/forward/2", "fixed-input keys and values must have equal lengths")
    names = input_names(transform)
    names === nothing && (names = _asdf_input_names(transform))
    fixed = map(keys, values) do key, value
        unwrapped = metadata_value(key)
        converted = unwrapped isa Integer ? Int(unwrapped) + 1 : _fixed_index(String(unwrapped), names, input_dimension(transform))
        converted => metadata_value(value)
    end
    try
        return fix_inputs(transform, fixed)
    catch error
        error isa ArgumentError || error isa DimensionMismatch || rethrow()
        _metadata_error(path, sprint(showerror, error))
    end
end

function _parse_remap(node, context)
    mapping = _tuple_ints(required_metadata(node, "mapping", context); zero_based=true)
    n_inputs = Int(get(node, "n_inputs", maximum(mapping)))
    return RemapAxes{n_inputs}(mapping)
end

function _parse_polynomial(node, context)
    coefficients = _coefficient_array(required_metadata(node, "coefficients", context))
    domain, window = _domain_window(node)
    ndims(coefficients) == 1 && return Polynomial1D(coefficients; domain, window)
    ndims(coefficients) == 2 && return Polynomial2D(coefficients; domain, window)
    _metadata_error(metadata_path(context), "only one- and two-dimensional polynomials are supported")
end

const _EULER_ORDERS = ("zxz", "xyx", "yzy", "zyz", "xzx", "yxy", "xyz", "yzx", "zxy", "xzy", "zyx", "yxz")

function _spherical_rotation(angles, axes_order)
    return compose(SphericalToCartesian(180), Rotation3D(collect(angles), axes_order), CartesianToSpherical(180))
end

function _parse_rotate3d(node, context)
    phi = _numeric_property(node, "phi", context)
    theta = _numeric_property(node, "theta", context)
    psi = _numeric_property(node, "psi", context)
    directions = (_EULER_ORDERS..., "native2celestial", "celestial2native")
    direction = _metadata_enum(node, "direction", directions, context)
    direction == "native2celestial" && return NativeToCelestial(phi, theta, psi)
    direction == "celestial2native" && return CelestialToNative(phi, theta, psi)
    return _spherical_rotation((phi, theta, psi), direction)
end

function _parse_rotation(node, context)
    rotation_type = _metadata_enum(node, "rotation_type", ("spherical", "cartesian"), context)
    angles = _numeric_vector(_metadata_property(node, "angles", context), _property_path(context, "angles"))
    isempty(angles) && _metadata_error(_property_path(context, "angles"), "rotation angles must not be empty")
    axes_order = _metadata_string(node, "axes_order", context)
    length(angles) == length(axes_order) || _metadata_error(metadata_path(context), "the number of angles must match axes_order")
    all(axis -> axis in ('x', 'y', 'z'), axes_order) || _metadata_error(_property_path(context, "axes_order"), "rotation axes must contain only x, y, or z")
    rotation_type == "cartesian" && return Rotation3D(angles, axes_order)
    return _spherical_rotation(angles, axes_order)
end

function _parse_spherical_cartesian(node, context)
    transform_type = String(metadata_value(required_metadata(node, "transform_type", context)))
    wrap_lon_at = Int(metadata_value(required_metadata(node, "wrap_lon_at", context)))
    transform_type == "spherical_to_cartesian" && return SphericalToCartesian(wrap_lon_at)
    transform_type == "cartesian_to_spherical" && return CartesianToSpherical(wrap_lon_at)
    _metadata_error(metadata_path(context), "unsupported spherical/cartesian transform type $transform_type")
end

function _parse_transform_annotation(transform, node, context)
    has_annotation = haskey(node, "name") || haskey(node, "inputs") || haskey(node, "outputs")
    has_annotation || return transform
    name = if haskey(node, "name")
        value = _metadata_string(node, "name", context)
        isempty(value) && _metadata_error(_property_path(context, "name"), "transform name must not be empty")
        value
    else
        transform_name(transform)
    end
    inputs = haskey(node, "inputs") ? _metadata_axis_names(node, "inputs", input_dimension(transform), "input", context) : input_names(transform)
    outputs = haskey(node, "outputs") ? _metadata_axis_names(node, "outputs", output_dimension(transform), "output", context) : output_names(transform)
    return annotate_transform(transform; name, inputs, outputs)
end

function parse_transform(node; path="transform", registry=default_transform_registry())
    tag, version = _split_metadata_tag(_tag(node), path)
    registration = _lookup_transform_registration(registry, tag, version, path)
    context = TransformParseContext(tag, version, String(path), registry)
    transform = registration.constructor(node, context)
    transform isa AbstractTransform || _metadata_error(path, "registered constructor for $tag returned $(typeof(transform)), expected AbstractTransform")
    if haskey(node, "inverse")
        backward = parse_transform(node["inverse"]; path="$path/inverse", registry)
        input_dimension(transform) == output_dimension(backward) && output_dimension(transform) == input_dimension(backward) || _metadata_error(path, "explicit inverse dimensions do not match the forward transform")
        transform = WithInverse(transform, backward)
    end
    return _parse_transform_annotation(transform, node, context)
end

function parse_transform(node, context::TransformParseContext; path=metadata_path(context))
    return parse_transform(node; path, registry=context.registry)
end

function _register_builtin_transform_parsers!()
    asdf_transform = "tag:stsci.edu:asdf/transform"
    register_transform!("$asdf_transform/compose", _parse_compose)
    register_transform!("$asdf_transform/concatenate", _parse_concatenate)
    register_transform!("$asdf_transform/add", _parse_add; max_version=v"1.5.0")
    register_transform!("$asdf_transform/subtract", _parse_subtract; max_version=v"1.5.0")
    register_transform!("$asdf_transform/multiply", _parse_multiply; max_version=v"1.5.0")
    register_transform!("$asdf_transform/divide", _parse_divide; max_version=v"1.5.0")
    register_transform!("$asdf_transform/power", _parse_power; max_version=v"1.5.0")
    register_transform!("$asdf_transform/shift", _parse_shift)
    register_transform!("$asdf_transform/scale", _parse_scale)
    register_transform!("$asdf_transform/identity", _parse_identity; max_version=v"1.5.0")
    register_transform!("$asdf_transform/affine", _parse_affine; max_version=v"1.6.0")
    register_transform!("$asdf_transform/gnomonic", _parse_gnomonic; max_version=v"1.5.0")
    register_transform!("$asdf_transform/rotate3d", _parse_rotate3d; max_version=v"1.6.0")
    register_transform!("$asdf_transform/constant", _parse_constant)
    register_transform!("$asdf_transform/fix_inputs", _parse_fixed_inputs; min_version=v"1.1.0")
    register_transform!("$asdf_transform/remap_axes", _parse_remap)
    register_transform!("$asdf_transform/polynomial", _parse_polynomial)
    register_transform!("$asdf_transform/rotate_sequence_3d", _parse_rotation; max_version=v"1.3.0")
    register_transform!("tag:stsci.edu:gwcs/spherical_cartesian", _parse_spherical_cartesian)
    return nothing
end

_register_builtin_transform_parsers!()

function _parse_reference_frame(value)
    value === nothing && return nothing
    try
        _tag(value)
    catch error
        error isa ArgumentError || rethrow()
        return _unwrap(value)
    end
    value isa AbstractDict && _tag_name(value) == "icrs" && return :ICRS
    return value
end

_metadata_sequence(value) = value isa AbstractVector || value isa Tuple ? value : (value,)

function _metadata_frame_strings(node, key, default)
    haskey(node, key) || return default
    return Tuple(begin
        unwrapped = _unwrap(value)
        unwrapped === nothing ? "" : String(unwrapped)
    end for value in _metadata_sequence(node[key]))
end

function _frame_order(node, count)
    values = haskey(node, "axes_order") ? _metadata_sequence(node["axes_order"]) : 0:(count - 1)
    return _tuple_ints(values; zero_based=true)
end

function _frame_axis_types(node, count)
    values = haskey(node, "axes_type") ? _metadata_sequence(node["axes_type"]) : ntuple(_ -> "UNDEFINED", count)
    return Tuple(_axis_type(String(_unwrap(value))) for value in values)
end

function _temporal_reference(node)
    node === nothing && return nothing, nothing, "", "", nothing
    reference = _parse_reference_frame(node)
    node isa AbstractDict || return reference, _unwrap(node), "", "", nothing
    epoch = haskey(node, "value") ? _unwrap(node["value"]) : nothing
    scale = haskey(node, "scale") ? String(_unwrap(node["scale"])) : ""
    format = haskey(node, "format") ? String(_unwrap(node["format"])) : ""
    location = get(node, "location", nothing)
    return reference, epoch, scale, format, location
end

function _parse_composite_frame(node, path)
    length(node) == 2 || _metadata_error(path, "composite frame must contain only name and frames")
    frame_name = String(_unwrap(_required(node, "name", path)))
    values = _required(node, "frames", path)
    frames = map(enumerate(values)) do (index, value)
        _parse_frame(value, "$path/frames/$index"; validate_order=false)
    end
    try
        return CompositeFrame(frames; name=frame_name)
    catch error
        error isa ArgumentError || error isa DimensionMismatch || rethrow()
        _metadata_error(path, sprint(showerror, error))
    end
end

function _parse_frame(node, path; validate_order=true)
    name = _tag_name(node)
    name == "composite_frame" && return _parse_composite_frame(node, path)
    name in ("frame", "frame2d", "celestial_frame", "spectral_frame", "temporal_frame", "stokes_frame") || _metadata_error(path, "unsupported coordinate frame tag $(_tag(node))")
    frame_name = String(_unwrap(_required(node, "name", path)))
    try
        frame = if name == "frame"
            N = Int(_unwrap(_required(node, "naxes", path)))
            axes_names = _metadata_frame_strings(node, "axes_names", ntuple(i -> "x$(i - 1)", N))
            axes_order = _frame_order(node, N)
            axis_types = _frame_axis_types(node, N)
            physical_types = haskey(node, "axis_physical_types") ? _metadata_frame_strings(node, "axis_physical_types", ()) : nothing
            units = _metadata_frame_strings(node, "unit", ntuple(_ -> "", N))
            reference_frame = _parse_reference_frame(get(node, "reference_frame", nothing))
            CoordinateFrame{N}(; name=frame_name, axes_names, axes_order, axis_types, axis_physical_types=physical_types, unit=units, reference_frame)
        elseif name == "frame2d"
            axes_names = _metadata_frame_strings(node, "axes_names", ("x", "y"))
            axes_order = _frame_order(node, 2)
            axis_types = haskey(node, "axes_type") ? _frame_axis_types(node, 2) : nothing
            physical_types = haskey(node, "axis_physical_types") ? _metadata_frame_strings(node, "axis_physical_types", ()) : nothing
            units = _metadata_frame_strings(node, "unit", ("", ""))
            reference_frame = _parse_reference_frame(get(node, "reference_frame", nothing))
            Frame2D(; name=frame_name, axes_names, axes_order, axis_types, axis_physical_types=physical_types, unit=units, reference_frame)
        elseif name == "celestial_frame"
            N = haskey(node, "axes_names") ? length(_metadata_sequence(node["axes_names"])) :
                haskey(node, "axes_order") ? length(_metadata_sequence(node["axes_order"])) :
                haskey(node, "unit") ? length(_metadata_sequence(node["unit"])) : 2
            axes_names = _metadata_frame_strings(node, "axes_names", ntuple(i -> i == 1 ? "lon" : i == 2 ? "lat" : "distance", N))
            axes_order = _frame_order(node, N)
            physical_types = haskey(node, "axis_physical_types") ? _metadata_frame_strings(node, "axis_physical_types", ()) : nothing
            units = _metadata_frame_strings(node, "unit", ntuple(i -> i <= 2 ? "deg" : "", N))
            reference_frame = _parse_reference_frame(get(node, "reference_frame", nothing))
            CelestialFrame{N}(; name=frame_name, axes_names, axes_order, axis_physical_types=physical_types, unit=units, reference_frame)
        elseif name == "spectral_frame"
            axes_names = _metadata_frame_strings(node, "axes_names", ("spectral",))
            axes_order = _frame_order(node, 1)
            physical_types = haskey(node, "axis_physical_types") ? _metadata_frame_strings(node, "axis_physical_types", ()) : nothing
            units = _metadata_frame_strings(node, "unit", ("",))
            reference_frame = _parse_reference_frame(get(node, "reference_frame", nothing))
            reference_position = String(_unwrap(get(node, "reference_position", "geocenter")))
            SpectralFrame(; name=frame_name, axes_names, axes_order, axis_physical_types=physical_types, unit=units, reference_frame, reference_position)
        elseif name == "temporal_frame"
            axes_names = _metadata_frame_strings(node, "axes_names", ("time",))
            axes_order = _frame_order(node, 1)
            physical_types = _metadata_frame_strings(node, "axis_physical_types", ("time",))
            units = _metadata_frame_strings(node, "unit", ("s",))
            reference, epoch, scale, format, location = _temporal_reference(get(node, "reference_frame", nothing))
            TemporalFrame(; name=frame_name, axes_names, axes_order, axis_physical_types=physical_types, unit=units,
                reference_frame=reference, reference_epoch=epoch, time_scale=scale, time_format=format, location)
        elseif name == "stokes_frame"
            axes_names = _metadata_frame_strings(node, "axes_names", ("stokes",))
            axes_order = _frame_order(node, 1)
            physical_types = _metadata_frame_strings(node, "axis_physical_types", ("phys.polarization.stokes",))
            StokesFrame(; name=frame_name, axes_names, axes_order, axis_physical_types=physical_types)
        end
        return validate_order ? _validate_frame_axes_order(frame) : frame
    catch error
        error isa ArgumentError || error isa DimensionMismatch || rethrow()
        _metadata_error(path, sprint(showerror, error))
    end
end

function _parse_step(node, path, registry)
    _tag_name(node) == "step" || _metadata_error(path, "expected a GWCS step")
    frame = _parse_frame(_required(node, "frame", path), "$path/frame")
    transform_node = get(node, "transform", nothing)
    return transform_node === nothing ? Step(frame) : Step(frame, parse_transform(transform_node; path="$path/transform", registry))
end

function _parse_bounding_box(node, path)
    _tag_name(node) == "bounding_box" || _metadata_error(path, "unsupported bounding-box tag $(_tag(node))")
    raw_intervals = _required(node, "intervals", path)
    keys_sorted = sort(collect(keys(raw_intervals)); by=key -> parse(Int, match(r"(\d+)$", String(key)).captures[1]))
    intervals = Tuple(begin
        values = raw_intervals[key]
        (float(_unwrap(values[1])), float(_unwrap(values[2])))
    end for key in keys_sorted)
    order = Symbol(String(_unwrap(get(node, "order", "F"))))
    return BoundingBox(intervals; order)
end

function _metadata_native_order(values, frame::AbstractCoordinateFrame{N}, path) where {N}
    length(values) == N || _metadata_error(path, "metadata value has length $(length(values)), expected $N")
    return ntuple(index -> values[frame.axes_order[index]], N)
end

function parse_wcs(node::AbstractDict; registry=default_transform_registry())
    _tag_name(node) == "wcs" || throw(ArgumentError("expected a GWCS WCS metadata node"))
    raw_steps = _required(node, "steps", "wcs")
    steps = map(enumerate(raw_steps)) do (index, value)
        _parse_step(value, "wcs/steps/$index", registry)
    end
    first_transform = get(first(raw_steps), "transform", nothing)
    bounding_box = first_transform === nothing ? nothing : get(first_transform, "bounding_box", nothing)
    parsed_box = bounding_box === nothing ? nothing : _parse_bounding_box(bounding_box, "wcs/steps/1/transform/bounding_box")
    input = first(steps).frame
    parsed_box === nothing || (parsed_box = BoundingBox(_metadata_native_order(parsed_box.intervals, input, "wcs/steps/1/transform/bounding_box")))
    name = String(_unwrap(get(node, "name", "")))
    pixel_shape = get(node, "pixel_shape", nothing)
    pixel_shape === nothing || (pixel_shape = _metadata_native_order(Tuple(Int(_unwrap(value)) for value in pixel_shape), input, "wcs/pixel_shape"))
    return WCS(steps...; name, pixel_shape, bounding_box=parsed_box)
end

parse_wcs(metadata::AbstractString; registry=default_transform_registry()) = parse_wcs(_load_tagged_yaml(metadata); registry)
