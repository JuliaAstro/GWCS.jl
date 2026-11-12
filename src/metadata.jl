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
# Whether a node retains a tag; the ASDF extension adds ASDF.jl's tagged nodes.
_has_tag(node) = false
_has_tag(::Union{TaggedMapping,TaggedSequence,TaggedScalar}) = true

_unwrap(node::TaggedScalar) = node.value
_unwrap(node) = node
_materialize_array(value) = value

"""
    metadata_value(value)

Return the plain value of a metadata scalar, such as the string inside a tagged scalar node, and any other value
unchanged.
"""
metadata_value(value) = _unwrap(value)

"""
    materialize_metadata(value)

Return the Julia array of a metadata array, reading a lazy ASDF block as needed. A matrix keeps the file's row and
column order rather than being transposed to an image's `data[x, y]` order, so `m[i + 1, j + 1]` is NumPy's `m[i, j]`.
Any other value is returned unchanged.
"""
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

function _metadata_error(path, message)
    throw(ArgumentError("$message at $path"))
end

struct _NoMetadataDefault end
const _NO_METADATA_DEFAULT = _NoMetadataDefault()

function _required(node, key, path)
    haskey(node, key) || _metadata_error(path, "missing required property $key")
    return node[key]
end

"""
    required_metadata(node, key, context::TransformParseContext)

Return the property `key` of a metadata node. A missing property throws an `ArgumentError` naming its path.
"""
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
        tag = try
            _tag(materialized)
        catch error
            error isa ArgumentError || rethrow()
            nothing
        end
        if tag !== nothing
            family, _ = _split_metadata_tag(tag, path)
        end
        if tag !== nothing && family == "tag:stsci.edu:asdf/core/ndarray"
            materialized, _ = _prepare_metadata(materialized, path, default_transform_registry())
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

function _metadata_matrix(value, path)
    materialized, declared_shape = _array_metadata(value, path)
    matrix = if materialized isa AbstractMatrix
        map(metadata_value, materialized)
    elseif materialized isa AbstractVector
        isempty(materialized) && _metadata_error(path, "expected a nonempty two-dimensional array")
        all(row -> row isa AbstractVector, materialized) || _metadata_error(path, "expected a two-dimensional array")
        width = length(first(materialized))
        width > 0 || _metadata_error(path, "expected a nonempty two-dimensional array")
        all(row -> length(row) == width, materialized) || _metadata_error(path, "array rows must have equal lengths")
        rows = map(row -> map(metadata_value, row), materialized)
        permutedims(reduce(hcat, rows))
    else
        _metadata_error(path, "expected a two-dimensional array")
    end
    declared_shape === nothing || declared_shape == size(matrix) || _metadata_error("$path/shape", "declared shape $declared_shape does not match array size $(size(matrix))")
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
    return compose((parse_transform(value; path="$path/forward/$index", registry=context.registry) for (index, value) in enumerate(values))...)
end

function _parse_concatenate(node, context)
    values = required_metadata(node, "forward", context)
    path = metadata_path(context)
    transforms = Any[parse_transform(value; path="$path/forward/$index", registry=context.registry) for (index, value) in enumerate(values)]
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
_parse_rotate2d(node, context) = Rotation2D(metadata_value(required_metadata(node, "angle", context)))

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

function _parse_projection(node, context, projection)
    direction = _metadata_enum(node, "direction", ("pix2sky", "sky2pix"), context; default="pix2sky")
    return direction == "pix2sky" ? IntermediateToNative(projection) : NativeToIntermediate(projection)
end

_parse_gnomonic(node, context) = _parse_projection(node, context, GnomonicProjection())

function _projection_parameter(node, key, default)
    return haskey(node, key) ? metadata_value(node[key]) : default
end

function _construct_projection(constructor, context)
    try
        return constructor()
    catch error
        error isa ArgumentError || rethrow()
        _metadata_error(metadata_path(context), sprint(showerror, error))
    end
end

function _parse_zenithal_perspective(node, context)
    projection = _construct_projection(context) do
        ZenithalPerspective(_projection_parameter(node, "mu", 0), _projection_parameter(node, "gamma", 0))
    end
    return _parse_projection(node, context, projection)
end

function _parse_slant_zenithal_perspective(node, context)
    projection = _construct_projection(context) do
        SlantZenithalPerspective(_projection_parameter(node, "mu", 0), _projection_parameter(node, "phi0", 0), _projection_parameter(node, "theta0", 90))
    end
    return _parse_projection(node, context, projection)
end

_parse_stereographic(node, context) = _parse_projection(node, context, StereographicProjection())

function _parse_slant_orthographic(node, context)
    projection = _construct_projection(context) do
        SlantOrthographic(_projection_parameter(node, "xi", 0), _projection_parameter(node, "eta", 0))
    end
    return _parse_projection(node, context, projection)
end

_parse_zenithal_equidistant(node, context) = _parse_projection(node, context, ZenithalEquidistant())
_parse_zenithal_equal_area(node, context) = _parse_projection(node, context, ZenithalEqualArea())
function _parse_airy(node, context)
    projection = _construct_projection(context) do
        AiryProjection(_projection_parameter(node, "theta_b", 90))
    end
    return _parse_projection(node, context, projection)
end

# The schemas name the cylinder radius `lambda`, but asdf-astropy reads and writes astropy's `lam` and ignores
# `lambda`. Either is read; the writer uses `lam`.
function _cylinder_radius(node, context)
    haskey(node, "lam") && haskey(node, "lambda") && _metadata_error(metadata_path(context), "the cylinder radius is given as both lam and lambda")
    return _projection_parameter(node, haskey(node, "lambda") ? "lambda" : "lam", 1)
end

function _parse_cylindrical_equal_area(node, context)
    projection = _construct_projection(context) do
        CylindricalEqualArea(_cylinder_radius(node, context))
    end
    return _parse_projection(node, context, projection)
end

function _parse_cylindrical_perspective(node, context)
    projection = _construct_projection(context) do
        CylindricalPerspective(_projection_parameter(node, "mu", 1), _cylinder_radius(node, context))
    end
    return _parse_projection(node, context, projection)
end

# The schemas give sigma no default. astropy's default of 90 puts both standard parallels at the pole, so a missing
# sigma is an error, as in WCSLIB.
function _conic_parser(constructor)
    return function (node, context)
        sigma = metadata_value(required_metadata(node, "sigma", context))
        projection = _construct_projection(context) do
            constructor(sigma, _projection_parameter(node, "delta", 0))
        end
        return _parse_projection(node, context, projection)
    end
end

function _parse_bonne_equal_area(node, context)
    projection = _construct_projection(context) do
        BonneProjection(_projection_parameter(node, "theta1", 0))
    end
    return _parse_projection(node, context, projection)
end

function _parse_healpix(node, context)
    H, K = _projection_parameter(node, "H", 4), _projection_parameter(node, "X", 3)
    isinteger(H) && isinteger(K) || _metadata_error(metadata_path(context), "HEALPix H and X must be integers, got $H and $K")
    projection = _construct_projection(context) do
        HEALPixProjection(Int(H), Int(K))
    end
    return _parse_projection(node, context, projection)
end

function _parse_constant(node, context)
    value = metadata_value(required_metadata(node, "value", context))
    dimensions = Int(metadata_value(required_metadata(node, "dimensions", context)))
    dimensions in (0, 1, 2) || _metadata_error(metadata_path(context), "constant dimensions must be 0, 1, or 2")
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

_asdf_output_names(transform::_Polynomial2DPair) = _asdf_parallel_names((transform.x, transform.y), _asdf_output_names)

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

_remap_axes(transform::RemapAxes) = transform
_remap_axes(transform::AnnotatedTransform) = _remap_axes(transform.transform)
_remap_axes(transform) = nothing

function _metadata_inputs(node, context; default=nothing)
    haskey(node, "inputs") || return default
    values = node["inputs"]
    values isa AbstractVector || _metadata_error(_property_path(context, "inputs"), "expected an array of strings")
    inputs = map(enumerate(values)) do (index, value)
        unwrapped = metadata_value(value)
        unwrapped isa AbstractString || _metadata_error("$(_property_path(context, "inputs"))/$index", "expected a string")
        String(unwrapped)
    end
    isempty(inputs) && _metadata_error(_property_path(context, "inputs"), "label-mapper inputs must not be empty")
    length(unique(inputs)) == length(inputs) || _metadata_error(_property_path(context, "inputs"), "label-mapper inputs must be unique")
    return Tuple(inputs)
end

function _parse_mapper_mapping(node, context, inputs, default_dimension)
    haskey(node, "inputs_mapping") || return inputs === nothing ? ntuple(identity, default_dimension) : ntuple(identity, min(length(inputs), default_dimension))
    path = _property_path(context, "inputs_mapping")
    transform = parse_transform(node["inputs_mapping"]; path, registry=context.registry)
    remap = _remap_axes(transform)
    remap === nothing && _metadata_error(path, "expected a remap_axes transform")
    inputs === nothing || input_dimension(remap) == length(inputs) || _metadata_error(path, "input mapping has input dimension $(input_dimension(remap)), expected $(length(inputs))")
    output_dimension(remap) == default_dimension || _metadata_error(path, "input mapping has output dimension $(output_dimension(remap)), expected $default_dimension")
    return remap.mapping
end

function _metadata_label_type(labels, path)
    isempty(labels) && _metadata_error(path, "selector labels must not be empty")
    types = unique(typeof.(labels))
    length(types) == 1 || _metadata_error(path, "selector labels must have one concrete type")
    type = only(types)
    isconcretetype(type) || _metadata_error(path, "selector label type must be concrete")
    return type
end

const _LABEL_MAPPER_PROPERTIES = ("name", "inputs", "mapper", "inputs_mapping", "no_label", "atol")

function _parse_label_mapper(node, path, registry; expected_label_type=nothing)
    node, context = _prepare_metadata(node, path, registry)
    metadata_tag(context) == "$_GWCS_TAG/label_mapper" || _metadata_error(path, "expected a GWCS label_mapper metadata node")
    mapper_value = _metadata_property(node, "mapper", context)
    inputs = _metadata_inputs(node, context)
    name = haskey(node, "name") ? _metadata_string(node, "name", context) : nothing
    has_no_label = haskey(node, "no_label")
    no_label = has_no_label ? metadata_value(node["no_label"]) : nothing
    path = _property_path(context, "mapper")
    materialized, _ = _array_metadata(mapper_value, path)

    if materialized isa AbstractMatrix || materialized isa AbstractVector
        labels = _metadata_matrix(mapper_value, path)
        if expected_label_type !== nothing && eltype(labels) != expected_label_type
            labels = try
                map(expected_label_type, labels)
            catch error
                error isa InexactError || error isa MethodError || error isa ArgumentError || rethrow()
                _metadata_error(path, "array labels cannot be converted from $(eltype(labels)) to $expected_label_type")
            end
        end
        ninputs = inputs === nothing ? 2 : length(inputs)
        mapping = _parse_mapper_mapping(node, context, inputs, 2)
        ninputs = max(ninputs, maximum(mapping))
        default_no_label = eltype(labels) <: Number ? zero(eltype(labels)) : eltype(labels) <: AbstractString ? "" : nothing
        has_no_label || (no_label = default_no_label)
        mapper = try
            ArrayLabelMapper{ninputs}(labels; inputs_mapping=mapping, no_label, inputs, name)
        catch error
            error isa ArgumentError || error isa DimensionMismatch || rethrow()
            _metadata_error(metadata_path(context), sprint(showerror, error))
        end
        return mapper
    end

    tagged = try
        _tag(mapper_value)
    catch error
        error isa ArgumentError || rethrow()
        nothing
    end
    if tagged !== nothing
        transform = parse_transform(mapper_value; path, registry=context.registry)
        inputs === nothing && (inputs = _asdf_input_names(transform))
        mapping = _parse_mapper_mapping(node, context, inputs, input_dimension(transform))
        if mapping != ntuple(identity, length(mapping)) || input_dimension(transform) != length(inputs)
            transform = compose(RemapAxes{length(inputs)}(mapping), transform)
        end
        L = expected_label_type === nothing ? Float64 : expected_label_type
        return TransformLabelMapper{L}(transform; no_label, inputs, name)
    end

    materialized isa AbstractDict || _metadata_error(path, "expected an array, transform, or label/model mapping")
    labels = _metadata_property(materialized, "labels", TransformParseContext(context.tag, context.version, path, context.registry))
    models = _metadata_property(materialized, "models", TransformParseContext(context.tag, context.version, path, context.registry))
    labels isa AbstractVector || _metadata_error("$path/labels", "expected an array")
    models isa AbstractVector || _metadata_error("$path/models", "expected an array")
    length(labels) == length(models) || _metadata_error(path, "label and model arrays must have equal lengths")
    inputs === nothing && _metadata_error(_property_path(context, "inputs"), "dictionary and range label mappers require inputs")
    parsed_models = map(enumerate(models)) do (index, model)
        parse_transform(model; path="$path/models/$index", registry=context.registry)
    end
    L = expected_label_type === nothing ? Float64 : expected_label_type
    first_label = metadata_value(first(labels))
    mapping = _parse_mapper_mapping(node, context, inputs, 1)
    mapper_entries = if first_label isa AbstractVector
        map(enumerate(labels), parsed_models) do index_label, model
            index, label = index_label
            length(label) == 2 || _metadata_error("$path/labels/$index", "range labels must contain two bounds")
            (_numeric_metadata(label[1], "$path/labels/$index/1"), _numeric_metadata(label[2], "$path/labels/$index/2")) => model
        end
    else
        map(enumerate(labels), parsed_models) do index_label, model
            index, label = index_label
            _numeric_metadata(label, "$path/labels/$index") => model
        end
    end
    try
        if first_label isa AbstractVector
            return RangeLabelMapper{length(inputs),L}(Tuple(mapper_entries); inputs_mapping=mapping, inputs, name)
        end
        atol = haskey(node, "atol") ? _numeric_property(node, "atol", context) : 1e-8
        return DictionaryLabelMapper{length(inputs),L}(Tuple(mapper_entries); inputs_mapping=mapping, atol, inputs, name)
    catch error
        error isa ArgumentError || error isa DimensionMismatch || rethrow()
        _metadata_error(metadata_path(context), sprint(showerror, error))
    end
end

function _parse_region_selector(node, context)
    path = metadata_path(context)
    selector = _metadata_property(node, "selector", context)
    selector isa AbstractDict || _metadata_error(_property_path(context, "selector"), "expected a mapping")
    labels = _metadata_property(selector, "labels", TransformParseContext(context.tag, context.version, "$path/selector", context.registry))
    transforms = _metadata_property(selector, "transforms", TransformParseContext(context.tag, context.version, "$path/selector", context.registry))
    labels isa AbstractVector || _metadata_error("$path/selector/labels", "expected an array")
    transforms isa AbstractVector || _metadata_error("$path/selector/transforms", "expected an array")
    length(labels) == length(transforms) || _metadata_error("$path/selector", "label and transform arrays must have equal lengths")
    label_values = map(metadata_value, labels)
    L = _metadata_label_type(label_values, "$path/selector/labels")
    mapper_node = _metadata_property(node, "label_mapper", context)
    mapper = _parse_label_mapper(mapper_node, "$path/label_mapper", context.registry; expected_label_type=L)
    parsed_transforms = map(enumerate(transforms)) do (index, transform)
        parse_transform(transform; path="$path/selector/transforms/$index", registry=context.registry)
    end
    entries = ntuple(index -> label_values[index] => parsed_transforms[index], length(label_values))
    undefined_value = haskey(node, "undefined_transform_value") ? _numeric_property(node, "undefined_transform_value", context) : NaN
    selector = try
        RegionSelector(mapper, entries; undefined_value)
    catch error
        error isa ArgumentError || error isa DimensionMismatch || rethrow()
        _metadata_error(path, sprint(showerror, error))
    end
    # Label mappers are not transforms, so the selector keeps its label mapper's unmodeled properties.
    pass_through = _add_pass_through!(_PassThrough(), "label_mapper", mapper_node, _LABEL_MAPPER_PROPERTIES)
    return isempty(pass_through) ? selector : AnnotatedTransform(selector; pass_through)
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

const _QUANTITY_TAG = "tag:stsci.edu:asdf/unit/quantity"

# Conversions from each angle unit to degrees and back, since the transform schemas define plain angles in degrees.
const _ANGLE_UNIT_CONVERSIONS = Dict{String,Tuple{Function,Function}}("deg" => (identity, identity), "rad" => (rad2deg, deg2rad),
    "arcmin" => (x -> x / 60, x -> x * 60), "arcsec" => (x -> x / 3600, x -> x * 3600), "mas" => (x -> x / 3.6e6, x -> x * 3.6e6),
    "hourangle" => (x -> 15x, x -> x / 15))

# Angles are converted by their units. astropy's projections pass parameter values to WCSLIB without converting them
# (`theta_b = 1.2 rad` acts as 1.2 degrees), so for a projection angle in another unit GWCS.jl deliberately differs
# from Python. A dimensionless parameter accepts only the empty unit: there is no standard set of other dimensionless
# units to convert, and astropy would read `mu = 200 %` as 200.
function _quantity_conversion(kind, unit, path)
    if kind === :angle
        haskey(_ANGLE_UNIT_CONVERSIONS, unit) || _metadata_error(path, "expected an angle quantity, got unit $(repr(unit))")
        return _ANGLE_UNIT_CONVERSIONS[unit]
    end
    isempty(unit) || _metadata_error(path, "expected a dimensionless quantity with an empty unit, got unit $(repr(unit))")
    return (identity, identity)
end

# Converts the quantity values of the registration's quantity properties to plain numbers before validation, and
# returns the node the constructor sees with the units to record. Other tags fail validation as before.
function _normalize_quantities(registration, node, context)
    units = _ParameterUnit[]
    isempty(registration.quantities) && return node, units
    replacements = Pair{Any,Any}[]
    for (property, kind) in registration.quantities
        value = get(node, property, nothing)
        value !== nothing && _has_tag(value) || continue
        path = _property_path(context, property)
        family, version = _split_metadata_tag(_tag(value), path)
        family == _QUANTITY_TAG || continue
        _lookup_metadata_schema(context.registry, family, version, path)
        unit_node = _required(value, "unit", path)
        unit = String(_unwrap(unit_node))
        to_degrees, _ = _quantity_conversion(kind, unit, "$path/unit")
        raw = _materialize_array(_unwrap(_required(value, "value", path)))
        raw isa Real || raw isa AbstractArray{<:Real} || _metadata_error("$path/value", "expected a number or numeric array")
        datatype = get(value, "datatype", nothing)
        push!(units, _ParameterUnit(registration.tag, property, kind, unit, _tag(value), _has_tag(unit_node) ? _tag(unit_node) : "$_UNIT_TAG-1.0.0", _unwrap(datatype) isa AbstractString ? String(_unwrap(datatype)) : nothing))
        push!(replacements, property => (raw isa Real ? to_degrees(raw) : to_degrees.(collect(Float64, raw))))
    end
    isempty(replacements) && return node, units
    values = Dict{Any,Any}(node)
    foreach(((key, item),) -> values[key] = item, replacements)
    return TaggedMapping(_tag(node), values), units
end

# The properties every transform has, which the annotation and `parse_transform` model.
const _TRANSFORM_PROPERTIES = ("name", "inputs", "outputs", "inverse")

# A frame2d is written for a two-axis spatial `frame`, so the two share their unmodeled properties.
_pass_through_family(tag) = replace(first(_split_metadata_tag(tag, "")), r"/frame2d$" => "/frame")

# Unmodeled values are kept as GWCS.jl's own metadata nodes, without ASDF.jl's lazy arrays, which hold the file open.
function _pass_through_value(value::AbstractDict)
    items = Dict{Any,Any}(_pass_through_value(key) => _pass_through_value(item) for (key, item) in value)
    return _has_tag(value) ? TaggedMapping(_tag(value), items) : items
end

function _pass_through_value(value::AbstractVector)
    items = Any[_pass_through_value(item) for item in value]
    return _has_tag(value) ? TaggedSequence(_tag(value), items) : items
end

_pass_through_value(value::AbstractString) = _has_tag(value) ? TaggedScalar(_tag(value), String(value)) : String(value)
_pass_through_value(value) = value

# Records the properties of `node` outside `modeled` under `relative`, the node's path below the owner's.
function _add_pass_through!(pass_through, relative, node, modeled)
    properties = Dict{Any,Any}(_unwrap(key) => _pass_through_value(value) for (key, value) in node if !(_unwrap(key) in modeled))
    isempty(properties) || (pass_through[relative] = (_pass_through_family(_tag(node)), properties))
    return pass_through
end

function _parse_transform_annotation(@nospecialize(transform), node, context, units=_ParameterUnit[], pass_through=_NO_PASS_THROUGH)
    has_annotation = haskey(node, "name") || haskey(node, "inputs") || haskey(node, "outputs") || !isempty(units) || !isempty(pass_through)
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
    return annotate_transform(transform; name, inputs, outputs, parameter_units=units, pass_through)
end

"""
    parse_transform(node; path="transform", registry=default_transform_registry()) -> AbstractTransform
    parse_transform(node, context::TransformParseContext; path=metadata_path(context)) -> AbstractTransform

Build a transform from a tagged metadata node. Its tag and version select the constructor registered in `registry`.
An `inverse` property becomes an explicit inverse ([`WithInverse`](@ref)), and `name`, `inputs`, and `outputs` become
annotations ([`AnnotatedTransform`](@ref)). Malformed metadata or an unsupported tag throws an `ArgumentError` naming
`path`. A registered constructor parses a child node with the second form, which keeps its context's registry.
"""
function parse_transform(node; path="transform", registry=default_transform_registry())
    tag, version = _split_metadata_tag(_tag(node), path)
    registration = _lookup_transform_registration(registry, tag, version, path)
    context = TransformParseContext(tag, version, String(path), registry)
    normalized, units = _normalize_quantities(registration, node, context)
    prepared = _prepare_transform_metadata(registration, normalized, context)
    transform = registration.constructor(prepared, context)
    transform isa AbstractTransform || _metadata_error(path, "registered constructor for $tag returned $(typeof(transform)), expected AbstractTransform")
    if haskey(node, "inverse")
        backward = parse_transform(node["inverse"]; path="$path/inverse", registry)
        input_dimension(transform) == output_dimension(backward) && output_dimension(transform) == input_dimension(backward) || _metadata_error(path, "explicit inverse dimensions do not match the forward transform")
        transform = WithInverse(transform, backward)
    end
    pass_through = _pass_through(transform)
    if registration.properties !== nothing
        pass_through = _add_pass_through!(copy(pass_through), "", node, (_TRANSFORM_PROPERTIES..., registration.properties...))
    end
    return _parse_transform_annotation(transform, node, context, units, pass_through)
end

function parse_transform(node, context::TransformParseContext; path=metadata_path(context))
    return parse_transform(node; path, registry=context.registry)
end

function _register_builtin_transform_parsers!()
    asdf_transform = "tag:stsci.edu:asdf/transform"
    register_transform!("$asdf_transform/compose", _parse_compose; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=_validate_compound_transform, properties=("forward",))
    register_transform!("$asdf_transform/concatenate", _parse_concatenate; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=_validate_compound_transform, properties=("forward",))
    for (name, parser) in (("add", _parse_add), ("subtract", _parse_subtract), ("multiply", _parse_multiply), ("divide", _parse_divide), ("power", _parse_power))
        register_transform!("$asdf_transform/$name", parser; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=_validate_binary_transform, properties=("forward",))
    end
    register_transform!("$asdf_transform/shift", _parse_shift; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=(node, context) -> _validate_numeric_transform(node, context, "offset"), properties=("offset",))
    register_transform!("$asdf_transform/scale", _parse_scale; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=(node, context) -> _validate_numeric_transform(node, context, "factor"), properties=("factor",))
    register_transform!("$asdf_transform/rotate2d", _parse_rotate2d; versions=_SCHEMA_VERSIONS_1_0_TO_1_5, validate=(node, context) -> _validate_numeric_transform(node, context, "angle"), angles=("angle",), properties=("angle",))
    register_transform!("$asdf_transform/identity", _parse_identity; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=_validate_identity_schema, properties=("n_dims",))
    register_transform!("$asdf_transform/affine", _parse_affine; versions=_SCHEMA_VERSIONS_1_0_TO_1_5, validate=_validate_affine_schema, properties=("matrix", "translation"))
    register_transform!("$asdf_transform/gnomonic", _parse_gnomonic; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=_validate_gnomonic_schema, properties=("direction",))
    register_transform!("$asdf_transform/zenithal_perspective", _parse_zenithal_perspective; versions=_SCHEMA_VERSIONS_1_0_TO_1_5, validate=(node, context) -> _validate_projection_schema(node, context, ("mu", "gamma")), angles=("gamma",), dimensionless=("mu",), properties=("direction", "mu", "gamma"))
    register_transform!("$asdf_transform/slant_zenithal_perspective", _parse_slant_zenithal_perspective; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=(node, context) -> _validate_projection_schema(node, context, ("mu", "phi0", "theta0")), angles=("phi0", "theta0"), dimensionless=("mu",), properties=("direction", "mu", "phi0", "theta0"))
    register_transform!("$asdf_transform/stereographic", _parse_stereographic; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=_validate_gnomonic_schema, properties=("direction",))
    register_transform!("$asdf_transform/slant_orthographic", _parse_slant_orthographic; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=(node, context) -> _validate_projection_schema(node, context, ("xi", "eta")), dimensionless=("xi", "eta"), properties=("direction", "xi", "eta"))
    register_transform!("$asdf_transform/zenithal_equidistant", _parse_zenithal_equidistant; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=_validate_gnomonic_schema, properties=("direction",))
    register_transform!("$asdf_transform/zenithal_equal_area", _parse_zenithal_equal_area; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=_validate_gnomonic_schema, properties=("direction",))
    register_transform!("$asdf_transform/airy", _parse_airy; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=(node, context) -> _validate_projection_schema(node, context, ("theta_b",)), angles=("theta_b",), properties=("direction", "theta_b"))
    for (name, projection) in (("plate_carree", PlateCarreeProjection()), ("mercator", MercatorProjection()), ("sanson_flamsteed", SansonFlamsteedProjection()),
                               ("parabolic", ParabolicProjection()), ("molleweide", MollweideProjection()), ("hammer_aitoff", HammerAitoffProjection()),
                               ("polyconic", PolyconicProjection()), ("tangential_spherical_cube", TangentialSphericalCube()), ("cobe_quad_spherical_cube", COBESphericalCube()),
                               ("quad_spherical_cube", QuadrilateralizedSphericalCube()), ("healpix_polar", HEALPixPolarProjection()))
        register_transform!("$asdf_transform/$name", (node, context) -> _parse_projection(node, context, projection); versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=_validate_gnomonic_schema, properties=("direction",))
    end
    register_transform!("$asdf_transform/cylindrical_equal_area", _parse_cylindrical_equal_area; versions=_SCHEMA_VERSIONS_1_0_TO_1_5, validate=(node, context) -> _validate_projection_schema(node, context, ("lam", "lambda")), dimensionless=("lam", "lambda"), properties=("direction", "lam", "lambda"))
    register_transform!("$asdf_transform/cylindrical_perspective", _parse_cylindrical_perspective; versions=_SCHEMA_VERSIONS_1_0_TO_1_5, validate=(node, context) -> _validate_projection_schema(node, context, ("mu", "lam", "lambda")), dimensionless=("mu", "lam", "lambda"), properties=("direction", "mu", "lam", "lambda"))
    for (name, constructor) in (("conic_perspective", ConicPerspective), ("conic_equidistant", ConicEquidistant), ("conic_equal_area", ConicEqualArea), ("conic_orthomorphic", ConicOrthomorphic))
        register_transform!("$asdf_transform/$name", _conic_parser(constructor); versions=_SCHEMA_VERSIONS_1_0_TO_1_5, validate=(node, context) -> _validate_projection_schema(node, context, ("sigma", "delta")), angles=("sigma", "delta"), properties=("direction", "sigma", "delta"))
    end
    register_transform!("$asdf_transform/bonne_equal_area", _parse_bonne_equal_area; versions=_SCHEMA_VERSIONS_1_0_TO_1_5, validate=(node, context) -> _validate_projection_schema(node, context, ("theta1",)), angles=("theta1",), properties=("direction", "theta1"))
    register_transform!("$asdf_transform/healpix", _parse_healpix; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=(node, context) -> _validate_projection_schema(node, context, ("H", "X")), properties=("direction", "H", "X"))
    register_transform!("$asdf_transform/rotate3d", _parse_rotate3d; versions=_SCHEMA_VERSIONS_1_0_TO_1_5, validate=_validate_rotate3d_schema, angles=("phi", "theta", "psi"), properties=("direction", "phi", "theta", "psi"))
    register_transform!("$asdf_transform/constant", _parse_constant; versions=_SCHEMA_VERSIONS_1_0_TO_1_6, validate=_validate_constant_schema, normalize=_normalize_constant_schema, properties=("value", "dimensions"))
    register_transform!("$asdf_transform/fix_inputs", _parse_fixed_inputs; versions=(v"1.1.0", v"1.2.0", v"1.3.0", v"1.4.0"), validate=_validate_fixed_inputs_schema, properties=("forward",))
    register_transform!("$asdf_transform/remap_axes", _parse_remap; versions=_SCHEMA_VERSIONS_1_0_TO_1_5, validate=_validate_remap_schema, properties=("mapping", "n_inputs"))
    register_transform!("$asdf_transform/polynomial", _parse_polynomial; versions=_SCHEMA_VERSIONS_1_0_TO_1_3, validate=_validate_polynomial_schema, properties=("coefficients", "domain", "window"))
    register_transform!("$asdf_transform/rotate_sequence_3d", _parse_rotation; versions=_SCHEMA_VERSIONS_1_0_TO_1_2, validate=_validate_rotation_sequence_schema, angles=("angles",), properties=("rotation_type", "angles", "axes_order"))
    register_transform!("$_GWCS_TAG/spherical_cartesian", _parse_spherical_cartesian; versions=_SCHEMA_VERSIONS_1_0_TO_1_3, validate=_validate_spherical_cartesian_schema, properties=("transform_type", "wrap_lon_at"))
    register_transform!("$_GWCS_TAG/regions_selector", _parse_region_selector; versions=_SCHEMA_VERSIONS_1_0_TO_1_3, validate=_validate_region_selector_schema, properties=("label_mapper", "selector", "undefined_transform_value"))
    return nothing
end

_register_builtin_transform_parsers!()
_register_builtin_metadata_schemas!()

function _parse_reference_frame(value, path, registry)
    value === nothing && return nothing
    tag = try
        _tag(value)
    catch error
        error isa ArgumentError || rethrow()
        return _unwrap(value)
    end
    family, version = _split_metadata_tag(tag, path)
    if family == "tag:astropy.org:astropy/coordinates/frames/icrs"
        _lookup_metadata_schema(registry, family, version, path)
        # An ICRS frame that carries more, such as the position data of a `wcs_from_fiducial` frame, is kept as stored,
        # like the other astropy frames, so writing reproduces it.
        all(key -> key == "frame_attributes" && isempty(value[key]), keys(value)) && return :ICRS
    end
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

function _temporal_reference(node, path, registry)
    node === nothing && return nothing, nothing, "", "", nothing
    reference = _parse_reference_frame(node, path, registry)
    node isa AbstractDict || return reference, _unwrap(node), "", "", nothing
    epoch = haskey(node, "value") ? _unwrap(node["value"]) : nothing
    scale = haskey(node, "scale") ? String(_unwrap(node["scale"])) : ""
    format = haskey(node, "format") ? String(_unwrap(node["format"])) : ""
    location = get(node, "location", nothing)
    return reference, epoch, scale, format, location
end

function _parse_composite_frame(node, path, registry)
    # Python GWCS rejects any other property here.
    length(node) == 2 || _metadata_error(path, "composite frame must contain only name and frames")
    frame_name = String(_unwrap(_required(node, "name", path)))
    values = _required(node, "frames", path)
    frames = map(enumerate(values)) do (index, value)
        _parse_frame(value, "$path/frames/$index"; validate_order=false, registry)
    end
    try
        return CompositeFrame(frames; name=frame_name)
    catch error
        error isa ArgumentError || error isa DimensionMismatch || rethrow()
        _metadata_error(path, sprint(showerror, error))
    end
end

function _parse_frame(node, path; validate_order=true, registry=default_transform_registry())
    node, context = _prepare_metadata(node, path, registry)
    name = split(metadata_tag(context), '/')[end]
    name == "composite_frame" && return _parse_composite_frame(node, path, registry)
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
            reference_frame = _parse_reference_frame(get(node, "reference_frame", nothing), "$path/reference_frame", registry)
            CoordinateFrame{N}(; name=frame_name, axes_names, axes_order, axis_types, axis_physical_types=physical_types, unit=units, reference_frame)
        elseif name == "frame2d"
            # The frame schemas do not define axes_type. Python GWCS's frame reader applies it only with naxes, which its
            # Frame2D does not accept, so a frame2d's axes are always spatial there; GWCS.jl ignores axes_type likewise
            # and writes it back unchanged.
            axes_names = _metadata_frame_strings(node, "axes_names", ("x", "y"))
            axes_order = _frame_order(node, 2)
            physical_types = haskey(node, "axis_physical_types") ? _metadata_frame_strings(node, "axis_physical_types", ()) : nothing
            units = _metadata_frame_strings(node, "unit", ("", ""))
            reference_frame = _parse_reference_frame(get(node, "reference_frame", nothing), "$path/reference_frame", registry)
            Frame2D(; name=frame_name, axes_names, axes_order, axis_physical_types=physical_types, unit=units, reference_frame)
        elseif name == "celestial_frame"
            N = haskey(node, "axes_names") ? length(_metadata_sequence(node["axes_names"])) :
                haskey(node, "axes_order") ? length(_metadata_sequence(node["axes_order"])) :
                haskey(node, "unit") ? length(_metadata_sequence(node["unit"])) : 2
            axes_names = _metadata_frame_strings(node, "axes_names", ntuple(i -> i == 1 ? "lon" : i == 2 ? "lat" : "distance", N))
            axes_order = _frame_order(node, N)
            physical_types = haskey(node, "axis_physical_types") ? _metadata_frame_strings(node, "axis_physical_types", ()) : nothing
            units = _metadata_frame_strings(node, "unit", ntuple(i -> i <= 2 ? "deg" : "", N))
            reference_frame = _parse_reference_frame(get(node, "reference_frame", nothing), "$path/reference_frame", registry)
            CelestialFrame{N}(; name=frame_name, axes_names, axes_order, axis_physical_types=physical_types, unit=units, reference_frame)
        elseif name == "spectral_frame"
            axes_names = _metadata_frame_strings(node, "axes_names", ("spectral",))
            axes_order = _frame_order(node, 1)
            physical_types = haskey(node, "axis_physical_types") ? _metadata_frame_strings(node, "axis_physical_types", ()) : nothing
            units = _metadata_frame_strings(node, "unit", ("",))
            reference_frame = _parse_reference_frame(get(node, "reference_frame", nothing), "$path/reference_frame", registry)
            reference_position = String(_unwrap(get(node, "reference_position", "geocenter")))
            SpectralFrame(; name=frame_name, axes_names, axes_order, axis_physical_types=physical_types, unit=units, reference_frame, reference_position)
        elseif name == "temporal_frame"
            axes_names = _metadata_frame_strings(node, "axes_names", ("time",))
            axes_order = _frame_order(node, 1)
            physical_types = _metadata_frame_strings(node, "axis_physical_types", ("time",))
            units = _metadata_frame_strings(node, "unit", ("s",))
            reference, epoch, scale, format, location = _temporal_reference(get(node, "reference_frame", nothing), "$path/reference_frame", registry)
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

# The properties each frame tag's reader and writer model.
const _FRAME_PROPERTIES = Dict("frame" => ("name", "naxes", "axes_type", "axes_names", "axes_order", "axis_physical_types", "unit", "reference_frame"),
    "frame2d" => ("name", "axes_names", "axes_order", "axis_physical_types", "unit", "reference_frame"),
    "celestial_frame" => ("name", "axes_names", "axes_order", "axis_physical_types", "unit", "reference_frame"),
    "spectral_frame" => ("name", "axes_names", "axes_order", "axis_physical_types", "unit", "reference_frame", "reference_position"),
    "temporal_frame" => ("name", "axes_names", "axes_order", "axis_physical_types", "unit", "reference_frame"),
    "stokes_frame" => ("name", "axes_names", "axes_order", "axis_physical_types"), "composite_frame" => ("name", "frames"))

function _frame_pass_through!(pass_through, node, relative)
    node isa AbstractDict || return pass_through
    name = split(first(_split_metadata_tag(_tag(node), relative)), '/')[end]
    _add_pass_through!(pass_through, relative, node, _FRAME_PROPERTIES[name])
    if name == "composite_frame"
        foreach(((index, frame),) -> _frame_pass_through!(pass_through, frame, "$relative/frames/$index"), enumerate(node["frames"]))
    end
    return pass_through
end

# A frame written as a bare string is an EmptyFrame, sized by the step's transform or, for the last step, by the
# previous step's outputs (`previous`).
function _parse_step(node, path, registry, previous=nothing)
    node, context = _prepare_metadata(node, path, registry)
    metadata_tag(context) == "$_GWCS_TAG/step" || _metadata_error(path, "expected a GWCS step")
    frame_node = _unwrap(_required(node, "frame", path))
    transform_node = get(node, "transform", nothing)
    transform = transform_node === nothing ? nothing : parse_transform(transform_node; path="$path/transform", registry)
    if frame_node isa AbstractString
        sizing = transform === nothing ? previous : transform
        sizing === nothing && _metadata_error("$path/frame", "a frame given only by name needs a transform to set its dimension")
        frame = EmptyFrame{transform === nothing ? output_dimension(sizing) : input_dimension(sizing)}(frame_node)
    else
        frame = _parse_frame(frame_node, "$path/frame"; registry)
    end
    pass_through = _frame_pass_through!(_add_pass_through!(_PassThrough(), "", node, ("frame", "transform")), frame_node, "frame")
    return Step(frame, transform, isempty(pass_through) ? _NO_PASS_THROUGH : pass_through)
end

function _bounding_box_order(node, path)
    value = String(_unwrap(get(node, "order", "C")))
    value in ("C", "F") || _metadata_error("$path/order", "expected one of C, F, got $value")
    return Symbol(value)
end

function _bounding_box_names(transform)
    names = input_names(transform)
    return names === nothing ? _asdf_input_names(transform) : names
end

function _bounding_box_axis(name, names, path)
    value = String(_unwrap(name))
    indices = findall(==(value), names)
    isempty(indices) && _metadata_error(path, "bounding-box input $value is not among $names")
    length(indices) == 1 || _metadata_error(path, "bounding-box input $value is ambiguous in $names")
    return only(indices)
end

function _bounding_box_name_axes(values, names, path)
    values isa AbstractVector || _metadata_error(path, "expected an array of input names")
    axes = map(enumerate(values)) do (index, value)
        _bounding_box_axis(value, names, "$path/$index")
    end
    length(unique(axes)) == length(axes) || _metadata_error(path, "bounding-box input names must be unique")
    return Tuple(axes)
end

function _parse_simple_bounding_box(node, names, path; implicit_ignored=())
    _bounding_box_order(node, path)
    raw_intervals = _required(node, "intervals", path)
    raw_intervals isa AbstractDict || _metadata_error("$path/intervals", "expected a mapping of input names to intervals")
    isempty(raw_intervals) && _metadata_error("$path/intervals", "bounding-box intervals must not be empty")
    axes = Int[]
    intervals = Tuple{Real,Real}[]
    for (key, raw_interval) in raw_intervals
        name = String(_unwrap(key))
        interval_path = "$path/intervals/$name"
        axis = _bounding_box_axis(key, names, interval_path)
        axis in axes && _metadata_error(interval_path, "bounding-box interval input names must be unique")
        raw_interval isa AbstractVector || raw_interval isa Tuple || _metadata_error(interval_path, "expected a two-element interval")
        length(raw_interval) == 2 || _metadata_error(interval_path, "expected a two-element interval")
        lower = _numeric_metadata(raw_interval[1], "$interval_path/1")
        upper = _numeric_metadata(raw_interval[2], "$interval_path/2")
        lower <= upper || _metadata_error(interval_path, "bounding-box lower limit must not exceed upper limit")
        push!(axes, axis)
        push!(intervals, (lower, upper))
    end
    ignored = haskey(node, "ignore") ? _bounding_box_name_axes(node["ignore"], names, "$path/ignore") : ()
    all_ignored = Tuple(unique((Tuple(implicit_ignored)..., ignored...)))
    overlap = intersect(axes, all_ignored)
    isempty(overlap) || _metadata_error(path, "bounding-box inputs cannot be both bounded and ignored")
    accounted = union(Set(axes), Set(all_ignored))
    accounted == Set(1:length(names)) || _metadata_error(path, "bounding-box intervals and ignored inputs must account for every transform input")
    try
        return BoundingBox{length(names)}(Tuple(intervals), Tuple(axes))
    catch error
        error isa ArgumentError || error isa DimensionMismatch || rethrow()
        _metadata_error(path, sprint(showerror, error))
    end
end

function _parse_selector_args(node, names, path)
    raw = _required(node, "selector_args", path)
    raw isa AbstractVector || _metadata_error("$path/selector_args", "expected an array of selector arguments")
    isempty(raw) && _metadata_error("$path/selector_args", "selector arguments must not be empty")
    selectors = map(enumerate(raw)) do (index, value)
        item_path = "$path/selector_args/$index"
        value isa AbstractDict || _metadata_error(item_path, "expected a selector-argument mapping")
        axis = _bounding_box_axis(_required(value, "argument", item_path), names, "$item_path/argument")
        ignore = _unwrap(_required(value, "ignore", item_path))
        ignore isa Bool || _metadata_error("$item_path/ignore", "expected a boolean")
        (axis, ignore)
    end
    axes = Tuple(first.(selectors))
    length(unique(axes)) == length(axes) || _metadata_error("$path/selector_args", "selector input names must be unique")
    return axes, Tuple(last.(selectors))
end

function _parse_compound_bounding_box(node, names, path, registry)
    _bounding_box_order(node, path)
    selectors, selector_flags = _parse_selector_args(node, names, path)
    globally_ignored = haskey(node, "ignore") ? _bounding_box_name_axes(node["ignore"], names, "$path/ignore") : ()
    isempty(intersect(selectors, globally_ignored)) || _metadata_error(path, "selector inputs cannot also be globally ignored")
    implicit_ignored = (globally_ignored..., Tuple(selectors[index] for index in eachindex(selectors) if selector_flags[index])...)
    raw_entries = _required(node, "cbbox", path)
    raw_entries isa AbstractVector || _metadata_error("$path/cbbox", "expected an array of compound bounding-box entries")
    isempty(raw_entries) && _metadata_error("$path/cbbox", "compound bounding-box entries must not be empty")
    entries = Pair[]
    seen = Set()
    for (index, value) in enumerate(raw_entries)
        item_path = "$path/cbbox/$index"
        value isa AbstractDict || _metadata_error(item_path, "expected a compound bounding-box entry mapping")
        raw_key = _required(value, "key", item_path)
        raw_key isa AbstractVector || raw_key isa Tuple || _metadata_error("$item_path/key", "expected an array of selector values")
        length(raw_key) == length(selectors) || _metadata_error("$item_path/key", "selector key has length $(length(raw_key)), expected $(length(selectors))")
        key = Tuple(_numeric_metadata(item, "$item_path/key/$key_index") for (key_index, item) in enumerate(raw_key))
        key in seen && _metadata_error("$item_path/key", "compound bounding-box selector keys must be unique")
        push!(seen, key)
        nested = _required(value, "bbox", item_path)
        nested, nested_context = _prepare_metadata(nested, "$item_path/bbox", registry)
        metadata_tag(nested_context) == "$_ASDF_TRANSFORM_TAG/property/bounding_box" || _metadata_error("$item_path/bbox", "expected a simple bounding box")
        box = _parse_simple_bounding_box(nested, names, "$item_path/bbox"; implicit_ignored)
        push!(entries, key => box)
    end
    try
        return CompoundBoundingBox{length(names)}(Dict(entries), selectors; selector_ignored=selector_flags, ignored=globally_ignored)
    catch error
        error isa ArgumentError || error isa DimensionMismatch || rethrow()
        _metadata_error(path, sprint(showerror, error))
    end
end

function _parse_bounding_box(node, transform, path, registry)
    node, context = _prepare_metadata(node, path, registry)
    names = _bounding_box_names(transform)
    metadata_tag(context) == "$_ASDF_TRANSFORM_TAG/property/bounding_box" && return _parse_simple_bounding_box(node, names, path)
    metadata_tag(context) == "$_ASDF_TRANSFORM_TAG/property/compound_bounding_box" && return _parse_compound_bounding_box(node, names, path, registry)
    _metadata_error(path, "unsupported bounding-box tag $(_tag(node))")
end

"""
    parse_wcs(node::AbstractDict; registry=default_transform_registry()) -> WCS
    parse_wcs(yaml::AbstractString; registry=default_transform_registry()) -> WCS

Build a [`WCS`](@ref) from a `tag:stsci.edu:gwcs/wcs` metadata node, such as `ASDF.load_file(filename)["wcs"]` with
ASDF.jl loaded, or from YAML text holding one. An astropy FITS WCS stored with asdf-astropy's `wcs` tag is read
through [`parse_fits_wcs`](@ref). Transforms are parsed by the tags registered in `registry`; see
[`register_transform!`](@ref).

Metadata that is malformed, uses an unsupported tag or schema version, or describes an invalid pipeline throws an
`ArgumentError` naming its path in the tree. Properties GWCS.jl does not model are kept and written back by
[`write_asdf`](@ref).
"""
function parse_wcs(node::AbstractDict; registry=default_transform_registry())
    node, context = _prepare_metadata(node, "wcs", registry)
    metadata_tag(context) == _ASTROPY_WCS_TAG && return _parse_astropy_wcs(node, registry)
    metadata_tag(context) == "$_GWCS_TAG/wcs" || _metadata_error("wcs", "expected a GWCS WCS or astropy FITS WCS metadata node")
    raw_steps = _required(node, "steps", "wcs")
    steps = Any[]
    for (index, value) in enumerate(raw_steps)
        push!(steps, _parse_step(value, "wcs/steps/$index", registry, isempty(steps) ? nothing : last(steps).transform))
    end
    first_transform = get(first(raw_steps), "transform", nothing)
    bounding_box = first_transform === nothing ? nothing : get(first_transform, "bounding_box", nothing)
    parsed_box = bounding_box === nothing ? nothing : _parse_bounding_box(bounding_box, first(steps).transform, "wcs/steps/1/transform/bounding_box", registry)
    name = String(_unwrap(get(node, "name", "")))
    pixel_shape = get(node, "pixel_shape", nothing)
    pixel_shape === nothing || (pixel_shape = Tuple(Int(_unwrap(value)) for value in pixel_shape))
    pass_through = _add_pass_through!(_PassThrough(), "", node, ("name", "steps", "pixel_shape"))
    return _with_pass_through(WCS(steps...; name, pixel_shape, bounding_box=parsed_box), pass_through)
end

parse_wcs(metadata::AbstractString; registry=default_transform_registry()) = parse_wcs(_load_tagged_yaml(metadata); registry)

# asdf-astropy does not record an alternate key, since astropy's WCS has no `key` attribute, so astropy cannot read back
# an alternate-key WCS it wrote. `WCS.to_fits` writes only the WCS's own representation, which WCSLIB always begins with
# WCSAXESa, so that card identifies the key.
function _astropy_wcs_key(header)
    found = sort!([isempty(m.captures[1]) ? ' ' : only(m.captures[1]) for m in (match(r"^WCSAXES([A-Z]?)$", name) for name in keys(header)) if m !== nothing])
    length(found) <= 1 || _metadata_error("wcs/hdulist", "the header holds WCS representations with keys $(join(repr.(found), ", ")) and attrs names none")
    return isempty(found) ? ' ' : only(found)
end

# asdf-astropy writes an astropy FITS WCS as the HDU list `WCS.to_fits` returns plus the attributes its header lacks, and
# reads it back as `WCS(header, fobj=hdulist; attrs...)`. The header goes to `parse_fits_wcs`, so it follows WCSLIB, and
# APE 14's `pixel_bounds`, which Python GWCS reports from its bounding box, becomes the bounding box. Lookup and -TAB
# arrays would come from the other HDUs, which are not read yet.
function _parse_astropy_wcs(node, registry)
    hdulist, _ = _prepare_metadata(_required(node, "hdulist", "wcs"), "wcs/hdulist", registry)
    attrs = _required(node, "attrs", "wcs")
    for name in ("colsel", "keysel")
        _unwrap(get(attrs, name, nothing)) === nothing || _metadata_error("wcs/attrs/$name", "astropy FITS WCS attribute $name is not supported")
    end
    length(hdulist) == 1 || _metadata_error("wcs/hdulist", "astropy FITS WCS extension HDUs, which hold lookup distortions and -TAB arrays, are not supported")
    header = _fits_card_dictionary((_unwrap(card[1]), length(card) > 1 ? _unwrap(card[2]) : "") for card in first(hdulist)["header"])
    # Like asdf-astropy, a nonzero attribute replaces the header's NAXIS, which `WCS.to_fits` writes as 0.
    naxis = _unwrap(get(attrs, "naxis", nothing))
    naxis === nothing || iszero(naxis) || (header["NAXIS"] = Int(naxis))
    key = _unwrap(get(attrs, "key", nothing))
    wcs = parse_fits_wcs(header; key=key === nothing ? _astropy_wcs_key(header) : only(String(key)))
    shape = _unwrap(get(attrs, "pixel_shape", nothing))
    shape === nothing || (wcs = with_pixel_shape(wcs, Tuple(Int(_unwrap(value)) for value in shape)))
    bounds = _unwrap(get(attrs, "pixel_bounds", nothing))
    bounds === nothing || (wcs = with_bounding_box(wcs, BoundingBox(Tuple((Float64(_unwrap(bound[1])), Float64(_unwrap(bound[2]))) for bound in bounds))))
    return wcs
end
