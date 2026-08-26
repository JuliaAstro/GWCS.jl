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

function _required(node, key, path)
    haskey(node, key) || _metadata_error(path, "missing required property $key")
    return node[key]
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

function _parse_compose(node, path)
    values = _required(node, "forward", path)
    transforms = map(enumerate(values)) do (index, value)
        parse_transform(value; path="$path/forward/$index")
    end
    return compose(transforms...)
end

function _parse_concatenate(node, path)
    values = _required(node, "forward", path)
    transforms = map(enumerate(values)) do (index, value)
        parse_transform(value; path="$path/forward/$index")
    end
    return parallel(transforms...)
end

_parse_shift(node, path) = Shift(_unwrap(_required(node, "offset", path)))
_parse_scale(node, path) = Scale(_unwrap(_required(node, "factor", path)))

function _parse_remap(node, path)
    mapping = _tuple_ints(_required(node, "mapping", path); zero_based=true)
    n_inputs = Int(get(node, "n_inputs", maximum(mapping)))
    return RemapAxes{n_inputs}(mapping)
end

function _parse_polynomial(node, path)
    coefficients = _coefficient_array(_required(node, "coefficients", path))
    domain, window = _domain_window(node)
    ndims(coefficients) == 1 && return Polynomial1D(coefficients; domain, window)
    ndims(coefficients) == 2 && return Polynomial2D(coefficients; domain, window)
    _metadata_error(path, "only one- and two-dimensional polynomials are supported")
end

function _parse_rotation(node, path)
    rotation_type = String(_unwrap(_required(node, "rotation_type", path)))
    rotation_type == "cartesian" || _metadata_error(path, "only cartesian rotation sequences are supported")
    angles = collect(_required(node, "angles", path))
    axes_order = String(_unwrap(_required(node, "axes_order", path)))
    return Rotation3D(angles, axes_order)
end

function _parse_spherical_cartesian(node, path)
    transform_type = String(_unwrap(_required(node, "transform_type", path)))
    wrap_lon_at = Int(_unwrap(_required(node, "wrap_lon_at", path)))
    transform_type == "spherical_to_cartesian" && return SphericalToCartesian(wrap_lon_at)
    transform_type == "cartesian_to_spherical" && return CartesianToSpherical(wrap_lon_at)
    _metadata_error(path, "unsupported spherical/cartesian transform type $transform_type")
end

const _TRANSFORM_PARSERS = Dict{String,Function}(
    "compose" => _parse_compose,
    "concatenate" => _parse_concatenate,
    "shift" => _parse_shift,
    "scale" => _parse_scale,
    "remap_axes" => _parse_remap,
    "polynomial" => _parse_polynomial,
    "rotate_sequence_3d" => _parse_rotation,
    "spherical_cartesian" => _parse_spherical_cartesian
)

function parse_transform(node; path="transform")
    name = _tag_name(node)
    parser = get(_TRANSFORM_PARSERS, name, nothing)
    parser === nothing && _metadata_error(path, "unsupported transform tag $(_tag(node))")
    transform = parser(node, path)
    if haskey(node, "inverse")
        backward = parse_transform(node["inverse"]; path="$path/inverse")
        input_dimension(transform) == output_dimension(backward) && output_dimension(transform) == input_dimension(backward) || _metadata_error(path, "explicit inverse dimensions do not match the forward transform")
        transform = WithInverse(transform, backward)
    end
    return transform
end

function _parse_reference_frame(value)
    value === nothing && return nothing
    value isa AbstractDict || return _unwrap(value)
    _tag_name(value) == "icrs" && return :ICRS
    return _tag(value)
end

function _parse_frame(node, path)
    name = _tag_name(node)
    frame_name = String(_unwrap(_required(node, "name", path)))
    axes_names = _tuple_strings(_required(node, "axes_names", path))
    N = length(axes_names)
    axes_order = _tuple_ints(get(node, "axes_order", 0:(N - 1)); zero_based=true)
    physical_types = _tuple_strings(get(node, "axis_physical_types", ntuple(i -> "custom:$(axes_names[i])", N)))
    units = _tuple_strings(get(node, "unit", ntuple(i -> "", N)))
    reference_frame = _parse_reference_frame(get(node, "reference_frame", nothing))
    name == "frame2d" && N == 2 && return Frame2D(; name=frame_name, axes_names, axes_order, axis_physical_types=physical_types, unit=units, reference_frame)
    name == "celestial_frame" && return CelestialFrame{N}(; name=frame_name, axes_names, axes_order, axis_physical_types=physical_types, unit=units, reference_frame)
    _metadata_error(path, "unsupported coordinate frame tag $(_tag(node))")
end

function _parse_step(node, path)
    _tag_name(node) == "step" || _metadata_error(path, "expected a GWCS step")
    frame = _parse_frame(_required(node, "frame", path), "$path/frame")
    transform_node = get(node, "transform", nothing)
    return transform_node === nothing ? Step(frame) : Step(frame, parse_transform(transform_node; path="$path/transform"))
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

function parse_wcs(node::AbstractDict)
    _tag_name(node) == "wcs" || throw(ArgumentError("expected a GWCS WCS metadata node"))
    raw_steps = _required(node, "steps", "wcs")
    steps = map(enumerate(raw_steps)) do (index, value)
        _parse_step(value, "wcs/steps/$index")
    end
    first_transform = get(first(raw_steps), "transform", nothing)
    bounding_box = first_transform === nothing ? nothing : get(first_transform, "bounding_box", nothing)
    parsed_box = bounding_box === nothing ? nothing : _parse_bounding_box(bounding_box, "wcs/steps/1/transform/bounding_box")
    name = String(_unwrap(get(node, "name", "")))
    pixel_shape = get(node, "pixel_shape", nothing)
    return WCS(steps...; name, pixel_shape, bounding_box=parsed_box)
end

parse_wcs(metadata::AbstractString) = parse_wcs(_load_tagged_yaml(metadata))
