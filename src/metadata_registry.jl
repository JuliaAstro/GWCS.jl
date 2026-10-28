"""A metadata constructor registered for one tag family and supported versions."""
struct TransformRegistration{F,V,N,S}
    tag::String
    min_version::VersionNumber
    max_version::VersionNumber
    constructor::F
    validator::V
    normalizer::N
    versions::S
    # Properties that may hold quantities, as `property => :angle` or `:dimensionless`.
    quantities::Vector{Pair{String,Symbol}}
    # The properties the constructor and serializer model, or `nothing` when none are declared.
    properties::Union{Nothing,Vector{String}}
end

"""Validation and normalization rules for a non-transform metadata tag."""
struct MetadataSchemaRegistration{V,N,S}
    tag::String
    versions::S
    validator::V
    normalizer::N
end

"""Canonical metadata serializer registered for one transform type."""
struct TransformSerializerRegistration{F}
    transform_type
    tag::String
    version::VersionNumber
    serializer::F
end

"""
    TransformRegistry()

Mutable construction-time registry for metadata transform constructors and
serializers.
"""
mutable struct TransformRegistry
    registrations::Dict{String,Vector{Any}}
    schemas::Dict{String,MetadataSchemaRegistration}
    serializers::Vector{Any}
    lock::ReentrantLock
end

function TransformRegistry()
    schemas = isdefined(@__MODULE__, :_DEFAULT_TRANSFORM_REGISTRY) ? copy(_DEFAULT_TRANSFORM_REGISTRY.schemas) : Dict{String,MetadataSchemaRegistration}()
    return TransformRegistry(Dict{String,Vector{Any}}(), schemas, Any[], ReentrantLock())
end

function Base.copy(registry::TransformRegistry)
    return lock(registry.lock) do
        registrations = Dict(tag => copy(values) for (tag, values) in registry.registrations)
        TransformRegistry(registrations, copy(registry.schemas), copy(registry.serializers), ReentrantLock())
    end
end

const _DEFAULT_TRANSFORM_REGISTRY = TransformRegistry()

"""Return the process-wide transform registry used by metadata parsers."""
default_transform_registry() = _DEFAULT_TRANSFORM_REGISTRY

"""Metadata tag, version, path, and registry supplied to a constructor."""
struct TransformParseContext
    tag::String
    version::VersionNumber
    path::String
    registry::TransformRegistry
end

metadata_tag(context::TransformParseContext) = context.tag
metadata_version(context::TransformParseContext) = context.version
metadata_path(context::TransformParseContext) = context.path

"""Metadata tag, version, path, and registry supplied to a serializer."""
struct TransformWriteContext
    tag::String
    version::VersionNumber
    path::String
    registry::TransformRegistry
end

metadata_tag(context::TransformWriteContext) = context.tag
metadata_version(context::TransformWriteContext) = context.version
metadata_path(context::TransformWriteContext) = context.path

function _version_range(registration::TransformRegistration)
    registration.versions === nothing && return "[$(registration.min_version), $(registration.max_version))"
    return "{" * join(registration.versions, ", ") * "}"
end

_no_metadata_validation(node, context) = nothing
_no_metadata_normalization(node, context) = node

function _metadata_versions(versions)
    values = map(versions) do version
        version isa VersionNumber ? version : VersionNumber(version)
    end
    isempty(values) && throw(ArgumentError("metadata versions must not be empty"))
    return Tuple(sort!(unique!(collect(values))))
end

function _registration_accepts(registration::TransformRegistration, version)
    registration.versions === nothing && return registration.min_version <= version < registration.max_version
    return version in registration.versions
end

function _registrations_overlap(left::TransformRegistration, right::TransformRegistration)
    if left.versions === nothing && right.versions === nothing
        return left.min_version < right.max_version && right.min_version < left.max_version
    elseif left.versions === nothing
        return any(version -> left.min_version <= version < left.max_version, right.versions)
    elseif right.versions === nothing
        return any(version -> right.min_version <= version < right.max_version, left.versions)
    end
    return !isempty(intersect(left.versions, right.versions))
end

"""
    register_transform!([registry], tag, constructor; min_version, max_version, versions, validate, normalize, angles, dimensionless, properties)

Register a native transform constructor for either an explicit collection of
schema `versions` or the half-open semantic version range
`[min_version, max_version)`. Existing registrations default to `[1.0.0, 2.0.0)`.
`tag` is the full metadata tag family without its trailing version.

`angles` and `dimensionless` name properties that may hold an ASDF quantity, which is replaced by a plain number
before validation, so the constructor sees what a unitless file holds. An angle quantity in any angle unit is converted
to degrees; a dimensionless one must have an empty unit. The parsed transform records each unit, and writing restores
it.

`properties` names the properties the constructor and the type's serializer handle, besides the `name`, `inputs`,
`outputs`, and `inverse` every transform has. The others in a parsed node are kept on the transform and written back
with it, when it is written with the same tag family. A registration that declares no `properties` keeps none.
"""
function register_transform!(registry::TransformRegistry, tag::AbstractString, constructor; min_version=nothing, max_version=nothing, versions=nothing, validate=_no_metadata_validation, normalize=_no_metadata_normalization, angles=(), dimensionless=(), properties=nothing)
    if versions === nothing
        minimum = min_version === nothing ? v"1.0.0" : min_version
        maximum = max_version === nothing ? v"2.0.0" : max_version
        minimum isa VersionNumber || throw(ArgumentError("transform registration minimum version must be a VersionNumber"))
        maximum isa VersionNumber || throw(ArgumentError("transform registration maximum version must be a VersionNumber"))
        minimum < maximum || throw(ArgumentError("transform registration minimum version must be less than its maximum version"))
        supported = nothing
    else
        min_version === nothing && max_version === nothing || throw(ArgumentError("transform registration versions cannot be combined with min_version or max_version"))
        supported = _metadata_versions(versions)
        minimum = first(supported)
        maximum = last(supported)
    end
    family = String(tag)
    isempty(family) && throw(ArgumentError("transform registration tag must not be empty"))
    occursin(r"-\d+\.\d+\.\d+$", family) && throw(ArgumentError("transform registration tag must omit the schema version: $family"))
    quantities = Pair{String,Symbol}[(String(name) => :angle for name in angles)..., (String(name) => :dimensionless for name in dimensionless)...]
    allunique(first.(quantities)) || throw(ArgumentError("transform registration for $family lists a quantity property twice"))
    modeled = properties === nothing ? nothing : String[String(name) for name in properties]
    registration = TransformRegistration(family, minimum, maximum, constructor, validate, normalize, supported, quantities, modeled)
    lock(registry.lock) do
        registrations = get!(registry.registrations, family, Any[])
        for existing in registrations
            _registrations_overlap(registration, existing) || continue
            throw(ArgumentError("transform registration for $family $(_version_range(registration)) overlaps existing range $(_version_range(existing))"))
        end
        push!(registrations, registration)
        sort!(registrations; by=value -> value.min_version)
    end
    return registry
end

function register_transform!(tag::AbstractString, constructor; registry=default_transform_registry(), kwargs...)
    return register_transform!(registry, tag, constructor; kwargs...)
end

"""
    register_transform_serializer!([registry], type, tag, serializer; version)

Register a canonical metadata serializer for an `AbstractTransform` subtype.
The serializer receives the transform and a `TransformWriteContext` and must
return an untagged mapping. The registry attaches `tag-version`.
"""
function register_transform_serializer!(registry::TransformRegistry, transform_type, tag::AbstractString, serializer; version)
    transform_type <: AbstractTransform || throw(ArgumentError("transform serializer type must subtype AbstractTransform"))
    version isa VersionNumber || throw(ArgumentError("transform serializer version must be a VersionNumber"))
    version >= v"0.0.0" || throw(ArgumentError("transform serializer version must be non-negative"))
    family = String(tag)
    isempty(family) && throw(ArgumentError("transform serializer tag must not be empty"))
    occursin(r"-\d+\.\d+\.\d+$", family) && throw(ArgumentError("transform serializer tag must omit the schema version: $family"))
    registration = TransformSerializerRegistration(transform_type, family, version, serializer)
    lock(registry.lock) do
        parsers = get(registry.registrations, family, Any[])
        any(parser -> _registration_accepts(parser, version), parsers) || throw(ArgumentError("transform serializer tag $family-$version has no matching parser registration"))
        any(existing -> existing.transform_type == transform_type, registry.serializers) && throw(ArgumentError("transform serializer for $transform_type is already registered"))
        push!(registry.serializers, registration)
    end
    return registry
end

function register_transform_serializer!(transform_type, tag::AbstractString, serializer; registry=default_transform_registry(), version)
    return register_transform_serializer!(registry, transform_type, tag, serializer; version)
end

function _lookup_transform_serializer(registry::TransformRegistry, transform::AbstractTransform, path)
    transform_type = typeof(transform)
    candidates = lock(registry.lock) do
        filter(registration -> transform_type <: registration.transform_type, registry.serializers)
    end
    isempty(candidates) && _metadata_error(path, "no metadata serializer is registered for $transform_type")
    selected = first(candidates)
    for candidate in Iterators.drop(candidates, 1)
        if candidate.transform_type <: selected.transform_type
            selected = candidate
        elseif !(selected.transform_type <: candidate.transform_type)
            _metadata_error(path, "metadata serializers for $(selected.transform_type) and $(candidate.transform_type) are ambiguous for $transform_type")
        end
    end
    return selected
end

function _split_metadata_tag(tag::AbstractString, path)
    matched = match(r"^(.*)-(\d+\.\d+\.\d+)$", tag)
    matched === nothing && _metadata_error(path, "unsupported metadata tag $tag")
    return matched.captures[1], VersionNumber(matched.captures[2])
end

function _lookup_transform_registration(registry::TransformRegistry, tag, version, path)
    registrations = lock(registry.lock) do
        copy(get(registry.registrations, tag, Any[]))
    end
    isempty(registrations) && _metadata_error(path, "unsupported transform tag $tag-$version")
    index = findfirst(registration -> _registration_accepts(registration, version), registrations)
    if index === nothing
        ranges = join(map(_version_range, registrations), ", ")
        description = all(registration -> registration.versions !== nothing, registrations) ? "supported versions are" : "registered ranges are"
        _metadata_error(path, "unsupported transform tag version $tag-$version; $description $ranges")
    end
    return registrations[index]
end

function _prepare_transform_metadata(registration::TransformRegistration, node, context)
    registration.validator(node, context)
    return registration.normalizer(node, context)
end

function _register_metadata_schema!(registry::TransformRegistry, tag::AbstractString; versions, validate=_no_metadata_validation, normalize=_no_metadata_normalization)
    family = String(tag)
    isempty(family) && throw(ArgumentError("metadata schema tag must not be empty"))
    occursin(r"-\d+\.\d+\.\d+$", family) && throw(ArgumentError("metadata schema tag must omit the schema version: $family"))
    registration = MetadataSchemaRegistration(family, _metadata_versions(versions), validate, normalize)
    lock(registry.lock) do
        haskey(registry.schemas, family) && throw(ArgumentError("metadata schema for $family is already registered"))
        registry.schemas[family] = registration
    end
    return registry
end

function _lookup_metadata_schema(registry::TransformRegistry, tag, version, path)
    registration = lock(registry.lock) do
        get(registry.schemas, tag, nothing)
    end
    registration === nothing && _metadata_error(path, "unsupported metadata tag $tag-$version")
    version in registration.versions || _metadata_error(path, "unsupported metadata tag version $tag-$version; supported versions are {$(join(registration.versions, ", "))}")
    return registration
end

function _prepare_metadata(node, path, registry::TransformRegistry)
    tag, version = _split_metadata_tag(_tag(node), path)
    context = TransformParseContext(tag, version, String(path), registry)
    registration = _lookup_metadata_schema(registry, tag, version, path)
    registration.validator(node, context)
    return registration.normalizer(node, context), context
end
