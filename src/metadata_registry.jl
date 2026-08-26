"""A metadata constructor registered for one tag family and version interval."""
struct TransformRegistration{F}
    tag::String
    min_version::VersionNumber
    max_version::VersionNumber
    constructor::F
end

"""
    TransformRegistry()

Mutable construction-time registry for metadata transform constructors.
"""
mutable struct TransformRegistry
    registrations::Dict{String,Vector{Any}}
    lock::ReentrantLock
end

TransformRegistry() = TransformRegistry(Dict{String,Vector{Any}}(), ReentrantLock())

function Base.copy(registry::TransformRegistry)
    return lock(registry.lock) do
        registrations = Dict(tag => copy(values) for (tag, values) in registry.registrations)
        TransformRegistry(registrations, ReentrantLock())
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

function _version_range(registration::TransformRegistration)
    return "[$(registration.min_version), $(registration.max_version))"
end

"""
    register_transform!([registry], tag, constructor; min_version, max_version)

Register a native transform constructor for the half-open semantic version range
`[min_version, max_version)`. `tag` is the full metadata tag family without its
trailing version.
"""
function register_transform!(registry::TransformRegistry, tag::AbstractString, constructor; min_version=v"1.0.0", max_version=v"2.0.0")
    min_version < max_version || throw(ArgumentError("transform registration minimum version must be less than its maximum version"))
    family = String(tag)
    isempty(family) && throw(ArgumentError("transform registration tag must not be empty"))
    occursin(r"-\d+\.\d+\.\d+$", family) && throw(ArgumentError("transform registration tag must omit the schema version: $family"))
    registration = TransformRegistration(family, min_version, max_version, constructor)
    lock(registry.lock) do
        registrations = get!(registry.registrations, family, Any[])
        for existing in registrations
            overlaps = min_version < existing.max_version && existing.min_version < max_version
            overlaps || continue
            throw(ArgumentError("transform registration for $family $(_version_range(registration)) overlaps existing range $(_version_range(existing))"))
        end
        push!(registrations, registration)
        sort!(registrations; by=value -> value.min_version)
    end
    return registry
end

function register_transform!(tag::AbstractString, constructor; registry=default_transform_registry(), min_version=v"1.0.0", max_version=v"2.0.0")
    return register_transform!(registry, tag, constructor; min_version, max_version)
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
    index = findfirst(registration -> registration.min_version <= version < registration.max_version, registrations)
    if index === nothing
        ranges = join(map(_version_range, registrations), ", ")
        _metadata_error(path, "unsupported transform tag version $tag-$version; registered ranges are $ranges")
    end
    return registrations[index]
end
