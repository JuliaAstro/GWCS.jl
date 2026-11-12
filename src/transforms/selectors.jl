function _region_entries(selector, ::Val{Nin}, ::Val{Nout}) where {Nin,Nout}
    entries = _selector_entries(selector)
    isempty(entries) && throw(ArgumentError("a region selector requires at least one transform"))
    labels = map(entry -> entry.first, entries)
    label_types = unique(typeof.(labels))
    length(label_types) == 1 || throw(ArgumentError("region labels must have one concrete type"))
    L = only(label_types)
    isconcretetype(L) || throw(ArgumentError("region label type must be concrete"))
    transforms = Tuple(entry.second for entry in entries)
    all(transform -> transform isa AbstractTransform{Nin,Nout}, transforms) || throw(DimensionMismatch("all selected transforms must have dimensions $Nin -> $Nout"))
    lookup = Dict{L,Int}()
    for (index, label) in enumerate(labels)
        haskey(lookup, label) && throw(ArgumentError("region labels must be unique"))
        lookup[label] = index
    end
    return L, collect(L, labels), transforms, lookup
end

"""
    RegionSelector(label_mapper, selector; undefined_value=NaN)

A transform that applies a different branch in each region, such as each slice of an integral field unit.
`label_mapper` assigns a region label to the input coordinate, and `selector` holds `label => transform` pairs with
labels of the mapper's label type. The selected branch evaluates the original coordinate. A coordinate in no region
gives `NaN`, which [`pixel_to_world`](@ref) reports as having no result.

It has an inverse when the label mapper has one, as a [`WithInverseLabelMapper`](@ref), and every branch has one.
`undefined_value` is kept only so ASDF metadata round-trips exactly.

Python equivalent: `gwcs.selector.RegionsSelector`, whose `undefined_transform_value` is the value it returns outside
every region.
"""
struct RegionSelector{Nin,Nout,L,M<:AbstractLabelMapper{Nin,L},Ts<:Tuple,D<:AbstractDict{L,Int},U<:Real,I<:Union{Nothing,Tuple}} <: AbstractTransform{Nin,Nout}
    label_mapper::M
    labels::Vector{L}
    transforms::Ts
    # Branch inverses, kept only when the label mapper has no inverse: the selector then has none either, and they
    # seed the numerical inverse. An invertible mapper gives an analytical inverse selector instead.
    inverses::I
    lookup::D
    # Python GWCS's `undefined_transform_value`, kept only so ASDF round trips are exact. Evaluation always
    # returns NaN for a coordinate in no region, so the API treats it as having no result (`fill_value` or `strict`).
    undefined_value::U
end

function RegionSelector(label_mapper::AbstractLabelMapper{Nin,L}, selector; undefined_value=NaN) where {Nin,L}
    undefined_value isa Real || throw(ArgumentError("undefined_value must be a real number"))
    entries = _selector_entries(selector)
    isempty(entries) && throw(ArgumentError("a region selector requires at least one transform"))
    Nout = output_dimension(first(entries).second)
    label_type_value, labels, transforms, lookup = _region_entries(entries, Val(Nin), Val(Nout))
    label_type_value == L || throw(ArgumentError("region labels have type $label_type_value, expected $L from the label mapper"))
    inverses = _inverse_label_mapper_checked(label_mapper) isa NoInverse ? map(_inverse_transform_checked, transforms) : nothing
    return RegionSelector{Nin,Nout,L,typeof(label_mapper),typeof(transforms),typeof(lookup),typeof(undefined_value),typeof(inverses)}(label_mapper, labels, transforms, inverses, lookup, undefined_value)
end

@inline function _undefined_region_output(::RegionSelector{Nin,Nout}, coords::StaticVector{Nin}) where {Nin,Nout}
    T = _evaluation_type(eltype(coords))
    return SVector{Nout,T}(ntuple(_ -> T(NaN), Nout))
end

@inline function (transform::RegionSelector{Nin,Nout})(coords::StaticVector{Nin}) where {Nin,Nout}
    label = map_label(transform.label_mapper, coords)
    label === nothing && return _undefined_region_output(transform, coords)
    index = get(transform.lookup, label, 0)
    iszero(index) && return _undefined_region_output(transform, coords)
    return _apply_transform_at(transform.transforms, index, coords)
end

"""
    selected_transform(selector::RegionSelector, label) -> AbstractTransform

Return the transform `selector` applies in the region with `label`. A label with no transform throws a `KeyError`.

Python equivalent: `gwcs.selector.RegionsSelector.selector[label]`.
"""
function selected_transform(transform::RegionSelector, label)
    index = get(transform.lookup, label, 0)
    iszero(index) && throw(KeyError(label))
    return transform.transforms[index]
end

function compile_transform(transform::RegionSelector)
    mapper = _compile_label_mapper_checked(transform.label_mapper)
    transforms = map(_compile_transform_checked, transform.transforms)
    selector = ntuple(index -> transform.labels[index] => transforms[index], length(transforms))
    return RegionSelector(mapper, selector; undefined_value=transform.undefined_value)
end

function inverse_transform(transform::RegionSelector)
    mapper = _inverse_label_mapper_checked(transform.label_mapper)
    mapper isa NoInverse && return mapper
    inverses = map(_inverse_transform_checked, transform.transforms)
    any(inverse -> inverse isa NoInverse, inverses) && return NoInverse()
    input_dimension(mapper) == output_dimension(transform) || throw(DimensionMismatch("inverse label mapper has input dimension $(input_dimension(mapper)), expected $(output_dimension(transform))"))
    selector = ntuple(index -> transform.labels[index] => inverses[index], length(inverses))
    # Like Python GWCS, the inverse selector does not inherit a custom undefined value.
    return RegionSelector(mapper, selector)
end

# Without an inverse label mapper, start Newton from the first branch inverse that lands back in its own region;
# Newton then verifies the candidate against the forward selector.
@inline _region_initialize(::Tuple{}, labels, mapper, coords, index) = (coords, false)

@inline function _region_initialize(inverses::Tuple, labels, mapper, coords, index)
    inverse = first(inverses)
    if !(inverse isa NoInverse)
        candidate = inverse(coords)
        all(isfinite, candidate) && map_label(mapper, candidate) == labels[index] && return (candidate, true)
    end
    return _region_initialize(Base.tail(inverses), labels, mapper, coords, index + 1)
end

@inline function _inverse_initialize(transform::RegionSelector{N,N,L,M,Ts,D,U,<:Tuple}, coords::StaticVector{N}) where {N,L,M,Ts,D,U}
    return _region_initialize(transform.inverses, transform.labels, transform.label_mapper, coords, 1)
end

function transform_jacobian(transform::RegionSelector{Nin,Nout}, coords::StaticVector{Nin}) where {Nin,Nout}
    # Outside every region the selector has no derivative.
    undefined = SMatrix{Nout,Nin}(ntuple(_ -> _evaluation_type(eltype(coords))(NaN), Nin * Nout))
    label = map_label(transform.label_mapper, coords)
    label === nothing && return undefined
    index = get(transform.lookup, label, 0)
    iszero(index) && return undefined
    return _transform_jacobian_at(transform.transforms, index, coords)
end

function transform_dependencies(transform::RegionSelector{Nin,Nout}) where {Nin,Nout}
    dependencies = MMatrix{Nout,Nin,Bool}(false for _ in 1:Nout, _ in 1:Nin)
    for branch in transform.transforms
        dependencies .|= _transform_dependencies_checked(branch)
    end
    selection = _label_dependencies_checked(transform.label_mapper)
    for row in 1:Nout, column in 1:Nin
        dependencies[row, column] |= selection[column]
    end
    return SMatrix(dependencies)
end
