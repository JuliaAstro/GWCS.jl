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

# Writes back the properties a node was read with that GWCS.jl does not model. A property the writer produced wins, and a
# node now written with another tag family, such as an Euler rotate3d written as a composition, gets none.
function _restore_pass_through!(node, pass_through, path, registry)
    for (relative, (family, properties)) in pass_through
        target = node
        for key in split(relative, '/'; keepempty=false)
            target = target isa AbstractDict ? get(target, key, nothing) : target isa AbstractVector ? get(target, parse(Int, key), nothing) : nothing
        end
        target isa TaggedMapping && _pass_through_family(target.tag) == family || continue
        node_path = isempty(relative) ? path : "$path/$relative"
        for (key, value) in properties
            haskey(target, key) || (target.value[key] = _to_metadata(value, "$node_path/$key", registry))
        end
    end
    return node
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
_to_metadata(value::Union{Bool,Number}, path, registry) = value
# Stored metadata, such as an astropy frame kept from a parsed file, keeps its tags: ASDF.jl's tagged nodes are also
# strings, mappings, and vectors.
_to_metadata(value::AbstractString, path, registry) = _has_tag(value) ? TaggedScalar(_tag(value), String(value)) : value
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
    properties = Dict{Any,Any}(_to_metadata(key, "$path/key", registry) => _to_metadata(item, "$path/$key", registry) for (key, item) in value)
    return _has_tag(value) ? TaggedMapping(_tag(value), properties) : properties
end

function _to_metadata(value::Tuple, path, registry)
    return [_to_metadata(item, "$path/$index", registry) for (index, item) in enumerate(value)]
end

function _to_metadata(value::AbstractArray, path, registry)
    items = map(CartesianIndices(value)) do index
        indices = join(Tuple(index), "/")
        _to_metadata(value[index], "$path/$indices", registry)
    end
    return _has_tag(value) ? TaggedSequence(_tag(value), items) : items
end

# Writes a parameter read as a quantity back in its unit; the serializer wrote it as a plain number in degrees.
function _quantity_metadata(node, unit::_ParameterUnit, property, path)
    haskey(node, property) || _metadata_error("$path/$property", "the serializer omitted parameter $property, which was read as a quantity in $(repr(unit.unit))")
    _, from_degrees = _quantity_conversion(unit.kind, unit.unit, "$path/$property")
    value = metadata_value(node[property])
    quantity = Dict{Any,Any}("value" => value isa Real ? from_degrees(value) : from_degrees.(collect(Float64, value)), "unit" => TaggedScalar(unit.unit_tag, unit.unit))
    unit.datatype === nothing || (quantity["datatype"] = unit.datatype)
    return TaggedMapping(unit.tag, quantity)
end

function _to_metadata(transform::AnnotatedTransform, path, registry)
    node = _metadata_mapping(_to_metadata(transform.transform, path, registry), path)
    properties = Pair{String,Any}[]
    transform.name === nothing || push!(properties, "name" => transform.name)
    transform.inputs === nothing || push!(properties, "inputs" => _metadata_vector(transform.inputs))
    transform.outputs === nothing || push!(properties, "outputs" => _metadata_vector(transform.outputs))
    # A transform written with another tag, such as an Euler rotate3d written as rotate_sequence_3d, has other
    # parameters; its plain values are in degrees, as the schemas define them.
    family = first(_split_metadata_tag(node.tag, path))
    for unit in transform.parameter_units
        unit.transform_tag == family || continue
        # A cylinder radius read from the schemas' `lambda` is written as `lam`, the key Python GWCS reads.
        property = unit.property == "lambda" && haskey(node, "lam") ? "lam" : unit.property
        push!(properties, property => _quantity_metadata(node, unit, property, path))
    end
    return _restore_pass_through!(_decorate_metadata(node, properties), transform.pass_through, path, registry)
end

function _to_metadata(transform::WithInverse, path, registry)
    node = _metadata_mapping(_to_metadata(transform.forward, path, registry), path)
    inverse = _to_metadata(transform.backward, "$path/inverse", registry)
    return _decorate_metadata(node, ("inverse" => inverse,))
end

# FITS distortions have no ASDF tags (asdf-astropy leaves SIP out as FITS-specific), so they are written as the standard
# polynomial transforms they equal, which Python GWCS evaluates. A reverse polynomial, such as SIP's AP and BP, only
# starts the numerical inverse and is not an inverse, so it is not written.
_to_metadata(transform::Union{SIPTransform,TPVTransform,_FITSTPDTransform,_FITSPriorDistortion}, path, registry) = _to_metadata(_polynomial_form(transform, path), path, registry)
_to_metadata(transform::_ApproximateInverseTransform, path, registry) = _to_metadata(transform.forward, path, registry)
_to_metadata(::FITSLookupDistortion, path, registry) = _metadata_error(path, "FITS lookup-table distortions (D2IM, CPDIS) have no ASDF form")

# Coefficients with one more in the linear term of `axis`, which adds that input to the polynomial.
function _with_linear_term(coefficients, axis, weight=1.0)
    padded = zeros(max(size(coefficients, 1), 2), max(size(coefficients, 2), 2))
    padded[axes(coefficients)...] .= coefficients
    padded[axis == 1 ? 2 : 1, axis == 1 ? 1 : 2] += weight
    return padded
end

# Both outputs of a pair of polynomials in the two inputs.
_polynomial_pair(first, second) = ComposedTransform((RemapAxes{2}((1, 2, 1, 2)), parallel(Polynomial2D(first), Polynomial2D(second))))

function _polynomial_form(transform::SIPTransform, path)
    (transform.a, transform.b) isa NTuple{2,Polynomial2D{<:Any,<:Any,Nothing,Nothing}} || _metadata_error(path, "SIP polynomials with a domain or window have no polynomial form")
    return _polynomial_pair(_with_linear_term(transform.a.coefficients, 1), _with_linear_term(transform.b.coefficients, 2))
end

function _polynomial_form(transform::TPVTransform, path)
    any(!iszero, (transform.x[2]..., transform.y[2]...)) && _metadata_error(path, "TPV radial terms have no polynomial form")
    # The y polynomial takes (y, x).
    return _polynomial_pair(transform.x[1], permutedims(transform.y[1]))
end

function _polynomial_form(transform::_FITSPriorDistortion, path)
    any(table -> table !== nothing, (transform.detector_to_image..., transform.cpdis...)) && _metadata_error(path, "FITS lookup-table distortions (D2IM, CPDIS) have no ASDF form")
    # Without lookup tables, the axes a TPD model leaves out keep the SIP output, so this is SIP followed by TPD.
    return ComposedTransform((transform.sip, transform.tpd))
end

# Each output axis takes the inputs its model reads, or its own input where it has no model.
function _polynomial_form(transform::_FITSTPDTransform{N}, path) where {N}
    selected = Int[]
    outputs = AbstractTransform[]
    for (axis, model) in enumerate(transform.models)
        append!(selected, model === nothing ? (axis,) : model.axes)
        push!(outputs, model === nothing ? IdentityTransform{1}() : _tpd_axis_form(model, axis, "$path/$axis"))
    end
    return ComposedTransform((RemapAxes{N}(Tuple(selected)), parallel(outputs...)))
end

# A TPD axis evaluates its polynomial at (u, v), an affine map of the inputs it reads (an offset, a scale, and an
# auxiliary affine), and with DOCORR adds its own input, which inverting that map writes as linear terms in (u, v).
function _tpd_axis_form(model::_FITSTPDAxis{M}, axis, path) where {M}
    coefficients, radial = model.polynomial
    any(!iszero, radial) && _metadata_error(path, "FITS TPD radial terms have no polynomial form")
    position = findfirst(==(axis), model.axes)
    model.docorr && position === nothing && _metadata_error(path, "a FITS TPD correction to an axis its polynomial does not read has no polynomial form")
    if M == 1
        model.auxiliary === nothing || _metadata_error(path, "a one-variable FITS TPD with an auxiliary affine has no polynomial form")
        offset, scale = Float64(model.offset[1]), Float64(model.scale[1])
        values = Vector{Float64}(coefficients[:, 1])
        if model.docorr
            values[1] += offset
            length(values) < 2 && push!(values, 0.0)
            values[2] += 1 / scale
        end
        return ComposedTransform((Shift(-offset), Scale(scale), Polynomial1D(values)))
    end
    matrix = SMatrix{2,2,Float64}(Diagonal(SVector{2,Float64}(model.scale)))
    translation = -(matrix * SVector{2,Float64}(model.offset))
    if model.auxiliary !== nothing
        auxiliary = SMatrix{2,3,Float64}(model.auxiliary)
        linear = auxiliary[:, SVector(2, 3)]
        matrix, translation = linear * matrix, auxiliary[:, 1] + linear * translation
    end
    polynomial = Matrix{Float64}(coefficients)
    if model.docorr
        iszero(det(matrix)) && _metadata_error(path, "a FITS TPD correction with a singular variable map has no polynomial form")
        # The input is an affine function of (u, v): row `position` of inv(matrix) * ((u, v) - translation).
        inverse = inv(matrix)
        polynomial = _with_linear_term(_with_linear_term(polynomial, 1, inverse[position, 1]), 2, inverse[position, 2])
        polynomial[1, 1] -= (inverse * translation)[position]
    end
    return ComposedTransform((AffineTransform(matrix, translation), Polynomial2D(polynomial)))
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
_serialize_rotate2d(transform, context) = Dict("angle" => transform.angle)

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

const _ParameterlessProjection = Union{GnomonicProjection,StereographicProjection,ZenithalEquidistant,ZenithalEqualArea,PlateCarreeProjection,MercatorProjection,
    SansonFlamsteedProjection,ParabolicProjection,MollweideProjection,HammerAitoffProjection,PolyconicProjection,TangentialSphericalCube,COBESphericalCube,
    QuadrilateralizedSphericalCube,HEALPixPolarProjection}

_projection_properties(::_ParameterlessProjection) = Dict{String,Any}()
_projection_properties(projection::ZenithalPerspective) = Dict{String,Any}("mu" => projection.mu, "gamma" => projection.gamma)
_projection_properties(projection::SlantZenithalPerspective) = Dict{String,Any}("mu" => projection.mu, "phi0" => projection.phi0, "theta0" => projection.theta0)
_projection_properties(projection::SlantOrthographic) = Dict{String,Any}("xi" => projection.xi, "eta" => projection.eta)
_projection_properties(projection::AiryProjection) = Dict{String,Any}("theta_b" => projection.theta_b)
# Written as astropy's `lam`, which Python GWCS reads; it ignores the schemas' `lambda`.
_projection_properties(projection::CylindricalEqualArea) = Dict{String,Any}("lam" => projection.lambda)
_projection_properties(projection::CylindricalPerspective) = Dict{String,Any}("mu" => projection.mu, "lam" => projection.lambda)
_projection_properties(projection::AbstractConicProjection) = Dict{String,Any}("sigma" => projection.sigma, "delta" => projection.delta)
_projection_properties(projection::BonneProjection) = Dict{String,Any}("theta1" => projection.theta1)
_projection_properties(projection::HEALPixProjection) = Dict{String,Any}("H" => Float64(projection.H), "X" => Float64(projection.K))

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
    register_transform_serializer!(registry, Rotation2D, "$_ASDF_TRANSFORM_TAG/rotate2d", _serialize_rotate2d; version=v"1.5.0")
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
        (AiryProjection, "airy", v"1.4.0"),
        (PlateCarreeProjection, "plate_carree", v"1.4.0"),
        (CylindricalEqualArea, "cylindrical_equal_area", v"1.5.0"),
        (CylindricalPerspective, "cylindrical_perspective", v"1.5.0"),
        (MercatorProjection, "mercator", v"1.4.0"),
        (SansonFlamsteedProjection, "sanson_flamsteed", v"1.4.0"),
        (ParabolicProjection, "parabolic", v"1.4.0"),
        (MollweideProjection, "molleweide", v"1.4.0"),
        (HammerAitoffProjection, "hammer_aitoff", v"1.4.0"),
        (PolyconicProjection, "polyconic", v"1.4.0"),
        (ConicPerspective, "conic_perspective", v"1.5.0"),
        (ConicEquidistant, "conic_equidistant", v"1.5.0"),
        (ConicEqualArea, "conic_equal_area", v"1.5.0"),
        (ConicOrthomorphic, "conic_orthomorphic", v"1.5.0"),
        (BonneProjection, "bonne_equal_area", v"1.5.0"),
        (TangentialSphericalCube, "tangential_spherical_cube", v"1.4.0"),
        (COBESphericalCube, "cobe_quad_spherical_cube", v"1.4.0"),
        (QuadrilateralizedSphericalCube, "quad_spherical_cube", v"1.4.0"),
        (HEALPixProjection, "healpix", v"1.4.0"),
        (HEALPixPolarProjection, "healpix_polar", v"1.4.0")
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

# The frame schemas require a unit tag for every axis; an axis without a unit is dimensionless, as Python GWCS writes it.
_unit_metadata(value) = TaggedScalar("$_UNIT_TAG-1.0.0", String(value))

function _reference_metadata(value, path, registry)
    value === nothing && return nothing
    # Metadata kept from a FITS header that the ASDF frame schemas cannot hold.
    value isa NamedTuple && _metadata_error(path, "FITS reference metadata $value has no ASDF representation")
    value === :ICRS && return _tagged_metadata(_ICRS_TAG, v"1.1.0", Dict("frame_attributes" => Dict()))
    return _to_metadata(value, path, registry)
end

function _frame_properties(frame, path, registry)
    # Per-axis lists are in the frame's own order, which Python GWCS reads, and writes since gwcs#687.
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
    # A two-axis spatial frame without a reference frame is what Python GWCS's Frame2D holds, and it reads frame2d as
    # one; its reader rejects naxes there, and ignores axes_type.
    if frame_dimension(frame) == 2 && all(==(SPATIAL), frame.axis_types) && frame.reference_frame === nothing
        return _tagged_metadata("$_GWCS_TAG/frame2d", v"1.2.0", properties)
    end
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
    frame.reference_frame isa NamedTuple && return _reference_metadata(frame.reference_frame, path, registry)
    frame.reference_frame === :ICRS && return _reference_metadata(:ICRS, path, registry)
    frame.reference_frame isa Union{TaggedMapping,TaggedSequence,TaggedScalar} && return _to_metadata(frame.reference_frame, path, registry)
    frame.reference_epoch === nothing && return _reference_metadata(frame.reference_frame, path, registry)
    properties = Dict{String,Any}("value" => _to_metadata(frame.reference_epoch, "$path/value", registry))
    isempty(frame.time_scale) || (properties["scale"] = frame.time_scale)
    isempty(frame.time_format) || (properties["format"] = frame.time_format)
    frame.location === nothing || (properties["location"] = _location_metadata(frame.location, "$path/location", registry))
    return _tagged_metadata(_TIME_TAG, v"1.4.0", properties)
end

# A FITS OBSGEO location is geocentric (x, y, z) in meters.
_location_metadata(location::NTuple{3,Real}, path, registry) = Dict(name => TaggedMapping("$_QUANTITY_TAG-1.3.0", Dict{Any,Any}("value" => value, "unit" => TaggedScalar("$_UNIT_TAG-1.0.0", "m"))) for (name, value) in zip(("x", "y", "z"), location))
_location_metadata(location, path, registry) = _to_metadata(location, path, registry)

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

function _serialize_bounding_box(box::BoundingBox, names, path, registry; include_ignored=true)
    intervals = Dict{Any,Any}()
    for (index, axis) in enumerate(box.axes)
        intervals[names[axis]] = _metadata_vector(box.intervals[index])
    end
    properties = Dict{String,Any}("intervals" => intervals, "order" => "F")
    if include_ignored
        ignored = [names[axis] for axis in ignored_axes(box)]
        isempty(ignored) || (properties["ignore"] = ignored)
    end
    return _tagged_metadata("$_ASDF_TRANSFORM_TAG/property/bounding_box", v"1.2.0", properties)
end

function _serialize_bounding_box(box::CompoundBoundingBox, names, path, registry)
    selector_args = [Dict("argument" => names[box.selector_axes[index]], "ignore" => box.selector_ignored[index]) for index in eachindex(box.selector_axes)]
    entries = Any[]
    for (index, entry) in enumerate(box.boxes)
        nested = _serialize_bounding_box(entry.second, names, "$path/cbbox/$index/bbox", registry; include_ignored=false)
        push!(entries, Dict("key" => _metadata_vector(entry.first), "bbox" => nested))
    end
    properties = Dict{String,Any}("selector_args" => selector_args, "cbbox" => entries, "order" => "F")
    ignored = [names[axis] for axis in box.globally_ignored]
    isempty(ignored) || (properties["ignore"] = ignored)
    return _tagged_metadata("$_ASDF_TRANSFORM_TAG/property/compound_bounding_box", v"1.2.0", properties)
end

function _to_metadata(box::AbstractBoundingBox, path, registry)
    names = ntuple(index -> "x$(index - 1)", input_dimension(box))
    return _serialize_bounding_box(box, names, path, registry)
end

# The first step's transform carries the WCS's bounding box, which replaces any it was read with, even when it is
# `nothing`.
function _serialize_step(step::Step, path, registry; bounding_box=_UNSPECIFIED_PIPELINE_METADATA)
    # The step schema requires an astropy reference frame, which a temporal frame's time is not; in a composite it is valid.
    frame = step.frame isa TemporalFrame ? CompositeFrame((step.frame,); name=step.frame.name) : step.frame
    properties = Dict{String,Any}("frame" => _to_metadata(frame, "$path/frame", registry))
    if step.transform !== nothing
        transform = _metadata_mapping(_to_metadata(step.transform, "$path/transform", registry), "$path/transform")
        if bounding_box !== _UNSPECIFIED_PIPELINE_METADATA
            names = _transform_axis_names(step.transform)
            box = bounding_box === nothing ? nothing : _serialize_bounding_box(bounding_box, names, "$path/transform/bounding_box", registry)
            transform = _decorate_metadata(transform, ("bounding_box" => box,))
        end
        properties["transform"] = transform
    end
    return _restore_pass_through!(_tagged_metadata("$_GWCS_TAG/step", v"1.3.0", properties), step.pass_through, path, registry)
end

function _to_metadata(step::Step, path, registry)
    return _serialize_step(step, path, registry)
end

function _to_metadata(wcs::WCS, path, registry)
    steps = Any[]
    for (index, step) in enumerate(wcs.steps)
        push!(steps, index == 1 ? _serialize_step(step, "$path/steps/$index", registry; wcs.bounding_box) : _serialize_step(step, "$path/steps/$index", registry))
    end
    properties = Dict{String,Any}("name" => wcs.name, "steps" => steps)
    wcs.pixel_shape === nothing || any(isnothing, wcs.pixel_shape) || (properties["pixel_shape"] = collect(wcs.pixel_shape))
    return _restore_pass_through!(_tagged_metadata("$_GWCS_TAG/wcs", v"1.4.0", properties), wcs.pass_through, path, registry)
end

function _to_metadata(value, path, registry)
    _metadata_error(path, "value of type $(typeof(value)) cannot be serialized canonically")
end

function write_asdf(filename, value; kwargs...)
    throw(ArgumentError("load ASDF.jl before calling write_asdf"))
end

_register_builtin_transform_serializers!()
register_transform_serializer!(default_transform_registry(), RegionSelector, "$_GWCS_TAG/regions_selector", _serialize_region_selector; version=v"1.3.0")
