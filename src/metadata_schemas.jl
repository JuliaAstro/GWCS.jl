const _ASDF_TRANSFORM_TAG = "tag:stsci.edu:asdf/transform"
const _GWCS_TAG = "tag:stsci.edu:gwcs"
const _SCHEMA_VERSIONS_1_0_TO_1_2 = (v"1.0.0", v"1.1.0", v"1.2.0")
const _SCHEMA_VERSIONS_1_0_TO_1_3 = (_SCHEMA_VERSIONS_1_0_TO_1_2..., v"1.3.0")
const _SCHEMA_VERSIONS_1_0_TO_1_4 = (_SCHEMA_VERSIONS_1_0_TO_1_3..., v"1.4.0")
const _SCHEMA_VERSIONS_1_0_TO_1_5 = (_SCHEMA_VERSIONS_1_0_TO_1_4..., v"1.5.0")
const _SCHEMA_VERSIONS_1_0_TO_1_6 = (_SCHEMA_VERSIONS_1_0_TO_1_5..., v"1.6.0")

function _schema_property(node, key, context; required=false)
    if !haskey(node, key)
        required && _metadata_error(metadata_path(context), "missing required property $key")
        return nothing
    end
    return metadata_value(node[key])
end

function _schema_mapping(value, path)
    value isa AbstractDict || _metadata_error(path, "expected a mapping")
    return value
end

function _schema_sequence(value, path)
    value isa AbstractVector || value isa Tuple || _metadata_error(path, "expected an array")
    return value
end

function _schema_array(value, path)
    value isa AbstractArray || value isa Tuple || _metadata_error(path, "expected an array")
    return value
end

function _schema_metadata_array(value, path)
    materialized, _ = _array_metadata(value, path)
    return _schema_array(materialized, path)
end

function _schema_string(value, path)
    value isa AbstractString || _metadata_error(path, "expected a string")
    return value
end

function _schema_integer(value, path)
    value isa Integer && !(value isa Bool) || _metadata_error(path, "expected an integer")
    return value
end

function _schema_number(value, path)
    value isa Real && !(value isa Bool) || _metadata_error(path, "expected a unitless numeric value")
    return value
end

function _schema_string_sequence(value, path)
    values = _schema_sequence(value, path)
    for (index, item) in enumerate(values)
        _schema_string(metadata_value(item), "$path/$index")
    end
    return values
end

function _schema_tag(value, family, path)
    value isa AbstractDict || _metadata_error(path, "expected tagged metadata")
    tag, version = _split_metadata_tag(_tag(value), path)
    tag == family || _metadata_error(path, "expected metadata tag $family, got $tag-$version")
    return version
end

function _validate_transform_annotations(node, context)
    path = metadata_path(context)
    haskey(node, "name") && _schema_string(metadata_value(node["name"]), "$path/name")
    haskey(node, "inputs") && _schema_string_sequence(node["inputs"], "$path/inputs")
    haskey(node, "outputs") && _schema_string_sequence(node["outputs"], "$path/outputs")
    haskey(node, "inverse") && _schema_mapping(node["inverse"], "$path/inverse")
    return nothing
end

function _validate_required_transform(node, context, properties)
    _schema_mapping(node, metadata_path(context))
    _validate_transform_annotations(node, context)
    for property in properties
        _schema_property(node, property, context; required=true)
    end
    return nothing
end

function _validate_compound_transform(node, context)
    _validate_required_transform(node, context, ("forward",))
    forward = _schema_sequence(node["forward"], "$(metadata_path(context))/forward")
    isempty(forward) && _metadata_error("$(metadata_path(context))/forward", "transform array must not be empty")
    for (index, transform) in enumerate(forward)
        _schema_mapping(transform, "$(metadata_path(context))/forward/$index")
    end
    return nothing
end

function _validate_binary_transform(node, context)
    _validate_compound_transform(node, context)
    length(node["forward"]) == 2 || _metadata_error("$(metadata_path(context))/forward", "arithmetic metadata must contain exactly two transforms")
    return nothing
end

function _validate_numeric_transform(node, context, property)
    _validate_required_transform(node, context, (property,))
    _schema_number(metadata_value(node[property]), "$(metadata_path(context))/$property")
    return nothing
end

function _validate_affine_schema(node, context)
    _validate_required_transform(node, context, ("matrix",))
    _schema_metadata_array(node["matrix"], "$(metadata_path(context))/matrix")
    haskey(node, "translation") && _schema_metadata_array(node["translation"], "$(metadata_path(context))/translation")
    return nothing
end

function _validate_gnomonic_schema(node, context)
    _validate_required_transform(node, context, ())
    haskey(node, "direction") && _schema_string(metadata_value(node["direction"]), "$(metadata_path(context))/direction")
    return nothing
end

function _validate_rotate3d_schema(node, context)
    _validate_required_transform(node, context, ("phi", "theta", "psi", "direction"))
    path = metadata_path(context)
    for property in ("phi", "theta", "psi")
        _schema_number(metadata_value(node[property]), "$path/$property")
    end
    _schema_string(metadata_value(node["direction"]), "$path/direction")
    return nothing
end

function _validate_fixed_inputs_schema(node, context)
    _validate_required_transform(node, context, ("forward",))
    path = metadata_path(context)
    forward = _schema_sequence(node["forward"], "$path/forward")
    length(forward) == 2 || _metadata_error("$path/forward", "fixed-input metadata must contain a transform and fixed values")
    _schema_mapping(forward[1], "$path/forward/1")
    specification = _schema_mapping(forward[2], "$path/forward/2")
    haskey(specification, "keys") || _metadata_error("$path/forward/2", "missing required property keys")
    haskey(specification, "values") || _metadata_error("$path/forward/2", "missing required property values")
    keys = _schema_sequence(specification["keys"], "$path/forward/2/keys")
    values = _schema_sequence(specification["values"], "$path/forward/2/values")
    length(keys) == length(values) || _metadata_error("$path/forward/2", "fixed-input keys and values must have equal lengths")
    return nothing
end

function _validate_remap_schema(node, context)
    _validate_required_transform(node, context, ("mapping",))
    path = metadata_path(context)
    mapping = _schema_sequence(node["mapping"], "$path/mapping")
    for (index, value) in enumerate(mapping)
        _schema_integer(metadata_value(value), "$path/mapping/$index")
    end
    haskey(node, "n_inputs") && _schema_integer(metadata_value(node["n_inputs"]), "$path/n_inputs")
    return nothing
end

function _validate_polynomial_schema(node, context)
    _validate_required_transform(node, context, ("coefficients",))
    _schema_metadata_array(node["coefficients"], "$(metadata_path(context))/coefficients")
    return nothing
end

function _validate_rotation_sequence_schema(node, context)
    _validate_required_transform(node, context, ("angles", "axes_order", "rotation_type"))
    path = metadata_path(context)
    angles = _schema_sequence(node["angles"], "$path/angles")
    for (index, angle) in enumerate(angles)
        _schema_number(metadata_value(angle), "$path/angles/$index")
    end
    _schema_string(metadata_value(node["axes_order"]), "$path/axes_order")
    _schema_string(metadata_value(node["rotation_type"]), "$path/rotation_type")
    return nothing
end

function _validate_identity_schema(node, context)
    _validate_required_transform(node, context, ())
    haskey(node, "n_dims") && _schema_integer(metadata_value(node["n_dims"]), "$(metadata_path(context))/n_dims")
    return nothing
end

function _validate_constant_schema(node, context)
    _validate_required_transform(node, context, ("value",))
    path = metadata_path(context)
    _schema_number(metadata_value(node["value"]), "$path/value")
    if metadata_version(context) < v"1.4.0"
        haskey(node, "dimensions") && _metadata_error("$path/dimensions", "property dimensions is not supported by $(metadata_tag(context))-$(metadata_version(context))")
    else
        dimensions = _schema_property(node, "dimensions", context; required=true)
        _schema_integer(dimensions, "$path/dimensions")
        dimensions in (1, 2) || _metadata_error("$path/dimensions", "constant dimensions must be 1 or 2")
    end
    return nothing
end

function _validate_ndarray_schema(node, context)
    path = metadata_path(context)
    _schema_mapping(node, path)
    haskey(node, "data") || _metadata_error(path, "missing required property data")
    if haskey(node, "shape")
        shape = _schema_sequence(node["shape"], "$path/shape")
        for (index, value) in enumerate(shape)
            _schema_integer(metadata_value(value), "$path/shape/$index")
        end
    end
    return nothing
end

function _normalize_constant_schema(node, context)
    metadata_version(context) >= v"1.4.0" && return node
    normalized = Dict{Any,Any}(node)
    normalized["dimensions"] = 0
    return normalized
end

function _validate_wcs_schema(node, context)
    path = metadata_path(context)
    _schema_mapping(node, path)
    _schema_string(_schema_property(node, "name", context; required=true), "$path/name")
    steps = _schema_sequence(_schema_property(node, "steps", context; required=true), "$path/steps")
    isempty(steps) && _metadata_error("$path/steps", "WCS steps must not be empty")
    for (index, step) in enumerate(steps)
        _schema_tag(step, "$_GWCS_TAG/step", "$path/steps/$index")
    end
    if haskey(node, "pixel_shape")
        metadata_version(context) >= v"1.2.0" || _metadata_error("$path/pixel_shape", "property pixel_shape is not supported by $(metadata_tag(context))-$(metadata_version(context))")
        value = metadata_value(node["pixel_shape"])
        value === nothing || _schema_sequence(value, "$path/pixel_shape")
    end
    return nothing
end

function _validate_step_schema(node, context)
    path = metadata_path(context)
    _schema_mapping(node, path)
    _schema_mapping(_schema_property(node, "frame", context; required=true), "$path/frame")
    haskey(node, "transform") && metadata_value(node["transform"]) !== nothing && _schema_mapping(node["transform"], "$path/transform")
    return nothing
end

function _validate_frame_schema(node, context)
    path = metadata_path(context)
    _schema_mapping(node, path)
    _schema_string(_schema_property(node, "name", context; required=true), "$path/name")
    endswith(metadata_tag(context), "/frame") && _schema_integer(_schema_property(node, "naxes", context; required=true), "$path/naxes")
    haskey(node, "axes_names") && _schema_sequence(node["axes_names"], "$path/axes_names")
    haskey(node, "axes_order") && _schema_sequence(node["axes_order"], "$path/axes_order")
    haskey(node, "axes_type") && _schema_sequence(node["axes_type"], "$path/axes_type")
    haskey(node, "axis_physical_types") && _schema_sequence(node["axis_physical_types"], "$path/axis_physical_types")
    haskey(node, "unit") && _schema_sequence(node["unit"], "$path/unit")
    return nothing
end

function _validate_composite_frame_schema(node, context)
    _validate_frame_schema(node, context)
    _schema_sequence(_schema_property(node, "frames", context; required=true), "$(metadata_path(context))/frames")
    return nothing
end

function _validate_label_mapper_schema(node, context)
    path = metadata_path(context)
    _validate_required_transform(node, context, ("mapper",))
    haskey(node, "inputs") && _schema_string_sequence(node["inputs"], "$path/inputs")
    haskey(node, "inputs_mapping") && _schema_mapping(node["inputs_mapping"], "$path/inputs_mapping")
    haskey(node, "atol") && _schema_number(metadata_value(node["atol"]), "$path/atol")
    if haskey(node, "no_label")
        value = metadata_value(node["no_label"])
        value isa Number || value isa AbstractString || _metadata_error("$path/no_label", "expected a number or string")
    end
    return nothing
end

function _validate_region_selector_schema(node, context)
    path = metadata_path(context)
    _validate_required_transform(node, context, ("label_mapper", "inputs", "outputs", "selector"))
    _schema_mapping(node["label_mapper"], "$path/label_mapper")
    _schema_tag(node["label_mapper"], "$_GWCS_TAG/label_mapper", "$path/label_mapper")
    _schema_string_sequence(node["inputs"], "$path/inputs")
    _schema_string_sequence(node["outputs"], "$path/outputs")
    selector = _schema_mapping(node["selector"], "$path/selector")
    haskey(selector, "labels") || _metadata_error("$path/selector", "missing required property labels")
    haskey(selector, "transforms") || _metadata_error("$path/selector", "missing required property transforms")
    labels = _schema_sequence(selector["labels"], "$path/selector/labels")
    transforms = _schema_sequence(selector["transforms"], "$path/selector/transforms")
    length(labels) == length(transforms) || _metadata_error("$path/selector", "label and transform arrays must have equal lengths")
    haskey(node, "undefined_transform_value") && _schema_number(metadata_value(node["undefined_transform_value"]), "$path/undefined_transform_value")
    return nothing
end

function _validate_spherical_cartesian_schema(node, context)
    path = metadata_path(context)
    _validate_required_transform(node, context, ("transform_type", "wrap_lon_at"))
    _schema_string(metadata_value(node["transform_type"]), "$path/transform_type")
    _schema_integer(metadata_value(node["wrap_lon_at"]), "$path/wrap_lon_at")
    return nothing
end

function _validate_bounding_box_schema(node, context)
    path = metadata_path(context)
    _schema_mapping(node, path)
    haskey(node, "order") && _schema_string(metadata_value(node["order"]), "$path/order")
    haskey(node, "ignore") && _schema_string_sequence(node["ignore"], "$path/ignore")
    if endswith(metadata_tag(context), "/bounding_box")
        _schema_mapping(_schema_property(node, "intervals", context; required=true), "$path/intervals")
    else
        _schema_sequence(_schema_property(node, "cbbox", context; required=true), "$path/cbbox")
        _schema_sequence(_schema_property(node, "selector_args", context; required=true), "$path/selector_args")
    end
    return nothing
end

function _register_builtin_metadata_schemas!()
    registry = default_transform_registry()
    _register_metadata_schema!(registry, "$_GWCS_TAG/wcs"; versions=_SCHEMA_VERSIONS_1_0_TO_1_4, validate=_validate_wcs_schema)
    _register_metadata_schema!(registry, "$_GWCS_TAG/step"; versions=_SCHEMA_VERSIONS_1_0_TO_1_3, validate=_validate_step_schema)
    for family in ("frame", "frame2d", "celestial_frame", "spectral_frame", "temporal_frame")
        _register_metadata_schema!(registry, "$_GWCS_TAG/$family"; versions=_SCHEMA_VERSIONS_1_0_TO_1_2, validate=_validate_frame_schema)
    end
    _register_metadata_schema!(registry, "$_GWCS_TAG/stokes_frame"; versions=(v"1.0.0",), validate=_validate_frame_schema)
    _register_metadata_schema!(registry, "$_GWCS_TAG/composite_frame"; versions=(v"1.0.0",), validate=_validate_composite_frame_schema)
    _register_metadata_schema!(registry, "$_GWCS_TAG/label_mapper"; versions=_SCHEMA_VERSIONS_1_0_TO_1_3, validate=_validate_label_mapper_schema)
    _register_metadata_schema!(registry, "$_ASDF_TRANSFORM_TAG/property/bounding_box"; versions=_SCHEMA_VERSIONS_1_0_TO_1_2, validate=_validate_bounding_box_schema)
    _register_metadata_schema!(registry, "$_ASDF_TRANSFORM_TAG/property/compound_bounding_box"; versions=_SCHEMA_VERSIONS_1_0_TO_1_2, validate=_validate_bounding_box_schema)
    _register_metadata_schema!(registry, "tag:astropy.org:astropy/coordinates/frames/icrs"; versions=_SCHEMA_VERSIONS_1_0_TO_1_3)
    _register_metadata_schema!(registry, "tag:stsci.edu:asdf/core/ndarray"; versions=_SCHEMA_VERSIONS_1_0_TO_1_2, validate=_validate_ndarray_schema)
    return nothing
end
