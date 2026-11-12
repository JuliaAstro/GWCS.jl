const _FITS_LONGITUDE_TYPES = ("RA", "GLON", "ELON", "HLON", "SLON", "HPLN")
const _FITS_LATITUDE_TYPES = ("DEC", "GLAT", "ELAT", "HLAT", "SLAT", "HPLT")

_fits_is_longitude_type(coordinate) = coordinate == "RA" || ncodeunits(coordinate) == 4 && (endswith(coordinate, "LON") || endswith(coordinate, "LN"))
_fits_is_latitude_type(coordinate) = coordinate == "DEC" || ncodeunits(coordinate) == 4 && (endswith(coordinate, "LAT") || endswith(coordinate, "LT"))
_fits_celestial_partner(coordinate) = coordinate == "RA" ? "DEC" : coordinate == "DEC" ? "RA" : endswith(coordinate, "LON") ? coordinate[1:1] * "LAT" : endswith(coordinate, "LAT") ? coordinate[1:1] * "LON" : endswith(coordinate, "LN") ? coordinate[1:2] * "LT" : coordinate[1:2] * "LN"
const _FITS_SPECTRAL_TYPES = ("FREQ", "AFRQ", "ENER", "WAVN", "WAVE", "AWAV", "VRAD", "VOPT", "ZOPT", "VELO", "BETA")
const _FITS_TIME_TYPES = ("TIME", "UTC", "TAI", "TT", "TDB", "TCG", "TCB", "GPS")

function _fits_header(source::AbstractDict)
    header = Dict{String,Any}()
    for (raw_key, value) in source
        key = uppercase(strip(String(raw_key)))
        isempty(key) && continue
        haskey(header, key) && header[key] != value && throw(ArgumentError("FITS header contains conflicting values for keyword $key"))
        header[key] = value
    end
    return header
end

function _fits_key(key::Char)
    key == ' ' || 'A' <= key <= 'Z' || throw(ArgumentError("key must be ' ' or 'A'-'Z', got $(repr(key))"))
    return key == ' ' ? "" : string(key)
end

function _fits_number(value, keyword)
    value isa Real && !(value isa Bool) || throw(ArgumentError("FITS keyword $keyword must be numeric"))
    result = Float64(value)
    isfinite(result) || throw(ArgumentError("FITS keyword $keyword must be finite"))
    return result
end

function _fits_integer(value, keyword; minimum=0)
    value isa Integer && !(value isa Bool) || throw(ArgumentError("FITS keyword $keyword must be an integer"))
    result = Int(value)
    result >= minimum || throw(ArgumentError("FITS keyword $keyword must be at least $minimum"))
    return result
end

function _fits_axis_indices(keyword::AbstractString, suffix::AbstractString)
    core = if isempty(suffix)
        keyword
    elseif endswith(keyword, suffix) && ncodeunits(keyword) > 1
        keyword[1:(end - 1)]
    else
        return ()
    end
    match_result = match(r"^(?:CRPIX|CRVAL|CDELT|CTYPE|CUNIT|CNAME|CRDER|CSYER)([1-9][0-9]*)$", core)
    match_result !== nothing && return (parse(Int, only(match_result.captures)),)
    match_result = match(r"^(?:PC|CD)([1-9][0-9]*)_([1-9][0-9]*)$", core)
    match_result !== nothing && return (parse(Int, match_result.captures[1]), parse(Int, match_result.captures[2]))
    match_result = match(r"^(?:PV|PS)([1-9][0-9]*)_([0-9]+)$", core)
    return match_result === nothing ? () : (parse(Int, match_result.captures[1]),)
end

function _fits_axis_count(header, key::Char)
    suffix = _fits_key(key)
    inferred = 0
    for keyword in keys(header), axis in _fits_axis_indices(keyword, suffix)
        inferred = max(inferred, axis)
    end
    wcsaxes = "WCSAXES$suffix"
    if haskey(header, wcsaxes)
        count = _fits_integer(header[wcsaxes], wcsaxes; minimum=1)
    elseif haskey(header, "NAXIS")
        count = max(_fits_integer(header["NAXIS"], "NAXIS"), inferred)
    else
        inferred > 0 || throw(ArgumentError("FITS header contains neither WCSAXES$suffix nor NAXIS"))
        count = inferred
    end
    count >= 1 || throw(ArgumentError("FITS WCS must have at least one axis"))
    for keyword in keys(header), axis in _fits_axis_indices(keyword, suffix)
        axis <= count || throw(ArgumentError("FITS keyword $keyword uses axis $axis but WCSAXES/NAXIS is $count"))
    end
    return count
end

function _fits_ctype(value)
    text = uppercase(strip(String(value)))
    isempty(text) && return ("", :linear, "", "")
    if ncodeunits(text) >= 8 && text[5] == '-'
        coordinate = rstrip(text[1:4], '-')
        algorithm = text[6:8]
        suffix = ncodeunits(text) > 8 ? text[9:end] : ""
    else
        coordinate, algorithm, suffix = text, "", ""
    end
    kind = _fits_is_longitude_type(coordinate) ? :longitude :
        _fits_is_latitude_type(coordinate) ? :latitude :
        coordinate in _FITS_SPECTRAL_TYPES ? :spectral :
        coordinate in _FITS_TIME_TYPES ? :time :
        coordinate == "STOKES" ? :stokes : :linear
    return coordinate, kind, algorithm, suffix
end

function _fits_reject_external_data(header, suffix; allow_lookup=false, allow_tab=false)
    axis_suffix = isempty(suffix) ? raw"$" : "$(suffix)\$"
    for (keyword, value) in header
        if occursin(Regex("^D2IMDIS[1-9][0-9]*$axis_suffix"), keyword) ||
                occursin(Regex("^CPDIS[1-9][0-9]*$axis_suffix"), keyword) && uppercase(strip(String(value))) == "LOOKUP"
            allow_lookup || throw(ArgumentError("FITS WCS references external distortion or table data; parse it from a linked FITS HDU"))
        end
        match_result = match(Regex("^CTYPE[1-9][0-9]*$axis_suffix"), keyword)
        match_result === nothing && continue
        value isa AbstractString && occursin("-TAB", uppercase(value)) && !allow_tab && throw(ArgumentError("FITS WCS references external distortion or table data; parse it from a linked FITS HDU"))
    end
end

struct _FITSTabularSpec{I}
    axis::Int
    table_axis::Int
    extname::String
    extver::Int
    extlev::Int
    coordinate_column::String
    index_column::I
end

function _fits_tabular_specs(header, key, parsed)
    suffix = _fits_key(key)
    specs = _FITSTabularSpec[]
    for axis in eachindex(parsed)
        parsed[axis][3] == "TAB" || continue
        extension_keyword = "PS$(axis)_0$suffix"
        coordinate_keyword = "PS$(axis)_1$suffix"
        haskey(header, extension_keyword) || throw(ArgumentError("FITS TAB axis CTYPE$axis requires $extension_keyword"))
        haskey(header, coordinate_keyword) || throw(ArgumentError("FITS TAB axis CTYPE$axis requires $coordinate_keyword"))
        table_axis_keyword = "PV$(axis)_3$suffix"
        extver_keyword = "PV$(axis)_1$suffix"
        extlev_keyword = "PV$(axis)_2$suffix"
        table_axis = _fits_integer(get(header, table_axis_keyword, 1), table_axis_keyword; minimum=1)
        extver = _fits_integer(get(header, extver_keyword, 1), extver_keyword; minimum=1)
        extlev = _fits_integer(get(header, extlev_keyword, 1), extlev_keyword; minimum=1)
        index_keyword = "PS$(axis)_2$suffix"
        index_column = haskey(header, index_keyword) ? String(header[index_keyword]) : nothing
        push!(specs, _FITSTabularSpec(axis, table_axis, String(header[extension_keyword]), extver, extlev, String(header[coordinate_keyword]), index_column))
    end
    return specs
end

function _fits_tabular_coordinates(data, dimensions)
    dimensions == 1 && return reshape(Float64.(vec(data)), 1, :)
    ndims(data) >= 2 || throw(ArgumentError("TAB coordinate array for $dimensions axes must have at least two dimensions"))
    size(data, 1) == dimensions && return Array{Float64}(data)
    size(data, ndims(data)) == dimensions || throw(ArgumentError("TAB coordinate array must have a component dimension of length $dimensions"))
    return Array{Float64}(permutedims(data, (ndims(data), (1:(ndims(data) - 1))...)))
end

function _fits_tabular_data(loader, header, key)
    parsed = _fits_ctype.(String.([get(header, "CTYPE$i$(_fits_key(key))", "") for i in 1:_fits_axis_count(header, key)]))
    specs = _fits_tabular_specs(header, key, parsed)
    isempty(specs) && return ()
    groups = Dict{Tuple{String,Int,Int,String},Vector{_FITSTabularSpec}}()
    for spec in specs
        group_key = (spec.extname, spec.extver, spec.extlev, spec.coordinate_column)
        push!(get!(groups, group_key, _FITSTabularSpec[]), spec)
    end
    tables = map(collect(values(groups))) do group
        ordered = sort(group; by=spec -> spec.table_axis)
        table_axes = getproperty.(ordered, :table_axis)
        table_axes == collect(1:length(ordered)) || throw(ArgumentError("TAB table axes must define every dimension from 1 through $(length(ordered))"))
        first_spec = first(ordered)
        coordinate_data = loader(first_spec.extname, first_spec.extver, first_spec.extlev, first_spec.coordinate_column)
        coordinates = _fits_tabular_coordinates(coordinate_data, length(ordered))
        ndims(coordinates) == length(ordered) + 1 || throw(ArgumentError("TAB coordinate array has $(ndims(coordinates) - 1) table axes, expected $(length(ordered))"))
        indices = ntuple(length(ordered)) do position
            spec = ordered[position]
            length = size(coordinates, position + 1)
            spec.index_column === nothing ? collect(1.0:length) : Float64.(vec(loader(spec.extname, spec.extver, spec.extlev, spec.index_column)))
        end
        FITSTabularTable(Tuple(getproperty.(ordered, :axis)), coordinates, indices)
    end
    return tuple(tables...)
end

struct _FITSLookupSpec{N}
    extname::String
    extver::Int
    axes::NTuple{N,Int}
    naxis::Int
end

function _fits_lookup_spec(header, axis, count, distortion_prefix, parameter_prefix, error_prefix, extname, suffix, minerr)
    distortion_keyword = "$distortion_prefix$axis$suffix"
    haskey(header, distortion_keyword) || return nothing
    distortion_type = uppercase(strip(String(header[distortion_keyword])))
    # Header-only CPDIS TPD and SIP records are parsed by `_parse_fits_wcs`; D2IMDIS supports only LOOKUP.
    distortion_prefix == "CPDIS" && distortion_type in ("TPD", "SIP") && return nothing
    distortion_type == "LOOKUP" || throw(ArgumentError("unsupported Paper IV distortion type $(repr(header[distortion_keyword])) in $distortion_keyword"))
    error_keyword = "$error_prefix$axis$suffix"
    error_estimate = _fits_number(get(header, error_keyword, 0.0), error_keyword)
    error_estimate < minerr && return nothing
    parameter = "$parameter_prefix$axis$suffix"
    naxes_keyword = "$parameter.NAXES"
    haskey(header, naxes_keyword) || throw(ArgumentError("Paper IV lookup distortion $distortion_keyword requires $naxes_keyword"))
    dimensions = _fits_integer(header[naxes_keyword], naxes_keyword; minimum=1)
    dimensions <= min(count, 2) || throw(ArgumentError("$naxes_keyword must be 1 or 2 and not exceed the $count WCS axes, got $dimensions"))
    axes = ntuple(dimensions) do table_axis
        axis_keyword = "$parameter.AXIS.$table_axis"
        _fits_integer(get(header, axis_keyword, table_axis), axis_keyword; minimum=1)
    end
    all(value -> value <= count, axes) || throw(ArgumentError("$parameter axis map contains an axis outside 1:$count"))
    allunique(axes) || throw(ArgumentError("$parameter axis map contains duplicate axes"))
    extver_keyword = "$parameter.EXTVER"
    extver = _fits_integer(get(header, extver_keyword, 1), extver_keyword; minimum=1)
    return _FITSLookupSpec(extname, extver, axes, count)
end

function _fits_lookup_data(loader, header, key, minerr)
    isfinite(minerr) && minerr >= 0 || throw(ArgumentError("minerr must be finite and nonnegative"))
    suffix = _fits_key(key)
    referenced = any(axis -> haskey(header, "D2IMDIS$axis$suffix") || haskey(header, "CPDIS$axis$suffix"), 1:2)
    referenced || return nothing
    count = _fits_axis_count(header, key)
    count >= 2 || throw(ArgumentError("FITS lookup distortions require at least two pixel axes"))
    detector_specs = ntuple(axis -> _fits_lookup_spec(header, axis, count, "D2IMDIS", "D2IM", "D2IMERR", "D2IMARR", suffix, minerr), 2)
    cpdis_specs = ntuple(axis -> _fits_lookup_spec(header, axis, count, "CPDIS", "DP", "CPERR", "WCSDVARR", suffix, minerr), 2)
    load(spec) = spec === nothing ? nothing : loader(spec)
    return (detector_to_image=map(load, detector_specs), cpdis=map(load, cpdis_specs))
end

function _fits_lookup_table(data::AbstractArray, header, axes::NTuple{N,Int}, naxis) where {N}
    ndims(data) == N || throw(ArgumentError("FITS WCS lookup extension declares $N independent axes but contains a $(ndims(data))-dimensional image"))
    return FITSLookupTable(copy(data), axes, Val(naxis);
        crpix=ntuple(i -> get(header, "CRPIX$i", 0.0), N), crval=ntuple(i -> get(header, "CRVAL$i", 0.0), N), cdelt=ntuple(i -> get(header, "CDELT$i", 1.0), N))
end

function _fits_tpd_index(xpower, ypower)
    for index in 0:59
        radial, xp, yp = _tpd_term_descriptor(index)
        !radial && xp == xpower && yp == ypower && return index
    end
    throw(ArgumentError("legacy FITS polynomial term x^$xpower*y^$ypower exceeds TPD degree 9"))
end

function _fits_set_tpd_axis!(header, axis, coefficients, auxiliary; docorr)
    prefix = "DQ$axis"
    header["CQDIS$axis"] = "TPD"
    header["$prefix.NAXES"] = 2
    header["$prefix.DOCORR"] = docorr ? 1 : 0
    for (position, coefficient) in pairs(coefficients)
        iszero(coefficient) || (header["$prefix.TPD.FWD.$(position - 1)"] = coefficient)
    end
    if auxiliary !== nothing
        for variable in 1:2, coefficient in 0:2
            header["$prefix.AUX.$variable.COEFF.$coefficient"] = auxiliary[variable, coefficient + 1]
        end
    end
    return header
end

function _fits_dss_coefficients(header, axis, scale)
    coefficients = zeros(Float64, 24)
    prefix = axis == 1 ? "AMDX" : "AMDY"
    for (keyword, value) in header
        matched = match(Regex("^$(prefix)([0-9]+)\$"), keyword)
        matched === nothing && continue
        parameter = parse(Int, only(matched.captures))
        coefficient = _fits_number(value, keyword)
        iszero(coefficient) && continue
        1 <= parameter <= 13 || throw(ArgumentError("DSS polynomial $keyword must use parameter numbers 1 through 13"))
        coefficient /= scale
        if parameter == 1
            coefficients[2] = coefficient
        elseif parameter == 2
            coefficients[3] = coefficient
        elseif parameter == 3
            coefficients[1] = coefficient
        elseif parameter == 4
            coefficients[5] += coefficient
        elseif parameter == 5
            coefficients[6] = coefficient
        elseif parameter == 6
            coefficients[7] += coefficient
        elseif parameter == 7
            coefficients[5] += coefficient
            coefficients[7] += coefficient
        elseif parameter == 8
            coefficients[8] += coefficient
        elseif parameter == 9
            coefficients[9] = coefficient
        elseif parameter == 10
            coefficients[10] += coefficient
        elseif parameter == 11
            coefficients[11] = coefficient
        elseif parameter == 12
            coefficients[8] += coefficient
            coefficients[10] += coefficient
        else
            coefficients[18] = coefficient
            coefficients[20] = 2coefficient
            coefficients[22] = coefficient
        end
    end
    return coefficients
end

function _fits_normalize_dss(header)
    any(keyword -> occursin(r"^AMDX[0-9]+$", keyword), keys(header)) || return header
    required = ("CNPIX1", "CNPIX2", "PPO3", "PPO6", "XPIXELSZ", "YPIXELSZ", "PLTRAH", "PLTRAM", "PLTRAS", "PLTDECSN", "PLTDECD", "PLTDECM", "PLTDECS", "AMDX1", "AMDX2", "AMDX3", "AMDY1", "AMDY2", "AMDY3")
    missing = filter(keyword -> !haskey(header, keyword), required)
    isempty(missing) || throw(ArgumentError("DSS WCS is missing required keyword $(first(missing))"))
    number(keyword) = _fits_number(header[keyword], keyword)
    A1, A2, A3 = number("AMDX1"), number("AMDX2"), number("AMDX3")
    B1, B2, B3 = number("AMDY1"), number("AMDY2"), number("AMDY3")
    determinant = A1 * B1 - A2 * B2
    scale = sqrt(abs(determinant))
    iszero(scale) && throw(ArgumentError("DSS coefficient scale is zero"))
    X0 = (A2 * B3 - A3 * B1) / determinant
    Y0 = (A3 * B2 - A1 * B3) / determinant
    Xc, Yc = number("PPO3") / 1000, number("PPO6") / 1000
    Rx, Ry = number("XPIXELSZ") / 1000, number("YPIXELSZ") / 1000
    !iszero(Rx) && !iszero(Ry) || throw(ArgumentError("DSS pixel sizes must be nonzero"))
    normalized = copy(header)
    normalized["CRPIX1"] = (Xc - X0) / Rx - (number("CNPIX1") - 0.5)
    normalized["CRPIX2"] = (Yc + Y0) / Ry - (number("CNPIX2") - 0.5)
    normalized["PC1_1"], normalized["PC1_2"] = A1 * Rx / scale, -A2 * Ry / scale
    normalized["PC2_1"], normalized["PC2_2"] = -B2 * Rx / scale, B1 * Ry / scale
    normalized["CDELT1"], normalized["CDELT2"] = -scale / 3600, scale / 3600
    normalized["CRVAL1"] = (number("PLTRAH") + (number("PLTRAM") + number("PLTRAS") / 60) / 60) * 15
    declination = number("PLTDECD") + (number("PLTDECM") + number("PLTDECS") / 60) / 60
    normalized["CRVAL2"] = startswith(strip(String(header["PLTDECSN"])), "-") ? -declination : declination
    normalized["CTYPE1"], normalized["CTYPE2"] = "RA---TAN", "DEC--TAN"
    normalized["WCSNAME"] = "DSS PLATEID $(first(String(get(header, "PLATEID", "    ")), min(4, length(String(get(header, "PLATEID", "    "))))))"
    for keyword in collect(keys(normalized))
        occursin(r"^CD[12]_[12]$", keyword) && delete!(normalized, keyword)
    end
    auxiliary1 = [X0 -B1 / scale -A2 / scale; Y0 B2 / scale A1 / scale]
    auxiliary2 = [Y0 B2 / scale A1 / scale; X0 -B1 / scale -A2 / scale]
    _fits_set_tpd_axis!(normalized, 1, _fits_dss_coefficients(header, 1, -scale), auxiliary1; docorr=false)
    _fits_set_tpd_axis!(normalized, 2, _fits_dss_coefficients(header, 2, scale), auxiliary2; docorr=false)
    return normalized
end

function _fits_wat_string(header, axis)
    cards = Pair{Int,String}[]
    pattern = Regex("^WAT$(axis)_([0-9]+)\$")
    for (keyword, value) in header
        matched = match(pattern, keyword)
        matched === nothing || push!(cards, parse(Int, only(matched.captures)) => String(value))
    end
    isempty(cards) && return nothing
    sort!(cards; by=first)
    return join(last.(cards))
end

function _fits_basis_coefficients(kind, degree)
    degree == 0 && return [1.0]
    degree == 1 && return [0.0, 1.0]
    previous, current = [1.0], [0.0, 1.0]
    for order in 2:degree
        next = zeros(Float64, order + 1)
        if kind == 1
            next[2:end] .+= 2 .* current
            next[1:length(previous)] .-= previous
        else
            next[2:end] .+= ((2order - 1) / order) .* current
            next[1:length(previous)] .-= ((order - 1) / order) .* previous
        end
        previous, current = current, next
    end
    return current
end

function _fits_wat_type(text)
    type_match = match(r"wtype=([^ ]+)", lowercase(text))
    return type_match === nothing ? nothing : only(type_match.captures)
end

function _fits_wat_axis(text, axis)
    wtype = _fits_wat_type(text)
    wtype === nothing && throw(ArgumentError("WAT$axis distortion has no wtype"))
    wtype in ("tnx", "zpx") || throw(ArgumentError("unsupported WAT$axis wtype $(repr(wtype))"))
    correction_match = match(axis == 1 ? r"lngcor\s*=\s*\"([^\"]*)\"" : r"latcor\s*=\s*\"([^\"]*)\"", lowercase(text))
    correction_match === nothing && throw(ArgumentError("WAT$axis $wtype distortion has no correction polynomial"))
    values = tryparse.(Float64, split(only(correction_match.captures)))
    any(isnothing, values) && throw(ArgumentError("WAT$axis correction polynomial contains a nonnumeric value"))
    numbers = something.(values)
    length(numbers) >= 8 || throw(ArgumentError("WAT$axis correction polynomial is incomplete"))
    kind, xorder, yorder, cross = Int.(numbers[1:4])
    kind in 1:3 || throw(ArgumentError("WAT$axis polynomial type must be 1, 2, or 3"))
    1 <= xorder <= 10 && 1 <= yorder <= 10 || throw(ArgumentError("WAT$axis polynomial orders must be between 1 and 10"))
    cross in 0:2 || throw(ArgumentError("WAT$axis cross-term type must be 0, 1, or 2"))
    xmin, xmax, ymin, ymax = numbers[5:8]
    coefficients = zeros(Float64, 60)
    cursor = 9
    maxorder = max(xorder, yorder)
    for n in 0:(yorder - 1), m in 0:(xorder - 1)
        cross == 0 && m != 0 && n != 0 && continue
        cross == 2 && m + n > maxorder - 1 && continue
        cursor <= length(numbers) || throw(ArgumentError("WAT$axis correction polynomial has too few coefficients"))
        coefficient = numbers[cursor]
        cursor += 1
        if kind == 3
            coefficients[_fits_tpd_index(m, n) + 1] = coefficient
        else
            mx, ny = _fits_basis_coefficients(kind, m), _fits_basis_coefficients(kind, n)
            for im in 0:m, jn in 0:n
                coefficients[_fits_tpd_index(im, jn) + 1] += coefficient * mx[im + 1] * ny[jn + 1]
            end
        end
    end
    cursor == length(numbers) + 1 || throw(ArgumentError("WAT$axis correction polynomial has too many coefficients"))
    auxiliary = if kind in (1, 2)
        dx, dy = (xmax - xmin) / 2, (ymax - ymin) / 2
        !iszero(dx) && !iszero(dy) || throw(ArgumentError("WAT$axis polynomial domain has zero span"))
        [-((xmax + xmin) / 2) / dx 1 / dx 0; -((ymax + ymin) / 2) / dy 0 1 / dy]
    else
        nothing
    end
    return wtype, coefficients, auxiliary
end

function _fits_normalize_wat(header)
    texts = (_fits_wat_string(header, 1), _fits_wat_string(header, 2))
    all(isnothing, texts) && return header
    # Like WCSLIB, ignore WAT strings that do not encode TNX/ZPX on both axes (for example IRAF wtype=linear),
    # unless CTYPE declares TNX or ZPX and therefore requires them.
    declared = any(axis -> _fits_ctype(get(header, "CTYPE$axis", ""))[3] in ("TNX", "ZPX"), 1:2)
    declared || all(text -> text !== nothing && _fits_wat_type(text) in ("tnx", "zpx"), texts) || return header
    all(!isnothing, texts) || throw(ArgumentError("WAT distortion requires both WAT1_n and WAT2_n strings"))
    first_axis = _fits_wat_axis(texts[1], 1)
    second_axis = _fits_wat_axis(texts[2], 2)
    first_axis[1] == second_axis[1] || throw(ArgumentError("WAT axes declare different projection types"))
    normalized = copy(header)
    code = first_axis[1] == "tnx" ? "TAN" : "ZPN"
    normalized["CTYPE1"], normalized["CTYPE2"] = "RA---$code", "DEC--$code"
    if code == "ZPN"
        for matched in eachmatch(r"projp([0-9]+)=([^ ]+)", lowercase(texts[2]))
            parameter = parse(Int, matched.captures[1])
            parameter < 30 || throw(ArgumentError("ZPX projection parameters must use indices 0 through 29"))
            value = tryparse(Float64, matched.captures[2])
            value === nothing && throw(ArgumentError("invalid ZPX projp$parameter value"))
            normalized["PV2_$parameter"] = value
        end
    end
    _fits_set_tpd_axis!(normalized, 1, first_axis[2], first_axis[3]; docorr=true)
    _fits_set_tpd_axis!(normalized, 2, second_axis[2], second_axis[3]; docorr=true)
    return normalized
end

_fits_normalize_legacy(header, suffix) = isempty(suffix) ? _fits_normalize_wat(_fits_normalize_dss(header)) : header

const _FITS_LEGACY_MODES = (:warn, :quiet, :error)

function _fits_report_repairs(repairs, key, legacy)
    isempty(repairs) && return nothing
    message = "FITS WCS key $(repr(key)) uses legacy conventions:\n" * join(("  - " * repair for repair in repairs), "\n")
    legacy == :error && throw(ArgumentError(message * "\nPass legacy=:warn or legacy=:quiet to repair them as WCSLIB and astropy do."))
    legacy == :warn && @warn message
    return nothing
end

# Copy-on-write header edits; `rename!` moves a legacy keyword to its standard spelling.
mutable struct _FITSHeaderEdit
    header::Dict{String,Any}
    copied::Bool
end

function _fits_edit!(edit::_FITSHeaderEdit, keyword, value)
    edit.copied || (edit.header = copy(edit.header); edit.copied = true)
    edit.header[keyword] = value
end

function _fits_rename!(edit::_FITSHeaderEdit, repairs, legacy, standard; conflict=:error, value=edit.header[legacy], note="")
    if haskey(edit.header, standard)
        edit.header[standard] == value || conflict != :error || throw(ArgumentError("FITS keywords $legacy and $standard give conflicting values"))
        return nothing
    end
    _fits_edit!(edit, standard, value)
    push!(repairs, "$legacy read as $standard$note")
end

const _FITS_DATE_SYNONYMS = ("DATE-REF" => "DATEREF", "MJD-REF" => "MJDREF", "MJD-REFI" => "MJDREFI", "MJD-REFF" => "MJDREFF",
    "JD-REF" => "JDREF", "JD-REFI" => "JDREFI", "JD-REFF" => "JDREFF")

# Keyword-level legacy forms accepted by WCSLIB's wcspih with all relaxations (astropy's default).
function _fits_normalize_keywords(header, key, repairs)
    suffix = _fits_key(key)
    edit = _FITSHeaderEdit(header, false)
    for keyword in collect(keys(header))
        matched = match(r"^(PC|CD)00([1-9])00([1-9])$", keyword)
        if matched !== nothing
            isempty(suffix) && _fits_rename!(edit, repairs, keyword, "$(matched[1])$(matched[2])_$(matched[3])")
            continue
        end
        matched = match(Regex("^(PC|CD|PV|PS)0*([1-9][0-9]*)_([0-9]+)$suffix\$"), keyword)
        matched === nothing && continue
        standard = "$(matched[1])$(matched[2])_$(parse(Int, matched[3]))$suffix"
        standard == keyword || _fits_rename!(edit, repairs, keyword, standard)
    end
    if isempty(suffix)
        projp = filter(keyword -> occursin(r"^PROJP[0-9]$", keyword), collect(keys(header)))
        latitude = nothing
        for (keyword, value) in header
            matched = match(r"^CTYPE([1-9][0-9]*)$", keyword)
            matched !== nothing && value isa AbstractString && _fits_ctype(value)[2] == :latitude && (latitude = parse(Int, matched[1]))
        end
        latitude === nothing || foreach(keyword -> _fits_rename!(edit, repairs, keyword, "PV$(latitude)_$(keyword[end])"), projp)
        haskey(header, "RADECSYS") && _fits_rename!(edit, repairs, "RADECSYS", "RADESYS"; conflict=:ignore)
    end
    # Like WCSLIB, the standard keyword wins when both spellings are present.
    haskey(header, "EPOCH$suffix") && _fits_rename!(edit, repairs, "EPOCH$suffix", "EQUINOX$suffix"; conflict=:ignore)
    if haskey(header, "VSOURCE$suffix")
        beta = _fits_number(header["VSOURCE$suffix"], "VSOURCE$suffix") / _FITS_SPEED_OF_LIGHT
        abs(beta) < 1 || throw(ArgumentError("FITS VSOURCE$suffix must be slower than light"))
        _fits_rename!(edit, repairs, "VSOURCE$suffix", "ZSOURCE$suffix"; conflict=:ignore, value=(1 + beta) / sqrt(1 - beta^2) - 1, note=" (relativistic Doppler redshift)")
    end
    for (legacy, standard) in _FITS_DATE_SYNONYMS
        haskey(header, legacy) && _fits_rename!(edit, repairs, legacy, standard; conflict=:ignore)
    end
    return edit.header
end

const _FITS_AIPS_FRAMES = ("LSRK", "BARYCENT", "TOPOCENT", "LSRD", "GEOCENTR", "SOURCE", "GALACTOC")
const _FITS_AIPS_SUFFIX_FRAMES = Dict("LSR" => "LSRK", "HEL" => "BARYCENT", "OBS" => "TOPOCENT")

# AIPS celestial (NCP, GLS) and spectral ({FREQ,VELO,FELO}-{LSR,HEL,OBS}, VELREF) types, translated as WCSLIB's wcsset does.
function _fits_normalize_aips(header, key, count, repairs)
    suffix = _fits_key(key)
    edit = _FITSHeaderEdit(header, false)
    ctype(axis) = uppercase(strip(String(get(header, "CTYPE$axis$suffix", ""))))
    kinds = [_fits_ctype(ctype(axis))[2] for axis in 1:count]
    longitude, latitude = findfirst(==(:longitude), kinds), findfirst(==(:latitude), kinds)
    if longitude !== nothing && latitude !== nothing
        lon_text, lat_text = ctype(longitude), ctype(latitude)
        code = ncodeunits(lon_text) >= 8 && ncodeunits(lat_text) >= 8 ? (lon_text[6:8], lat_text[6:8]) : ("", "")
        if code == ("NCP", "NCP") || code == ("GLS", "GLS")
            unit = _fits_unit_translation(String(get(header, "CUNIT$latitude$suffix", "")))
            reference = _fits_number(get(header, "CRVAL$latitude$suffix", 0.0), "CRVAL$latitude$suffix") * _fits_angle_scale(unit, "CUNIT$latitude$suffix")
            replacement = first(code) == "NCP" ? "SIN" : "SFL"
            _fits_edit!(edit, "CTYPE$longitude$suffix", lon_text[1:5] * replacement * lon_text[9:end])
            _fits_edit!(edit, "CTYPE$latitude$suffix", lat_text[1:5] * replacement * lat_text[9:end])
            if replacement == "SIN"
                iszero(reference) && throw(ArgumentError("FITS NCP projection is undefined for a zero reference latitude (CRVAL$latitude$suffix), as in WCSLIB"))
                _fits_edit!(edit, "PV$(latitude)_1$suffix", 0.0)
                _fits_edit!(edit, "PV$(latitude)_2$suffix", cosd(reference) / sind(reference))
                push!(repairs, "NCP projection read as SIN with PV$(latitude)_1$suffix=0 and PV$(latitude)_2$suffix=cot(CRVAL$latitude$suffix)")
            else
                # AIPS GLS translates the map to the reference latitude rather than tilting the graticule.
                _fits_edit!(edit, "PV$(longitude)_0$suffix", 1.0)
                _fits_edit!(edit, "PV$(longitude)_1$suffix", 0.0)
                _fits_edit!(edit, "PV$(longitude)_2$suffix", reference)
                push!(repairs, "GLS projection read as SFL with its fiducial point offset to CRVAL$latitude$suffix")
            end
        end
    end
    velref_keyword = "VELREF$suffix"
    velref = haskey(header, velref_keyword) ? _fits_integer(header[velref_keyword], velref_keyword) : 0
    frame_code, velocity_type = mod(velref, 256), div(velref, 256)
    frame_code <= length(_FITS_AIPS_FRAMES) || throw(ArgumentError("FITS $velref_keyword=$velref has an invalid AIPS reference frame code $frame_code"))
    for axis in 1:count
        text = ctype(axis)
        ncodeunits(text) >= 4 && text[1:4] in ("FREQ", "VELO", "FELO") || continue
        code = ncodeunits(text) == 4 ? "" : ncodeunits(text) == 8 && text[5] == '-' ? text[6:8] : continue
        suffix_frame = get(_FITS_AIPS_SUFFIX_FRAMES, code, nothing)
        suffix_frame === nothing && !(code in ("", "LOG", "TAB")) && continue
        frame = frame_code > 0 ? _FITS_AIPS_FRAMES[frame_code] : suffix_frame
        coordinate = text[1:4]
        if coordinate == "FELO"
            code in ("LOG", "TAB") && throw(ArgumentError("FITS CTYPE$axis$suffix=$(repr(text)) combines the AIPS FELO type with $code"))
            translated = "VOPT-F2W"
        elseif frame === nothing
            continue
        elseif coordinate == "VELO"
            velocity_type in (0, 1) || throw(ArgumentError("FITS $velref_keyword=$velref has an invalid AIPS velocity type $velocity_type"))
            translated = (velocity_type == 0 ? "VOPT" : "VRAD") * (code in ("LOG", "TAB") ? "-$code" : "")
        else
            translated = suffix_frame === nothing ? text : "FREQ"
        end
        translated == text || _fits_edit!(edit, "CTYPE$axis$suffix", translated)
        specsys = frame !== nothing && !haskey(header, "SPECSYS$suffix")
        specsys && _fits_edit!(edit, "SPECSYS$suffix", frame)
        (translated != text || specsys) && push!(repairs, "AIPS spectral CTYPE$axis$suffix=$(repr(text))$(velref == 0 ? "" : " with $velref_keyword=$velref") read as $(repr(translated))$(specsys ? " with SPECSYS$suffix=$(repr(frame))" : "")")
    end
    return edit.header
end

function _fits_projection(code, header, latitude_axis, suffix)
    parameter(index, default) = _fits_number(get(header, "PV$(latitude_axis)_$(index)$suffix", default), "PV$(latitude_axis)_$(index)$suffix")
    function required(index)
        keyword = "PV$(latitude_axis)_$(index)$suffix"
        haskey(header, keyword) || throw(ArgumentError("FITS $code projection requires $keyword, as in WCSLIB"))
        return parameter(index, 0.0)
    end
    code == "TAN" && return GnomonicProjection()
    code == "STG" && return StereographicProjection()
    code == "ARC" && return ZenithalEquidistant()
    code == "ZEA" && return ZenithalEqualArea()
    code == "SIN" && return SlantOrthographic(parameter(1, 0.0), parameter(2, 0.0))
    code == "AZP" && return ZenithalPerspective(parameter(1, 0.0), parameter(2, 0.0))
    code == "SZP" && return SlantZenithalPerspective(parameter(1, 0.0), parameter(2, 0.0), parameter(3, 90.0))
    code == "AIR" && return AiryProjection(parameter(1, 90.0))
    if code == "ZPN"
        final_index = maximum((index for index in 0:30 if haskey(header, "PV$(latitude_axis)_$(index)$suffix")); init=1)
        # WCSLIB defaults every ZPN coefficient to zero and rejects degenerate polynomials.
        coefficients = ntuple(index -> parameter(index - 1, 0.0), final_index + 1)
        return ZenithalPolynomial(coefficients)
    end
    code == "CAR" && return PlateCarreeProjection()
    code == "CEA" && return CylindricalEqualArea(parameter(1, 1.0))
    code == "CYP" && return CylindricalPerspective(parameter(1, 1.0), parameter(2, 1.0))
    code == "MER" && return MercatorProjection()
    code == "SFL" && return SansonFlamsteedProjection()
    code == "PAR" && return ParabolicProjection()
    code == "MOL" && return MollweideProjection()
    code == "PCO" && return PolyconicProjection()
    code == "AIT" && return HammerAitoffProjection()
    code == "COP" && return ConicPerspective(required(1), parameter(2, 0.0))
    code == "COD" && return ConicEquidistant(required(1), parameter(2, 0.0))
    code == "COE" && return ConicEqualArea(required(1), parameter(2, 0.0))
    code == "COO" && return ConicOrthomorphic(required(1), parameter(2, 0.0))
    code == "BON" && return BonneProjection(required(1))
    code == "TSC" && return TangentialSphericalCube()
    code == "CSC" && return COBESphericalCube()
    code == "QSC" && return QuadrilateralizedSphericalCube()
    if code == "HPX"
        H, K = parameter(1, 4.0), parameter(2, 3.0)
        isinteger(H) && isinteger(K) || throw(ArgumentError("HPX H and K parameters must be integers"))
        return HEALPixProjection(Int(H), Int(K))
    end
    code == "XPH" && return HEALPixPolarProjection()
    code == "TPV" && return GnomonicProjection()
    code == "TPD" && throw(ArgumentError("TPD is a Paper IV distortion function, not a celestial projection code"))
    throw(ArgumentError("unsupported FITS spherical projection code $(repr(code))"))
end

function _fits_tpv_coefficients(header, axis, suffix)
    pattern = Regex("^PV$(axis)_([0-9]+)$(suffix)\$")
    indices = Int[]
    for keyword in keys(header)
        matched = match(pattern, keyword)
        matched === nothing && continue
        index = parse(Int, only(matched.captures))
        index <= 39 || throw(ArgumentError("TPV coefficients must use parameter numbers 0 through 39"))
        push!(indices, index)
    end
    isempty(indices) && return ()
    final_index = maximum(indices)
    return ntuple(final_index + 1) do position
        keyword = "PV$(axis)_$(position - 1)$suffix"
        _fits_number(get(header, keyword, 0.0), keyword)
    end
end

function _fits_tpd_coefficients(header, prefix, direction)
    pattern = Regex("^$(prefix)\\.TPD\\.$(direction)\\.([0-9]+)\$")
    found = Pair{Int,Float64}[]
    for keyword in keys(header)
        matched = match(pattern, keyword)
        matched === nothing && continue
        index = parse(Int, only(matched.captures))
        index <= 59 || throw(ArgumentError("TPD coefficients must use parameter numbers 0 through 59"))
        push!(found, index => _fits_number(header[keyword], keyword))
    end
    isempty(found) && return ()
    final_index = maximum(first, found)
    values = zeros(Float64, final_index + 1)
    for (index, value) in found
        values[index + 1] = value
    end
    return Tuple(values)
end

function _fits_tpd_auxiliary(header, prefix)
    has_auxiliary = any(haskey(header, "$prefix.AUX.$variable.COEFF.$coefficient") for variable in 1:2 for coefficient in 0:2)
    has_auxiliary || return nothing
    values = zeros(Float64, 2, 3)
    values[1, 2] = 1.0
    values[2, 3] = 1.0
    for variable in 1:2, coefficient in 0:2
        keyword = "$prefix.AUX.$variable.COEFF.$coefficient"
        haskey(header, keyword) && (values[variable, coefficient + 1] = _fits_number(header[keyword], keyword))
    end
    return SMatrix{2,3}(values)
end

function _fits_valid_tpd_field(field)
    field in ("DOCORR", "NAXES") && return true
    occursin(r"^(AXIS|OFFSET|SCALE)\.[1-9][0-9]*$", field) && return true
    occursin(r"^TPD\.(FWD|REV)\.([0-9]|[1-5][0-9])$", field) && return true
    return occursin(r"^AUX\.[12]\.COEFF\.[0-2]$", field)
end

function _fits_sip_record_terms(header, prefix, direction)
    pattern = Regex("^$(prefix)\\.SIP\\.$(direction)\\.([0-9]+)_([0-9]+)\$")
    terms = Tuple{Float64,Bool,Int,Int}[]
    for (keyword, value) in header
        matched = match(pattern, keyword)
        matched === nothing && continue
        xpower, ypower = parse.(Int, matched.captures)
        xpower <= 9 && ypower <= 9 && xpower + ypower <= 9 || throw(ArgumentError("invalid SIP powers ($xpower, $ypower) in $keyword"))
        coefficient = _fits_number(value, keyword)
        iszero(coefficient) || push!(terms, (Float64(coefficient), false, xpower, ypower))
    end
    return Tuple(terms)
end

function _fits_sip_record_axis(header, axis, count, suffix, minerr)
    distortion_keyword = "CPDIS$axis$suffix"
    haskey(header, distortion_keyword) || return nothing
    uppercase(strip(String(header[distortion_keyword]))) == "SIP" || return nothing
    error_keyword = "CPERR$axis$suffix"
    _fits_number(get(header, error_keyword, 0.0), error_keyword) < minerr && return nothing
    prefix = "DP$axis$suffix"
    naxes_keyword = "$prefix.NAXES"
    haskey(header, naxes_keyword) || throw(ArgumentError("SIP distortion $distortion_keyword requires $naxes_keyword"))
    _fits_integer(header[naxes_keyword], naxes_keyword) == 2 || throw(ArgumentError("$naxes_keyword must be 2 for SIP"))
    axes = ntuple(2) do variable
        keyword = "$prefix.AXIS.$variable"
        _fits_integer(get(header, keyword, variable), keyword; minimum=1)
    end
    all(value -> value <= count, axes) || throw(ArgumentError("$prefix axis map contains an axis outside 1:$count"))
    allunique(axes) || throw(ArgumentError("$prefix axis map contains duplicate axes"))
    offset = SVector{2}(ntuple(variable -> _fits_number(get(header, "$prefix.OFFSET.$variable", 0.0), "$prefix.OFFSET.$variable"), 2))
    scale = SVector{2}(ntuple(variable -> _fits_number(get(header, "$prefix.SCALE.$variable", 1.0), "$prefix.SCALE.$variable"), 2))
    dotted_prefix = "$prefix."
    for keyword in keys(header)
        startswith(keyword, dotted_prefix) || continue
        field = keyword[(length(dotted_prefix) + 1):end]
        # Like WCSLIB's sipset, reject DOCORR: SIP is always an additive correction.
        field == "NAXES" || occursin(r"^(AXIS|OFFSET|SCALE)\.[12]$", field) || occursin(r"^SIP\.(FWD|REV)\.[0-9]+_[0-9]+$", field) || throw(ArgumentError("unrecognized SIP parameter $keyword"))
    end
    forward = _fits_sip_record_terms(header, prefix, "FWD")
    reverse = _fits_sip_record_terms(header, prefix, "REV")
    return _FITSTPDAxis(axes, offset, scale, nothing, forward, reverse, true)
end

_fits_sip_record_models(header, count, suffix, minerr) = ntuple(axis -> _fits_sip_record_axis(header, axis, count, suffix, minerr), count)

function _fits_tpd_axis(header, axis, count, suffix, distortion_prefix, parameter_prefix, error_prefix, minerr)
    distortion_keyword = "$distortion_prefix$axis$suffix"
    haskey(header, distortion_keyword) || return nothing
    distortion_type = uppercase(strip(String(header[distortion_keyword])))
    distortion_type in ("LOOKUP", "SIP") && distortion_prefix == "CPDIS" && return nothing
    distortion_type == "LOOKUP" && throw(ArgumentError("sequent LOOKUP distortion $distortion_keyword is not supported"))
    distortion_type == "TPD" || throw(ArgumentError("unsupported Paper IV distortion type $(repr(header[distortion_keyword])) in $distortion_keyword"))

    error_keyword = "$error_prefix$axis$suffix"
    _fits_number(get(header, error_keyword, 0.0), error_keyword) < minerr && return nothing
    prefix = "$parameter_prefix$axis$suffix"
    naxes_keyword = "$prefix.NAXES"
    haskey(header, naxes_keyword) || throw(ArgumentError("TPD distortion $distortion_keyword requires $naxes_keyword"))
    dimensions = _fits_integer(header[naxes_keyword], naxes_keyword; minimum=1)
    dimensions in (1, 2) || throw(ArgumentError("$naxes_keyword must be 1 or 2, got $dimensions"))

    dotted_prefix = "$prefix."
    for keyword in keys(header)
        startswith(keyword, dotted_prefix) || continue
        field = keyword[(length(dotted_prefix) + 1):end]
        _fits_valid_tpd_field(field) || throw(ArgumentError("unrecognized TPD parameter $keyword"))
        mapped = match(r"^(AXIS|OFFSET|SCALE)\.([1-9][0-9]*)$", field)
        mapped !== nothing && parse(Int, mapped.captures[2]) > dimensions && throw(ArgumentError("$keyword exceeds the $dimensions independent axes declared by $naxes_keyword"))
    end

    axes = ntuple(dimensions) do variable
        keyword = "$prefix.AXIS.$variable"
        _fits_integer(get(header, keyword, variable), keyword; minimum=1)
    end
    all(value -> value <= count, axes) || throw(ArgumentError("$prefix axis map contains an axis outside 1:$count"))
    allunique(axes) || throw(ArgumentError("$prefix axis map contains duplicate axes"))
    offset = SVector{dimensions}(ntuple(variable -> _fits_number(get(header, "$prefix.OFFSET.$variable", 0.0), "$prefix.OFFSET.$variable"), dimensions))
    scale = SVector{dimensions}(ntuple(variable -> _fits_number(get(header, "$prefix.SCALE.$variable", 1.0), "$prefix.SCALE.$variable"), dimensions))
    forward = _fits_tpd_coefficients(header, prefix, "FWD")
    reverse = _fits_tpd_coefficients(header, prefix, "REV")
    auxiliary = _fits_tpd_auxiliary(header, prefix)
    docorr_keyword = "$prefix.DOCORR"
    docorr_value = get(header, docorr_keyword, 1)
    docorr_value isa Integer && !(docorr_value isa Bool) || throw(ArgumentError("FITS keyword $docorr_keyword must be an integer"))
    terms = _fits_tpd_compiled_terms(forward, Val(dimensions))
    reverse_terms = _fits_tpd_compiled_terms(reverse, Val(dimensions))
    return _FITSTPDAxis(axes, offset, scale, auxiliary, terms, reverse_terms, !iszero(docorr_value))
end

function _fits_tpd_models(header, count, suffix, distortion_prefix, parameter_prefix, error_prefix, minerr)
    return ntuple(axis -> _fits_tpd_axis(header, axis, count, suffix, distortion_prefix, parameter_prefix, error_prefix, minerr), count)
end

function _fits_tpv_models(header, count, suffix, longitude_axis, latitude_axis)
    xcoeff = _fits_tpv_coefficients(header, longitude_axis, suffix)
    ycoeff = _fits_tpv_coefficients(header, latitude_axis, suffix)
    offset = SVector(0.0, 0.0)
    scale = SVector(1.0, 1.0)
    # Unlike WCSLIB's legacy translation, bind TPV variables to the celestial
    # intermediate axes so unrelated axes cannot enter a mixed-dimensional WCS.
    xmodel = _FITSTPDAxis((longitude_axis, latitude_axis), offset, scale, nothing, _fits_tpd_compiled_terms(xcoeff, Val(2)), false)
    ymodel = _FITSTPDAxis((latitude_axis, longitude_axis), offset, scale, nothing, _fits_tpd_compiled_terms(ycoeff, Val(2)), false)
    return ntuple(axis -> axis == longitude_axis ? xmodel : axis == latitude_axis ? ymodel : nothing, count)
end

# Returns the projection code and whether SCAMP PV distortion terms must be ignored.
function _fits_projection_code(header, parsed, longitude_axis, latitude_axis, suffix)
    code = parsed[longitude_axis][3]
    code == "TAN" || return code, false
    any(haskey(header, "PV$(axis)_$(parameter)$suffix") for axis in (longitude_axis, latitude_axis) for parameter in 5:59) || return code, false
    # Like astropy (Shupe et al. 2012, Sec. 7): -TAN-SIP prefers SIP over SCAMP PV terms; plain -TAN means TPV.
    parsed[longitude_axis][4] == "-SIP" && parsed[latitude_axis][4] == "-SIP" && return code, true
    return "TPV", false
end

function _fits_sip_polynomial(header, prefix, suffix)
    order_keyword = "$(prefix)_ORDER$suffix"
    haskey(header, order_keyword) || return nothing
    order = _fits_integer(header[order_keyword], order_keyword)
    order <= 9 || throw(ArgumentError("FITS SIP $order_keyword must not exceed 9"))
    coefficients = zeros(Float64, order + 1, order + 1)
    for x_power in 0:order, y_power in 0:(order - x_power)
        keyword = "$(prefix)_$(x_power)_$(y_power)$suffix"
        haskey(header, keyword) && (coefficients[x_power + 1, y_power + 1] = _fits_number(header[keyword], keyword))
    end
    return Polynomial2D(coefficients)
end

function _fits_sip_transform(header, suffix, crpix, count, ctype_suffixes, repairs; enabled=true)
    enabled || return IdentityTransform{count}()
    has_sip = any(==("-SIP"), ctype_suffixes) || haskey(header, "A_ORDER$suffix") || haskey(header, "B_ORDER$suffix")
    has_sip || return IdentityTransform{count}()
    count >= 2 || throw(ArgumentError("FITS SIP distortion requires at least two pixel axes"))
    any(==("-SIP"), ctype_suffixes[3:end]) && throw(ArgumentError("FITS SIP must be declared on the first two CTYPE axes"))
    # Like astropy, apply SIP coefficients even when CTYPE omits the -SIP suffix.
    all(==("-SIP"), ctype_suffixes[1:2]) || push!(repairs, "SIP coefficients applied although CTYPE lacks the -SIP suffix")
    a = _fits_sip_polynomial(header, "A", suffix)
    b = _fits_sip_polynomial(header, "B", suffix)
    (a === nothing) == (b === nothing) || throw(ArgumentError("FITS SIP requires both A_ORDER and B_ORDER"))
    a === nothing && throw(ArgumentError("FITS CTYPE declares SIP but A_ORDER and B_ORDER are absent"))
    forward = SIPTransform(a, b)
    ap = _fits_sip_polynomial(header, "AP", suffix)
    bp = _fits_sip_polynomial(header, "BP", suffix)
    (ap === nothing) == (bp === nothing) || throw(ArgumentError("FITS SIP inverse requires both AP_ORDER and BP_ORDER"))
    identity = Matrix{Float64}(I, 2, 2)
    center = AffineTransform(identity, -crpix[1:2])
    restore = AffineTransform(identity, crpix[1:2])
    correction = center | forward | restore
    if ap !== nothing
        reverse = center | SIPTransform(ap, bp) | restore
        correction = _ApproximateInverseTransform(correction, reverse)
    end
    return count == 2 ? correction : parallel(correction, IdentityTransform{count - 2}())
end

_fits_has_matrix(header, prefix, count, suffix) = any(haskey(header, "$prefix$(i)_$(j)$suffix") for i in 1:count for j in 1:count)

# Like WCSLIB, PCi_j (with CDELTi) takes precedence when a nonconforming header also has CDi_j.
_fits_uses_cd(header, count, suffix) = _fits_has_matrix(header, "CD", count, suffix) && !_fits_has_matrix(header, "PC", count, suffix)

function _fits_cd_matrix(header, count, cdelt, suffix, longitude_axis, latitude_axis, repairs)
    has_pc = _fits_has_matrix(header, "PC", count, suffix)
    has_pc && _fits_has_matrix(header, "CD", count, suffix) && push!(repairs, "both PCi_j and CDi_j are present; using PCi_j$suffix and CDELTi$suffix and ignoring CDi_j$suffix")
    matrix = zeros(Float64, count, count)
    if _fits_uses_cd(header, count, suffix)
        for i in 1:count, j in 1:count
            keyword = "CD$(i)_$(j)$suffix"
            matrix[i, j] = haskey(header, keyword) ? _fits_number(header[keyword], keyword) : 0.0
        end
        # WCSLIB's cdfix (applied by astropy): an axis with no CDi_j terms gets a unit diagonal; CDELTi is ignored with CDi_j.
        for i in 1:count
            if all(iszero, view(matrix, i, :)) && all(iszero, view(matrix, :, i))
                matrix[i, i] = 1.0
                push!(repairs, "CD$(i)_$(i)$suffix set to 1 because axis $i has no CDi_j terms (CDELT$i$suffix is ignored with CDi_j)")
            end
        end
        return _fits_nonsingular(matrix, "CDi_j$suffix")
    end
    pc = Matrix{Float64}(I, count, count)
    if has_pc
        for i in 1:count, j in 1:count
            keyword = "PC$(i)_$(j)$suffix"
            haskey(header, keyword) && (pc[i, j] = _fits_number(header[keyword], keyword))
        end
    elseif longitude_axis != 0
        # Like WCSLIB, CROTA applies only to a celestial pair and uses the latitude axis's value.
        i, j = longitude_axis, latitude_axis
        crota = "CROTA$j$suffix"
        longitude_crota = "CROTA$i$suffix"
        !haskey(header, crota) && haskey(header, longitude_crota) && !iszero(_fits_number(header[longitude_crota], longitude_crota)) &&
            push!(repairs, "$longitude_crota on the longitude axis is ignored; only the latitude axis CROTA applies")
        if haskey(header, crota)
            iszero(cdelt[i]) && throw(ArgumentError("FITS CROTA$j$suffix gives a singular transformation matrix because CDELT$i$suffix is zero"))
            angle = deg2rad(_fits_number(header[crota], crota))
            cosine, sine = cos(angle), sin(angle)
            pc[i, i] = pc[j, j] = cosine
            pc[i, j] = -cdelt[j] * sine / cdelt[i]
            pc[j, i] = cdelt[i] * sine / cdelt[j]
        end
    end
    for i in 1:count, j in 1:count
        matrix[i, j] = cdelt[i] * pc[i, j]
    end
    return _fits_nonsingular(matrix, has_pc ? "PCi_j$suffix and CDELTi$suffix" : "CDELTi$suffix")
end

function _fits_nonsingular(matrix, description)
    iszero(det(matrix)) && throw(ArgumentError("FITS linear transformation matrix from $description is singular"))
    return matrix
end

# WCSLIB's wcsutrn safe translations of non-standard unit tokens; astropy applies the same list by default.
# The ambiguous D, H, and S (debye, henry, siemens versus day, hour, second) are deliberately absent.
const _FITS_UNIT_ALIASES = Dict(
    "Angstroms" => "Angstrom", "angstrom" => "Angstrom", "angstroms" => "Angstrom",
    "arcmins" => "arcmin", "ARCMIN" => "arcmin", "ARCMINS" => "arcmin", "arcsecs" => "arcsec", "ARCSEC" => "arcsec", "ARCSECS" => "arcsec",
    "BEAM" => "beam", "Byte" => "byte", "day" => "d", "days" => "d", "DAY" => "d", "DAYS" => "d",
    "degree" => "deg", "degrees" => "deg", "Deg" => "deg", "Degree" => "deg", "Degrees" => "deg", "DEG" => "deg", "DEGREE" => "deg", "DEGREES" => "deg",
    "GHZ" => "GHz", "hr" => "h", "HR" => "h", "hz" => "Hz", "HZ" => "Hz", "KHZ" => "kHz", "JY" => "Jy", "KM" => "km",
    "metre" => "m", "metres" => "m", "meter" => "m", "meters" => "m", "M" => "m", "METRE" => "m", "METRES" => "m", "METER" => "m", "METERS" => "m",
    "MIN" => "min", "MHZ" => "MHz", "Ohm" => "ohm", "pixels" => "pixel", "PIXEL" => "pixel", "PIXELS" => "pixel",
    "radian" => "rad", "radians" => "rad", "RAD" => "rad", "RADIAN" => "rad", "RADIANS" => "rad",
    "sec" => "s", "second" => "s", "seconds" => "s", "SEC" => "s", "SECOND" => "s", "SECONDS" => "s",
    "year" => "yr", "years" => "yr", "YR" => "yr", "YEAR" => "yr", "YEARS" => "yr",
)

# Translate each alphabetic token as WCSLIB's unitfix does and collapse whitespace.
_fits_unit_translation(unit) = join(split(replace(strip(unit), r"[A-Za-z]+" => token -> get(_FITS_UNIT_ALIASES, token, token))), " ")

const _FITS_ANGLE_UNITS = Dict("" => 1.0, "deg" => 1.0, "rad" => 180 / pi, "arcmin" => 1 / 60, "arcsec" => 1 / 3600, "mas" => 1 / 3.6e6)
const _FITS_VELOCITY_UNITS = Dict("m/s" => 1.0, "m s-1" => 1.0, "m.s-1" => 1.0, "km/s" => 1e3, "km s-1" => 1e3, "km.s-1" => 1e3)
const _FITS_SPECTRAL_UNITS = Dict(
    "FREQ" => Dict("Hz" => 1.0, "kHz" => 1e3, "MHz" => 1e6, "GHz" => 1e9, "THz" => 1e12),
    "AFRQ" => Dict("rad/s" => 1.0, "rad s-1" => 1.0, "rad.s-1" => 1.0),
    "WAVE" => Dict("m" => 1.0, "cm" => 1e-2, "mm" => 1e-3, "um" => 1e-6, "nm" => 1e-9, "Angstrom" => 1e-10),
    "ENER" => Dict("J" => 1.0, "erg" => 1e-7, "eV" => 1.602176634e-19, "keV" => 1.602176634e-16, "MeV" => 1.602176634e-13),
    "WAVN" => Dict("1/m" => 1.0, "m-1" => 1.0, "m^-1" => 1.0, "/m" => 1.0, "1/cm" => 100.0, "cm-1" => 100.0, "cm^-1" => 100.0, "/cm" => 100.0),
    "VRAD" => _FITS_VELOCITY_UNITS, "VOPT" => _FITS_VELOCITY_UNITS, "VELO" => _FITS_VELOCITY_UNITS,
)
# "a"/"yr" and "cy" are the Julian year and century of the FITS time standard, as in WCSLIB.
const _FITS_TIME_UNITS = Dict("" => 1.0, "s" => 1.0, "min" => 60.0, "h" => 3600.0, "d" => 86400.0, "a" => 31557600.0, "yr" => 31557600.0,
    "cy" => 3155760000.0, "ms" => 1e-3, "us" => 1e-6, "ns" => 1e-9)

# Units are case-sensitive FITS strings that have already passed through _fits_unit_translation.
function _fits_angle_scale(unit, keyword)
    haskey(_FITS_ANGLE_UNITS, unit) || throw(ArgumentError("unsupported celestial $keyword=$(repr(unit))"))
    return _FITS_ANGLE_UNITS[unit]
end

function _fits_spectral_scale(coordinate, unit, keyword)
    isempty(unit) && return 1.0
    coordinate in ("ZOPT", "BETA") && throw(ArgumentError("spectral $keyword must be empty for dimensionless $coordinate, got $(repr(unit))"))
    factors = get(_FITS_SPECTRAL_UNITS, coordinate == "AWAV" ? "WAVE" : coordinate, nothing)
    factors === nothing && return 1.0
    haskey(factors, unit) || throw(ArgumentError("unsupported spectral $keyword=$(repr(unit))"))
    return factors[unit]
end

function _fits_time_scale(unit, keyword)
    unit in ("S", "H", "D") && throw(ArgumentError("ambiguous time $keyword=$(repr(unit)); FITS reads it as siemens, henry, or debye, so use s, h, or d"))
    haskey(_FITS_TIME_UNITS, unit) || throw(ArgumentError("unsupported time $keyword=$(repr(unit))"))
    return _FITS_TIME_UNITS[unit]
end

function _fits_time_reference(header, suffix)
    value(keyword, default=nothing) = get(header, "$keyword$suffix", get(header, keyword, default))
    reference = value("MJDREF")
    if reference === nothing
        integer = value("MJDREFI")
        fraction = value("MJDREFF")
        if integer !== nothing || fraction !== nothing
            reference = _fits_number(something(integer, 0.0), "MJDREFI$suffix") + _fits_number(something(fraction, 0.0), "MJDREFF$suffix")
        end
    else
        reference = _fits_number(reference, "MJDREF$suffix")
    end
    if reference === nothing
        # JDREF/JDREFI/JDREFF are the Julian-date forms of the same reference epoch.
        jd, integer, fraction = value("JDREF"), value("JDREFI"), value("JDREFF")
        if jd !== nothing
            reference = _fits_number(jd, "JDREF$suffix") - 2400000.5
        elseif integer !== nothing || fraction !== nothing
            reference = (_fits_number(something(integer, 0.0), "JDREFI$suffix") - 2400000.5) + _fits_number(something(fraction, 0.0), "JDREFF$suffix")
        end
    end
    dateref = value("DATEREF")
    if dateref !== nothing
        # As in WCSLIB's datfix, MJDREF or JDREF takes precedence, and DATEREF must agree with it to 0.001 day.
        from_date = _fits_date_mjd(dateref, "DATEREF$suffix")
        reference === nothing && (reference = from_date)
        abs(reference - from_date) <= 0.001 || throw(ArgumentError("FITS DATEREF$suffix=$(repr(dateref)) is inconsistent with the MJDREF or JDREF reference time $reference"))
    end
    # WCSLIB defaults MJDREF to zero when no reference date is given.
    reference === nothing && (reference = 0.0)
    location = if all(keyword -> haskey(header, keyword), ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z"))
        ntuple(i -> _fits_number(header[("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z")[i]], ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z")[i]), 3)
    else
        nothing
    end
    return (reference_epoch=reference, time_scale=String(value("TIMESYS", "")),
        reference_position=String(value("TREFPOS", "")), reference_direction=value("TREFDIR"), location=location)
end

# DATEREF in the FITS standard's ISO 8601 form YYYY-MM-DD[Thh:mm:ss[.s]], as an MJD by WCSLIB's integer algorithm.
# WCSLIB's datfix also accepts legacy forms such as DD/MM/YY; those are rejected rather than guessed.
function _fits_date_mjd(value, keyword)
    text = strip(String(value))
    matched = match(r"^(\d{4})-(\d{2})-(\d{2})(?:T(\d{2}):(\d{2}):(\d{2}(?:\.\d*)?))?$", text)
    matched === nothing && throw(ArgumentError("FITS $keyword=$(repr(text)) is not an ISO 8601 date"))
    year, month, day = parse.(Int, matched.captures[1:3])
    days = (1461 * (year - (12 - month) ÷ 10 + 4712)) ÷ 4 + (306 * ((month + 9) % 12) + 5) ÷ 10 - (3 * ((year - (12 - month) ÷ 10 + 4900) ÷ 100)) ÷ 4 + day - 2399904
    matched.captures[4] === nothing && return Float64(days)
    hour, minute = parse.(Int, matched.captures[4:5])
    return days + (hour + (minute + parse(Float64, matched.captures[6]) / 60) / 60) / 24
end

# A reference frame as asdf-astropy writes it, so Python reads it as the astropy frame.
_fits_astropy_frame(name, attributes...) = _tagged_metadata("tag:astropy.org:astropy/coordinates/frames/$name", v"1.2.0", Dict("frame_attributes" => Dict{Any,Any}(attributes...)))
_fits_epoch(value, format) = _tagged_metadata(_TIME_TAG, v"1.4.0", Dict("value" => value, "format" => format, "scale" => "tt"))

# Like astropy's wcs_to_celestial_frame, RA/DEC take ICRS, FK5, FK4, or FK4-NO-E from RADESYS and EQUINOX, with FK4's
# obstime the equinox, and GLON/GLAT is Galactic. astropy also applies RADESYS to other pairs, so that it labels ecliptic
# coordinates equatorial; here only RA/DEC are. Systems without an astropy frame tag keep their FITS metadata, which has
# no ASDF representation.
function _fits_celestial_reference(longitude, latitude, radesys, equinox)
    if longitude == "RA"
        radesys == "ICRS" && return :ICRS
        radesys == "FK5" && return _fits_astropy_frame("fk5", "equinox" => _fits_epoch(equinox, "jyear"))
        if radesys in ("FK4", "FK4-NO-E")
            epoch = _fits_epoch(equinox, "byear")
            return _fits_astropy_frame(radesys == "FK4" ? "fk4" : "fk4noeterms", "equinox" => epoch, "obstime" => epoch)
        end
    elseif longitude == "GLON"
        return _fits_astropy_frame("galactic")
    end
    return (system=(longitude, latitude), radesys, equinox)
end

# FITS spectral reference metadata, which the ASDF spectral frame cannot represent.
function _fits_spectral_reference(header, suffix, restfrq, restwav)
    text(keyword) = String(get(header, "$keyword$suffix", ""))
    number(keyword) = haskey(header, "$keyword$suffix") ? _fits_number(header["$keyword$suffix"], "$keyword$suffix") : nothing
    reference = (restfrq, restwav, specsys=text("SPECSYS"), ssysobs=text("SSYSOBS"), velosys=number("VELOSYS"), zsource=number("ZSOURCE"), ssyssrc=text("SSYSSRC"))
    return all(value -> value === nothing || value == "", reference) ? nothing : reference
end

const _FITS_DEPRECATED_TIME_SCALES = Dict("TDT" => "tt", "ET" => "tt", "GMT" => "utc", "UT" => "utc", "IAT" => "tai")
const _ASDF_TIME_SCALES = ("utc", "tai", "tcb", "tcg", "tdb", "tt", "ut1")

# The reference time of a FITS time axis as the ASDF time a TemporalFrame holds, following astropy's FITS WCS: the scale
# is the CTYPE, or TIMESYS for TIME; GPS is TAI shifted by 19 s; TREFPOS TOPOCENTER takes OBSGEO and GEOCENTER the
# geocenter. Metadata the ASDF time cannot hold, such as a sub-scale, another TREFPOS, or TREFDIR, is kept in
# `reference_frame` and has no ASDF representation.
function _fits_temporal_frame(header, suffix, coordinate; kwargs...)
    reference = _fits_time_reference(header, suffix)
    written = uppercase(coordinate == "TIME" ? (isempty(reference.time_scale) ? "UTC" : reference.time_scale) : coordinate)
    scale, subscale = occursin('(', written) ? (written[1:(findfirst('(', written) - 1)], written[(findfirst('(', written) + 1):(end - 1)]) : (written, "")
    epoch = reference.reference_epoch
    if scale == "GPS"
        scale, epoch = "TAI", epoch + 19 / 86400
    end
    scale = get(_FITS_DEPRECATED_TIME_SCALES, scale, lowercase(scale))
    position = uppercase(reference.reference_position)
    location = startswith(position, "TOPOCENT") ? reference.location : position == "GEOCENTER" ? (0.0, 0.0, 0.0) : nothing
    leftover = merge(scale in _ASDF_TIME_SCALES ? (;) : (time_scale=scale,), isempty(subscale) ? (;) : (; subscale),
        isempty(position) || location !== nothing ? (;) : (reference_position=reference.reference_position,),
        reference.reference_direction === nothing ? (;) : (reference_direction=reference.reference_direction,))
    return TemporalFrame(; kwargs..., reference_frame=isempty(leftover) ? nothing : leftover, reference_epoch=epoch, time_scale=scale, time_format="mjd", location)
end

# A CelestialFrame for the celestial pair, a SpectralFrame, TemporalFrame, or StokesFrame for each such axis, and one
# CoordinateFrame for the linear axes, in a CompositeFrame when there is more than one. axes_order places each frame's
# axes at their header positions, so the world coordinates stay in header order.
function _fits_world_frame(header, suffix, parsed, names, descriptions, units, longitude_axis, latitude_axis, restfrq, restwav)
    kinds = map(item -> item[2], parsed)
    linear = Tuple(axis for axis in eachindex(kinds) if kinds[axis] == :linear)
    single = (longitude_axis != 0) + count(in((:spectral, :time, :stokes)), kinds) + !isempty(linear) == 1
    name(kind) = single ? "world" : kind
    frames = AbstractCoordinateFrame[]
    if longitude_axis != 0
        radesys, equinox = _fits_reference_system(header, suffix, parsed)
        axes = (longitude_axis, latitude_axis)
        reference = _fits_celestial_reference(parsed[longitude_axis][1], parsed[latitude_axis][1], radesys, equinox)
        push!(frames, CelestialFrame{2}(; name=name("celestial"), axes_names=map(i -> names[i], axes), axes_order=axes,
            axis_physical_types=map(i -> descriptions[i][2], axes), unit=map(i -> units[i], axes), reference_frame=reference))
    end
    for (axis, kind) in enumerate(kinds)
        properties = (axes_names=(names[axis],), axes_order=(axis,), axis_physical_types=(descriptions[axis][2],))
        if kind == :spectral
            push!(frames, SpectralFrame(; name=name("spectral"), properties..., unit=(units[axis],), reference_frame=_fits_spectral_reference(header, suffix, restfrq, restwav)))
        elseif kind == :time
            push!(frames, _fits_temporal_frame(header, suffix, parsed[axis][1]; name=name("temporal"), properties..., unit=(units[axis],)))
        elseif kind == :stokes
            push!(frames, StokesFrame(; name=name("stokes"), properties...))
        end
    end
    if !isempty(linear)
        push!(frames, CoordinateFrame{length(linear)}(; name=name("linear"), axes_names=map(i -> names[i], linear), axes_order=linear,
            axis_types=map(i -> descriptions[i][1], linear), axis_physical_types=map(i -> descriptions[i][2], linear), unit=map(i -> units[i], linear)))
    end
    return single ? only(frames) : CompositeFrame(Tuple(frames); name="world")
end

function _fits_reference_system(header, suffix, parsed)
    applies = any(item -> item[1] in ("RA", "ELON", "HLON"), parsed)
    applies || return "", nothing
    equinox_value = get(header, "EQUINOX$suffix", nothing)
    equinox = equinox_value === nothing ? nothing : _fits_number(equinox_value, "EQUINOX$suffix")
    radesys = uppercase(strip(String(get(header, "RADESYS$suffix", ""))))
    if isempty(radesys)
        radesys = equinox === nothing ? "ICRS" : equinox < 1984 ? "FK4" : "FK5"
    elseif radesys in ("ICRS", "GAPPT")
        equinox = nothing
    elseif equinox === nothing
        radesys == "FK5" && (equinox = 2000.0)
        radesys in ("FK4", "FK4-NO-E") && (equinox = 1950.0)
    end
    return radesys, equinox
end

function _fits_axis_scale(parsed_ctype, unit, keyword)
    coordinate, kind, _, _ = parsed_ctype
    kind in (:longitude, :latitude) && return _fits_angle_scale(unit, keyword)
    kind == :spectral && return _fits_spectral_scale(coordinate, unit, keyword)
    kind == :time && return _fits_time_scale(unit, keyword)
    return 1.0
end

function _fits_axis_metadata(parsed_ctype, unit, index)
    coordinate, kind, _, _ = parsed_ctype
    if kind == :longitude
        physical = coordinate == "RA" ? "pos.eq.ra" : coordinate == "GLON" ? "pos.galactic.lon" : "custom:$(lowercase(coordinate))"
        return SPATIAL, physical, "deg"
    elseif kind == :latitude
        physical = coordinate == "DEC" ? "pos.eq.dec" : coordinate == "GLAT" ? "pos.galactic.lat" : "custom:$(lowercase(coordinate))"
        return SPATIAL, physical, "deg"
    elseif kind == :spectral
        canonical = coordinate == "FREQ" ? "Hz" : coordinate == "AFRQ" ? "rad/s" : coordinate in ("WAVE", "AWAV") ? "m" : coordinate == "ENER" ? "J" : coordinate == "WAVN" ? "m-1" : coordinate in ("ZOPT", "BETA") ? "" : "m/s"
        physical = coordinate == "FREQ" ? "em.freq" : coordinate == "ENER" ? "em.energy" : coordinate == "WAVN" ? "em.wavenumber" :
            coordinate == "WAVE" ? "em.wl" : coordinate == "AWAV" ? "em.wl;obs.atmos" : coordinate == "VRAD" ? "spect.dopplerVeloc.radio" :
            coordinate == "VOPT" ? "spect.dopplerVeloc.opt" : coordinate == "ZOPT" ? "src.redshift" : coordinate == "VELO" ? "spect.dopplerVeloc" :
            coordinate == "BETA" ? "custom:spect.dopplerVeloc.beta" : "custom:spectral"
        return SPECTRAL, physical, canonical
    elseif kind == :time
        return TIME, "time", "s"
    elseif kind == :stokes
        return STOKES, "phys.polarization.stokes", ""
    end
    name = isempty(coordinate) ? "axis$index" : lowercase(coordinate)
    return UNDEFINED, "custom:$name", unit
end

function _fits_selected_transform(count, first_axis, second_axis, selected)
    remaining = Tuple(filter(axis -> axis != first_axis && axis != second_axis, 1:count))
    order = (first_axis, second_axis, remaining...)
    arrange = RemapAxes{count}(order)
    restore = inverse_transform(arrange)
    count == 2 && return arrange | selected | restore
    return arrange | parallel(selected, IdentityTransform{count - 2}()) | restore
end

function _fits_selected_axis_transform(count, axis, selected)
    remaining = Tuple(filter(index -> index != axis, 1:count))
    arrange = RemapAxes{count}((axis, remaining...))
    restore = inverse_transform(arrange)
    count == 1 && return selected
    return arrange | parallel(selected, IdentityTransform{count - 1}()) | restore
end

const _FITS_SPECTRAL_PARENTS = Dict(
    "FREQ" => :F, "AFRQ" => :F, "ENER" => :F, "WAVN" => :F, "VRAD" => :F,
    "WAVE" => :W, "VOPT" => :W, "ZOPT" => :W, "AWAV" => :A, "VELO" => :V, "BETA" => :V,
)
const _FITS_SPECTRAL_ALGORITHMS = Dict(
    "F2W" => (:F, :W), "W2F" => (:W, :F), "F2V" => (:F, :V), "V2F" => (:V, :F),
    "W2V" => (:W, :V), "V2W" => (:V, :W), "F2A" => (:F, :A), "A2F" => (:A, :F),
    "W2A" => (:W, :A), "A2W" => (:A, :W), "V2A" => (:V, :A), "A2V" => (:A, :V),
)

function _fits_spectral_transform(parsed_ctype, reference, rest_frequency, rest_wavelength, output_scale)
    coordinate, _, algorithm, _ = parsed_ctype
    parent = _FITS_SPECTRAL_PARENTS[coordinate]
    haskey(_FITS_SPECTRAL_ALGORITHMS, algorithm) || throw(ArgumentError("unsupported FITS nonlinear spectral algorithm $(repr(algorithm))"))
    x_type, algorithm_parent = _FITS_SPECTRAL_ALGORITHMS[algorithm]
    algorithm_parent == parent || throw(ArgumentError("FITS spectral coordinate $coordinate is incompatible with algorithm $algorithm"))
    !isfinite(rest_frequency) && isfinite(rest_wavelength) && (rest_frequency = _FITS_SPEED_OF_LIGHT / rest_wavelength)
    !isfinite(rest_wavelength) && isfinite(rest_frequency) && (rest_wavelength = _FITS_SPEED_OF_LIGHT / rest_frequency)
    requires_rest = coordinate in ("VRAD", "VOPT", "ZOPT", "VELO", "BETA") || x_type == :V || parent == :V
    requires_rest && !isfinite(rest_frequency) && throw(ArgumentError("FITS spectral algorithm $algorithm requires RESTFRQ or RESTWAV"))
    transform = FITSSpectralTransform(Val(x_type), Val(parent), Val(Symbol(coordinate)), Val(Symbol(algorithm)), reference, rest_frequency, rest_wavelength; output_scale)
    derivative = ForwardDiff.derivative(reference) do value
        parent_value = _fits_spectral_to_parent(value, Val(Symbol(coordinate)), rest_frequency, rest_wavelength)
        _fits_parent_to_x(parent_value, Val(x_type), Val(parent), rest_frequency, rest_wavelength)
    end
    isfinite(derivative) && !iszero(derivative) || throw(ArgumentError("FITS spectral algorithm $algorithm has a singular reference value"))
    return transform, derivative
end

function _fits_grism_transform(header, suffix, axis, parsed_ctype, reference, rest_frequency, rest_wavelength, output_scale)
    coordinate, _, algorithm, _ = parsed_ctype
    parent = _FITS_SPECTRAL_PARENTS[coordinate]
    x_type = algorithm == "GRI" ? :W : :A
    requires_rest = coordinate in ("VRAD", "VOPT", "ZOPT", "VELO", "BETA")
    requires_rest && !isfinite(rest_frequency) && !isfinite(rest_wavelength) && throw(ArgumentError("FITS spectral algorithm $algorithm requires RESTFRQ or RESTWAV"))
    defaults = (0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0)
    parameters = ntuple(7) do position
        keyword = "PV$(axis)_$(position - 1)$suffix"
        _fits_number(get(header, keyword, defaults[position]), keyword)
    end
    transform = FITSGrismTransform(Val(x_type), Val(parent), Val(Symbol(coordinate)), Val(Symbol(algorithm)),
        reference, rest_frequency, rest_wavelength, parameters; output_scale)
    return transform, transform.linear_scale
end

function _fits_world_transform(count, crval, longitude_axis, latitude_axis, projection, lonpole, latpole, native_reference, celestial_offset, output_scaling)
    offset = copy(crval)
    if projection !== nothing
        offset[longitude_axis] = 0.0
        offset[latitude_axis] = 0.0
    end
    result = AffineTransform(Matrix{Float64}(I, count, count), offset)
    if projection !== nothing
        celestial = celestial_transform(projection, (crval[longitude_axis], crval[latitude_axis]); lonpole, latpole, native_reference, offset=celestial_offset)
        result = result | _fits_selected_transform(count, longitude_axis, latitude_axis, celestial)
    end
    all(isone, output_scaling) || (result = result | AffineTransform(Diagonal(output_scaling)))
    return result
end

function _fits_pixel_shape(header, count)
    haskey(header, "NAXIS") || return nothing
    image_axes = _fits_integer(header["NAXIS"], "NAXIS")
    return ntuple(count) do axis
        keyword = "NAXIS$axis"
        axis <= image_axes && haskey(header, keyword) ? _fits_integer(header[keyword], keyword) : nothing
    end
end

function _parse_fits_wcs(header::Dict{String,Any}; key::Char=' ', preserve_units::Bool=false, minerr::Real=0.0, legacy::Symbol=:warn, lookup_data=nothing, tabular_data=())
    isfinite(minerr) && minerr >= 0 || throw(ArgumentError("minerr must be finite and nonnegative"))
    legacy in _FITS_LEGACY_MODES || throw(ArgumentError("legacy must be :warn, :quiet, or :error, got $(repr(legacy))"))
    suffix = _fits_key(key)
    repairs = String[]
    legacy_header = _fits_normalize_legacy(header, suffix)
    legacy_header === header || (header = _fits_header(legacy_header))
    header = _fits_normalize_keywords(header, key, repairs)
    _fits_reject_external_data(header, suffix; allow_lookup=lookup_data !== nothing, allow_tab=!isempty(tabular_data))
    count = _fits_axis_count(header, key)
    header = _fits_normalize_aips(header, key, count, repairs)
    prior_tpd_models = _fits_tpd_models(header, count, suffix, "CPDIS", "DP", "CPERR", minerr)
    prior_sip_models = _fits_sip_record_models(header, count, suffix, minerr)
    prior_polynomial_models = ntuple(axis -> prior_sip_models[axis] === nothing ? prior_tpd_models[axis] : prior_sip_models[axis], count)
    sequent_tpd_models = _fits_tpd_models(header, count, suffix, "CQDIS", "DQ", "CQERR", minerr)
    keywords(prefix, default) = [get(header, "$prefix$i$suffix", default) for i in 1:count]
    ctype = String.(keywords("CTYPE", ""))
    parsed = _fits_ctype.(ctype)
    for (axis, (_, kind, algorithm, extra)) in enumerate(parsed)
        isempty(extra) || extra == "-SIP" || throw(ArgumentError("FITS CTYPE$axis=$(repr(ctype[axis])) uses unsupported distortion suffix $(repr(extra))"))
        if kind == :spectral && !isempty(algorithm) && algorithm != "TAB" && algorithm != "LOG" && algorithm ∉ ("GRI", "GRA") && !haskey(_FITS_SPECTRAL_ALGORITHMS, algorithm)
            throw(ArgumentError("FITS CTYPE$axis=$(repr(ctype[axis])) uses unsupported nonlinear algorithm $(repr(algorithm))"))
        elseif kind in (:time, :linear) && !isempty(algorithm) && algorithm ∉ ("TAB", "LOG")
            throw(ArgumentError("FITS CTYPE$axis=$(repr(ctype[axis])) uses unsupported nonlinear algorithm $(repr(algorithm))"))
        end
    end
    longitude_axes = findall(value -> value[2] == :longitude, parsed)
    latitude_axes = findall(value -> value[2] == :latitude, parsed)
    length(longitude_axes) <= 1 || throw(ArgumentError("FITS header has multiple celestial longitude axes"))
    length(latitude_axes) <= 1 || throw(ArgumentError("FITS header has multiple celestial latitude axes"))
    isempty(longitude_axes) == isempty(latitude_axes) || throw(ArgumentError("FITS celestial longitude and latitude axes must be paired"))
    longitude_axis = isempty(longitude_axes) ? 0 : only(longitude_axes)
    latitude_axis = isempty(latitude_axes) ? 0 : only(latitude_axes)
    if longitude_axis != 0
        longitude_type = parsed[longitude_axis][1]
        latitude_type = parsed[latitude_axis][1]
        latitude_type == _fits_celestial_partner(longitude_type) || throw(ArgumentError("FITS celestial axis types $(repr(longitude_type)) and $(repr(latitude_type)) are not a matching longitude/latitude pair"))
    end

    crpix = [_fits_number(value, "CRPIX$i$suffix") for (i, value) in enumerate(keywords("CRPIX", 0.0))]
    raw_crval = [_fits_number(value, "CRVAL$i$suffix") for (i, value) in enumerate(keywords("CRVAL", 0.0))]
    cdelt = [_fits_number(value, "CDELT$i$suffix") for (i, value) in enumerate(keywords("CDELT", 1.0))]
    cunit = String.(keywords("CUNIT", ""))
    timeunit = String(get(header, "TIMEUNIT$suffix", get(header, "TIMEUNIT", "")))
    unit_keywords = ["CUNIT$i$suffix" for i in 1:count]
    for i in eachindex(cunit)
        parsed[i][2] == :time && isempty(strip(cunit[i])) && (cunit[i] = timeunit)
        translated = _fits_unit_translation(cunit[i])
        if translated != join(split(cunit[i]), " ")
            push!(repairs, "CUNIT$i$suffix=$(repr(cunit[i])) read as $(repr(translated))")
            unit_keywords[i] = "CUNIT$i$suffix (written $(repr(cunit[i])))"
        end
        cunit[i] = translated
    end
    scaling = [_fits_axis_scale(parsed[i], cunit[i], unit_keywords[i]) for i in 1:count]
    restfrq_value = get(header, "RESTFRQ$suffix", get(header, "RESTFREQ$suffix", nothing))
    restwav_value = get(header, "RESTWAV$suffix", nothing)
    restfrq = restfrq_value === nothing ? nothing : _fits_number(restfrq_value, "RESTFRQ$suffix")
    restwav = restwav_value === nothing ? nothing : _fits_number(restwav_value, "RESTWAV$suffix")
    spectral_transforms = Any[nothing for _ in 1:count]
    spectral_derivatives = ones(Float64, count)
    for axis in 1:count
        algorithm = parsed[axis][3]
        algorithm == "LOG" && begin
            output_scale = preserve_units ? scaling[axis] : 1.0
            spectral_transforms[axis] = FITSLogTransform(raw_crval[axis] * scaling[axis]; output_scale)
            continue
        end
        parsed[axis][2] == :spectral || continue
        isempty(algorithm) || algorithm == "TAB" || begin
            output_scale = preserve_units ? scaling[axis] : 1.0
            if algorithm in ("GRI", "GRA")
                transform, derivative = _fits_grism_transform(header, suffix, axis, parsed[axis], raw_crval[axis] * scaling[axis], something(restfrq, NaN), something(restwav, NaN), output_scale)
            else
                transform, derivative = _fits_spectral_transform(parsed[axis], raw_crval[axis] * scaling[axis], something(restfrq, NaN), something(restwav, NaN), output_scale)
            end
            spectral_transforms[axis] = transform
            spectral_derivatives[axis] = derivative
        end
    end
    linear_scaling = [parsed[i][3] == "TAB" ? 1.0 : scaling[i] * spectral_derivatives[i] for i in 1:count]
    crval = raw_crval .* linear_scaling
    for axis in eachindex(spectral_transforms)
        spectral_transforms[axis] === nothing || (crval[axis] = 0.0)
    end
    raw_cd = _fits_cd_matrix(header, count, cdelt, suffix, longitude_axis, latitude_axis, repairs)
    cd = copy(raw_cd)
    for i in 1:count
        cd[i, :] .*= linear_scaling[i]
    end

    projection = nothing
    projection_code = ""
    lonpole = nothing
    native_reference = nothing
    celestial_offset = false
    latpole = _fits_number(get(header, "LATPOLE$suffix", 90.0), "LATPOLE$suffix")
    if longitude_axis != 0
        lon_code, ignore_scamp_pv = _fits_projection_code(header, parsed, longitude_axis, latitude_axis, suffix)
        ignore_scamp_pv && push!(repairs, "Ignoring SCAMP PVi_m distortion keywords on the celestial axes because CTYPE declares SIP")
        lat_code = parsed[latitude_axis][3]
        lon_code == "TPV" && lat_code == "TAN" && (lat_code = lon_code)
        lon_code == lat_code || throw(ArgumentError("FITS celestial axes use different projection codes $(repr(lon_code)) and $(repr(lat_code))"))
        # Like WCSLIB, celestial axes without a projection code are evaluated linearly.
        isempty(lon_code) && push!(repairs, "celestial CTYPE$longitude_axis$suffix and CTYPE$latitude_axis$suffix have no projection code and are evaluated linearly")
        projection_code = lon_code
        if !isempty(projection_code) && projection_code != "TAB"
            projection = _fits_projection(projection_code, header, latitude_axis, suffix)
            default_native = native_fiducial(projection)
            use_reference_pv = projection_code != "TPV" && !ignore_scamp_pv
            celestial_offset = use_reference_pv && !iszero(_fits_number(get(header, "PV$(longitude_axis)_0$suffix", 0.0), "PV$(longitude_axis)_0$suffix"))
            native_reference = use_reference_pv ? (
                _fits_number(get(header, "PV$(longitude_axis)_1$suffix", default_native[1]), "PV$(longitude_axis)_1$suffix"),
                _fits_number(get(header, "PV$(longitude_axis)_2$suffix", default_native[2]), "PV$(longitude_axis)_2$suffix")
            ) : default_native
            default_pole = default_lonpole(crval[latitude_axis], native_reference[2], native_reference[1])
            lonpole_keyword = use_reference_pv && haskey(header, "PV$(longitude_axis)_3$suffix") ? "PV$(longitude_axis)_3$suffix" : "LONPOLE$suffix"
            latpole_keyword = use_reference_pv && haskey(header, "PV$(longitude_axis)_4$suffix") ? "PV$(longitude_axis)_4$suffix" : "LATPOLE$suffix"
            lonpole = _fits_number(get(header, lonpole_keyword, default_pole), lonpole_keyword)
            latpole = _fits_number(get(header, latpole_keyword, latpole), latpole_keyword)
        end
    end

    detector = CoordinateFrame{count}(; name="detector", axes_names=ntuple(i -> "pixel$i", count), axis_types=ntuple(_ -> PIXEL, count), unit=ntuple(_ -> "pix", count))
    intermediate = CoordinateFrame{count}(; name="intermediate", axes_names=ntuple(i -> "intermediate$i", count), axis_types=ntuple(_ -> SPATIAL, count), unit=Tuple(map(item -> item[3], _fits_axis_metadata.(parsed, cunit, 1:count))))
    descriptions = [_fits_axis_metadata(parsed[i], cunit[i], i) for i in 1:count]
    names = ntuple(i -> String(get(header, "CNAME$i$suffix", isempty(parsed[i][1]) ? "world$i" : parsed[i][1])), count)
    units = ntuple(i -> preserve_units && !isempty(strip(cunit[i])) ? cunit[i] : descriptions[i][3], count)
    world = _fits_world_frame(header, suffix, parsed, names, descriptions, units, longitude_axis, latitude_axis, restfrq, restwav)
    tpv_models = projection_code == "TPV" ? _fits_tpv_models(header, count, suffix, longitude_axis, latitude_axis) : ntuple(_ -> nothing, count)
    projection_code == "TPV" && any(model -> model !== nothing, sequent_tpd_models) && throw(ArgumentError("TPV cannot be combined with a CQDIS sequent distortion"))
    sequent_models = projection_code == "TPV" ? tpv_models : sequent_tpd_models
    has_sequent = any(model -> model !== nothing, sequent_models)
    if has_sequent
        if _fits_uses_cd(header, count, suffix)
            pre_matrix = raw_cd
            post_scaling = linear_scaling
        else
            all(!iszero, cdelt) || throw(ArgumentError("FITS CDELT values must be nonzero with a sequent distortion"))
            pre_matrix = Diagonal(inv.(cdelt)) * raw_cd
            post_scaling = cdelt .* linear_scaling
        end
        sequent = _FITSTPDTransform(sequent_models)
        pixel_to_intermediate = AffineTransform(pre_matrix, -(pre_matrix * crpix)) | sequent | AffineTransform(Diagonal(post_scaling))
    else
        pixel_to_intermediate = AffineTransform(cd, -(cd * crpix))
    end
    sip = _fits_sip_transform(header, suffix, crpix, count, getindex.(parsed, 4), repairs; enabled=projection_code != "TPV")
    has_lookup = lookup_data !== nothing && any(table -> table !== nothing, (lookup_data.detector_to_image..., lookup_data.cpdis...))
    has_prior_polynomial = any(model -> model !== nothing, prior_polynomial_models)
    if has_prior_polynomial
        detector_to_image = lookup_data === nothing ? (nothing, nothing) : lookup_data.detector_to_image
        cpdis = lookup_data === nothing ? (nothing, nothing) : lookup_data.cpdis
        distortion = _FITSPriorDistortion(Val(count), detector_to_image, sip, cpdis, _FITSTPDTransform(prior_polynomial_models))
        pixel_to_intermediate = distortion | pixel_to_intermediate
    elseif !has_lookup
        pixel_to_intermediate = sip | pixel_to_intermediate
    else
        distortion = FITSLookupDistortion(Val(count), lookup_data.detector_to_image, sip, lookup_data.cpdis)
        pixel_to_intermediate = distortion | pixel_to_intermediate
    end
    output_scaling = [parsed[i][3] == "TAB" || spectral_transforms[i] !== nothing ? 1.0 : preserve_units ? inv(scaling[i]) : 1.0 for i in 1:count]
    intermediate_to_world = _fits_world_transform(count, crval, longitude_axis, latitude_axis, projection, lonpole, latpole, native_reference, celestial_offset, output_scaling)
    for axis in eachindex(spectral_transforms)
        transform = spectral_transforms[axis]
        transform === nothing || (intermediate_to_world = intermediate_to_world | _fits_selected_axis_transform(count, axis, transform))
    end
    if !isempty(tabular_data)
        tabular_scaling = ntuple(i -> parsed[i][3] == "TAB" && !preserve_units ? scaling[i] : 1.0, count)
        intermediate_to_world = intermediate_to_world | FITSTabularTransform(Val(count), tabular_data, tabular_scaling)
    end
    name = String(get(header, "WCSNAME$suffix", ""))
    _fits_report_repairs(repairs, key, legacy)
    # GWCS pixel coordinates are zero-based, like Python GWCS and astropy's APE 14 API, while FITS numbers the first
    # pixel center 1. Shifting at the boundary keeps CRPIX, SIP, and the lookup and prior distortions in FITS pixels.
    fits_pixels = AffineTransform(Matrix{Float64}(I, count, count), ones(count))
    return WCS(Step(detector, fits_pixels | pixel_to_intermediate), Step(intermediate, intermediate_to_world), Step(world); name, pixel_shape=_fits_pixel_shape(header, count))
end

"""
    parse_fits_wcs(header::AbstractDict; key=' ', preserve_units=false, minerr=0.0, legacy=:warn) -> WCS
    parse_fits_wcs(header::FITSIO.FITSHeader; kwargs...) -> WCS
    parse_fits_wcs(hdu::FITSIO.HDU; kwargs...) -> WCS
    parse_fits_wcs(cards::FITSFiles.Cards; kwargs...) -> WCS
    parse_fits_wcs(hdu::FITSFiles.HDU; kwargs...) -> WCS
    parse_fits_wcs(hdus::AbstractVector{<:FITSFiles.HDU}, index; kwargs...) -> WCS

Build a [`WCS`](@ref) from a FITS image header. `key=' '` selects the primary WCS, and a letter selects that alternate
WCS. A header that is not a valid WCS throws an `ArgumentError`.

The header is a dictionary of keywords to values, or, with FITSIO.jl or FITSFiles.jl loaded, one of the objects above.
To read the lookup-table and `-TAB` arrays a header refers to, pass a FITSIO.jl HDU, or FITSFiles.jl HDUs with the
`index` of the HDU holding the WCS.

Pixel coordinates are zero-based, as everywhere in GWCS.jl, so FITS pixel `p` is GWCS pixel `p .- 1`; `CRPIX` needs no
adjustment. World coordinates keep the header's axis order. They are in degrees for celestial axes, SI units for
spectral axes, and seconds from the reference epoch for time axes, unless `preserve_units=true` keeps the header's
`CUNIT` values.

A distortion whose `D2IMERR` or `CPERR` error estimate is below `minerr` is dropped, as astropy does. Legacy conventions
that WCSLIB and astropy accept, such as `NCP`, AIPS spectral types, `EPOCH`, and zero-padded keyword indices, are
repaired with one warning per parse when `legacy=:warn`, silently when `legacy=:quiet`, and rejected when
`legacy=:error`.

Python equivalent: `astropy.wcs.WCS`, which evaluates the header with WCSLIB rather than building a GWCS pipeline.
"""
parse_fits_wcs(source::AbstractDict; key::Char=' ', preserve_units::Bool=false, minerr::Real=0.0, legacy::Symbol=:warn) = _parse_fits_wcs(_fits_header(source); key, preserve_units, minerr, legacy)

function _fits_alternate_keys(header)
    keys_found = Set{Char}()
    (haskey(header, "WCSAXES") || any(keyword -> !isempty(_fits_axis_indices(keyword, "")), keys(header))) && push!(keys_found, ' ')
    for keyword in keys(header), key in 'A':'Z'
        (keyword == "WCSAXES$key" || !isempty(_fits_axis_indices(keyword, string(key)))) && push!(keys_found, key)
    end
    isempty(keys_found) && haskey(header, "NAXIS") && push!(keys_found, ' ')
    return sort!(collect(keys_found))
end

function _fits_malformed_error(error, key)
    error isa Union{ArgumentError,DimensionMismatch,DomainError,InexactError} || throw(error)
    return ArgumentError("FITS WCS key $(repr(key)) is malformed: $(sprint(showerror, error))")
end

function _find_all_fits_wcs(parse, header)
    result = Dict{Char,WCS}()
    for key in _fits_alternate_keys(header)
        try
            result[key] = parse(key)
        catch error
            error isa InterruptException && rethrow()
            throw(_fits_malformed_error(error, key))
        end
    end
    return result
end

"""
    find_all_fits_wcs(header; preserve_units=false, minerr=0.0, legacy=:warn) -> Dict{Char,WCS}

Build a [`WCS`](@ref) from every primary and alternate WCS representation in a FITS header, keyed by its letter
(`' '` for the primary). The header forms and keywords are those of [`parse_fits_wcs`](@ref). A malformed
representation throws an `ArgumentError` naming its key.

Python equivalent: `astropy.wcs.find_all_wcs`.
"""
function find_all_fits_wcs(source::AbstractDict; preserve_units::Bool=false, minerr::Real=0.0, legacy::Symbol=:warn)
    header = _fits_header(source)
    return _find_all_fits_wcs(key -> _parse_fits_wcs(header; key, preserve_units, minerr, legacy), header)
end

function _fits_record_card(value)
    value isa AbstractString || return nothing
    parts = split(String(value), ':'; limit=2)
    length(parts) == 2 || return nothing
    name, raw = uppercase(strip(parts[1])), strip(parts[2])
    integer = tryparse(Int, raw)
    integer === nothing || return name, integer
    parsed = tryparse(Float64, raw)
    return name, parsed === nothing ? raw : parsed
end

# Container adapters pass `(keyword, value)` cards; Paper IV record-valued cards such as
# `DP1 = 'AXIS.1: 1'` repeat a keyword, so expand them to `DP1.AXIS.1` before they collide.
function _fits_card_dictionary(cards)
    result = Dict{String,Any}()
    for (raw_key, value) in cards
        key = uppercase(strip(String(raw_key)))
        result[key] = value
        parsed = _fits_record_card(value)
        parsed === nothing && continue
        name, parsed_value = parsed
        if name in ("EXTVER", "NAXES", "DOCORR") || any(prefix -> startswith(name, prefix), ("AXIS.", "OFFSET.", "SCALE.", "TPD.", "SIP.", "AUX."))
            result["$key.$name"] = parsed_value
        end
    end
    return result
end

"""Return the unique `item` from `item => header` candidates whose header matches `EXTNAME`, `EXTVER`, and `EXTLEVEL`."""
function _fits_find_extension(candidates, extname, extver, extlev=1)
    target = uppercase(strip(String(extname)))
    matches = [item for (item, header) in candidates if haskey(header, "EXTNAME") && uppercase(strip(String(header["EXTNAME"]))) == target &&
        Int(get(header, "EXTVER", 1)) == extver && Int(get(header, "EXTLEVEL", 1)) == extlev]
    isempty(matches) && throw(ArgumentError("FITS WCS extension $(repr(extname)) with EXTVER=$extver and EXTLEVEL=$extlev was not found"))
    length(matches) == 1 || throw(ArgumentError("FITS WCS extension $(repr(extname)) with EXTVER=$extver and EXTLEVEL=$extlev is ambiguous ($(length(matches)) matching HDUs)"))
    return only(matches)
end

# Construction-time loaders supplied by container extensions: `image(extname, extver)` returns
# `(data, header)` and `column(extname, extver, extlev, name)` returns one table column.
function _parse_linked_fits_wcs(image, column, header; key::Char=' ', preserve_units::Bool=false, minerr::Real=0.0, legacy::Symbol=:warn)
    lookup = _fits_lookup_data(header, key, minerr) do spec
        data, table_header = image(spec.extname, spec.extver)
        data isa AbstractArray || throw(ArgumentError("FITS WCS lookup extension $(repr(spec.extname)) must contain an image"))
        _fits_lookup_table(data, _fits_header(table_header), spec.axes, spec.naxis)
    end
    tabular = _fits_tabular_data(column, header, key)
    return _parse_fits_wcs(header; key, preserve_units, minerr, legacy, lookup_data=lookup, tabular_data=tabular)
end
