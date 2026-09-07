const _FITS_LONGITUDE_TYPES = ("RA", "GLON", "ELON", "HLON", "SLON", "HPLN")
const _FITS_LATITUDE_TYPES = ("DEC", "GLAT", "ELAT", "HLAT", "SLAT", "HPLT")
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
    kind = coordinate in _FITS_LONGITUDE_TYPES ? :longitude :
        coordinate in _FITS_LATITUDE_TYPES ? :latitude :
        coordinate in _FITS_SPECTRAL_TYPES ? :spectral :
        coordinate in _FITS_TIME_TYPES ? :time :
        coordinate == "STOKES" ? :stokes : :linear
    return coordinate, kind, algorithm, suffix
end

function _fits_reject_external_data(header, suffix; allow_lookup=false, allow_tab=false)
    axis_suffix = isempty(suffix) ? raw"$" : "$(suffix)\$"
    for (keyword, value) in header
        if occursin(Regex("^(?:CPDIS|D2IMDIS)[1-9][0-9]*$axis_suffix"), keyword)
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

function _fits_tabular_data(loader, header, key, parsed)
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

struct _FITSLookupSpec
    extname::String
    extver::Int
    transpose::Bool
end

function _fits_lookup_spec(header, axis, distortion_prefix, parameter_prefix, error_prefix, extname, suffix, minerr)
    distortion_keyword = "$distortion_prefix$axis$suffix"
    haskey(header, distortion_keyword) || return nothing
    uppercase(strip(String(header[distortion_keyword]))) == "LOOKUP" || throw(ArgumentError("unsupported Paper IV distortion type $(repr(header[distortion_keyword])) in $distortion_keyword"))
    error_keyword = "$error_prefix$axis$suffix"
    error_estimate = _fits_number(get(header, error_keyword, 0.0), error_keyword)
    error_estimate < minerr && return nothing
    parameter = "$parameter_prefix$axis$suffix"
    axis_keyword = "$parameter.AXIS.$axis"
    haskey(header, axis_keyword) || throw(ArgumentError("Paper IV lookup distortion $distortion_keyword requires $axis_keyword"))
    table_axis = _fits_integer(header[axis_keyword], axis_keyword; minimum=1)
    table_axis in (1, 2) || throw(ArgumentError("$axis_keyword must be 1 or 2, got $table_axis"))
    extver_keyword = "$parameter.EXTVER"
    extver = _fits_integer(get(header, extver_keyword, 1), extver_keyword; minimum=1)
    return _FITSLookupSpec(extname, extver, table_axis != axis)
end

function _fits_lookup_data(loader, header, key, minerr)
    isfinite(minerr) && minerr >= 0 || throw(ArgumentError("minerr must be finite and nonnegative"))
    suffix = _fits_key(key)
    referenced = any(axis -> haskey(header, "D2IMDIS$axis$suffix") || haskey(header, "CPDIS$axis$suffix"), 1:2)
    referenced || return nothing
    detector_specs = ntuple(axis -> _fits_lookup_spec(header, axis, "D2IMDIS", "D2IM", "D2IMERR", "D2IMARR", suffix, minerr), 2)
    cpdis_specs = ntuple(axis -> _fits_lookup_spec(header, axis, "CPDIS", "DP", "CPERR", "WCSDVARR", suffix, minerr), 2)
    load(spec) = spec === nothing ? nothing : loader(spec)
    return (detector_to_image=map(load, detector_specs), cpdis=map(load, cpdis_specs))
end

function _fits_lookup_table(data::AbstractMatrix, header, transpose)
    values = transpose ? permutedims(data) : copy(data)
    return BilinearLookupTable(values; crpix=(get(header, "CRPIX1", 0.0), get(header, "CRPIX2", 0.0)),
        crval=(get(header, "CRVAL1", 0.0), get(header, "CRVAL2", 0.0)), cdelt=(get(header, "CDELT1", 1.0), get(header, "CDELT2", 1.0)))
end

function _fits_projection(code, header, latitude_axis, suffix)
    parameter(index, default) = _fits_number(get(header, "PV$(latitude_axis)_$(index)$suffix", default), "PV$(latitude_axis)_$(index)$suffix")
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
        coefficients = ntuple(index -> parameter(index - 1, index == 2 ? 1.0 : 0.0), final_index + 1)
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
    code == "COP" && return ConicPerspective(parameter(1, 45.0), parameter(2, 0.0))
    code == "COD" && return ConicEquidistant(parameter(1, 45.0), parameter(2, 0.0))
    code == "COE" && return ConicEqualArea(parameter(1, 45.0), parameter(2, 0.0))
    code == "COO" && return ConicOrthomorphic(parameter(1, 45.0), parameter(2, 0.0))
    code == "BON" && return BonneProjection(parameter(1, 45.0))
    code == "TSC" && return TangentialSphericalCube()
    code == "CSC" && return COBESphericalCube()
    code == "QSC" && return QuadrilateralizedSphericalCube()
    if code == "HPX"
        H, K = parameter(1, 4.0), parameter(2, 3.0)
        isinteger(H) && isinteger(K) || throw(ArgumentError("HPX H and K parameters must be integers"))
        return HEALPixProjection(Int(H), Int(K))
    end
    code == "XPH" && return HEALPixPolarProjection()
    code in ("TPV", "TPD") && return GnomonicProjection()
    throw(ArgumentError("unsupported FITS spherical projection code $(repr(code))"))
end

function _fits_tpv_coefficients(header, axis, suffix, defaults)
    indices = Tuple(index for index in 0:59 if haskey(header, "PV$(axis)_$(index)$suffix"))
    isempty(indices) && return defaults
    final_index = maximum(indices)
    return ntuple(final_index + 1) do position
        keyword = "PV$(axis)_$(position - 1)$suffix"
        _fits_number(get(header, keyword, 0.0), keyword)
    end
end

function _fits_tpv_transform(header, suffix, count, longitude_axis, latitude_axis, code)
    code in ("TPV", "TPD") || return IdentityTransform{count}()
    xcoeff = _fits_tpv_coefficients(header, longitude_axis, suffix, (0.0, 1.0))
    ycoeff = _fits_tpv_coefficients(header, latitude_axis, suffix, (0.0, 1.0))
    return _fits_selected_transform(count, longitude_axis, latitude_axis, TPVTransform(xcoeff, ycoeff))
end

function _fits_projection_code(header, parsed, longitude_axis, latitude_axis, suffix)
    code = parsed[longitude_axis][3]
    code == "TAN" || return code
    for axis in (longitude_axis, latitude_axis), parameter in 5:59
        haskey(header, "PV$(axis)_$(parameter)$suffix") && return "TPV"
    end
    return code
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

function _fits_sip_transform(header, suffix, crpix, count, ctype_suffixes; enabled=true)
    enabled || return IdentityTransform{count}()
    has_sip = any(==("-SIP"), ctype_suffixes) || haskey(header, "A_ORDER$suffix") || haskey(header, "B_ORDER$suffix")
    has_sip || return IdentityTransform{count}()
    count >= 2 || throw(ArgumentError("FITS SIP distortion requires at least two pixel axes"))
    all(value -> value == "-SIP", ctype_suffixes[1:2]) || throw(ArgumentError("FITS SIP must be declared on the first two CTYPE axes"))
    a = _fits_sip_polynomial(header, "A", suffix)
    b = _fits_sip_polynomial(header, "B", suffix)
    (a === nothing) == (b === nothing) || throw(ArgumentError("FITS SIP requires both A_ORDER and B_ORDER"))
    a === nothing && throw(ArgumentError("FITS CTYPE declares SIP but A_ORDER and B_ORDER are absent"))
    forward = SIPTransform(a, b)
    ap = _fits_sip_polynomial(header, "AP", suffix)
    bp = _fits_sip_polynomial(header, "BP", suffix)
    (ap === nothing) == (bp === nothing) || throw(ArgumentError("FITS SIP inverse requires both AP_ORDER and BP_ORDER"))
    correction = ap === nothing ? forward : WithInverse(forward, SIPTransform(ap, bp))
    identity = Matrix{Float64}(I, 2, 2)
    correction = AffineTransform(identity, -crpix[1:2]) | correction | AffineTransform(identity, crpix[1:2])
    return count == 2 ? correction : parallel(correction, IdentityTransform{count - 2}())
end

function _fits_cd_matrix(header, count, cdelt, suffix)
    has_cd = any(haskey(header, "CD$(i)_$(j)$suffix") for i in 1:count for j in 1:count)
    has_pc = any(haskey(header, "PC$(i)_$(j)$suffix") for i in 1:count for j in 1:count)
    has_cd && has_pc && throw(ArgumentError("FITS CDi_j and PCi_j keywords must not be mixed"))
    matrix = zeros(Float64, count, count)
    if has_cd
        for i in 1:count, j in 1:count
            keyword = "CD$(i)_$(j)$suffix"
            matrix[i, j] = haskey(header, keyword) ? _fits_number(header[keyword], keyword) : 0.0
        end
        return matrix
    end
    pc = Matrix{Float64}(I, count, count)
    if has_pc
        for i in 1:count, j in 1:count
            keyword = "PC$(i)_$(j)$suffix"
            haskey(header, keyword) && (pc[i, j] = _fits_number(header[keyword], keyword))
        end
    elseif count >= 2 && haskey(header, "CROTA2$suffix")
        angle = deg2rad(_fits_number(header["CROTA2$suffix"], "CROTA2$suffix"))
        cosine, sine = cos(angle), sin(angle)
        pc[1, 1] = cosine
        pc[1, 2] = -cdelt[2] * sine / cdelt[1]
        pc[2, 1] = cdelt[1] * sine / cdelt[2]
        pc[2, 2] = cosine
    end
    for i in 1:count, j in 1:count
        matrix[i, j] = cdelt[i] * pc[i, j]
    end
    return matrix
end

function _fits_angle_scale(unit, keyword)
    normalized = lowercase(strip(unit))
    normalized in ("", "deg", "degree", "degrees") && return 1.0
    normalized in ("rad", "radian", "radians") && return 180 / pi
    normalized in ("arcmin", "amin") && return 1 / 60
    normalized in ("arcsec", "asec") && return 1 / 3600
    normalized == "mas" && return 1 / 3.6e6
    throw(ArgumentError("unsupported celestial $keyword=$(repr(unit))"))
end

function _fits_spectral_scale(coordinate, unit, keyword)
    normalized = lowercase(replace(strip(unit), ' ' => ""))
    isempty(normalized) && return 1.0
    if coordinate == "FREQ"
        factors = Dict("hz" => 1.0, "khz" => 1e3, "mhz" => 1e6, "ghz" => 1e9, "thz" => 1e12)
    elseif coordinate == "AFRQ"
        factors = Dict("rad/s" => 1.0, "rads-1" => 1.0)
    elseif coordinate in ("WAVE", "AWAV")
        factors = Dict("m" => 1.0, "cm" => 1e-2, "mm" => 1e-3, "um" => 1e-6, "micron" => 1e-6, "nm" => 1e-9, "angstrom" => 1e-10, "aa" => 1e-10)
    elseif coordinate in ("VRAD", "VOPT", "VELO")
        factors = Dict("m/s" => 1.0, "ms-1" => 1.0, "km/s" => 1e3, "kms-1" => 1e3)
    elseif coordinate == "ENER"
        factors = Dict("j" => 1.0, "erg" => 1e-7, "ev" => 1.602176634e-19, "kev" => 1.602176634e-16, "mev" => 1.602176634e-13)
    elseif coordinate == "WAVN"
        factors = Dict("1/m" => 1.0, "m-1" => 1.0, "1/cm" => 100.0, "cm-1" => 100.0)
    elseif coordinate in ("ZOPT", "BETA")
        factors = Dict("1" => 1.0)
    else
        return 1.0
    end
    haskey(factors, normalized) || throw(ArgumentError("unsupported spectral $keyword=$(repr(unit))"))
    return factors[normalized]
end

function _fits_time_scale(unit, keyword)
    normalized = lowercase(strip(unit))
    factors = Dict("" => 1.0, "s" => 1.0, "sec" => 1.0, "min" => 60.0, "h" => 3600.0, "hr" => 3600.0, "d" => 86400.0, "day" => 86400.0, "ms" => 1e-3, "us" => 1e-6, "ns" => 1e-9)
    haskey(factors, normalized) || throw(ArgumentError("unsupported time $keyword=$(repr(unit))"))
    return factors[normalized]
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
    location = if all(keyword -> haskey(header, keyword), ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z"))
        ntuple(i -> _fits_number(header[("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z")[i]], ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z")[i]), 3)
    else
        nothing
    end
    return (reference_epoch=reference, time_scale=String(value("TIMESYS", "")), time_format=reference === nothing ? "" : "mjd",
        reference_position=String(value("TREFPOS", "")), reference_direction=value("TREFDIR"), location=location)
end

function _fits_reference_system(header, suffix, parsed)
    has_equatorial = any(item -> item[1] in ("RA", "DEC"), parsed)
    equinox_value = get(header, "EQUINOX$suffix", nothing)
    equinox = equinox_value === nothing ? nothing : _fits_number(equinox_value, "EQUINOX$suffix")
    radesys = String(get(header, "RADESYS$suffix", ""))
    if isempty(strip(radesys)) && has_equatorial
        radesys = equinox === nothing ? "ICRS" : equinox < 1984 ? "FK4" : "FK5"
    end
    return uppercase(strip(radesys)), equinox
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
    if algorithm == "LOG"
        x_type = parent
    else
        haskey(_FITS_SPECTRAL_ALGORITHMS, algorithm) || throw(ArgumentError("unsupported FITS nonlinear spectral algorithm $(repr(algorithm))"))
        x_type, algorithm_parent = _FITS_SPECTRAL_ALGORITHMS[algorithm]
        algorithm_parent == parent || throw(ArgumentError("FITS spectral coordinate $coordinate is incompatible with algorithm $algorithm"))
    end
    !isfinite(rest_frequency) && isfinite(rest_wavelength) && (rest_frequency = _FITS_SPEED_OF_LIGHT / rest_wavelength)
    !isfinite(rest_wavelength) && isfinite(rest_frequency) && (rest_wavelength = _FITS_SPEED_OF_LIGHT / rest_frequency)
    requires_rest = coordinate in ("VRAD", "VOPT", "ZOPT", "VELO", "BETA") || x_type == :V || parent == :V
    requires_rest && !isfinite(rest_frequency) && throw(ArgumentError("FITS spectral algorithm $algorithm requires RESTFRQ or RESTWAV"))
    transform = FITSSpectralTransform(Val(x_type), Val(parent), Val(Symbol(coordinate)), Val(Symbol(algorithm)), reference, rest_frequency, rest_wavelength; output_scale)
    derivative = algorithm == "LOG" ? 1.0 : ForwardDiff.derivative(reference) do value
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

function _fits_world_transform(count, crval, longitude_axis, latitude_axis, projection, lonpole, latpole, output_scaling)
    offset = copy(crval)
    if projection !== nothing
        offset[longitude_axis] = 0.0
        offset[latitude_axis] = 0.0
    end
    result = AffineTransform(Matrix{Float64}(I, count, count), offset)
    if projection !== nothing
        celestial = celestial_transform(projection, (crval[longitude_axis], crval[latitude_axis]); lonpole, latpole)
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

function _parse_fits_wcs(header::Dict{String,Any}; key::Char=' ', preserve_units::Bool=false, minerr::Real=0.0, lookup_data=nothing, tabular_data=())
    isfinite(minerr) && minerr >= 0 || throw(ArgumentError("minerr must be finite and nonnegative"))
    suffix = _fits_key(key)
    _fits_reject_external_data(header, suffix; allow_lookup=lookup_data !== nothing, allow_tab=!isempty(tabular_data))
    count = _fits_axis_count(header, key)
    keywords(prefix, default) = [get(header, "$prefix$i$suffix", default) for i in 1:count]
    ctype = String.(keywords("CTYPE", ""))
    parsed = _fits_ctype.(ctype)
    for (axis, (_, kind, algorithm, extra)) in enumerate(parsed)
        isempty(extra) || extra == "-SIP" || throw(ArgumentError("FITS CTYPE$axis=$(repr(ctype[axis])) uses unsupported distortion suffix $(repr(extra))"))
        if kind == :spectral && !isempty(algorithm) && algorithm != "TAB" && algorithm != "LOG" && algorithm ∉ ("GRI", "GRA") && !haskey(_FITS_SPECTRAL_ALGORITHMS, algorithm)
            throw(ArgumentError("FITS CTYPE$axis=$(repr(ctype[axis])) uses unsupported nonlinear algorithm $(repr(algorithm))"))
        elseif kind in (:time, :linear) && !isempty(algorithm) && algorithm != "TAB"
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

    crpix = [_fits_number(value, "CRPIX$i$suffix") for (i, value) in enumerate(keywords("CRPIX", 0.0))]
    raw_crval = [_fits_number(value, "CRVAL$i$suffix") for (i, value) in enumerate(keywords("CRVAL", 0.0))]
    cdelt = [_fits_number(value, "CDELT$i$suffix") for (i, value) in enumerate(keywords("CDELT", 1.0))]
    cunit = String.(keywords("CUNIT", ""))
    timeunit = String(get(header, "TIMEUNIT$suffix", get(header, "TIMEUNIT", "")))
    for i in eachindex(cunit)
        parsed[i][2] == :time && isempty(strip(cunit[i])) && (cunit[i] = timeunit)
    end
    scaling = [_fits_axis_scale(parsed[i], cunit[i], "CUNIT$i$suffix") for i in 1:count]
    restfrq_value = get(header, "RESTFRQ$suffix", get(header, "RESTFREQ$suffix", nothing))
    restwav_value = get(header, "RESTWAV$suffix", nothing)
    restfrq = restfrq_value === nothing ? nothing : _fits_number(restfrq_value, "RESTFRQ$suffix")
    restwav = restwav_value === nothing ? nothing : _fits_number(restwav_value, "RESTWAV$suffix")
    spectral_transforms = Any[nothing for _ in 1:count]
    spectral_derivatives = ones(Float64, count)
    for axis in 1:count
        parsed[axis][2] == :spectral || continue
        algorithm = parsed[axis][3]
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
    raw_cd = _fits_cd_matrix(header, count, cdelt, suffix)
    cd = copy(raw_cd)
    for i in 1:count
        cd[i, :] .*= linear_scaling[i]
    end

    projection = nothing
    projection_code = ""
    lonpole = nothing
    latpole = _fits_number(get(header, "LATPOLE$suffix", 90.0), "LATPOLE$suffix")
    if longitude_axis != 0
        lon_code = _fits_projection_code(header, parsed, longitude_axis, latitude_axis, suffix)
        lat_code = parsed[latitude_axis][3]
        lon_code in ("TPV", "TPD") && lat_code == "TAN" && (lat_code = lon_code)
        isempty(lon_code) && throw(ArgumentError("FITS celestial CTYPE$(longitude_axis) has no projection code"))
        lon_code == lat_code || throw(ArgumentError("FITS celestial axes use different projection codes $(repr(lon_code)) and $(repr(lat_code))"))
        projection_code = lon_code
        if projection_code != "TAB"
            projection = _fits_projection(projection_code, header, latitude_axis, suffix)
            default_pole = default_lonpole(crval[latitude_axis], native_fiducial(projection)[2])
            lonpole = _fits_number(get(header, "LONPOLE$suffix", default_pole), "LONPOLE$suffix")
        end
    end

    detector = CoordinateFrame{count}(; name="detector", axes_names=ntuple(i -> "pixel$i", count), axis_types=ntuple(_ -> PIXEL, count), unit=ntuple(_ -> "pix", count))
    intermediate = CoordinateFrame{count}(; name="intermediate", axes_names=ntuple(i -> "intermediate$i", count), axis_types=ntuple(_ -> SPATIAL, count), unit=Tuple(map(item -> item[3], _fits_axis_metadata.(parsed, cunit, 1:count))))
    radesys, equinox = _fits_reference_system(header, suffix, parsed)
    velosys_value, zsource_value = get(header, "VELOSYS$suffix", nothing), get(header, "ZSOURCE$suffix", nothing)
    spectral_metadata = (specsys=String(get(header, "SPECSYS$suffix", "")), ssysobs=String(get(header, "SSYSOBS$suffix", "")),
        velosys=velosys_value === nothing ? nothing : _fits_number(velosys_value, "VELOSYS$suffix"),
        zsource=zsource_value === nothing ? nothing : _fits_number(zsource_value, "ZSOURCE$suffix"), ssyssrc=String(get(header, "SSYSSRC$suffix", "")))
    metadata = (radesys=radesys, equinox=equinox, restfrq=restfrq, restwav=restwav, spectral=spectral_metadata, time=_fits_time_reference(header, suffix), key=key)
    descriptions = [_fits_axis_metadata(parsed[i], cunit[i], i) for i in 1:count]
    names = ntuple(i -> String(get(header, "CNAME$i$suffix", isempty(parsed[i][1]) ? "world$i" : parsed[i][1])), count)
    units = ntuple(i -> preserve_units && !isempty(strip(cunit[i])) ? cunit[i] : descriptions[i][3], count)
    world = CoordinateFrame{count}(; name="world", axes_names=names, axis_types=Tuple(first.(descriptions)), axis_physical_types=Tuple(getindex.(descriptions, 2)), unit=units, reference_frame=metadata)
    if projection_code in ("TPV", "TPD")
        has_cd = any(haskey(header, "CD$(i)_$(j)$suffix") for i in 1:count for j in 1:count)
        if has_cd
            pre_matrix = raw_cd
            post_scaling = scaling
        else
            all(!iszero, cdelt) || throw(ArgumentError("FITS CDELT values must be nonzero for TPV/TPD"))
            pre_matrix = Diagonal(inv.(cdelt)) * raw_cd
            post_scaling = cdelt .* scaling
        end
        tpv = _fits_tpv_transform(header, suffix, count, longitude_axis, latitude_axis, projection_code)
        pixel_to_intermediate = AffineTransform(pre_matrix, -(pre_matrix * crpix)) | tpv | AffineTransform(Diagonal(post_scaling))
    else
        pixel_to_intermediate = AffineTransform(cd, -(cd * crpix))
    end
    sip = _fits_sip_transform(header, suffix, crpix, count, getindex.(parsed, 4); enabled=!(projection_code in ("TPV", "TPD")))
    has_lookup = lookup_data !== nothing && any(table -> table !== nothing, (lookup_data.detector_to_image..., lookup_data.cpdis...))
    if !has_lookup
        pixel_to_intermediate = sip | pixel_to_intermediate
    else
        distortion = FITSLookupDistortion(Val(count), lookup_data.detector_to_image, sip, lookup_data.cpdis)
        pixel_to_intermediate = distortion | pixel_to_intermediate
    end
    output_scaling = [parsed[i][3] == "TAB" || spectral_transforms[i] !== nothing ? 1.0 : preserve_units ? inv(scaling[i]) : 1.0 for i in 1:count]
    intermediate_to_world = _fits_world_transform(count, crval, longitude_axis, latitude_axis, projection, lonpole, latpole, output_scaling)
    for axis in eachindex(spectral_transforms)
        transform = spectral_transforms[axis]
        transform === nothing || (intermediate_to_world = intermediate_to_world | _fits_selected_axis_transform(count, axis, transform))
    end
    if !isempty(tabular_data)
        tabular_scaling = ntuple(i -> parsed[i][3] == "TAB" && !preserve_units ? scaling[i] : 1.0, count)
        intermediate_to_world = intermediate_to_world | FITSTabularTransform(Val(count), tabular_data, tabular_scaling)
    end
    name = String(get(header, "WCSNAME$suffix", ""))
    return WCS(Step(detector, pixel_to_intermediate), Step(intermediate, intermediate_to_world), Step(world); name, pixel_shape=_fits_pixel_shape(header, count))
end

"""
    parse_fits_wcs(source; key=' ', preserve_units=false, minerr=0.0)

Parse one primary or alternate FITS image-header WCS into a compiled `WCS` pipeline.
Pixel coordinates use the FITS 1-based convention encoded by `CRPIX`.
"""
parse_fits_wcs(source::AbstractDict; key::Char=' ', preserve_units::Bool=false, minerr::Real=0.0) = _parse_fits_wcs(_fits_header(source); key, preserve_units, minerr)

function _fits_alternate_keys(header)
    keys_found = Set{Char}()
    any(keyword -> !isempty(_fits_axis_indices(keyword, "")), keys(header)) && push!(keys_found, ' ')
    for keyword in keys(header), key in 'A':'Z'
        isempty(_fits_axis_indices(keyword, string(key))) || push!(keys_found, key)
    end
    isempty(keys_found) && haskey(header, "NAXIS") && push!(keys_found, ' ')
    return sort!(collect(keys_found))
end

function _fits_malformed_error(error, key)
    error isa Union{ArgumentError,DimensionMismatch,DomainError,InexactError} || throw(error)
    return ArgumentError("FITS WCS key $(repr(key)) is malformed: $(sprint(showerror, error))")
end

"""Parse every primary and alternate FITS WCS representation in `source`."""
function find_all_fits_wcs(source::AbstractDict; preserve_units::Bool=false, minerr::Real=0.0)
    header = _fits_header(source)
    result = Dict{Char,WCS}()
    for key in _fits_alternate_keys(header)
        try
            result[key] = _parse_fits_wcs(header; key, preserve_units, minerr)
        catch error
            error isa InterruptException && rethrow()
            throw(_fits_malformed_error(error, key))
        end
    end
    return result
end
