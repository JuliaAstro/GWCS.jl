module GWCSFITSFilesExt

import FITSFiles
import GWCS
import GWCS: _fits_header, _fits_lookup_data, _fits_lookup_table, _fits_malformed_error, _fits_tabular_data, _parse_fits_wcs, find_all_fits_wcs, parse_fits_wcs

function _parameter_card(value)
    value isa AbstractString || return nothing
    parts = split(String(value), ':'; limit=2)
    length(parts) == 2 || return nothing
    name, raw = uppercase(strip(parts[1])), strip(parts[2])
    integer = tryparse(Int, raw)
    integer === nothing || return name, integer
    parsed = tryparse(Float64, raw)
    return name, parsed === nothing ? raw : parsed
end

function _set_header_card!(result, key, value)
    result[key] = value
    parsed = _parameter_card(value)
    parsed === nothing && return result
    name, parsed_value = parsed
    (name == "EXTVER" || name == "NAXES" || startswith(name, "AXIS.")) && (result["$key.$name"] = parsed_value)
    return result
end

function _header_dictionary(cards::FITSFiles.Cards)
    result = Dict{String,Any}()
    for card in cards
        _set_header_card!(result, uppercase(strip(String(card.key))), card.value)
    end
    return result
end

parse_fits_wcs(cards::FITSFiles.Cards; kwargs...) = parse_fits_wcs(_header_dictionary(cards); kwargs...)
find_all_fits_wcs(cards::FITSFiles.Cards; kwargs...) = find_all_fits_wcs(_header_dictionary(cards); kwargs...)
parse_fits_wcs(hdu::FITSFiles.HDU; kwargs...) = parse_fits_wcs(hdu.cards; kwargs...)
find_all_fits_wcs(hdu::FITSFiles.HDU; kwargs...) = find_all_fits_wcs(hdu.cards; kwargs...)

function _tabular_column(data, extname, column)
    target = uppercase(strip(column))
    if data isa NamedTuple || data isa AbstractDict
        for name in keys(data)
            uppercase(strip(String(name))) == target && return data[name]
        end
    end
    throw(ArgumentError("FITS TAB extension $(repr(extname)) has no column $(repr(column))"))
end

function parse_fits_wcs(hdus::AbstractVector{<:FITSFiles.HDU}, index::Integer; key::Char=' ', preserve_units::Bool=false, minerr::Real=0.0)
    checkbounds(hdus, index)
    header = _fits_header(_header_dictionary(hdus[index].cards))
    lookup = _fits_lookup_data(header, key, minerr) do spec
        table_hdu = _lookup_hdu(hdus, spec.extname, spec.extver)
        table_hdu.data isa AbstractMatrix || throw(ArgumentError("FITS WCS lookup extension $(repr(spec.extname)) must contain a two-dimensional image"))
        table_header = _fits_header(_header_dictionary(table_hdu.cards))
        _fits_lookup_table(table_hdu.data, table_header, spec.transpose)
    end
    parsed = GWCS._fits_ctype.(String.([get(header, "CTYPE$i$(GWCS._fits_key(key))", "") for i in 1:GWCS._fits_axis_count(header, key)]))
    tabular = _fits_tabular_data(header, key, parsed) do extname, extver, extlev, column
        table_hdu = _lookup_hdu(hdus, extname, extver, extlev)
        return _tabular_column(table_hdu.data, extname, column)
    end
    return _parse_fits_wcs(header; key, preserve_units, minerr, lookup_data=lookup, tabular_data=tabular)
end

function find_all_fits_wcs(hdus::AbstractVector{<:FITSFiles.HDU}, index::Integer; preserve_units::Bool=false, minerr::Real=0.0)
    checkbounds(hdus, index)
    header = _fits_header(_header_dictionary(hdus[index].cards))
    result = Dict{Char,GWCS.WCS}()
    for key in GWCS._fits_alternate_keys(header)
        try
            result[key] = parse_fits_wcs(hdus, index; key, preserve_units, minerr)
        catch error
            error isa InterruptException && rethrow()
            throw(_fits_malformed_error(error, key))
        end
    end
    return result
end

function _lookup_hdu(hdus, extname, extver, extlev=1)
    for hdu in hdus
        cards = hdu.cards
        haskey(cards, "EXTNAME") || continue
        uppercase(strip(String(cards["EXTNAME"]))) == uppercase(extname) || continue
        Int(get(cards, "EXTVER", 1)) == extver || continue
        Int(get(cards, "EXTLEVEL", 1)) == extlev && return hdu
    end
    throw(ArgumentError("FITS WCS extension $(repr(extname)) with EXTVER=$extver and EXTLEVEL=$extlev was not found"))
end

end
