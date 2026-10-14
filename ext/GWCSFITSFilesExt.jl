module GWCSFITSFilesExt

import FITSFiles
import GWCS: _find_all_fits_wcs, _fits_card_dictionary, _fits_find_extension, _fits_header, _parse_linked_fits_wcs, find_all_fits_wcs, parse_fits_wcs

# FITSFiles parses short real literals as Float32; Float64(x) would add binary noise (150.1 -> 150.10000610...),
# so recover the decimal literal from the shortest Float32 representation.
# shortcut: literals with 8-9 significant digits are already truncated inside FITSFiles; fixing that needs FITSFiles to parse reals as Float64.
_card_value(value::Float32) = parse(Float64, replace(string(value), 'f' => 'e'))
_card_value(value) = value

_header_dictionary(cards::FITSFiles.Cards) = _fits_card_dictionary((card.key, _card_value(card.value)) for card in cards)

parse_fits_wcs(cards::FITSFiles.Cards; kwargs...) = parse_fits_wcs(_header_dictionary(cards); kwargs...)
find_all_fits_wcs(cards::FITSFiles.Cards; kwargs...) = find_all_fits_wcs(_header_dictionary(cards); kwargs...)
parse_fits_wcs(hdu::FITSFiles.HDU; kwargs...) = parse_fits_wcs(hdu.cards; kwargs...)
find_all_fits_wcs(hdu::FITSFiles.HDU; kwargs...) = find_all_fits_wcs(hdu.cards; kwargs...)

function _tabular_column(hdu, extname, column)
    data = hdu.data
    target = uppercase(strip(column))
    if data isa NamedTuple || data isa AbstractDict
        for name in keys(data)
            uppercase(strip(String(name))) == target || continue
            values = data[name]
            column_index = findfirst(axis -> uppercase(strip(String(hdu.cards["TTYPE$axis"]))) == target, 1:Int(hdu.cards["TFIELDS"]))
            column_index === nothing && return values
            raw_dimensions = strip(strip(String(get(hdu.cards, "TDIM$column_index", ""))), ['(', ')'])
            isempty(raw_dimensions) && return values
            dimensions = tryparse.(Int, split(raw_dimensions, ','))
            any(isnothing, dimensions) && return values
            return reshape(values, Tuple(something.(dimensions)))
        end
    end
    throw(ArgumentError("FITS TAB extension $(repr(extname)) has no column $(repr(column))"))
end

_lookup_hdu(hdus, extname, extver, extlev=1) = _fits_find_extension((hdu => _header_dictionary(hdu.cards) for hdu in hdus), extname, extver, extlev)

function _parse_hdus(hdus, header; key::Char=' ', preserve_units::Bool=false, minerr::Real=0.0)
    function image(extname, extver)
        table_hdu = _lookup_hdu(hdus, extname, extver)
        return table_hdu.data, _header_dictionary(table_hdu.cards)
    end
    column(extname, extver, extlev, name) = _tabular_column(_lookup_hdu(hdus, extname, extver, extlev), extname, name)
    return _parse_linked_fits_wcs(image, column, header; key, preserve_units, minerr)
end

_hdu_header(hdus, index) = (checkbounds(hdus, index); _fits_header(_header_dictionary(hdus[index].cards)))

function parse_fits_wcs(hdus::AbstractVector{<:FITSFiles.HDU}, index::Integer; key::Char=' ', preserve_units::Bool=false, minerr::Real=0.0)
    return _parse_hdus(hdus, _hdu_header(hdus, index); key, preserve_units, minerr)
end

function find_all_fits_wcs(hdus::AbstractVector{<:FITSFiles.HDU}, index::Integer; preserve_units::Bool=false, minerr::Real=0.0)
    header = _hdu_header(hdus, index)
    return _find_all_fits_wcs(key -> _parse_hdus(hdus, header; key, preserve_units, minerr), header)
end

end
