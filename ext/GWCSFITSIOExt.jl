module GWCSFITSIOExt

import FITSIO
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

function _header_dictionary(header::FITSIO.FITSHeader)
    result = Dict{String,Any}()
    for (key, value) in zip(keys(header), values(header))
        _set_header_card!(result, uppercase(strip(String(key))), value)
    end
    return result
end

parse_fits_wcs(header::FITSIO.FITSHeader; kwargs...) = parse_fits_wcs(_header_dictionary(header); kwargs...)
find_all_fits_wcs(header::FITSIO.FITSHeader; kwargs...) = find_all_fits_wcs(_header_dictionary(header); kwargs...)

function _hdu_at(fitsfile, index)
    FITSIO.fits_movabs_hdu(fitsfile, index)
    type = FITSIO.fits_get_hdu_type(fitsfile)
    return type == :image_hdu ? FITSIO.ImageHDU(fitsfile, index) : type == :binary_table ? FITSIO.TableHDU(fitsfile, index) : FITSIO.ASCIITableHDU(fitsfile, index)
end

function _named_hdu(hdu::FITSIO.HDU, extname, extver, extlev=1)
    fitsfile = getfield(hdu, :fitsfile)
    for index in 1:Int(FITSIO.fits_get_num_hdus(fitsfile))
        candidate = _hdu_at(fitsfile, index)
        header = _header_dictionary(FITSIO.read_header(candidate))
        uppercase(strip(String(get(header, "EXTNAME", "")))) == uppercase(strip(extname)) || continue
        Int(get(header, "EXTVER", 1)) == extver || continue
        Int(get(header, "EXTLEVEL", 1)) == extlev && return candidate
    end
    throw(ArgumentError("FITS WCS extension $(repr(extname)) with EXTVER=$extver and EXTLEVEL=$extlev was not found"))
end

function _parse_hdu(hdu::FITSIO.HDU; key::Char=' ', preserve_units::Bool=false, minerr::Real=0.0)
    header = _fits_header(_header_dictionary(FITSIO.read_header(hdu)))
    lookup = _fits_lookup_data(header, key, minerr) do spec
        table_hdu = _named_hdu(hdu, spec.extname, spec.extver)
        data = FITSIO.read(table_hdu)
        data isa AbstractMatrix || throw(ArgumentError("FITS WCS lookup extension $(repr(spec.extname)) must contain a two-dimensional image"))
        table_header = _fits_header(_header_dictionary(FITSIO.read_header(table_hdu)))
        _fits_lookup_table(data, table_header, spec.transpose)
    end
    parsed = GWCS._fits_ctype.(String.([get(header, "CTYPE$i$(GWCS._fits_key(key))", "") for i in 1:GWCS._fits_axis_count(header, key)]))
    tabular = _fits_tabular_data(header, key, parsed) do extname, extver, extlev, column
        table_hdu = _named_hdu(hdu, extname, extver, extlev)
        table_hdu isa Union{FITSIO.TableHDU,FITSIO.ASCIITableHDU} || throw(ArgumentError("FITS TAB extension $(repr(extname)) must be a table HDU"))
        FITSIO.read(table_hdu, column)
    end
    return _parse_fits_wcs(header; key, preserve_units, minerr, lookup_data=lookup, tabular_data=tabular)
end

parse_fits_wcs(hdu::FITSIO.HDU; kwargs...) = _parse_hdu(hdu; kwargs...)

function find_all_fits_wcs(hdu::FITSIO.HDU; preserve_units::Bool=false, minerr::Real=0.0)
    header = _fits_header(_header_dictionary(FITSIO.read_header(hdu)))
    result = Dict{Char,GWCS.WCS}()
    for key in GWCS._fits_alternate_keys(header)
        try
            result[key] = _parse_hdu(hdu; key, preserve_units, minerr)
        catch error
            error isa InterruptException && rethrow()
            throw(_fits_malformed_error(error, key))
        end
    end
    return result
end

end
