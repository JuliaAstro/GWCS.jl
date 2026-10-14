module GWCSFITSIOExt

import FITSIO
import GWCS: _find_all_fits_wcs, _fits_card_dictionary, _fits_find_extension, _fits_header, _parse_linked_fits_wcs, find_all_fits_wcs, parse_fits_wcs

_header_dictionary(header::FITSIO.FITSHeader) = _fits_card_dictionary(zip(keys(header), values(header)))

parse_fits_wcs(header::FITSIO.FITSHeader; kwargs...) = parse_fits_wcs(_header_dictionary(header); kwargs...)
find_all_fits_wcs(header::FITSIO.FITSHeader; kwargs...) = find_all_fits_wcs(_header_dictionary(header); kwargs...)

function _hdu_at(fitsfile, index)
    FITSIO.fits_movabs_hdu(fitsfile, index)
    type = FITSIO.fits_get_hdu_type(fitsfile)
    return type == :image_hdu ? FITSIO.ImageHDU(fitsfile, index) : type == :binary_table ? FITSIO.TableHDU(fitsfile, index) : FITSIO.ASCIITableHDU(fitsfile, index)
end

function _named_hdu(hdu::FITSIO.HDU, extname, extver, extlev=1)
    fitsfile = getfield(hdu, :fitsfile)
    candidates = (candidate => _header_dictionary(FITSIO.read_header(candidate)) for candidate in (_hdu_at(fitsfile, index) for index in 1:Int(FITSIO.fits_get_num_hdus(fitsfile))))
    return _fits_find_extension(candidates, extname, extver, extlev)
end

function _parse_hdu(hdu::FITSIO.HDU, header; key::Char=' ', preserve_units::Bool=false, minerr::Real=0.0, legacy::Symbol=:warn)
    function image(extname, extver)
        table_hdu = _named_hdu(hdu, extname, extver)
        return FITSIO.read(table_hdu), _header_dictionary(FITSIO.read_header(table_hdu))
    end
    function column(extname, extver, extlev, name)
        table_hdu = _named_hdu(hdu, extname, extver, extlev)
        table_hdu isa Union{FITSIO.TableHDU,FITSIO.ASCIITableHDU} || throw(ArgumentError("FITS TAB extension $(repr(extname)) must be a table HDU"))
        data = FITSIO.read(table_hdu, name)
        return size(data, ndims(data)) == 1 ? dropdims(data; dims=ndims(data)) : data
    end
    return _parse_linked_fits_wcs(image, column, header; key, preserve_units, minerr, legacy)
end

_hdu_header(hdu::FITSIO.HDU) = _fits_header(_header_dictionary(FITSIO.read_header(hdu)))

function parse_fits_wcs(hdu::FITSIO.HDU; key::Char=' ', preserve_units::Bool=false, minerr::Real=0.0, legacy::Symbol=:warn)
    return _parse_hdu(hdu, _hdu_header(hdu); key, preserve_units, minerr, legacy)
end

function find_all_fits_wcs(hdu::FITSIO.HDU; preserve_units::Bool=false, minerr::Real=0.0, legacy::Symbol=:warn)
    header = _hdu_header(hdu)
    return _find_all_fits_wcs(key -> _parse_hdu(hdu, header; key, preserve_units, minerr, legacy), header)
end

end
