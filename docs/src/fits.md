# FITS WCS

```@meta
CurrentModule = GWCS
```

[`parse_fits_wcs`](@ref) parses one primary or alternate FITS image-header WCS into an ordinary compiled
[`WCS`](@ref), and [`find_all_fits_wcs`](@ref) returns every one in a header as a `Dict{Char,WCS}`. A header is a
dictionary of keywords to values, or, with FITSIO.jl or FITSFiles.jl loaded, a header or HDU from those packages.

```@example fits
using GWCS

header = Dict(
    "NAXIS" => 2, "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
    "CRPIX1" => 512.0, "CRPIX2" => 512.0,
    "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
    "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
)
wcs = parse_fits_wcs(header)
pixel_to_world(wcs, 511.0, 511.0)  # FITS pixel CRPIX = (512, 512) is GWCS pixel (511, 511)
```

**FITS numbers the first pixel center 1, and GWCS.jl numbers it 0.** The parser adds one to pixel coordinates before the
FITS stages (`CRPIX`, SIP, lookup tables, and prior distortions), so FITS pixel `p` is GWCS pixel `p .- 1`. Results at
the same physical pixel match WCSLIB, and results at the same pixel coordinates match astropy's
`pixel_to_world_values`.

## Reading FITS files

With FITSIO.jl loaded, `parse_fits_wcs` reads a header or HDU from a FITS file. This writes the header above to a file
and reads it back both ways:

```@example fits
using FITSIO

path = joinpath(mktempdir(), "image.fits")
cards = filter(!=("NAXIS"), collect(keys(header)))  # FITSIO writes NAXIS from the image
FITS(path, "w") do file
    write(file, zeros(Float32, 1024, 1024); header=FITSHeader(cards, Any[header[key] for key in cards], fill("", length(cards))))
end

FITS(path) do file
    from_hdu = parse_fits_wcs(file[1])
    from_header = parse_fits_wcs(read_header(file[1]))
    pixel_to_world(from_hdu, 511.0, 511.0), from_hdu.pixel_shape, pixel_to_world(from_header, 511.0, 511.0)
end
```

With FITSFiles.jl loaded, it reads an HDU, or the HDUs of a file with the index of the one holding the WCS:

```@example fits
using FITSFiles

hdus = FITSFiles.fits(path)
pixel_to_world(parse_fits_wcs(hdus, 1), 511.0, 511.0)
```

An HDU in a file, as a FITSIO.jl HDU or FITSFiles.jl HDUs with an index, also reads the lookup-table and `-TAB` arrays
in other HDUs that the header refers to; a header alone cannot.

## Output order and units

World coordinates keep the FITS axis order, including interleaved celestial, spectral, time, Stokes, and linear axes.
They are in degrees for celestial axes, SI units for supported linear spectral axes, and seconds for time axes;
`preserve_units=true` keeps supported `CUNIT` values instead. An image shape missing an `NAXISn` gives `nothing` for
that axis.

Time axes are normalized deliberately: values in `CUNITi`, or in `TIMEUNIT` when `CUNITi` is absent, are converted to
seconds from the header's reference epoch, so downstream code never needs to interpret `TIMEUNIT`; WCSLIB and astropy
instead return time values in the header's unit. The supported time units are `s`, `min`, `h`, `d`, `a` or `yr` (the
Julian year of 365.25 d), `cy` (the Julian century), `ms`, `us`, and `ns`. The single letters `S`, `H`, and `D` are
rejected, because FITS reads them as siemens, henry, and debye.

## Supported conventions

The parser supports:

- PC, CD, and legacy CROTA linear forms (as in WCSLIB, CROTA applies only to a celestial axis pair);
- all 28 FITS spherical projections with their PV parameters, and the pole parameters;
- SIP, TPV, pre-2012 SCAMP TAN-plus-PV headers, Paper IV prior and sequent TPD polynomials, and legacy DSS and
  WAT-encoded TNX/ZPX polynomials, normalized to TPD at construction;
- Paper IV D2IM/CPDIS lookup arrays and multilinear `-TAB` coordinate arrays;
- Paper III spectral algorithms and GRI/GRA grism models.

SIP AP/BP and Paper IV TPD.REV coefficients only start the numerical inverse: `world_to_pixel` seeds the whole-WCS
inverse with them and always refines against the forward WCS. Without them, the numerical inverse starts, like WCSLIB,
from the solution with distortions treated as identity. CSC (COBE quadrilateralized spherical cube) is evaluated in
double precision, whereas WCSLIB uses single precision, so the two differ by up to about 2e-5 degrees. Linked HDU
inputs copy referenced images and table columns into the transforms during construction. Separable monotonic 1D tables
have direct inverses, while coupled tables use the WCS numerical inverse.

## Frames and reference metadata

The world frame has a [`CelestialFrame`](@ref) for the celestial pair, a [`SpectralFrame`](@ref),
[`TemporalFrame`](@ref), or [`StokesFrame`](@ref) for each such axis, and one [`CoordinateFrame`](@ref) for the linear
axes, combined in a [`CompositeFrame`](@ref) named `world` when there is more than one; each frame's `axes_order` places
its axes at their header positions. As in astropy's `wcs_to_celestial_frame`, RA/DEC take ICRS, FK5, FK4, or FK4-NO-E
from `RADESYS` and `EQUINOX` (FK4's `obstime` is its equinox) and GLON/GLAT are Galactic, and these reference frames
are written to ASDF in the form astropy reads. Unlike astropy, `RADESYS` makes only RA/DEC equatorial, so an ecliptic
header is not labeled FK5 or ICRS.

A time axis's reference time follows astropy's FITS WCS: `MJDREF`, or `JDREF`, or an ISO 8601 `DATEREF` that must agree
with them to 0.001 day (legacy forms such as `DD/MM/YY`, which WCSLIB also reads, are rejected), defaulting to zero as
in WCSLIB. The scale is the CTYPE, or `TIMESYS` for `TIME`, with GPS read as TAI 19 s later, and `TREFPOS` TOPOCENTER
(with `OBSGEO-X/Y/Z`) or GEOCENTER gives the location.

Metadata the ASDF frame schemas cannot hold stays on the frame's `reference_frame` as a named tuple, and writing that
frame to ASDF raises an error. That metadata is other celestial systems with their `RADESYS` and `EQUINOX`; `RESTFRQ`,
`RESTWAV`, `SPECSYS`, `SSYSOBS`, `VELOSYS`, `ZSOURCE`, and `SSYSSRC`; and time sub-scales, the `LOCAL` scale, other
`TREFPOS` values, and `TREFDIR`.

## Legacy conventions

Headers that use deprecated or informal conventions are repaired the way WCSLIB (with all of its header-parser
relaxations) and astropy (which also applies WCSLIB's `wcsfix` repairs) read them. With the default `legacy=:warn`,
each parse that needs repairs emits one warning listing them; `legacy=:quiet` repairs silently, and `legacy=:error`
rejects any header that needs a repair. A legacy keyword is used only when its standard spelling is absent; when both
are present the standard keyword wins, except that conflicting zero-padded matrix keywords are rejected.

| Legacy form | Interpretation |
|---|---|
| `NCP` projection | `SIN` with `PV2_1 = 0` and `PV2_2 = cot(CRVAL2)` on the latitude axis; rejected for a zero reference latitude |
| `GLS` projection | `SFL` with its fiducial point offset to the reference latitude (`PV1_0 = 1`, `PV1_2 = CRVAL2` on the longitude axis) |
| `FREQ`, `VELO`, `FELO` with `-LSR`, `-HEL`, `-OBS`, and `VELREF` | `SPECSYS` from the suffix, or from `VELREF mod 256` when present; `VELO` becomes `VOPT`, or `VRAD` when `VELREF ÷ 256 = 1`; `FELO` becomes `VOPT-F2W` |
| Celestial `RA`/`DEC` (or other pairs) without a projection code | Evaluated linearly |
| `PC00i00j`, `CD00i00j`, and zero-padded `PC0i_0j`, `CD0i_0j`, `PV0i_0m`, `PS0i_0m` | The standard keyword |
| `PROJPn` | `PVi_n` on the latitude axis |
| `EPOCH`, `RADECSYS`, `VSOURCE` | `EQUINOX`, `RADESYS`, and the relativistic `ZSOURCE` |
| `DATE-REF`, `MJD-REF`, `MJD-REFI`, `MJD-REFF`, `JD-REF`, `JD-REFI`, `JD-REFF` | The hyphen-free standard keywords; `JDREF*` give the MJD reference epoch |
| Both `PCi_j` and `CDi_j` | `PCi_j` with `CDELTi`; `CDi_j` is ignored |
| `CDi_j` with no terms for an axis | Unit diagonal for that axis (WCSLIB's `cdfix`); `CDELTi` is ignored with `CDi_j` |
| SIP coefficients without a `-SIP` CTYPE suffix | SIP is applied (astropy) |
| `-TAN-SIP` with SCAMP `PVi_m` terms (m ≥ 5) | SIP is applied and the celestial `PVi_m` are ignored (astropy) |
| Non-standard unit spellings such as `DEG`, `HZ`, `MHZ`, `Angstroms`, `SEC` | WCSLIB's safe unit translations (astropy's default); units are otherwise case-sensitive |
| `CROTAi` on the longitude axis only | Ignored, as in WCSLIB, with a warning |

Singular linear transformation matrices are rejected. As in WCSLIB, the conic (COP, COD, COE, COO) and Bonne projections
require `PV2_1` rather than defaulting it, and IRAF `WATi_n` strings are ignored unless they encode TNX or ZPX on both
axes (or CTYPE declares TNX/ZPX). WCSLIB's `cylfix` (re-centering cylindrical all-sky maps whose pixels extend beyond
native longitude ±180°) is not applied, matching astropy's default, and `DATE-OBS`/`OBSGEO-L,B,H` repairs are not
modeled.
