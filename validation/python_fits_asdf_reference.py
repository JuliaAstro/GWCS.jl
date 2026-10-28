"""Read GWCS.jl's ASDF rewrite of FITS headers and compare its frames and coordinates with astropy's FITS WCS, which
applies the headers' distortions."""

import json
import sys
import warnings

import asdf
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import wcs_to_celestial_frame
from gwcs import coordinate_frames as cf


def fits_header(cards):
    """FITS cards from a keyword mapping, with Paper IV keywords such as DP1.NAXES as record-valued cards."""
    lines = []
    for key, value in cards.items():
        keyword, record = key.split(".", 1) if "." in key else (key, None)
        if record is not None:
            value = f"{record}: {value!r}"
        lines.append(f"{keyword:8}= " + (f"'{value}'" if isinstance(value, str) else repr(value)))
    return fits.Header.fromstring("\n".join(lines), sep="\n")


def frames(frame):
    return frame.frames if isinstance(frame, cf.CompositeFrame) else [frame]


def main():
    warnings.simplefilter("ignore")
    headers = json.load(open(sys.argv[2]))
    records = []
    with asdf.open(sys.argv[1]) as af:
        for name, gwcs in af["wcs"].items():
            pixel = np.arange(gwcs.pixel_n_dim, dtype=float) * 37.0 + 120.5
            values = np.array(gwcs.pixel_to_world_values(*pixel), dtype=float).ravel()
            record = {"name": name, "pixel": pixel.tolist(), "values": values.tolist()}
            try:
                fitswcs = WCS(fits_header(headers[name]), relax=True)
            except ValueError as error:
                # astropy reads CPDISja itself for lookup tables, so WCSLIB never sees a prior TPD.
                record["astropy_error"] = str(error).splitlines()[-1]
                records.append(record)
                continue
            expected = np.array(fitswcs.pixel_to_world_values(*pixel), dtype=float).ravel()
            record["astropy_error_value"] = float(np.max(np.abs(values - expected)))
            celestial = [frame for frame in frames(gwcs.output_frame) if isinstance(frame, cf.CelestialFrame)]
            temporal = [frame for frame in frames(gwcs.output_frame) if isinstance(frame, cf.TemporalFrame)]
            if celestial:
                reference, astropy_reference = celestial[0].reference_frame, wcs_to_celestial_frame(fitswcs)
                record["frame"], record["astropy_frame"] = repr(reference), repr(astropy_reference)
                record["frame_matches"] = reference.is_equivalent_frame(astropy_reference)
            elif not temporal:
                record["frame"], record["astropy_frame"], record["frame_matches"] = "no celestial or temporal frame", "", True
            else:
                # astropy's FITS WCS gives the time at a pixel as a Time from its reference time and offset.
                reference = temporal[0].reference_frame
                computed, astropy_time = (next(value for value in np.atleast_1d(world) if isinstance(value, Time)) if not isinstance(world, Time) else world
                    for world in (gwcs.pixel_to_world(*pixel), fitswcs.pixel_to_world(*pixel)))
                record["frame"] = f"Time {reference.scale} {reference.format} location={reference.location is not None}"
                record["astropy_frame"] = f"Time {astropy_time.scale} location={astropy_time.location is not None}"
                same_location = (reference.location is None) == (astropy_time.location is None) and (
                    reference.location is None or np.allclose(u.Quantity(reference.location.geocentric).to_value(u.m), u.Quantity(astropy_time.location.geocentric).to_value(u.m)))
                record["frame_matches"] = bool(computed.scale == astropy_time.scale and same_location
                    and abs((computed - astropy_time).sec) < 1e-5)
            records.append(record)
    print(json.dumps(records))


if __name__ == "__main__":
    main()
