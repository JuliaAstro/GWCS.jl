import json
import sys
import warnings

import asdf
import numpy as np
from astropy import units as u
from astropy.coordinates import FK5, ICRS, Galactic, SkyCoord, Supergalactic
from astropy.io import fits
from astropy.wcs import WCS as FITSWCS

FRAMES = {"icrs": lambda equinox: ICRS(), "galactic": lambda equinox: Galactic(), "supergalactic": lambda equinox: Supergalactic(),
          "fk5": lambda equinox: FK5(equinox=f"J{equinox}")}


def frame(spec):
    return FRAMES[spec["system"]](spec["equinox"])


def original_values(handle, name, pixels):
    """The original's world coordinates as (lon, lat), and whether its latitude comes first."""
    if name in handle["headers"]:
        header = fits.Header()
        for keyword, value in handle["headers"][name].items():
            header[keyword] = value
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            first, second = FITSWCS(header).pixel_to_world_values(*pixels.T)
        latitude_first = header["CTYPE1"][:4] in ("DEC-", "GLAT", "SLAT")
    else:
        first, second = handle["originals"][name](*pixels.T)
        latitude_first = False
    return ((second, first) if latitude_first else (first, second)), latitude_first


def main():
    with asdf.open(sys.argv[1]) as handle:
        pixels = np.asarray(handle["pixels"], dtype=float)
        records = []
        for name, wcs in handle["wcs"].items():
            (lon, lat), latitude_first = original_values(handle, name, pixels)
            sky = SkyCoord(lon * u.deg, lat * u.deg, frame=frame(handle["sources"][name])).transform_to(frame(handle["targets"][name]))
            lon, lat = sky.spherical.lon.wrap_at(360 * u.deg).deg, sky.spherical.lat.deg
            astropy = np.column_stack((lat, lon) if latitude_first else (lon, lat))
            gwcs = np.column_stack(wcs(*pixels.T, with_bounding_box=False))
            # Python GWCS inverts analytically through `backward_transform`, which raises when there is no inverse.
            try:
                backward = wcs.backward_transform
            except NotImplementedError:
                backward = None
            inverse = np.column_stack(backward(*gwcs.T)) if backward is not None else None
            records.append({"wcs": name, "frame": wcs.output_frame.reference_frame.name, "gwcs": gwcs.tolist(), "astropy": astropy.tolist(),
                            "inverse": None if inverse is None else inverse.tolist()})
    print(json.dumps(records))


if __name__ == "__main__":
    main()
