import json
import sys
import warnings

import asdf
import erfa
import numpy as np
from astropy import units as u
from astropy.coordinates import FK5, ICRS, FK4NoETerms, Galactic, SkyCoord, Supergalactic
from astropy.coordinates.matrix_utilities import rotation_matrix
from astropy.io import fits
from astropy.wcs import WCS as FITSWCS


def iau1976(equinox):
    return erfa.pmat76(2451545.0 + (equinox - 2000.0) * 365.25, 0.0)


# Each system as an astropy frame and the rotation from that frame's coordinates to the system's: astropy precesses FK5
# with IAU 2006, so another FK5 equinox is FK5 J2000 precessed by ERFA's IAU 1976 pmat76, and an FK5 ecliptic is that
# rotated by ERFA's IAU 1980 obliquity.
def frame(spec):
    system, equinox = spec["system"], spec["equinox"]
    if system == "fk5":
        return FK5(equinox="J2000"), iau1976(equinox)
    if system == "ecliptic_fk5":
        obliquity = np.degrees(erfa.obl80(2451545.0 + (equinox - 2000.0) * 365.25, 0.0))
        return FK5(equinox="J2000"), rotation_matrix(obliquity, "x") @ iau1976(equinox)
    if system == "fk4noeterms":
        return FK4NoETerms(equinox=f"B{equinox}", obstime=f"B{spec['obstime']}"), np.eye(3)
    return {"icrs": ICRS(), "galactic": Galactic(), "supergalactic": Supergalactic()}[system], np.eye(3)


def unit_vectors(lon, lat):
    lon, lat = np.radians(lon), np.radians(lat)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def lon_lat(vectors):
    lon = np.mod(np.degrees(np.arctan2(vectors[1], vectors[0])), 360.0)
    lon[lon >= 360.0] -= 360.0
    return lon, np.degrees(np.arctan2(vectors[2], np.hypot(vectors[0], vectors[1])))


def original_values(handle, name, pixels):
    """The original's world coordinates as (lon, lat), and whether its latitude comes first."""
    if name in handle["headers"]:
        header = fits.Header()
        for keyword, value in handle["headers"][name].items():
            header[keyword] = value
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            first, second = FITSWCS(header).pixel_to_world_values(*pixels.T)
        latitude_first = header["CTYPE1"][:4] in ("DEC-", "GLAT", "SLAT", "ELAT")
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
            (source_frame, source_rotation), (target_frame, target_rotation) = frame(handle["sources"][name]), frame(handle["targets"][name])
            base_lon, base_lat = lon_lat(source_rotation.T @ unit_vectors(lon, lat))
            # GWCS.jl converts from and to FK4 through FK5 J2000, so Murray's obstime correction always applies; a frame,
            # unlike a SkyCoord, does not carry the source's obstime into the target.
            sky = SkyCoord(base_lon * u.deg, base_lat * u.deg, frame=source_frame).frame
            if isinstance(source_frame, FK4NoETerms) or isinstance(target_frame, FK4NoETerms):
                sky = sky.transform_to(FK5(equinox="J2000"))
            lon, lat = lon_lat(target_rotation @ sky.transform_to(target_frame).cartesian.xyz.value)
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
