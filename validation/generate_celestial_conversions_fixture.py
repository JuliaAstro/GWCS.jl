"""Positions converted between celestial systems by astropy and ERFA, for `with_celestial_system`.

`points` holds longitudes and latitudes in degrees, including points near the poles and on either side of longitude
zero. `systems` names each system with its equinox (Julian for FK5, Besselian for FK4) and FK4 obstime (Besselian),
zero where they do not apply, and `converted[source][target]` holds the conversion of `points`, read as coordinates in
`source`, to `target`, with longitudes in [0, 360).

Each system is an astropy frame, rotated where astropy has no equivalent: an FK5 equinox other than J2000 is FK5 J2000
precessed by ERFA's IAU 1976 `pmat76`, since astropy precesses FK5 with IAU 2006; an FK5 ecliptic is that, rotated by
ERFA's IAU 1980 obliquity `obl80`; an FK4 ecliptic is FK4-NO-E at its equinox (and obstime), rotated by Newcomb's
obliquity as SLALIB's sla_ETRMS gives it; and the ICRS ecliptic is astropy's BarycentricMeanEcliptic at J2000. A
conversion from or to an FK4 system passes explicitly through FK5 J2000, so Murray's obstime correction always applies,
and an FK4 source first through the FK4-NO-E frame at its own equinox and obstime: astropy's frame graph converts FK4
to Galactic coordinates through the B1950 pole without the obstime, and gives intermediate frames a target's equinox.
`default[source][target]` holds astropy's own `transform_to` result where both systems are astropy frames.
"""

from pathlib import Path

import asdf
import erfa
import numpy as np
from astropy import units as u
from astropy.coordinates import FK4, FK5, ICRS, BarycentricMeanEcliptic, FK4NoETerms, Galactic, SkyCoord, Supergalactic
from astropy.coordinates.matrix_utilities import rotation_matrix

IDENTITY = np.eye(3)


def iau1976(equinox):
    return erfa.pmat76(2451545.0 + (equinox - 2000.0) * 365.25, 0.0)


def iau1980_obliquity(equinox):
    return np.degrees(erfa.obl80(2451545.0 + (equinox - 2000.0) * 365.25, 0.0))


def newcomb_obliquity(equinox):
    t = (equinox - 1950.0) * 1.00002135903e-2
    return (84404.836 - (46.8495 + (0.00319 + 0.00181 * t) * t) * t) / 3600.0


# name: (system, equinox, obstime, astropy frame, rotation from that frame's coordinates to the system's)
SYSTEMS = {
    "icrs": ("icrs", 0.0, 0.0, ICRS(), IDENTITY),
    "fk5_2000": ("fk5", 2000.0, 0.0, FK5(equinox="J2000"), IDENTITY),
    "fk5_1975": ("fk5", 1975.0, 0.0, FK5(equinox="J2000"), iau1976(1975.0)),
    "fk5_2025.5": ("fk5", 2025.5, 0.0, FK5(equinox="J2000"), iau1976(2025.5)),
    "fk4_1950": ("fk4", 1950.0, 1950.0, FK4(equinox="B1950", obstime="B1950"), IDENTITY),
    "fk4_1975": ("fk4", 1975.0, 1975.0, FK4(equinox="B1975", obstime="B1975"), IDENTITY),
    "fk4_1950_obs1990": ("fk4", 1950.0, 1990.0, FK4(equinox="B1950", obstime="B1990"), IDENTITY),
    "fk4noeterms_1950": ("fk4noeterms", 1950.0, 1950.0, FK4NoETerms(equinox="B1950", obstime="B1950"), IDENTITY),
    "fk4noeterms_1975_obs2000": ("fk4noeterms", 1975.0, 2000.0, FK4NoETerms(equinox="B1975", obstime="B2000"), IDENTITY),
    "galactic": ("galactic", 0.0, 0.0, Galactic(), IDENTITY),
    "supergalactic": ("supergalactic", 0.0, 0.0, Supergalactic(), IDENTITY),
    "ecliptic_fk5_2000": ("ecliptic_fk5", 2000.0, 0.0, FK5(equinox="J2000"), rotation_matrix(iau1980_obliquity(2000.0), "x")),
    "ecliptic_fk5_1975": ("ecliptic_fk5", 1975.0, 0.0, FK5(equinox="J2000"), rotation_matrix(iau1980_obliquity(1975.0), "x") @ iau1976(1975.0)),
    "ecliptic_fk4_1950": ("ecliptic_fk4", 1950.0, 1950.0, FK4NoETerms(equinox="B1950", obstime="B1950"), rotation_matrix(newcomb_obliquity(1950.0), "x")),
    "ecliptic_fk4_2010": ("ecliptic_fk4", 2010.0, 2010.0, FK4NoETerms(equinox="B2010", obstime="B2010"), rotation_matrix(newcomb_obliquity(2010.0), "x")),
    "ecliptic_fk4_1950_obs1985": ("ecliptic_fk4", 1950.0, 1985.0, FK4NoETerms(equinox="B1950", obstime="B1985"), rotation_matrix(newcomb_obliquity(1950.0), "x")),
    "ecliptic_icrs": ("ecliptic_icrs", 0.0, 0.0, BarycentricMeanEcliptic(equinox="J2000"), IDENTITY),
}


def unit_vectors(lon, lat):
    lon, lat = np.radians(lon), np.radians(lat)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def lon_lat(vectors):
    lon = np.mod(np.degrees(np.arctan2(vectors[1], vectors[0])), 360.0)
    lon[lon >= 360.0] -= 360.0
    return lon, np.degrees(np.arctan2(vectors[2], np.hypot(vectors[0], vectors[1])))


def convert(source, target, lon, lat):
    _, _, _, source_frame, source_rotation = SYSTEMS[source]
    _, _, _, target_frame, target_rotation = SYSTEMS[target]
    base_lon, base_lat = lon_lat(source_rotation.T @ unit_vectors(lon, lat))
    # A frame, unlike a SkyCoord, does not carry the source's obstime into a target frame that sets its own.
    coords = SkyCoord(base_lon * u.deg, base_lat * u.deg, frame=source_frame).frame
    if isinstance(source_frame, FK4):
        coords = coords.transform_to(FK4NoETerms(equinox=source_frame.equinox, obstime=source_frame.obstime))
    if isinstance(source_frame, (FK4, FK4NoETerms)) or isinstance(target_frame, (FK4, FK4NoETerms)):
        coords = coords.transform_to(FK5(equinox="J2000"))
    return lon_lat(target_rotation @ coords.transform_to(target_frame).cartesian.xyz.value)


def main():
    lon = np.array([0.0, 359.9, 0.1, 83.8221, 150.0, 266.4, 180.0, 10.0, 300.0, 45.0])
    lat = np.array([0.0, -10.0, 10.0, -5.3911, 30.0, -28.9, 89.99, -89.99, 60.0, -45.0])
    converted, default = {}, {}
    for source, (_, _, _, source_frame, source_rotation) in SYSTEMS.items():
        converted[source], default[source] = {}, {}
        for target, (_, _, _, target_frame, target_rotation) in SYSTEMS.items():
            target_lon, target_lat = convert(source, target, lon, lat)
            converted[source][target] = {"lon": target_lon, "lat": target_lat}
            if source_rotation is IDENTITY and target_rotation is IDENTITY:
                result = SkyCoord(lon * u.deg, lat * u.deg, frame=source_frame).transform_to(target_frame).spherical
                default[source][target] = {"lon": result.lon.wrap_at(360 * u.deg).deg, "lat": result.lat.deg}
    systems = {name: {"system": system, "equinox": equinox, "obstime": obstime} for name, (system, equinox, obstime, _, _) in SYSTEMS.items()}
    tree = {"points": {"lon": lon, "lat": lat}, "systems": systems, "converted": converted, "default": default}
    path = Path(__file__).parents[1] / "test" / "fixtures" / "celestial_conversions.asdf"
    asdf.AsdfFile(tree).write_to(path, all_array_storage="inline")
    print(path)


if __name__ == "__main__":
    main()
