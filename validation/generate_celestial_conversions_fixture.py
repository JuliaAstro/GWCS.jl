"""Positions converted between celestial systems by astropy, for `with_celestial_system`.

`points` holds longitudes and latitudes in degrees, including points near the poles and on either side of longitude
zero. `systems` names each system and its FK5 equinox (zero for the others), and `converted[source][target]` holds
astropy's conversion of `points`, read as coordinates in `source`, to `target`, with longitudes in [0, 360).
"""

from pathlib import Path

import asdf
import numpy as np
from astropy import units as u
from astropy.coordinates import FK5, ICRS, Galactic, SkyCoord, Supergalactic

SYSTEMS = {
    "icrs": ("icrs", 0.0, ICRS()),
    "fk5_2000": ("fk5", 2000.0, FK5(equinox="J2000")),
    "fk5_1975": ("fk5", 1975.0, FK5(equinox="J1975")),
    "fk5_2025.5": ("fk5", 2025.5, FK5(equinox="J2025.5")),
    "galactic": ("galactic", 0.0, Galactic()),
    "supergalactic": ("supergalactic", 0.0, Supergalactic()),
}


def main():
    lon = np.array([0.0, 359.9, 0.1, 83.8221, 150.0, 266.4, 180.0, 10.0, 300.0, 45.0])
    lat = np.array([0.0, -10.0, 10.0, -5.3911, 30.0, -28.9, 89.99, -89.99, 60.0, -45.0])
    converted = {}
    for source, (_, _, source_frame) in SYSTEMS.items():
        coords = SkyCoord(lon * u.deg, lat * u.deg, frame=source_frame)
        converted[source] = {}
        for target, (_, _, target_frame) in SYSTEMS.items():
            result = coords.transform_to(target_frame).spherical
            converted[source][target] = {"lon": result.lon.wrap_at(360 * u.deg).deg, "lat": result.lat.deg}
    systems = {name: {"system": system, "equinox": equinox} for name, (system, equinox, _) in SYSTEMS.items()}
    tree = {"points": {"lon": lon, "lat": lat}, "systems": systems, "converted": converted}
    path = Path(__file__).parents[1] / "test" / "fixtures" / "celestial_conversions.asdf"
    asdf.AsdfFile(tree).write_to(path, all_array_storage="inline")
    print(path)


if __name__ == "__main__":
    main()
