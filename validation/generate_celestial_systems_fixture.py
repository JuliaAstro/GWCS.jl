"""Python GWCS celestial frames in each astropy reference frame, with the system and epochs astropy reports.

Each WCS maps two inputs to a `CelestialFrame` with an identity transform. `references` holds astropy's frame name, the
equinox as a Julian epoch for FK5 and PrecessedGeocentric or a Besselian one for FK4, and FK4's obstime as a Besselian
epoch. The obstime of `fk4noeterms_utc` is a UTC date, which GWCS.jl does not read, so it has no reference.
"""

from pathlib import Path

import asdf
from astropy import units as u
from astropy.coordinates import FK4, FK5, GCRS, ICRS, FK4NoETerms, Galactic, PrecessedGeocentric
from astropy.modeling import models
from gwcs import WCS
from gwcs import coordinate_frames as cf

FRAMES = {
    "icrs": ICRS(),
    "fk5_default": FK5(),
    "fk5_1975": FK5(equinox="J1975"),
    "fk5_besselian": FK5(equinox="B1950"),
    "fk4_default": FK4(),
    "fk4_1975": FK4(equinox="B1975"),
    "fk4noeterms": FK4NoETerms(equinox="B1950", obstime="J2000"),
    "fk4noeterms_utc": FK4NoETerms(obstime="2000-01-01"),
    "galactic": Galactic(),
    "gcrs": GCRS(obstime="2020-01-01"),
    "precessedgeocentric": PrecessedGeocentric(equinox="J2010"),
}


def reference(frame):
    result = {"system": frame.name}
    if isinstance(frame, (FK5, PrecessedGeocentric)):
        result["equinox"] = float(frame.equinox.jyear)
    elif isinstance(frame, (FK4, FK4NoETerms)):
        result["equinox"] = float(frame.equinox.byear)
        result["obstime"] = float(frame.obstime.byear)
    return result


def main():
    tree = {"wcs": {}, "references": {}}
    for name, frame in FRAMES.items():
        world = cf.CelestialFrame(reference_frame=frame, unit=(u.deg, u.deg), name="world")
        tree["wcs"][name] = WCS(models.Identity(2), input_frame="detector", output_frame=world)
        if name != "fk4noeterms_utc":
            tree["references"][name] = reference(frame)
    path = Path(__file__).parents[1] / "test" / "fixtures" / "celestial_systems.asdf"
    asdf.AsdfFile(tree).write_to(path)
    print(path)


if __name__ == "__main__":
    main()
