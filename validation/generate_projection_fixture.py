"""Python projection models for every ASDF projection tag beyond the zenithal family, in both directions.

Each projection is evaluated at plane coordinates (pix2sky) or native coordinates (sky2pix) inside its domain. COE
gives its angles as degree quantities. `wcs_from_fiducial` pipelines with CAR and Bonne projections cover the tags
inside a Python GWCS WCS. The references come from evaluating each in Python.
"""

from pathlib import Path

import asdf
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.modeling import models
from gwcs.wcstools import wcs_from_fiducial

PROJECTIONS = {
    "car": ("PlateCarree", {}),
    "cea": ("CylindricalEqualArea", {"lam": 0.7}),
    "cyp": ("CylindricalPerspective", {"mu": 1.5, "lam": 0.8}),
    "mer": ("Mercator", {}),
    "sfl": ("SansonFlamsteed", {}),
    "par": ("Parabolic", {}),
    "mol": ("Molleweide", {}),
    "ait": ("HammerAitoff", {}),
    "pco": ("Polyconic", {}),
    "cop": ("ConicPerspective", {"sigma": 45.0, "delta": 10.0}),
    "cod": ("ConicEquidistant", {"sigma": -30.0, "delta": 15.0}),
    "coe": ("ConicEqualArea", {"sigma": 45 * u.deg, "delta": 10 * u.deg}),
    "coo": ("ConicOrthomorphic", {"sigma": 45.0, "delta": 10.0}),
    "bon": ("BonneEqualArea", {"theta1": 30.0}),
    "tsc": ("TangentialSphericalCube", {}),
    "csc": ("COBEQuadSphericalCube", {}),
    "qsc": ("QuadSphericalCube", {}),
    "hpx": ("HEALPix", {"H": 6.0, "X": 4.0}),
    "xph": ("HEALPixPolar", {}),
}

PLANE = np.array([[10.0, 20.0], [-30.0, 5.0], [0.5, -0.25]])
NATIVE = np.array([[10.0, 20.0], [-30.0, 45.0], [100.0, -30.0]])


def main():
    transforms, inputs = {}, {}
    for name, (model, parameters) in PROJECTIONS.items():
        for direction, prefix, points in (("pix2sky", "Pix2Sky", PLANE), ("sky2pix", "Sky2Pix", NATIVE)):
            transforms[f"{name}_{direction}"] = getattr(models, f"{prefix}_{model}")(**parameters)
            inputs[f"{name}_{direction}"] = points
    fiducial = SkyCoord(30 * u.deg, 45 * u.deg)
    wcss = {"car": wcs_from_fiducial(fiducial, projection=models.Pix2Sky_CAR()), "bon": wcs_from_fiducial(fiducial, projection=models.Pix2Sky_BonneEqualArea(30.0))}
    wcs_inputs = {name: np.array([[0.0, 0.0], [1.5, -2.0]]) for name in wcss}
    tree = {"wcs": wcss, "transforms": transforms, "wcs_inputs": wcs_inputs, "transform_inputs": inputs,
        "wcs_outputs": {name: np.column_stack(w(*wcs_inputs[name].T)) for name, w in wcss.items()},
        "transform_outputs": {name: np.column_stack(t(*inputs[name].T)) for name, t in transforms.items()}}
    for group in ("wcs_outputs", "transform_outputs"):
        for name, values in tree[group].items():
            assert np.all(np.isfinite(values)), f"{group}/{name} has non-finite references"
    path = Path(__file__).parents[1] / "test" / "fixtures" / "projections.asdf"
    asdf.AsdfFile(tree).write_to(path, all_array_storage="inline")
    print(path)


if __name__ == "__main__":
    main()
