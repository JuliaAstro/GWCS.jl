"""Python GWCS pipelines and transforms whose angle and dimensionless parameters are ASDF quantities.

`wcs_from_fiducial` writes its rotation angles as degree quantities and keeps the fiducial position in its ICRS
frame. `mixed` uses rad, arcmin, deg, and dimensionless quantities; `plain` is the same pipeline with plain numbers in
degrees, so GWCS.jl can check that units do not change the WCS type. `galactic` keeps a non-ICRS astropy frame. Each
standalone transform covers the remaining quantity parameters. The references come from evaluating each in Python.
"""

from pathlib import Path

import asdf
import numpy as np
from astropy import units as u
from astropy.coordinates import ICRS, Galactic, SkyCoord
from astropy.modeling import models
from gwcs import WCS
from gwcs import coordinate_frames as cf
from gwcs.wcstools import wcs_from_fiducial


def sky(frame=None):
    return cf.CelestialFrame(reference_frame=ICRS() if frame is None else frame, unit=(u.deg, u.deg), name="world")


def chain(angle, mu, gamma, lon, lat, lon_pole):
    # A named input frame keeps Python GWCS from attaching pixel units, which the unitless shifts cannot take.
    return (models.Shift(-10.0) & models.Shift(-20.0) | models.Rotation2D(angle) | models.Scale(1e-3) & models.Scale(1e-3)
            | models.Pix2Sky_AZP(mu=mu, gamma=gamma) | models.RotateNative2Celestial(lon, lat, lon_pole))


def main():
    wcss = {
        "fiducial": wcs_from_fiducial(SkyCoord(30 * u.deg, 45 * u.deg), projection=models.Pix2Sky_TAN()),
        "mixed": WCS(chain(0.3 * u.rad, 2 * u.one, 30 * u.deg, 1800 * u.arcmin, 45 * u.deg, np.pi * u.rad), input_frame="detector", output_frame=sky()),
        "plain": WCS(chain(np.rad2deg(0.3), 2.0, 30.0, 30.0, 45.0, 180.0), input_frame="detector", output_frame=sky()),
        "galactic": WCS(models.Scale(1e-3) & models.Scale(1e-3) | models.Pix2Sky_TAN() | models.RotateNative2Celestial(30.0, 45.0, 180.0),
            input_frame="detector", output_frame=sky(Galactic())),
    }
    transforms = {
        "szp": models.Pix2Sky_SZP(mu=0.5 * u.one, phi0=170 * u.deg, theta0=60 * u.deg),
        "sin": models.Pix2Sky_SIN(xi=0.1 * u.one, eta=-0.05 * u.one),
        # astropy passes projection angles to WCSLIB without converting their units, so Python's references are right
        # only in degrees.
        "airy": models.Pix2Sky_Airy(theta_b=60 * u.deg),
        "rotation": models.Rotation2D(2700 * u.arcsec),
        "euler": models.EulerAngleRotation(0.1 * u.rad, 20 * u.deg, 30 * u.arcmin, "zxz"),
    }
    wcs_inputs = {"fiducial": np.array([[0.0, 0.0], [0.01, -0.02]]), "mixed": np.array([[100.0, 200.0], [-50.0, 30.0]])}
    wcs_inputs["plain"] = wcs_inputs["mixed"]
    wcs_inputs["galactic"] = np.array([[10.0, -20.0]])
    transform_inputs = {"szp": np.array([[1.0, 2.0], [-3.0, 0.5]]), "sin": np.array([[1.0, 2.0], [-3.0, 0.5]]), "airy": np.array([[1.0, 2.0], [-3.0, 0.5]]),
        "rotation": np.array([[1.5, -2.0]]), "euler": np.array([[10.0, 20.0], [-30.0, 40.0]])}
    tree = {"wcs": wcss, "transforms": transforms, "wcs_inputs": wcs_inputs, "transform_inputs": transform_inputs,
        "wcs_outputs": {name: np.column_stack(w(*wcs_inputs[name].T)) for name, w in wcss.items()},
        "transform_outputs": {name: np.column_stack(t(*transform_inputs[name].T)) for name, t in transforms.items()}}
    path = Path(__file__).parents[1] / "test" / "fixtures" / "quantities.asdf"
    asdf.AsdfFile(tree).write_to(path, all_array_storage="inline")
    print(path)


if __name__ == "__main__":
    main()
