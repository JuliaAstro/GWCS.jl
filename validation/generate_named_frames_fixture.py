"""Python GWCS pipelines whose step frames are bare names, which Python reads as EmptyFrame.

`wcs_from_fiducial` names its input frame "detector" this way, but stores its rotation angles as quantities, so the
sky WCS here uses plain numbers; a WCS built from frame names alone has no frame objects at all. The references come
from evaluating each WCS in Python.
"""

from pathlib import Path

import asdf
import numpy as np
from astropy import units as u
from astropy.coordinates import ICRS
from astropy.modeling import models
from gwcs import WCS
from gwcs import coordinate_frames as cf


def main():
    sky = models.Scale(1e-3) & models.Scale(1e-3) | models.Pix2Sky_TAN() | models.RotateNative2Celestial(30.0, 45.0, 180.0)
    celestial = WCS(sky, input_frame="detector", output_frame=cf.CelestialFrame(reference_frame=ICRS(), unit=(u.deg, u.deg), name="world"))
    named = WCS(models.Shift(1.0) & models.Scale(2.0), input_frame="detector", output_frame="focal")
    celestial_inputs = np.array([[0.0, 0.0], [10.0, -20.0], [-30.0, 40.0]])
    named_inputs = np.array([[0.0, 0.0], [3.0, -4.0]])
    tree = {"celestial": celestial, "celestial_inputs": celestial_inputs, "celestial_outputs": np.column_stack(celestial(*celestial_inputs.T)),
        "named": named, "named_inputs": named_inputs, "named_outputs": np.column_stack(named(*named_inputs.T))}
    path = Path(__file__).parents[1] / "test" / "fixtures" / "named_frames.asdf"
    asdf.AsdfFile(tree).write_to(path, all_array_storage="inline")
    print(path)


if __name__ == "__main__":
    main()
