"""A GWCS whose only non-invertible step is a sky-space correction after the projection, near RA 0/360.

The numerical inverse cannot solve such a pipeline in intermediate coordinates near the seam, because the corrected
longitude can sit 360 degrees from the inverted tail's longitude; it must fall back to the world-space solve.
"""

from pathlib import Path

import asdf
import numpy as np
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling import models
from gwcs import WCS
from gwcs import coordinate_frames as cf
from gwcs.geometry import CartesianToSpherical, SphericalToCartesian


def sky_correction_wcs():
    linear = models.AffineTransformation2D(matrix=[[-2.7778e-4, 0.0], [0.0, 2.7778e-4]], translation=[512 * 2.7778e-4, -512 * 2.7778e-4])
    sky = models.Pix2Sky_TAN() | models.RotateNative2Celestial(0.0, 10.0, 180.0)
    correction = models.Mapping((0, 1, 0, 1)) | (models.Polynomial2D(2, c1_0=1.0, c0_2=1e-4) & models.Polynomial2D(2, c0_1=1.0, c0_2=2e-6))
    to_other = SphericalToCartesian(wrap_lon_at=360) | models.RotationSequence3D([20.0, -30.0, 40.0], axes_order="zyx") | CartesianToSpherical(wrap_lon_at=360)
    detector = cf.Frame2D(name="detector", axes_names=("x", "y"), unit=(u.pix, u.pix))
    world = cf.CelestialFrame(name="world", reference_frame=coord.ICRS(), unit=(u.deg, u.deg))
    return WCS(linear | sky | correction | to_other, input_frame=detector, output_frame=world, name="sky_correction")


def main():
    wcs = sky_correction_wcs()
    # Columns straddle the RA 0/360 seam at x = 512.
    x, y = np.meshgrid(np.arange(505.0, 521.0), [0.0, 512.0, 1024.0])
    pixels = np.column_stack([x.ravel(), y.ravel()])
    world = np.column_stack(wcs(pixels[:, 0], pixels[:, 1], with_bounding_box=False))
    # Python GWCS 1.0.3 cannot serve as the inverse reference here: without a bounding box its fixed-point
    # iteration starts at pixel (1, 1) and diverges for 9 of these points, and with one it fails to build its
    # approximate inverse under astropy 8. The fixture therefore records the forward map and the true pixels.
    path = Path(__file__).parents[1] / "test" / "fixtures" / "sky_correction.asdf"
    asdf.AsdfFile({"wcs": wcs, "pixels": pixels, "world": world}).write_to(path, all_array_storage="inline")
    print(path)


if __name__ == "__main__":
    main()
