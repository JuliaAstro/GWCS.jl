"""astropy FITS WCSs as asdf-astropy writes them (tag:astropy.org:astropy/wcs/wcs), with astropy's own results.

The tag stores the HDU list of `WCS.to_fits` plus the attributes the header does not hold. The references come from
astropy's APE 14 `pixel_to_world_values`, which is zero-based like GWCS.jl.
"""

from pathlib import Path

import asdf
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def sip_wcs():
    header = fits.Header({"NAXIS": 2, "CTYPE1": "RA---TAN-SIP", "CTYPE2": "DEC--TAN-SIP", "CRPIX1": 512.0, "CRPIX2": 512.0, "CRVAL1": 83.8,
        "CRVAL2": -5.4, "CD1_1": -2.7778e-4, "CD1_2": 1e-6, "CD2_1": 2e-6, "CD2_2": 2.7778e-4, "A_ORDER": 2, "A_2_0": 5e-6, "A_1_1": -2e-6,
        "B_ORDER": 2, "B_0_2": 2e-6, "B_2_0": 1e-6})
    wcs = WCS(header)
    wcs.pixel_shape = (1024, 1024)
    return wcs


def bounded_wcs():
    header = fits.Header({"NAXIS": 2, "CTYPE1": "GLON-CAR", "CTYPE2": "GLAT-CAR", "CRPIX1": 100.5, "CRPIX2": 50.5, "CRVAL1": 30.0, "CRVAL2": 0.0,
        "CDELT1": -0.05, "CDELT2": 0.05})
    wcs = WCS(header)
    wcs.pixel_bounds = [(-0.5, 199.5), (-0.5, 99.5)]
    return wcs


def alternate_wcs():
    # asdf-astropy writes the alternate representation's keywords but not its key, and astropy itself cannot read
    # the file back, so the reference comes from this in-memory WCS.
    header = fits.Header({"NAXIS": 2, "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN", "CRPIX1": 512.0, "CRPIX2": 512.0, "CRVAL1": 83.8,
        "CRVAL2": -5.4, "CDELT1": -2.7778e-4, "CDELT2": 2.7778e-4, "WCSNAMEA": "FOCAL", "CTYPE1A": "XLIN", "CTYPE2A": "YLIN",
        "CRPIX1A": 1.0, "CRPIX2A": 1.0, "CRVAL1A": 0.0, "CRVAL2A": 0.0, "CDELT1A": 0.5, "CDELT2A": 0.5, "CUNIT1A": "mm", "CUNIT2A": "mm"})
    return WCS(header, key="A")


def main():
    sip, bounded, alternate = sip_wcs(), bounded_wcs(), alternate_wcs()
    alternate_pixels = np.array([[0.0, 0.0], [9.0, 19.0]])
    sip_pixels = np.array([[0.0, 0.0], [99.0, 899.0], [511.0, 511.0], [1023.0, 1023.0]])
    bounded_pixels = np.array([[0.0, 0.0], [99.5, 49.5], [199.0, 99.0]])
    tree = {"sip": sip, "sip_pixels": sip_pixels, "sip_world": np.column_stack(sip.pixel_to_world_values(sip_pixels[:, 0], sip_pixels[:, 1])),
        "bounded": bounded, "bounded_pixels": bounded_pixels,
        "bounded_world": np.column_stack(bounded.pixel_to_world_values(bounded_pixels[:, 0], bounded_pixels[:, 1])),
        "alternate": alternate, "alternate_pixels": alternate_pixels,
        "alternate_world": np.column_stack(alternate.pixel_to_world_values(alternate_pixels[:, 0], alternate_pixels[:, 1]))}
    path = Path(__file__).parents[1] / "test" / "fixtures" / "astropy_wcs.asdf"
    asdf.AsdfFile(tree).write_to(path, all_array_storage="inline")
    print(path)


if __name__ == "__main__":
    main()
