from pathlib import Path

import asdf
import astropy.units as u
from astropy.modeling import bind_compound_bounding_box, models
from astropy.modeling.bounding_box import ModelBoundingBox
from gwcs import WCS
from gwcs.coordinate_frames import CoordinateFrame


def frames():
    kwargs = {
        "naxes": 3,
        "axes_type": ("PIXEL", "PIXEL", "PIXEL"),
        "axes_order": (0, 1, 2),
        "unit": (u.pix, u.pix, u.pix),
        "axes_names": ("x", "y", "slice"),
    }
    return CoordinateFrame(name="detector", **kwargs), CoordinateFrame(name="output", **kwargs)


def simple_wcs():
    input_frame, output_frame = frames()
    transform = models.Identity(3)
    transform.inputs = ("x", "y", "slice")
    transform.bounding_box = ModelBoundingBox(
        {"x": (0.0, 5.0), "y": (10.0, 15.0)},
        model=transform,
        ignored=["slice"],
        order="F",
    )
    return WCS(transform, input_frame=input_frame, output_frame=output_frame, name="simple")


def compound_wcs():
    input_frame, output_frame = frames()
    transform = models.Identity(3)
    transform.inputs = ("x", "y", "slice")
    bind_compound_bounding_box(
        transform,
        {
            (1,): {"x": (0.0, 5.0), "y": (0.0, 5.0)},
            (2,): {"x": (10.0, 15.0), "y": (10.0, 15.0)},
        },
        [("slice",)],
        order="F",
    )
    return WCS(transform, input_frame=input_frame, output_frame=output_frame, name="compound")


def main():
    path = Path(__file__).parents[1] / "test" / "fixtures" / "bounding_boxes.asdf"
    asdf.AsdfFile({"simple": simple_wcs(), "compound": compound_wcs()}).write_to(path, all_array_storage="inline")
    print(path)


if __name__ == "__main__":
    main()
