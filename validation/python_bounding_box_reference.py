import copy
import json

import numpy as np
from astropy.modeling import bind_compound_bounding_box, models
from astropy.modeling.bounding_box import ModelBoundingBox
from gwcs import coordinate_frames as cf
from gwcs import wcs
from gwcs.wcstools import grid_from_bounding_box


def identity_wcs(transform, attach_bounding_box=False):
    input_frame = cf.CoordinateFrame(
        naxes=3,
        axes_type=("PIXEL", "PIXEL", "PIXEL"),
        axes_order=(0, 1, 2),
        axes_names=("x", "y", "slice"),
        name="input",
    )
    output_frame = cf.CoordinateFrame(
        naxes=3,
        axes_type=("CUSTOM", "CUSTOM", "CUSTOM"),
        axes_order=(0, 1, 2),
        name="output",
    )
    result = wcs.WCS([(input_frame, copy.deepcopy(transform)), (output_frame, None)])
    if attach_bounding_box:
        result.bounding_box = copy.deepcopy(transform.bounding_box)
    return result


def in_image_result(transform, *world):
    try:
        return bool(transform.in_image(*world))
    except Exception as error:
        return type(error).__name__


def main():
    simple = models.Identity(3)
    simple.inputs = ("x", "y", "slice")
    simple.bounding_box = ModelBoundingBox(
        {"x": (0.0, 5.0), "y": (10.0, 15.0)},
        model=simple,
        ignored=["slice"],
        order="F",
    )

    compound = models.Identity(3)
    compound.inputs = ("x", "y", "slice")
    bind_compound_bounding_box(
        compound,
        {
            (1,): {"x": (0.0, 5.0), "y": (0.0, 5.0)},
            (2,): {"x": (10.0, 15.0), "y": (10.0, 15.0)},
        },
        [("slice",)],
        order="F",
    )

    dense = models.Scale(2) & models.Scale(3) & models.Scale(1)
    dense.inputs = ("x", "y", "slice")
    dense.bounding_box = ModelBoundingBox(
        {"x": (0.0, 5.0), "y": (10.0, 15.0), "slice": (0.0, 200.0)},
        model=dense,
        order="F",
    )

    compound_mask = models.Scale(2) & models.Scale(3) & models.Scale(1)
    compound_mask.inputs = ("x", "y", "slice")
    bind_compound_bounding_box(
        compound_mask,
        {
            (1,): {"x": (0.0, 5.0), "y": (0.0, 5.0)},
            (2,): {"x": (10.0, 15.0), "y": (10.0, 15.0)},
        },
        [("slice",)],
        order="F",
    )

    simple_grid = grid_from_bounding_box(((-1.0, 2.9), (6.0, 7.5)), step=(1.0, 0.5), center=False)
    compound_grid = grid_from_bounding_box(compound.bounding_box, selector=(1,))
    dense_wcs = identity_wcs(dense, attach_bounding_box=True)
    simple_wcs = identity_wcs(simple)
    compound_wcs = identity_wcs(compound_mask, attach_bounding_box=True)

    try:
        compound(2.0, 3.0, 3.0, with_bounding_box=True)
    except RuntimeError as error:
        missing_selector = type(error).__name__
    else:
        missing_selector = "no error"

    def serializable(values):
        return [None if np.isnan(value) else float(value) for value in values]

    result = {
        "simple_inside": list(simple(2.0, 12.0, 100.0, with_bounding_box=True)),
        "simple_outside": serializable(simple(7.0, 12.0, 100.0, with_bounding_box=True)),
        "simple_fill": list(simple(7.0, 12.0, 100.0, with_bounding_box=True, fill_value=-1.0)),
        "compound_inside": list(compound(2.0, 3.0, 1.0, with_bounding_box=True)),
        "compound_outside": serializable(compound(2.0, 3.0, 2.0, with_bounding_box=True)),
        "missing_selector": missing_selector,
        "dense_in_image": [
            bool(dense_wcs.in_image(4.0, 36.0, 100.0)),
            bool(dense_wcs.in_image(14.0, 36.0, 100.0)),
            bool(dense_wcs.in_image(14.0, 36.0, 100.0, with_bounding_box=False)),
        ],
        "dense_in_image_batch": np.asarray(dense_wcs.in_image([4.0, 14.0, 4.0], [36.0, 36.0, 36.0], [100.0, 100.0, np.nan])).tolist(),
        "ignored_in_image_outside": bool(simple_wcs.in_image(7.0, 12.0, 100.0)),
        "compound_in_image": [
            in_image_result(compound_wcs, 4.0, 9.0, 1.0),
            in_image_result(compound_wcs, 4.0, 9.0, 2.0),
        ],
        "compound_missing_in_image": in_image_result(compound_wcs, 4.0, 9.0, 3.0),
        "simple_grid": np.asarray(simple_grid).reshape(2, -1).tolist(),
        "compound_grid": np.asarray(compound_grid).reshape(2, -1).tolist(),
    }
    print(json.dumps(result))


if __name__ == "__main__":
    main()
