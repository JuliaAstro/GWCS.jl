import json

import numpy as np
from astropy.modeling import bind_compound_bounding_box, models
from astropy.modeling.bounding_box import ModelBoundingBox
from gwcs.wcstools import grid_from_bounding_box


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

    simple_grid = grid_from_bounding_box(((-1.0, 2.9), (6.0, 7.5)), step=(1.0, 0.5), center=False)
    compound_grid = grid_from_bounding_box(compound.bounding_box, selector=(1,))

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
        "simple_grid": np.asarray(simple_grid).reshape(2, -1).tolist(),
        "compound_grid": np.asarray(compound_grid).reshape(2, -1).tolist(),
    }
    print(json.dumps(result))


if __name__ == "__main__":
    main()
