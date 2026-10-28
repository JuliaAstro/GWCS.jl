"""Python GWCS pipelines whose frames reorder their axes with axes_order.

In Python GWCS, coordinates are in transform order and a frame's axes_order says which transform position each of its
own axes takes; the frame reports its axis names, physical types, and units in transform order. The file stores each
frame's lists in the frame's own order, as Python GWCS reads them and, since gwcs#687, writes them. `celestial` has a
permuted input frame, a celestial frame with latitude first, a bounding box, and a pixel shape. `composite` puts a
spectral axis between longitude and latitude. The references come from evaluating each in Python.
"""

from pathlib import Path

import asdf
import numpy as np
from astropy import units as u
from astropy.coordinates import ICRS
from astropy.modeling import models
from gwcs import WCS
from gwcs import coordinate_frames as cf
from gwcs.converters.wcs import FrameConverter

# Releases before gwcs#687 write the lists sorted into transform order, which their reader takes as the frame's order.
if not hasattr(cf.CoordinateFrame, "raw_properties"):
    _to_yaml_tree = FrameConverter._to_yaml_tree

    def _frame_order_tree(self, frame, tag, ctx):
        node = _to_yaml_tree(self, frame, tag, ctx)
        for key in ("axes_type", "axes_names", "unit", "axis_physical_types"):
            if key in node:
                node[key] = list(getattr(frame._prop, key))
        return node

    FrameConverter._to_yaml_tree = _frame_order_tree


def axis_metadata(w):
    return {"pixel_axis_names": list(w.input_frame.axes_names), "world_axis_names": list(w.output_frame.axes_names),
        "world_axis_physical_types": list(w.world_axis_physical_types), "world_axis_units": [str(unit) for unit in w.output_frame.unit]}


def main():
    detector = cf.Frame2D(name="detector", axes_names=("x", "y"), axes_order=(1, 0))
    sky = cf.CelestialFrame(reference_frame=ICRS(), unit=(u.deg, u.deg), axes_order=(1, 0), name="world")
    celestial = WCS(models.Shift(10.0) & models.Shift(20.0) | models.Scale(1e-3) & models.Scale(2e-3), input_frame=detector, output_frame=sky)
    celestial.bounding_box = ((-0.5, 99.5), (-0.5, 49.5))
    celestial.pixel_shape = (100, 50)

    pixels = cf.CoordinateFrame(naxes=3, axes_type=("PIXEL",) * 3, axes_order=(0, 1, 2), name="pixels", axes_names=("i", "j", "k"), unit=(u.pix,) * 3)
    world = cf.CompositeFrame([cf.CelestialFrame(reference_frame=ICRS(), unit=(u.deg, u.deg), axes_order=(0, 2), name="sky"),
        cf.SpectralFrame(unit=(u.um,), axes_order=(1,), axes_names=("wavelength",), name="spectral")], name="world")
    composite = WCS(models.Shift(1.0) & models.Shift(2.0) & models.Shift(3.0), input_frame=pixels, output_frame=world)

    wcss = {"celestial": celestial, "composite": composite}
    wcs_inputs = {"celestial": np.array([[0.0, 0.0], [20.0, 40.0], [99.0, 3.0]]), "composite": np.array([[0.0, 0.0, 0.0], [1.5, -2.0, 4.0]])}
    tree = {"wcs": wcss, "transforms": {}, "wcs_inputs": wcs_inputs, "transform_inputs": {}, "transform_outputs": {},
        "wcs_outputs": {name: np.column_stack(w(*wcs_inputs[name].T)) for name, w in wcss.items()},
        # Inside and outside the bounding box, which applies to the transform inputs whatever the frame order.
        "outside": {"celestial": np.column_stack(celestial(*np.array([[99.0, 49.0], [49.0, 99.0]]).T))},
        "axes": {name: axis_metadata(w) for name, w in wcss.items()}}
    for name, values in tree["wcs_outputs"].items():
        assert np.all(np.isfinite(values)), f"wcs_outputs/{name} has non-finite references"
    path = Path(__file__).parents[1] / "test" / "fixtures" / "axes_order.asdf"
    asdf.AsdfFile(tree).write_to(path, all_array_storage="inline")
    with asdf.open(path) as af:
        for name, w in af["wcs"].items():
            assert axis_metadata(w) == tree["axes"][name], f"Python reads wcs/{name} back with different axis metadata"
    print(path)


if __name__ == "__main__":
    main()
