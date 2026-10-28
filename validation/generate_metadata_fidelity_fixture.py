"""A Python GWCS pipeline with metadata GWCS.jl does not model, which a GWCS.jl rewrite must keep.

`constrained` has parameter constraints (`fixed`, `bounds`) and `input_units_equivalencies`, which Python GWCS reads and
writes, a frame2d input frame, a two-axis spatial `frame`, and a composite output frame. Python GWCS writes no other
properties, so after writing, the raw tree gets unknown properties on the WCS, a step, frames, a composite constituent,
and transforms, one of them an array stored in a binary block, plus an `axes_type` on the frame2d, which Python GWCS
ignores there. The composite frame itself gets none, since Python GWCS rejects a composite frame with other properties.
The references come from evaluating the pipeline in Python.
"""

import io
from pathlib import Path

import asdf
import numpy as np
from astropy import units as u
from astropy.modeling import models
from gwcs import WCS
from gwcs import coordinate_frames as cf


def constraints(model):
    return [{"fixed": {k: v for k, v in leaf.fixed.items() if v}, "bounds": {k: list(v) for k, v in leaf.bounds.items() if any(b is not None for b in v)}}
        for leaf in model.traverse_postorder() if leaf.n_submodels == 1]


def main():
    detector = cf.Frame2D(name="detector")
    focal = cf.CoordinateFrame(naxes=2, axes_type=("SPATIAL", "SPATIAL"), axes_order=(0, 1), name="focal", axes_names=("u", "v"), unit=(u.pix, u.pix))
    world = cf.CompositeFrame([cf.SpectralFrame(unit=(u.um,), axes_order=(0,), axes_names=("wavelength",), name="spectral"),
        cf.CoordinateFrame(naxes=1, axes_type=("SPATIAL",), axes_order=(1,), axes_names=("slit",), unit=(u.arcsec,), name="slit")], name="world")
    offset = models.Shift(1.0, fixed={"offset": True}, bounds={"offset": (0.0, 10.0)})
    factor = models.Scale(2.0, bounds={"factor": (1.0, None)})
    factor.input_units_equivalencies = {"x": u.spectral()}
    constrained = WCS([(detector, offset & models.Shift(2.0)), (focal, factor & models.Scale(3.0)), (world, None)])
    constrained.bounding_box = ((-0.5, 99.5), (-0.5, 49.5))

    buffer = io.BytesIO()
    asdf.AsdfFile({"wcs": {"constrained": constrained}}).write_to(buffer, all_array_storage="inline")
    buffer.seek(0)
    with asdf.open(buffer, _force_raw_types=True, lazy_load=False, memmap=False) as af:
        raw = af.tree["wcs"]["constrained"]
        steps = raw["steps"]
        raw["note"] = "kept by GWCS.jl"
        steps[0]["note"] = {"step": 1}
        steps[0]["frame"]["axes_type"] = ["SPECTRAL", "SPATIAL"]
        steps[0]["frame"]["note"] = [1, 2]
        steps[1]["frame"]["note"] = "focal"
        steps[2]["frame"]["frames"][1]["note"] = "slit"
        steps[0]["transform"]["note"] = "concatenate"
        steps[0]["transform"]["forward"][0]["note"] = np.arange(6.0).reshape(2, 3)
        tree = {"wcs": {"constrained": raw}}

    wcs_inputs = {"constrained": np.array([[0.0, 0.0], [20.0, 40.0], [99.0, 3.0]])}
    tree.update({"transforms": {}, "wcs_inputs": wcs_inputs, "transform_inputs": {}, "transform_outputs": {},
        "wcs_outputs": {"constrained": np.column_stack(constrained(*wcs_inputs["constrained"].T))},
        "constraints": {"constrained": constraints(constrained.forward_transform)}})
    path = Path(__file__).parents[1] / "test" / "fixtures" / "metadata_fidelity.asdf"
    output = asdf.AsdfFile(tree)
    for values in (*wcs_inputs.values(), *tree["wcs_outputs"].values()):
        output.set_array_storage(values, "inline")
    output.write_to(path)
    with asdf.open(path) as af:
        w = af["wcs"]["constrained"]
        assert constraints(w.forward_transform) == tree["constraints"]["constrained"], "Python reads the constraints back differently"
        assert type(w.input_frame) is cf.Frame2D and w.input_frame.axes_type == ("SPATIAL", "SPATIAL")
        assert np.allclose(np.column_stack(w(*wcs_inputs["constrained"].T)), tree["wcs_outputs"]["constrained"])
    print(path)


if __name__ == "__main__":
    main()
