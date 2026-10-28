"""Read a GWCS.jl rewrite of a test fixture and evaluate it against the fixture's Python results."""

import json
import sys
import warnings

import asdf
import numpy as np


def parameter_units(model):
    return {name: str(getattr(model, name).unit) for name in model.param_names if getattr(model, name).unit is not None}


def axis_metadata(wcs):
    # Some fixture WCSs have no input frame object.
    return tuple(getattr(frame, "axes_names", None) for frame in (wcs.input_frame, wcs.output_frame)) + (wcs.world_axis_physical_types,)


def main():
    warnings.simplefilter("ignore")
    rewrite, fixture = asdf.open(sys.argv[1]), asdf.open(sys.argv[2])
    records = []
    for group, inputs, outputs in (("wcs", "wcs_inputs", "wcs_outputs"), ("transforms", "transform_inputs", "transform_outputs")):
        for name, value in rewrite[group].items():
            computed = np.column_stack(value(*np.asarray(fixture[inputs][name]).T))
            model = value.forward_transform if group == "wcs" else value
            record = {"name": f"{group}/{name}", "error": float(np.max(np.abs(computed - fixture[outputs][name]))), "units": parameter_units(model)}
            if group == "wcs":
                record["axes_match"] = axis_metadata(value) == axis_metadata(fixture[group][name])
            if group == "wcs" and getattr(value.output_frame, "reference_frame", None) is not None:
                frame = value.output_frame.reference_frame
                record["reference_frame"] = type(frame).__name__
                record["frame_data"] = frame.has_data
            records.append(record)
    print(json.dumps(records))


if __name__ == "__main__":
    main()
