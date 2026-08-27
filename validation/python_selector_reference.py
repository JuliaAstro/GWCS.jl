import json
import warnings
import numpy as np

from astropy.modeling import models
from gwcs import selector


def constant_polynomial(value):
    return models.Polynomial2D(1, c0_0=value, c1_0=0, c0_1=0)


def serializable(value):
    array = np.asarray(value)
    if array.ndim == 0:
        scalar = array.item()
        return None if isinstance(scalar, float) and np.isnan(scalar) else scalar
    return [None if np.isnan(item) else float(item) for item in array]


def main():
    warnings.filterwarnings("ignore", message="The input positions are not inside any region.")
    warnings.filterwarnings("ignore", message="All data is outside the valid range")
    labels = np.array([[0, 2], [1, 0]])
    array_mapper = selector.LabelMapperArray(labels)
    branches = {
        1: models.Shift(1) & models.Scale(1),
        2: models.Shift(2) & models.Scale(2),
    }
    regions = selector.RegionsSelector(
        inputs=("x", "y"),
        outputs=("u", "v"),
        selector=branches,
        label_mapper=array_mapper,
        undefined_transform_value=-100,
    )

    dictionary_mapper = selector.LabelMapperDict(
        ("x", "y"),
        {
            -2.0: constant_polynomial(1.0),
            3.0: constant_polynomial(2.0),
        },
        inputs_mapping=models.Mapping((0,), n_inputs=2),
        atol=0.001,
    )
    range_mapper = selector.LabelMapperRange(
        ("x", "y"),
        {
            (0.0, 2.0): constant_polynomial(1.0),
            (2.0, 4.0): constant_polynomial(2.0),
        },
        inputs_mapping=models.Mapping((0,), n_inputs=2),
    )

    result = {
        "array_labels": [
            serializable(array_mapper(0.0, 1.0)),
            serializable(array_mapper(1.0, 0.0)),
            serializable(array_mapper(0.0, 0.0)),
        ],
        "region_outputs": [
            serializable(regions(0.0, 1.0)),
            serializable(regions(1.0, 0.0)),
            serializable(regions(0.0, 0.0)),
        ],
        "dictionary_labels": [
            serializable(dictionary_mapper(-2.0005, 8.0)),
            serializable(dictionary_mapper(3.0, 8.0)),
            serializable(dictionary_mapper(0.0, 8.0)),
        ],
        "range_labels": [
            serializable(range_mapper(1.0, 8.0)),
            serializable(range_mapper(3.0, 8.0)),
            serializable(range_mapper(2.0, 8.0)),
        ],
    }
    print(json.dumps(result))


if __name__ == "__main__":
    main()
