import json
import sys

import asdf
import numpy as np


def main():
    with asdf.open(sys.argv[1]) as handle:
        records = []
        for name, transform in handle["transforms"].items():
            inputs = np.asarray(handle["inputs"][name], dtype=float)
            outputs = transform(*inputs.T)
            outputs = np.atleast_2d(np.column_stack(outputs if isinstance(outputs, tuple) else (outputs,)))
            records.append({"transform": name, "type": type(transform).__name__, "outputs": outputs.tolist()})
    print(json.dumps(records))


if __name__ == "__main__":
    main()
