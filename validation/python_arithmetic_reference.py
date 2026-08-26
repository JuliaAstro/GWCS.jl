import json

from astropy.modeling import models


def main():
    left = models.Shift(1.0)
    right = models.Scale(2.0)
    transforms = {
        "add": left + right,
        "subtract": left - right,
        "multiply": left * right,
        "divide": left / right,
        "power": left**right,
        "nested": (left + right) * (left - right),
    }
    records = []
    for name, transform in transforms.items():
        for coordinate in (0.5, 3.0):
            records.append(
                {
                    "transform": name,
                    "input": coordinate,
                    "output": float(transform(coordinate)),
                }
            )
    print(json.dumps(records))


if __name__ == "__main__":
    main()
