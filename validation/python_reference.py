#!/usr/bin/env python3

import argparse
import json

import asdf


POINTS = (
    (0.0, 0.0),
    (2043.5, 2043.5),
    (4087.0, 4087.0),
    (1234.5, 2345.5),
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()

    records = []
    for filename in args.files:
        with asdf.open(filename, lazy_load=True, memmap=False) as asdf_file:
            wcs = asdf_file.tree["roman"]["meta"]["wcs"]
            for pixel in POINTS:
                world = wcs(*pixel)
                roundtrip = wcs.invert(*world)
                records.append(
                    {
                        "file": filename,
                        "pixel": list(pixel),
                        "world": [float(value) for value in world],
                        "roundtrip": [float(value) for value in roundtrip],
                    }
                )

    print(json.dumps(records))


if __name__ == "__main__":
    main()
