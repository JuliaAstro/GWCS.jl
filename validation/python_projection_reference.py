#!/usr/bin/env python3

import json

from astropy.modeling import models


PROJECTIONS = (
    ("TAN", models.Pix2Sky_Gnomonic(), models.Sky2Pix_Gnomonic()),
    ("STG", models.Pix2Sky_Stereographic(), models.Sky2Pix_Stereographic()),
    ("ARC", models.Pix2Sky_ZenithalEquidistant(), models.Sky2Pix_ZenithalEquidistant()),
    ("ZEA", models.Pix2Sky_ZenithalEqualArea(), models.Sky2Pix_ZenithalEqualArea()),
    ("SIN", models.Pix2Sky_SlantOrthographic(xi=0.1, eta=0.2), models.Sky2Pix_SlantOrthographic(xi=0.1, eta=0.2)),
    ("AZP", models.Pix2Sky_ZenithalPerspective(mu=1.5, gamma=15.0), models.Sky2Pix_ZenithalPerspective(mu=1.5, gamma=15.0)),
    ("SZP", models.Pix2Sky_SlantZenithalPerspective(mu=1.5, phi0=15.0, theta0=80.0), models.Sky2Pix_SlantZenithalPerspective(mu=1.5, phi0=15.0, theta0=80.0)),
    ("AIR", models.Pix2Sky_Airy(theta_b=75.8), models.Sky2Pix_Airy(theta_b=75.8)),
)

INTERMEDIATE_POINTS = ((0.0, 0.0), (1.0, 2.0), (-10.0, 5.0), (20.0, -15.0))
NATIVE_POINTS = ((0.0, 90.0), (30.0, 60.0), (-120.0, 45.0), (75.0, 80.0))


def values(result):
    return [float(value) for value in result]


def main():
    records = []
    for name, pix2sky, sky2pix in PROJECTIONS:
        records.append(
            {
                "name": name,
                "intermediate": [
                    {"input": point, "output": values(pix2sky(*point))}
                    for point in INTERMEDIATE_POINTS
                ],
                "native": [
                    {"input": point, "output": values(sky2pix(*point))}
                    for point in NATIVE_POINTS
                ],
            }
        )
    print(json.dumps(records))


if __name__ == "__main__":
    main()
