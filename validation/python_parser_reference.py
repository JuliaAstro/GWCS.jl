import json

from astropy.modeling import models


def main():
    transforms = {
        "identity": models.Identity(2),
        "affine": models.AffineTransformation2D(
            matrix=[[1.0, 0.2], [-0.1, 2.0]], translation=[3.0, 4.0]
        ),
        "rotate2d": models.Rotation2D(30.0),
        "gnomonic": models.Pix2Sky_Gnomonic(),
        "native2celestial": models.RotateNative2Celestial(10.0, 20.0, 180.0),
        "celestial2native": models.RotateCelestial2Native(10.0, 20.0, 180.0),
        "euler": models.EulerAngleRotation(1.0, 2.0, 3.0, "xyz"),
        "spherical_sequence": models.SphericalRotationSequence(
            [1.0, 2.0, 3.0], "zyx"
        ),
    }
    points = [(12.0, -30.0), (0.1, 45.0)]
    records = []
    for name, transform in transforms.items():
        for point in points:
            result = tuple(float(value) for value in transform(*point))
            roundtrip = tuple(float(value) for value in transform.inverse(*result))
            records.append(
                {
                    "transform": name,
                    "input": point,
                    "output": result,
                    "roundtrip": roundtrip,
                }
            )
    print(json.dumps(records))


if __name__ == "__main__":
    main()
