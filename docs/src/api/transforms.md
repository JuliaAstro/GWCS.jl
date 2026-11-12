# Transform API

```@meta
CurrentModule = GWCS
```

## The transform protocol

```@docs
AbstractTransform
inverse_transform
has_inverse
NoInverse
transform_jacobian
transform_dependencies
compile_transform
input_dimension
output_dimension
transform_name
input_names
output_names
```

## Composition and annotation

```@docs
compose
parallel
elementwise
AnnotatedTransform
annotate_transform
WithInverse
```

## Elementary transforms

```@docs
IdentityTransform
Shift
Scale
AffineTransform
RemapAxes
ConstantTransform
Rotation2D
Rotation3D
Polynomial1D
Polynomial2D
degree
CallableTransform
```

## Spherical coordinates

```@docs
SphericalToCartesian
CartesianToSpherical
NativeToCelestial
CelestialToNative
IntermediateToNative
NativeToIntermediate
celestial_transform
```

## Projections

```@docs
GnomonicProjection
StereographicProjection
ZenithalEquidistant
ZenithalEqualArea
ZenithalPerspective
SlantZenithalPerspective
SlantOrthographic
AiryProjection
ZenithalPolynomial
PlateCarreeProjection
CylindricalEqualArea
CylindricalPerspective
MercatorProjection
SansonFlamsteedProjection
ParabolicProjection
MollweideProjection
HammerAitoffProjection
ConicPerspective
ConicEquidistant
ConicEqualArea
ConicOrthomorphic
PolyconicProjection
BonneProjection
TangentialSphericalCube
COBESphericalCube
QuadrilateralizedSphericalCube
HEALPixProjection
HEALPixPolarProjection
```
