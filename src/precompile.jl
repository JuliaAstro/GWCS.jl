# A Python GWCS imaging pipeline: a polynomial distortion without an inverse, then a TAN projection and sky rotation.
# GWCSASDFExt parses it too, as ASDF.jl nodes.
const _WORKLOAD_ASDF = """
%YAML 1.1
%TAG ! tag:stsci.edu:asdf/
---
!<tag:stsci.edu:gwcs/wcs-1.4.0>
name: imaging
steps:
  - !<tag:stsci.edu:gwcs/step-1.3.0>
    frame: !<tag:stsci.edu:gwcs/frame2d-1.2.0>
      axes_names: [x, y]
      axes_order: [0, 1]
      name: detector
      unit: [!unit/unit-1.0.0 pixel, !unit/unit-1.0.0 pixel]
    transform: !transform/compose-1.4.0
      bounding_box: !transform/property/bounding_box-1.2.0
        ignore: []
        intervals: {x0: [-0.5, 1023.5], x1: [-0.5, 1023.5]}
        order: F
      forward:
        - !transform/concatenate-1.4.0
          forward: [!transform/shift-1.4.0 {offset: -512.0}, !transform/shift-1.4.0 {offset: -512.0}]
        - !transform/remap_axes-1.5.0 {mapping: [0, 1, 0, 1]}
        - !transform/concatenate-1.4.0
          forward:
            - !transform/polynomial-1.3.0
              coefficients: [[0.0, 0.0, 1e-9], [2.7778e-4, 1e-9, 0.0], [2e-9, 0.0, 0.0]]
            - !transform/polynomial-1.3.0
              coefficients: [[0.0, 2.7778e-4, -1e-9], [0.0, 3e-9, 0.0], [1e-9, 0.0, 0.0]]
        - !transform/gnomonic-1.4.0 {direction: pix2sky}
        - !transform/rotate3d-1.5.0 {direction: native2celestial, phi: 83.8, psi: 180.0, theta: -5.4}
  - !<tag:stsci.edu:gwcs/step-1.3.0>
    frame: !<tag:stsci.edu:gwcs/celestial_frame-1.2.0>
      axes_names: [lon, lat]
      axes_order: [0, 1]
      axis_physical_types: [pos.eq.ra, pos.eq.dec]
      name: world
      reference_frame: !<tag:astropy.org:astropy/coordinates/frames/icrs-1.1.0> {frame_attributes: {}}
      unit: [!unit/unit-1.0.0 deg, !unit/unit-1.0.0 deg]
    transform: null
...
"""

# Caches the type-independent construction code (FITS keyword and ASDF tree interpretation, canonicalization, and the
# solver scaffolding) and the kernels of these common pipelines; other pipeline types still compile their own kernels.
@setup_workload begin
    sky = Dict{String,Any}("NAXIS" => 2, "NAXIS1" => 1024, "NAXIS2" => 1024, "CRPIX1" => 512.0, "CRPIX2" => 512.0, "CRVAL1" => 83.8, "CRVAL2" => -5.4,
        "CD1_1" => -2.7778e-4, "CD1_2" => 3e-6, "CD2_1" => 2e-6, "CD2_2" => 2.7778e-4, "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN")
    sip = merge(sky, Dict{String,Any}("CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "DEC--TAN-SIP", "A_ORDER" => 3, "B_ORDER" => 3, "A_2_0" => 2e-5, "A_1_1" => -1e-5, "A_0_2" => 1.5e-5,
        "A_3_0" => 3e-8, "B_2_0" => -1e-5, "B_1_1" => 2e-5, "B_0_3" => 4e-8))
    tpv = merge(sky, Dict{String,Any}("CTYPE1" => "RA---TPV", "CTYPE2" => "DEC--TPV", "PV1_1" => 1.0, "PV1_4" => 2e-3, "PV2_1" => 1.0, "PV2_6" => 2e-3))
    cube = merge(sky, Dict{String,Any}("NAXIS" => 3, "NAXIS3" => 100, "CTYPE3" => "WAVE", "CRPIX3" => 1.0, "CRVAL3" => 5e-7, "CD3_3" => 1e-10, "CUNIT3" => "m"))
    @compile_workload begin
        for wcs in (parse_fits_wcs(sky), parse_fits_wcs(sip), parse_fits_wcs(tpv), parse_fits_wcs(cube), parse_wcs(_load_tagged_yaml(_WORKLOAD_ASDF)))
            N = input_dimension(wcs.forward)
            pixel = SVector(ntuple(i -> 100.0 * i, N))
            world = pixel_to_world(wcs, pixel)
            pixel_to_world(wcs, Tuple(pixel)...)
            world_to_pixel(wcs, world)
            world_to_pixel(wcs, Tuple(world)...)
            in_image(wcs, world)
            worlds = pixel_to_world(wcs, Matrix(hcat(pixel, pixel)))
            world_to_pixel(wcs, worlds)
            in_image(wcs, worlds)
        end
    end
end
