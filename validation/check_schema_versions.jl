using GWCS
using GWCS: default_transform_registry

include("python.jl")

function resource_path(package)
    python = validation_python()
    code = "import importlib.resources; root = importlib.resources.files('$package'); print(root if hasattr(root, '__fspath__') else next(iter(root.iterdir())).parent)"
    return strip(read(`$python -c $code`, String))
end

function available_versions(root, family)
    prefix = family * "-"
    versions = VersionNumber[]
    for path in readdir(root)
        startswith(path, prefix) && endswith(path, ".yaml") || continue
        push!(versions, VersionNumber(path[(length(prefix) + 1):(end - 5)]))
    end
    return Tuple(sort!(versions))
end

function available_versions(roots::Tuple, family)
    versions = VersionNumber[]
    for root in roots
        append!(versions, available_versions(root, family))
    end
    return Tuple(sort!(unique!(versions)))
end

function registered_versions(registry, family)
    registrations = registry.registrations[family]
    length(registrations) == 1 || error("expected one built-in registration for $family")
    versions = only(registrations).versions
    versions === nothing && error("built-in registration for $family uses a version interval")
    return versions
end

registry = default_transform_registry()
transform_root = joinpath(resource_path("asdf_transform_schemas.resources"), "stsci.edu", "schemas")
gwcs_root = joinpath(resource_path("asdf_wcs_schemas.resources"), "schemas", "stsci.edu", "gwcs")

for family in ("compose", "concatenate", "add", "subtract", "multiply", "divide", "power", "shift", "scale", "identity", "affine", "gnomonic",
               "zenithal_perspective", "slant_zenithal_perspective", "stereographic", "slant_orthographic", "zenithal_equidistant", "zenithal_equal_area", "airy",
               "plate_carree", "cylindrical_equal_area", "cylindrical_perspective", "mercator", "sanson_flamsteed", "parabolic", "molleweide", "hammer_aitoff",
               "polyconic", "conic_perspective", "conic_equidistant", "conic_equal_area", "conic_orthomorphic", "bonne_equal_area", "tangential_spherical_cube",
               "cobe_quad_spherical_cube", "quad_spherical_cube", "healpix", "healpix_polar", "rotate3d", "rotate2d", "constant", "fix_inputs", "remap_axes", "polynomial", "rotate_sequence_3d")
    registered = registered_versions(registry, "tag:stsci.edu:asdf/transform/$family")
    available = available_versions(transform_root, family)
    registered == available || error("$family versions differ: GWCS.jl has $registered, schema package has $available")
end

for family in ("wcs", "step", "frame", "frame2d", "celestial_frame", "spectral_frame", "temporal_frame", "stokes_frame", "composite_frame", "label_mapper")
    registered = registry.schemas["tag:stsci.edu:gwcs/$family"].versions
    available = available_versions(gwcs_root, family)
    registered == available || error("$family versions differ: GWCS.jl has $registered, schema package has $available")
end

for family in ("spherical_cartesian", "regions_selector")
    registered = registered_versions(registry, "tag:stsci.edu:gwcs/$family")
    available = available_versions(gwcs_root, family)
    registered == available || error("$family versions differ: GWCS.jl has $registered, schema package has $available")
end

property_root = joinpath(transform_root, "property")
for family in ("bounding_box", "compound_bounding_box")
    registered = registry.schemas["tag:stsci.edu:asdf/transform/property/$family"].versions
    available = available_versions(property_root, family)
    registered == available || error("$family versions differ: GWCS.jl has $registered, schema package has $available")
end

coordinate_root = joinpath(resource_path("asdf_coordinates_schemas.resources"), "schemas", "frames")
registered_icrs = registry.schemas["tag:astropy.org:astropy/coordinates/frames/icrs"].versions
available_icrs = available_versions(coordinate_root, "icrs")
registered_icrs == available_icrs || error("ICRS versions differ: GWCS.jl has $registered_icrs, schema package has $available_icrs")

astropy_root = joinpath(resource_path("asdf_astropy.resources"), "schemas")
for family in ("wcs", "fits")
    registered = registry.schemas["tag:astropy.org:astropy/$family/$family"].versions
    available = available_versions(joinpath(astropy_root, family), family)
    registered == available || error("astropy $family versions differ: GWCS.jl has $registered, schema package has $available")
end

asdf_root = resource_path("asdf_standard.resources")
ndarray_roots = (joinpath(asdf_root, "stable", "schemas", "stsci.edu", "asdf", "core"), joinpath(asdf_root, "unstable", "schemas", "stsci.edu", "asdf", "core"))
registered_ndarray = registry.schemas["tag:stsci.edu:asdf/core/ndarray"].versions
available_ndarray = available_versions(ndarray_roots, "ndarray")
registered_ndarray == available_ndarray || error("ndarray versions differ: GWCS.jl has $registered_ndarray, schema package has $available_ndarray")

unit_root = joinpath(asdf_root, "stable", "schemas", "stsci.edu", "asdf", "unit")
registered_quantity = registry.schemas["tag:stsci.edu:asdf/unit/quantity"].versions
available_quantity = available_versions(unit_root, "quantity")
registered_quantity == available_quantity || error("quantity versions differ: GWCS.jl has $registered_quantity, schema package has $available_quantity")

println("GWCS.jl metadata schema versions match the installed ASDF schema packages")
