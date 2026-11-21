# Rotations between celestial systems, defined as astropy's coordinate frames define them, so a converted WCS gives
# astropy's coordinates. Each matrix maps ICRS unit vectors to the system's; `_axis_rotation` is astropy's
# `rotation_matrix` convention.

# ICRS to FK5 J2000, the frame bias of USNO Circular 179 that astropy uses.
const _ICRS_TO_FK5_J2000 = _axis_rotation(19.9 / 3.6e6, 'x') * _axis_rotation(9.1 / 3.6e6, 'y') * _axis_rotation(-22.9 / 3.6e6, 'z')

# FK5 J2000 to Galactic, from astropy's J2000 north Galactic pole and longitude of the ascending node.
const _FK5_J2000_TO_GALACTIC = _axis_rotation(180 - 122.9319185680026, 'z') * _axis_rotation(90 - 27.12825118085622, 'y') * _axis_rotation(192.8594812065348, 'z')

# Galactic to Supergalactic, from astropy's north supergalactic pole at l = 47.37, b = 6.32.
const _GALACTIC_TO_SUPERGALACTIC = _axis_rotation(90.0, 'z') * _axis_rotation(90 - 6.32, 'y') * _axis_rotation(47.37, 'z')

# The IAU 2006 precession and frame bias from the mean J2000 equator and equinox to those of a Julian epoch, from the
# Fukushima-Williams angles of ERFA's pfw06 and fw2m.
function _fukushima_williams(epoch)
    t = (epoch - 2000) / 100
    gamb = @evalpoly(t, -0.052928, 10.556378, 0.4932044, -0.00031238, -0.000002788, 0.0000000260) / 3600
    phib = @evalpoly(t, 84381.412819, -46.811016, 0.0511268, 0.00053289, -0.000000440, -0.0000000176) / 3600
    psib = @evalpoly(t, -0.041775, 5038.481484, 1.5584175, -0.00018522, -0.000026452, -0.0000000148) / 3600
    epsa = @evalpoly(t, 84381.406, -46.836769, -0.0001831, 0.00200340, -0.000000576, -0.0000000434) / 3600
    return _axis_rotation(-epsa, 'x') * _axis_rotation(-psib, 'z') * _axis_rotation(phib, 'x') * _axis_rotation(gamb, 'z')
end

# The precession from the mean J2000 equinox to that of `epoch`, without frame bias, as astropy's FK5 takes it from
# ERFA's bp06.
_fk5_precession(epoch) = _fukushima_williams(epoch) * transpose(_fukushima_williams(2000.0))

const _CONVERTIBLE_SYSTEMS = (:icrs, :fk5, :galactic, :supergalactic)

function _icrs_to_system(system, equinox, what)
    system === :icrs && return SMatrix{3,3,Float64}(I)
    system === :fk5 && return _fk5_precession(equinox) * _ICRS_TO_FK5_J2000
    system === :galactic && return _FK5_J2000_TO_GALACTIC * _ICRS_TO_FK5_J2000
    system === :supergalactic && return _GALACTIC_TO_SUPERGALACTIC * _FK5_J2000_TO_GALACTIC * _ICRS_TO_FK5_J2000
    throw(ArgumentError("$what $(repr(system)) cannot be converted; supported systems are $(join(map(repr, _CONVERTIBLE_SYSTEMS), ", "))"))
end

function _system_frame(system, equinox, name, axes_order, unit)
    if system === :galactic
        return CelestialFrame{2}(; name, axes_names=("l", "b"), axes_order, axis_physical_types=("pos.galactic.lon", "pos.galactic.lat"), unit, reference_frame=_fits_astropy_frame("galactic"))
    elseif system === :supergalactic
        # asdf-astropy cannot write a supergalactic frame, so writing this one to ASDF throws.
        return CelestialFrame{2}(; name, axes_names=("sgl", "sgb"), axes_order, axis_physical_types=("pos.supergalactic.lon", "pos.supergalactic.lat"), unit, reference_frame=:Supergalactic)
    end
    reference = system === :icrs ? :ICRS : _fits_astropy_frame("fk5", "equinox" => _fits_epoch(equinox, "jyear"))
    return CelestialFrame{2}(; name, axes_names=("ra", "dec"), axes_order, axis_physical_types=("pos.eq.ra", "pos.eq.dec"), unit, reference_frame=reference)
end

# The matrix that takes native unit vectors to celestial ones: its columns are the celestial images of the native
# x, y, and z axes.
function _native_to_celestial_matrix(transform::NativeToCelestial)
    x, y, z = map(native -> SphericalToCartesian()(transform(SVector{2,Float64}(native))), ((0.0, 0.0), (90.0, 0.0), (0.0, 90.0)))
    return SMatrix{3,3,Float64}(x..., y..., z...)
end

# The native pole is at the celestial position of the native z axis, and the celestial pole at the native position
# given by the matrix's third row.
function _native_to_celestial(matrix)
    alpha_p = mod(atand(matrix[2, 3], matrix[1, 3]), 360)
    delta_p = atand(matrix[3, 3], hypot(matrix[1, 3], matrix[2, 3]))
    return NativeToCelestial(alpha_p, delta_p, atand(matrix[3, 2], matrix[3, 1]))
end

_bare_transform(transform::AnnotatedTransform) = _bare_transform(transform.transform)
_bare_transform(transform::AbstractTransform) = transform

# `inner` with the annotations of `original`. Unmodeled properties, such as Python's parameter bounds, are dropped from a
# node whose parameters the fusion changed, since they describe the old ones.
function _reannotate(original::AnnotatedTransform, inner)
    changed = original.transform isa Union{NativeToCelestial,Rotation3D,CartesianToSpherical}
    return AnnotatedTransform(inner; name=original.name, inputs=original.inputs, outputs=original.outputs, parameter_units=original.parameter_units,
        pass_through=changed ? _NO_PASS_THROUGH : original.pass_through)
end

# `transform`, whose outputs are a longitude and latitude in degrees in `order`, followed by the celestial `rotation`,
# folded into the sphere rotation that ends it: a `NativeToCelestial`, or a `Rotation3D` before the final
# `CartesianToSpherical`. `nothing` when it ends in neither, or when the pole angles of a fused `NativeToCelestial`,
# which are ill-conditioned near a pole, do not reproduce the rotation.
_fuse_sky_rotation(transform::AbstractTransform, rotation, order) = nothing

function _fuse_sky_rotation(transform::AnnotatedTransform, rotation, order)
    fused = _fuse_sky_rotation(transform.transform, rotation, order)
    return fused === nothing ? nothing : _reannotate(transform, fused)
end

function _fuse_sky_rotation(transform::NativeToCelestial, rotation, order)
    order == (1, 2) || return nothing
    matrix = rotation * _native_to_celestial_matrix(transform)
    fused = _native_to_celestial(matrix)
    return isapprox(_native_to_celestial_matrix(fused), matrix; atol=1e-13, rtol=0) ? fused : nothing
end

function _fuse_sky_rotation(transform::ComposedTransform, rotation, order)
    parts = transform.transforms
    final = _bare_transform(last(parts))
    length(parts) > 1 || return _fuse_sky_rotation(last(parts), rotation, order)
    if final isa RemapAxes{2,2} && all(i -> final.mapping[order[i]] == i, 1:2)
        # A final remapping that puts the frame's (lon, lat) at their positions: fuse before it.
        fused = _fuse_sky_rotation(ComposedTransform(Base.front(parts)), rotation, (1, 2))
    elseif final isa CartesianToSpherical && order == (1, 2)
        rotated = _rotate_cartesian(parts[end - 1], rotation)
        rotated === nothing && return nothing
        # Longitudes are in [0, 360), as an unfused conversion gives them.
        spherical = last(parts) isa AnnotatedTransform ? _reannotate(last(parts), CartesianToSpherical(360)) : CartesianToSpherical(360)
        return ComposedTransform((parts[1:(end - 2)]..., rotated, spherical))
    else
        fused = _fuse_sky_rotation(last(parts), rotation, order)
        return fused === nothing ? nothing : ComposedTransform((Base.front(parts)..., fused))
    end
    return fused === nothing ? nothing : _join_composed(fused, last(parts))
end

_join_composed(first::ComposedTransform, second) = ComposedTransform((first.transforms..., second))
_join_composed(first, second) = ComposedTransform((first, second))

# `transform`, whose outputs are a Cartesian vector, followed by `rotation`, when it ends in a `Rotation3D`.
_rotate_cartesian(transform::AbstractTransform, rotation) = nothing
_rotate_cartesian(transform::Rotation3D, rotation) = Rotation3D(rotation * transform.matrix)

function _rotate_cartesian(transform::AnnotatedTransform, rotation)
    rotated = _rotate_cartesian(transform.transform, rotation)
    return rotated === nothing ? nothing : _reannotate(transform, rotated)
end

function _rotate_cartesian(transform::ComposedTransform, rotation)
    rotated = _rotate_cartesian(last(transform.transforms), rotation)
    return rotated === nothing ? nothing : ComposedTransform((Base.front(transform.transforms)..., rotated))
end

"""
    with_celestial_system(wcs, system; equinox=nothing, keep_frame=false, name) -> WCS

Return a copy of `wcs` whose world coordinates are in the celestial `system`, with its output frame replaced by one in
that system.

`system` is `:icrs`, `:fk5`, `:galactic`, or `:supergalactic`, as [`celestial_system`](@ref) names them, and `equinox`
is the Julian epoch of an FK5 system, J2000 by default; the WCS's own system must be one of these too. The rotations
are astropy's: the frame bias of USNO Circular 179 between ICRS and FK5, the IAU 2006 precession between FK5 equinoxes,
and astropy's Galactic and Supergalactic poles.

By default the new output frame takes the old one's name, and the rotation is folded into the transform into it: into
its final `NativeToCelestial`, as for a FITS image, or into its final `Rotation3D`, as for a Python GWCS V2/V3-to-sky
rotation, so evaluation costs no more than before. Otherwise it is composed after that transform. With
`keep_frame=true`, the old output frame remains as an intermediate frame, so [`GWCS.transform`](@ref transform) still reaches the
original system, and the new frame, named `String(system)` by default, follows it through a separate rotation step.

The new celestial frame keeps the positions and units of the old one, so the world coordinates stay in the same order
and units, and the other axes of a composite output frame pass through unchanged. Each celestial axis may be in `deg`,
`rad`, `arcmin`, `arcsec`, `mas`, or `hourangle`, and longitudes are in `[0, 360)` degrees in that unit. The rotation has
an analytical inverse, so [`world_to_pixel`](@ref) takes coordinates in the new system. A separate rotation step's
inverse gives the original pipeline longitudes in `[0, 360)`, which a sky projection accepts in any range; a pipeline
whose longitudes are linear in the pixels, such as a FITS header without a projection code, inverts correctly only for
longitudes in that range. Writing the result to ASDF records
the new system as an astropy frame, except for a Supergalactic one, which asdf-astropy cannot write; nor can it write a
kept frame in a system it has no frame for.

A WCS without a celestial frame, a celestial axis in another unit, an unsupported system, an `equinox` for a system
other than FK5, and a `name` another frame has throw an `ArgumentError`.

Python equivalent: `astropy.coordinates.SkyCoord.transform_to`, which converts coordinate objects rather than a WCS.
"""
function with_celestial_system(wcs::WCS, system::Symbol; equinox=nothing, keep_frame::Bool=false, name=nothing)
    frame = celestial_frame(wcs)
    world = output_frame(wcs)
    frame === nothing && throw(ArgumentError("output frame $(world.name) has no celestial frame"))
    name = name === nothing ? (keep_frame ? String(system) : world.name) : String(name)
    name in available_frames(wcs) && (keep_frame || name != world.name) && throw(ArgumentError("the WCS already has a frame named $(repr(name))"))
    frame_dimension(frame) == 2 || throw(ArgumentError("celestial frame $(frame.name) has $(frame_dimension(frame)) axes; only longitude and latitude can be converted"))
    # The rotation works in degrees, so other angle units are scaled to degrees before it and back after.
    scales = map(frame.unit) do unit
        haskey(_ANGLE_UNIT_CONVERSIONS, unit) || throw(ArgumentError("celestial frame $(frame.name) has unit $(repr(unit)); its axes must be in an angle unit ($(join(sort!(collect(keys(_ANGLE_UNIT_CONVERSIONS))), ", ")))"))
        first(_ANGLE_UNIT_CONVERSIONS[unit])(1.0)
    end
    system === :fk5 || equinox === nothing || throw(ArgumentError("equinox applies only to an FK5 system, not $(repr(system))"))
    target_equinox = system === :fk5 ? Float64(something(equinox, 2000.0)) : nothing
    target_equinox === nothing || isfinite(target_equinox) || throw(ArgumentError("equinox must be finite, got $target_equinox"))
    source = celestial_system(frame)
    matrix = _icrs_to_system(system, target_equinox, "celestial system") * transpose(_icrs_to_system(source.system, source.equinox, "the WCS's celestial system"))
    celestial = _system_frame(system, target_equinox, world === frame ? name : String(system), frame.axes_order, frame.unit)
    new_world = world === frame ? celestial : CompositeFrame(map(leaf -> leaf === frame ? celestial : leaf, world.frames); name)
    incoming = wcs.steps[end - 1].transform
    fused = keep_frame || world !== frame || scales != (1.0, 1.0) ? nothing : _fuse_sky_rotation(incoming, matrix, frame.axes_order)
    if fused === nothing
        rotation = compose(SphericalToCartesian(), Rotation3D(matrix), CartesianToSpherical())
        if scales != (1.0, 1.0)
            rotation = compose(parallel(Scale(scales[1]), Scale(scales[2])), rotation, parallel(Scale(inv(scales[1])), Scale(inv(scales[2]))))
        end
        rotation = annotate_transform(rotation; name="$(source.system)_to_$system")
        N = frame_dimension(world)
        order = (frame.axes_order..., Tuple(axis for axis in 1:N if !(axis in frame.axes_order))...)
        step = N == 2 ? rotation : parallel(rotation, IdentityTransform{N - 2}())
        if order != ntuple(identity, N)
            step = compose(RemapAxes{N}(order), step, RemapAxes{N}(Tuple(invperm(collect(order)))))
        end
        keep_frame && return append_frame(wcs, step, new_world)
        fused = compose(incoming, step)
    end
    # The output step keeps its unmodeled metadata, except what described the replaced celestial frame.
    replaced = world === frame ? "frame" : "frame/frames/$(findfirst(leaf -> leaf === frame, world.frames))"
    pass_through = filter(((path, _),) -> path != replaced && !startswith(path, "$replaced/"), last(wcs.steps).pass_through)
    output = Step(new_world, nothing, isempty(pass_through) ? _NO_PASS_THROUGH : pass_through)
    steps = _replace_step(_replace_step(wcs.steps, length(wcs.steps) - 1, _with_transform(wcs.steps[end - 1], fused)), length(wcs.steps), output)
    return _rebuild_wcs(wcs, steps)
end
