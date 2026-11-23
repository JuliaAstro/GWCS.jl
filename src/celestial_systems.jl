# Rotations between celestial systems, mostly as astropy's coordinate frames define them; the pipelines guide's
# Conversion details lists the models and the differences from astropy. Each matrix maps ICRS unit vectors to the
# system's; `_axis_rotation` is astropy's `rotation_matrix` convention.

# ICRS to FK5 J2000, the frame bias of USNO Circular 179 that astropy uses.
const _ICRS_TO_FK5_J2000 = _axis_rotation(19.9 / 3.6e6, 'x') * _axis_rotation(9.1 / 3.6e6, 'y') * _axis_rotation(-22.9 / 3.6e6, 'z')

# FK5 J2000 to Galactic, from astropy's J2000 north Galactic pole and longitude of the ascending node.
const _FK5_J2000_TO_GALACTIC = _axis_rotation(180 - 122.9319185680026, 'z') * _axis_rotation(90 - 27.12825118085622, 'y') * _axis_rotation(192.8594812065348, 'z')

# Galactic to Supergalactic, from astropy's north supergalactic pole at l = 47.37, b = 6.32.
const _GALACTIC_TO_SUPERGALACTIC = _axis_rotation(90.0, 'z') * _axis_rotation(90 - 6.32, 'y') * _axis_rotation(47.37, 'z')

# The IAU 1976 precession of Lieske et al. (1977), the FK5 catalog's, from the mean J2000 equator and equinox to those
# of a Julian epoch: ERFA's prec76 and pmat76, as SLALIB's sla_PREC.
function _fk5_precession(epoch)
    t = (epoch - 2000) / 100
    zeta = (2306.2181 + (0.30188 + 0.017998 * t) * t) * t / 3600
    z = (2306.2181 + (1.09468 + 0.018203 * t) * t) * t / 3600
    theta = (2004.3109 + (-0.42665 - 0.041833 * t) * t) * t / 3600
    return _axis_rotation(-z, 'z') * _axis_rotation(theta, 'y') * _axis_rotation(-zeta, 'z')
end

# The mean obliquity of the ecliptic of a Julian epoch in the IAU 1980 theory (SLALIB's sla_ECMAT, ERFA's obl80), and of
# a Besselian epoch in Newcomb's, as re-expressed in Julian centuries from B1950 by SLALIB's sla_ETRMS after Smith et
# al. (1989), in degrees.
_iau1980_obliquity(epoch) = @evalpoly((epoch - 2000) / 100, 84381.448, -46.8150, -0.00059, 0.001813) / 3600
_newcomb_obliquity(epoch) = @evalpoly((epoch - 1950) * 1.00002135903e-2, 84404.836, -46.8495, -0.00319, -0.00181) / 3600

# ICRS to the IAU 2006 mean ecliptic and equinox of J2000, ERFA's ecm06 and astropy's BarycentricMeanEcliptic: the
# Fukushima-Williams angles of ERFA's pfw06 at J2000, whose final rotation by the obliquity ecm06 undoes.
const _ICRS_TO_ECLIPTIC = let
    gamb, phib, psib = -0.052928 / 3600, 84381.412819 / 3600, -0.041775 / 3600
    _axis_rotation(-psib, 'z') * _axis_rotation(phib, 'x') * _axis_rotation(gamb, 'z')
end

# FK4-NO-E B1950 to FK5 J2000 for a Besselian FK4 obstime: the matrix of Murray (1989) eqn 28 plus his eqn 29
# correction for FK4's rotation, whose time argument is the obstime as a Julian epoch, as astropy applies them.
const _B1950_TO_J2000 = SMatrix{3,3,Float64}([0.9999256794956877 -0.0111814832204662 -0.0048590038153592;
    0.0111814832391717 0.9999374848933135 -0.0000271625947142; 0.0048590037723143 -0.0000271702937440 0.9999881946023742])
const _FK4_CORRECTION = SMatrix{3,3,Float64}([-0.0026455262 -1.1539918689 2.1111346190; 1.1540628161 -0.0129042997 0.0236021478;
    -2.1112979048 -0.0056024448 0.0102587734]) * 1e-6
function _fk4_to_fk5_j2000(obstime)
    matrix = _B1950_TO_J2000 + _FK4_CORRECTION * ((_convert_epoch(obstime, "byear", :julian) - 1950) / 100)
    # The linear correction leaves the matrix about 1e-12 from orthogonal; its nearest rotation, by one step of the
    # Newton polar iteration, keeps the conversion a rotation with an exact inverse and a faithful ASDF form.
    return (matrix + transpose(inv(matrix))) / 2
end

# Newcomb's precession between Besselian epochs, as astropy's FK4 frames compute it.
function _newcomb_precession(from, to)
    t = (from - 1850) / 1000
    dt = (to - 1850) / 1000 - t
    zeta = ((17.995 * dt + (-0.27 * t + 30.240)) * dt + ((0.060 * t + 139.720) * t + 23035.545)) * dt / 3600
    z = ((18.325 * dt + (109.480 + 0.39 * t)) * dt + ((0.060 * t + 139.720) * t + 23035.545)) * dt / 3600
    theta = ((-41.8 * dt + (-0.37 * t - 42.65)) * dt + ((-0.37 * t - 85.29) * t + 20051.12)) * dt / 3600
    return _axis_rotation(-z, 'z') * _axis_rotation(theta, 'y') * _axis_rotation(-zeta, 'z')
end


# The E-terms of aberration that FK4 positions include, for a Besselian equinox, as astropy's `fk4_e_terms` gives them:
# the constant of aberration, the eccentricity and perigee of Earth's orbit, and ERFA's obl80 obliquity.
function _fk4_eterms(equinox)
    jd = 2415020.31352 + (equinox - 1900) * 365.242198781
    t = (jd - (2415020.31352 + 50 * 365.242198781)) / 36525
    eccentricity = deg2rad(0.0056932) * ((-0.000000126 * t - 0.00004193) * t + 0.01673011)
    perigee = deg2rad((((0.012 * t + 1.65) * t + 6190.67) * t + 1015489.951) / 3600)
    obliquity = deg2rad(@evalpoly((jd - 2451545) / 36525, 84381.448, -46.8150, -0.00059, 0.001813) / 3600)
    sin_g, cos_g = sincos(perigee)
    return SVector(eccentricity * sin_g, -eccentricity * cos_g * cos(obliquity), -eccentricity * cos_g * sin(obliquity))
end

"""The removal (or, with `remove=false`, the addition) of FK4's E-terms of aberration `eterms` on unit vectors."""
struct _FK4ETerms <: AbstractTransform{3,3}
    eterms::SVector{3,Float64}
    remove::Bool
end

@inline function (transform::_FK4ETerms)(coords::StaticVector{3})
    e = transform.eterms
    # SLALIB's sla_SUBET and sla_ADDET, after the Explanatory Supplement (1961), sect. 2D. astropy's fixed-point addition
    # gives the same direction, since each of its iterates is a multiple of `coords + e`.
    vector = transform.remove ? coords - e + dot(e, coords) * coords : coords + e
    return vector / norm(vector)
end

inverse_transform(transform::_FK4ETerms) = _FK4ETerms(transform.eterms, !transform.remove)

const _CONVERTIBLE_SYSTEMS = (:icrs, :fk5, :fk4, :fk4noeterms, :galactic, :supergalactic, :ecliptic_fk5, :ecliptic_fk4, :ecliptic_icrs)

function _check_convertible(system, what)
    system in _CONVERTIBLE_SYSTEMS || throw(ArgumentError("$what $(repr(system)) cannot be converted; supported systems are $(join(map(repr, _CONVERTIBLE_SYSTEMS), ", "))"))
end

# The rotations between a system and ICRS. An FK4 system, here without its E-terms, reaches ICRS through FK5 J2000:
# Newcomb's precession to B1950 then Murray's rotation, as astropy defines it; the reverse is its exact inverse.
function _icrs_to(spec)
    spec.system === :icrs && return SMatrix{3,3,Float64}(I)
    spec.system === :fk5 && return _fk5_precession(spec.equinox) * _ICRS_TO_FK5_J2000
    spec.system === :galactic && return _FK5_J2000_TO_GALACTIC * _ICRS_TO_FK5_J2000
    spec.system === :supergalactic && return _GALACTIC_TO_SUPERGALACTIC * _FK5_J2000_TO_GALACTIC * _ICRS_TO_FK5_J2000
    return transpose(_to_icrs(spec))
end

function _to_icrs(spec)
    spec.system in (:fk4, :fk4noeterms) || return transpose(_icrs_to(spec))
    return transpose(_ICRS_TO_FK5_J2000) * _fk4_to_fk5_j2000(spec.obstime) * _newcomb_precession(spec.equinox, 1950.0)
end


# The rotation from `source` to `target`, each `(; system, equinox, obstime)`, through their equatorial systems.
function _celestial_rotation(source, target)
    (equatorial_source, from), (equatorial_target, to) = _equatorial_base(source), _equatorial_base(target)
    return to * _equatorial_rotation(equatorial_source, equatorial_target) * transpose(from)
end

# A FITS ecliptic system as its equatorial system and the rotation from that system to it: FK5 by the IAU 1980
# obliquity, FK4 without its E-terms by Newcomb's, with the obstime its equinox, and ICRS by astropy's J2000 ecliptic.
function _equatorial_base(spec)
    spec.system === :ecliptic_fk5 && return ((; system=:fk5, equinox=spec.equinox, obstime=nothing), _axis_rotation(_iau1980_obliquity(spec.equinox), 'x'))
    spec.system === :ecliptic_fk4 && return ((; system=:fk4noeterms, equinox=spec.equinox, obstime=something(spec.obstime, spec.equinox)), _axis_rotation(_newcomb_obliquity(spec.equinox), 'x'))
    spec.system === :ecliptic_icrs && return ((; system=:icrs, equinox=nothing, obstime=nothing), _ICRS_TO_ECLIPTIC)
    return (spec, SMatrix{3,3,Float64}(I))
end

# Every system reaches every other through ICRS, so an FK4 system's obstime enters each conversion through Murray's
# correction, whatever the route; astropy instead converts FK4 to Galactic coordinates through the B1950 Galactic pole,
# and between FK4 systems by precession alone, without the obstime.
_equatorial_rotation(source, target) = _icrs_to(target) * _to_icrs(source)

function _system_frame(target, name, axes_order, unit)
    system = target.system
    if system in (:ecliptic_fk5, :ecliptic_fk4, :ecliptic_icrs)
        # astropy has no frame for a FITS ecliptic, so it is kept as the FITS metadata, which has no ASDF form.
        radesys = system === :ecliptic_fk5 ? "FK5" : system === :ecliptic_fk4 ? "FK4" : "ICRS"
        return CelestialFrame{2}(; name, axes_names=("elon", "elat"), axes_order, axis_physical_types=("pos.ecliptic.lon", "pos.ecliptic.lat"), unit,
            reference_frame=(system=("ELON", "ELAT"), radesys, equinox=target.equinox, obstime=target.obstime))
    elseif system === :galactic
        return CelestialFrame{2}(; name, axes_names=("l", "b"), axes_order, axis_physical_types=("pos.galactic.lon", "pos.galactic.lat"), unit, reference_frame=_fits_astropy_frame("galactic"))
    elseif system === :supergalactic
        # asdf-astropy cannot write a supergalactic frame, so writing this one to ASDF throws.
        return CelestialFrame{2}(; name, axes_names=("sgl", "sgb"), axes_order, axis_physical_types=("pos.supergalactic.lon", "pos.supergalactic.lat"), unit, reference_frame=:Supergalactic)
    end
    reference = system === :icrs ? :ICRS : system === :fk5 ? _fits_astropy_frame("fk5", "equinox" => _fits_epoch(target.equinox, "jyear")) :
        _fits_astropy_frame(String(system), "equinox" => _fits_epoch(target.equinox, "byear"), "obstime" => _fits_epoch(target.obstime, "byear"))
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
    with_celestial_system(wcs, system; equinox=nothing, obstime=nothing, keep_frame=false, name) -> WCS

Return a copy of `wcs` whose world coordinates are in the celestial `system`, with its output frame replaced by one in
that system.

`system` is `:icrs`, `:fk5`, `:fk4`, `:fk4noeterms`, `:galactic`, `:supergalactic`, or one of the FITS ecliptic systems
`:ecliptic_fk5`, `:ecliptic_fk4`, and `:ecliptic_icrs`, as [`celestial_system`](@ref) names them; the WCS's own system
must be one of these too. `equinox` is the Julian epoch of an FK5 system or ecliptic, J2000 by default, or the
Besselian epoch of an FK4 or FK4-NO-E system or FK4 ecliptic, B1950 by default, and `obstime` is the Besselian epoch
of observation of an FK4 or FK4-NO-E system or FK4 ecliptic, its equinox by default, which Murray's correction for
FK4's rotation uses. ICRS defines no ecliptic, so a conversion from
or to `:ecliptic_icrs` warns that it takes astropy's J2000 mean ecliptic. [Conversion details](@ref) gives the models
used and where the results differ from astropy's.

By default the new output frame takes the old one's name, and the rotation is folded into the transform into it: into
its final `NativeToCelestial`, as for a FITS image, or into its final `Rotation3D`, as for a Python GWCS V2/V3-to-sky
rotation, so evaluation costs no more than before. Otherwise it is composed after that transform; so is a conversion
from or to FK4, whose E-terms are not a rotation and have no ASDF form, so writing the result throws. With
`keep_frame=true`, the old output frame remains as an intermediate frame, so [`GWCS.transform`](@ref transform) still reaches the
original system, and the new frame, named `String(system)` by default, follows it through a separate rotation step.

The new celestial frame keeps the positions and units of the old one, so the world coordinates stay in the same order
and units, and the other axes of a composite output frame pass through unchanged. Each celestial axis may be in `deg`,
`rad`, `arcmin`, `arcsec`, `mas`, or `hourangle`, and longitudes are in `[0, 360)` degrees in that unit. The rotation has
an analytical inverse, so [`world_to_pixel`](@ref) takes coordinates in the new system. A separate rotation step's
inverse gives the original pipeline longitudes in `[0, 360)`, which a sky projection accepts in any range; a pipeline
whose longitudes are linear in the pixels, such as a FITS header without a projection code, inverts correctly only for
longitudes in that range. Writing the result to ASDF records the new system as an astropy frame, except for a
Supergalactic or ecliptic one, which asdf-astropy cannot write; nor can it write a kept frame in a system it has no
frame for.

A WCS without a celestial frame, a celestial axis in another unit, an unsupported system, an `equinox` for a system
without one, an `obstime` for a system other than FK4, FK4-NO-E, and FK4 ecliptics, and a `name` another frame has
throw an `ArgumentError`.

Python equivalent: `astropy.coordinates.SkyCoord.transform_to`, which converts coordinate objects rather than a WCS.
"""
function with_celestial_system(wcs::WCS, system::Symbol; equinox=nothing, obstime=nothing, keep_frame::Bool=false, name=nothing)
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
    _check_convertible(system, "celestial system")
    fk4 = system in (:fk4, :fk4noeterms, :ecliptic_fk4)
    julian, besselian = system in (:fk5, :ecliptic_fk5), system in (:fk4, :fk4noeterms, :ecliptic_fk4)
    julian || besselian || equinox === nothing || throw(ArgumentError("equinox applies only to FK5 and FK4 systems, not $(repr(system))"))
    fk4 || obstime === nothing || throw(ArgumentError("obstime applies only to FK4, FK4-NO-E, and FK4 ecliptic systems, not $(repr(system))"))
    target_equinox = julian ? Float64(something(equinox, 2000.0)) : besselian ? Float64(something(equinox, 1950.0)) : nothing
    target_obstime = fk4 ? Float64(something(obstime, target_equinox)) : nothing
    for (keyword, value) in (("equinox", target_equinox), ("obstime", target_obstime))
        value === nothing || isfinite(value) || throw(ArgumentError("$keyword must be finite, got $value"))
    end
    target = (; system, equinox=target_equinox, obstime=target_obstime)
    source = celestial_system(frame)
    _check_convertible(source.system, "the WCS's celestial system")
    :ecliptic_icrs in (source.system, system) && @warn "ICRS defines no ecliptic; the ICRS ecliptic is taken as astropy's BarycentricMeanEcliptic at J2000, the IAU 2006 mean ecliptic and equinox"
    matrix = _celestial_rotation(source, target)
    # FK4's E-terms are not a rotation, so a conversion from or to FK4 is a separate step around the rotation.
    removal = source.system === :fk4 ? (_FK4ETerms(_fk4_eterms(source.equinox), true),) : ()
    addition = system === :fk4 ? (_FK4ETerms(_fk4_eterms(target_equinox), false),) : ()
    celestial = _system_frame(target, world === frame ? name : String(system), frame.axes_order, frame.unit)
    new_world = world === frame ? celestial : CompositeFrame(map(leaf -> leaf === frame ? celestial : leaf, world.frames); name)
    incoming = wcs.steps[end - 1].transform
    pure = isempty(removal) && isempty(addition)
    fused = keep_frame || !pure || world !== frame || scales != (1.0, 1.0) ? nothing : _fuse_sky_rotation(incoming, matrix, frame.axes_order)
    if fused === nothing
        rotation = compose(SphericalToCartesian(), removal..., Rotation3D(matrix), addition..., CartesianToSpherical())
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
