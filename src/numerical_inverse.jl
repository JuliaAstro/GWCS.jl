"""Numerical inverse failed to reach the requested tolerance."""
struct NoConvergence{S,A} <: Exception
    message::String
    best_solution::S
    accuracy::A
    niter::Int
end

function Base.showerror(io::IO, error::NoConvergence)
    print(io, error.message)
end

function _box_axis_guess(wcs::WCS, box::BoundingBox, axis, ::Type{T}) where {T}
    interval = bounds(box, axis)
    interval !== nothing && return T(sum(interval) / 2)
    wcs.pixel_shape !== nothing && wcs.pixel_shape[axis] !== nothing && return T((wcs.pixel_shape[axis] - 1) / 2)
    return one(T)
end

function _box_initial_guess(wcs::WCS{N}, box::BoundingBox{N}, ::Type{T}) where {N,T}
    return SVector{N,T}(ntuple(axis -> _box_axis_guess(wcs, box, axis, T), N))
end

function _box_initial_guess(wcs::WCS{N}, box::CompoundBoundingBox{N,K}, ::Type{T}) where {N,K,T}
    length(box.boxes) == 1 || throw(ArgumentError("numerical inversion with a compound bounding box requires an initial_guess unless the box has one selector entry"))
    key, selected = only(box.boxes)
    return SVector{N,T}(ntuple(N) do axis
        position = findfirst(==(axis), box.selector_axes)
        position === nothing ? _box_axis_guess(wcs, selected, axis, T) : T(key[position])
    end)
end

# Splits a forward pipeline without an inverse into a non-invertible prefix and the inverse of its
# longest analytically invertible tail, such as pixels -> SIP -> linear and linear -> TAN -> rotation.
function _numerical_inverse_split(forward::AbstractTransform, bounding_box, pixel_shape)
    transforms = forward isa ComposedTransform ? forward.transforms : (forward,)
    first_invertible = length(transforms) + 1
    while first_invertible > 1 && !(_inverse_transform_checked(transforms[first_invertible - 1]) isa NoInverse)
        first_invertible -= 1
    end
    prefix = compose(transforms[1:(first_invertible - 1)]...)
    # Newton evaluates `prefix`; only the initial guess uses `initializer`, where the WithInverse leaves would cost
    # the solve itself a little.
    initializer = _with_initializer_inverses(prefix)
    tail_inverse = first_invertible > length(transforms) ? IdentityTransform{output_dimension(forward)}() : _inverse_transform_checked(compose(transforms[first_invertible:end]...))
    return (; prefix, initializer, tail_inverse, approximate=_approximate_inverse(prefix, initializer, bounding_box, pixel_shape))
end

# Like Python GWCS's approximate inverse, a polynomial fit of pixels as a function of the prefix output, sampled over
# the image, gives Newton a starting point inside it. Fitting the prefix output needs no projection, unlike Python's
# fit in a TAN plane, and the fit must reach 0.01 px rather than Python's 5 px: coarser guesses save no iterations.
_approximate_inverse(prefix, initializer, bounding_box, pixel_shape) = nothing

# The fit's degree and whether a fit exists depend on the coefficients, so neither may change the WCS type: the
# polynomials use Matrix storage, and without a fit an empty box, which no guess lies in, disables the fit.
function _approximate_inverse_entry(coefficients, domain, low, high)
    P = Polynomial2D{Float64,Matrix{Float64},NTuple{2,NTuple{2,Float64}},NTuple{2,NTuple{2,Float64}}}
    window = ((-1.0, 1.0), (-1.0, 1.0))
    return (; fit=_Polynomial2DPair(P(coefficients[1], domain, window), P(coefficients[2], domain, window)), low, high)
end

const _NO_APPROXIMATE_INVERSE = _approximate_inverse_entry((zeros(1, 1), zeros(1, 1)), ((-1.0, 1.0), (-1.0, 1.0)), SVector(Inf, Inf), SVector(-Inf, -Inf))

function _approximate_inverse(@nospecialize(prefix::AbstractTransform{2,2}), @nospecialize(initializer), bounding_box, pixel_shape)
    domain = ntuple(axis -> _approximate_inverse_interval(bounding_box, pixel_shape, axis), 2)
    any(isnothing, domain) && return _NO_APPROXIMATE_INVERSE
    low, high = SVector{2,Float64}(first.(domain)), SVector{2,Float64}(last.(domain))
    all(isfinite, low) && all(isfinite, high) && all(low .< high) || return _NO_APPROXIMATE_INVERSE
    # An explicit approximate reverse, such as SIP AP/BP, already starts Newton near the root.
    last(_inverse_initialize(initializer, prefix((low + high) / 2))) && return _NO_APPROXIMATE_INVERSE
    fit = _fit_approximate_inverse(_approximate_inverse_samples(prefix, low, high, 16, false)..., _approximate_inverse_samples(prefix, low, high, 16, true)...)
    return fit === nothing ? _NO_APPROXIMATE_INVERSE : _approximate_inverse_entry(fit..., low, high)
end

function _approximate_inverse_interval(bounding_box, pixel_shape, axis)
    interval = bounding_box isa BoundingBox ? bounds(bounding_box, axis) : nothing
    interval !== nothing && return Float64.(interval)
    (pixel_shape === nothing || pixel_shape[axis] === nothing) && return nothing
    # Zero-based pixel centers 0:n-1, extended to the outer pixel edges.
    return (-0.5, pixel_shape[axis] - 0.5)
end

# An `n` x `n` lattice spanning the domain, or with `stagger`, the `n - 1` x `n - 1` lattice of its cell centers.
function _approximate_inverse_samples(@nospecialize(prefix), low, high, n, stagger)
    step = (high - low) / (n - 1)
    start, stop, m = stagger ? (low + step / 2, high - step / 2, n - 1) : (low, high, n)
    pixels = [SVector(x, y) for y in range(start[2], stop[2]; length=m) for x in range(start[1], stop[1]; length=m)]
    return pixels, map(pixel -> SVector{2,Float64}(prefix(pixel)), pixels)
end

# Least squares in monomials of the normalized prefix output at odd degrees 3 to 9, which reduce the error far more
# than even ones, keeping the first fit whose largest error on the staggered lattice is within 0.01 px. A 16 x 16
# lattice tracks a dense check to within a factor of 1.3; divergent or folded samples never pass.
function _fit_approximate_inverse(pixels, intermediate, check_pixels, check_intermediate; tolerance=0.01, degrees=3:2:9)
    all(point -> all(isfinite, point), intermediate) && all(point -> all(isfinite, point), check_intermediate) || return nothing
    domain = (extrema(first, intermediate), extrema(last, intermediate))
    all(interval -> interval[1] < interval[2], domain) || return nothing
    normalize(point, axis) = (point[axis] - sum(domain[axis]) / 2) * 2 / (domain[axis][2] - domain[axis][1])
    targets, check_targets = permutedims(reduce(hcat, pixels)), permutedims(reduce(hcat, check_pixels))
    design(points, degree) = [normalize(point, 1)^i * normalize(point, 2)^j for point in points, (i, j) in [(i, j) for i in 0:degree for j in 0:(degree - i)]]
    for degree in degrees
        coefficients = design(intermediate, degree) \ targets
        maximum(abs, design(check_intermediate, degree) * coefficients - check_targets) <= tolerance || continue
        matrices = ntuple(2) do axis
            matrix = zeros(degree + 1, degree + 1)
            index = 0
            for i in 0:degree, j in 0:(degree - i)
                matrix[i + 1, j + 1] = coefficients[index += 1, axis]
            end
            matrix
        end
        return (matrices, domain)
    end
    return nothing
end

# The initial guess inverts each invertible prefix leaf on every solve; storing those inverses as WithInverse
# leaves, which evaluate as the leaf itself, keeps it from building transforms, such as an inverse region selector,
# at runtime. Leaves without an inverse keep their own initializers.
_with_initializer_inverses(transform::ComposedTransform) = ComposedTransform(map(_with_initializer_inverses, transform.transforms))
_with_initializer_inverses(transform::ParallelTransform) = ParallelTransform(map(_with_initializer_inverses, transform.transforms))
_with_initializer_inverses(transform::WithInverse) = transform

function _with_initializer_inverses(transform::AbstractTransform)
    inverse = _inverse_transform_checked(transform)
    return inverse isa NoInverse ? transform : WithInverse(transform, inverse)
end

function _initial_guess(wcs::WCS{N}, intermediate::StaticVector, ::Type{T}, initial_guess, solve, offset=false) where {N,T}
    if initial_guess !== nothing
        length(initial_guess) == N || throw(DimensionMismatch("initial guess has length $(length(initial_guess)), expected $N"))
        return SVector{N,T}(initial_guess) .- offset
    end
    approximate = wcs.inverse_split.approximate
    if approximate !== nothing
        # Outside the fitted image the polynomial extrapolates, so the usual guesses apply there; without a fit the box is empty.
        guess = SVector{N,T}(approximate.fit(intermediate))
        all(approximate.low .<= guess .<= approximate.high) && return guess
    end
    initialized, used = _inverse_initialize(wcs.inverse_split.initializer, intermediate)
    usable = initialized !== nothing && all(isfinite, initialized)
    if used && usable
        return initialized
    elseif wcs.bounding_box !== nothing
        return _box_initial_guess(wcs, wcs.bounding_box, T)
    elseif usable && all(isfinite, solve(SVector{N,T}(initialized)))
        # Like WCSLIB, start from the solution with distortions treated as identity when the map being solved is defined there.
        return SVector{N,T}(initialized)
    end
    return SVector{N,T}(ntuple(i -> one(T), N))
end

function _longitude_axes(frame::AbstractCoordinateFrame{N}) where {N}
    types = axis_physical_types(frame)
    return ntuple(i -> endswith(types[i], ".ra") || endswith(types[i], ".lon"), N)
end

@inline function _inverse_residual(solve::AbstractTransform{N,N}, pixels::StaticVector{N}, goal::StaticVector{N}, longitude_axes) where {N}
    predicted = solve(pixels)
    return SVector{N}(ntuple(i -> begin
        difference = predicted[i] - goal[i]
        longitude_axes[i] ? mod(difference + oftype(difference, 180), oftype(difference, 360)) - oftype(difference, 180) : difference
    end, N))
end

# Isolated coordinate singularities, such as the native pole of zenithal projections, have no
# chain-rule Jacobian even where the composed map is smooth; Newton then uses a nearby point's Jacobian.
function _inverse_jacobian(solve::AbstractTransform, estimate::StaticVector)
    jacobian = _transform_jacobian_checked(solve, estimate)
    all(isfinite, jacobian) && return jacobian
    T = eltype(estimate)
    nudged = estimate .+ sqrt(eps(T)) .* max.(one(T), abs.(estimate))
    nudged_jacobian = _transform_jacobian_checked(solve, nudged)
    return all(isfinite, nudged_jacobian) ? nudged_jacobian : jacobian
end

# The solvers return the raw solution; `_finish_inverse` applies validity and bounding-box handling.
function _numerical_solve(wcs::WCS, world::StaticVector; quiet=false, offset=false, kwargs...)
    converged, solution, accuracy, iterations = _numerical_result(wcs, world; offset, kwargs...)
    (converged || quiet) && return solution
    message = "numerical inverse failed to converge to the requested accuracy after $iterations iterations"
    throw(NoConvergence(message, solution .+ offset, accuracy, iterations))
end

# Returns `(converged, solution, accuracy, iterations)`, so callers that only need the outcome, like `in_image`,
# avoid constructing a `NoConvergence` per point. A target without a preimage converges to NaN.
function _numerical_result(wcs::WCS{0,0}, world::StaticVector{0}; tolerance=1e-5, maxiter=30, initial_guess=nothing, max_backtracks=12, offset=false)
    tolerance > 0 || throw(ArgumentError("tolerance must be positive"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    max_backtracks >= 0 || throw(ArgumentError("max_backtracks must be non-negative"))
    initial_guess === nothing || isempty(initial_guess) || throw(DimensionMismatch("initial guess has length $(length(initial_guess)), expected 0"))
    empty = SVector{0,_coordinate_float_type(eltype(world))}()
    return (true, empty, empty, 0)
end

function _numerical_result(wcs::WCS{N,N}, world::StaticVector{N}; tolerance=1e-5, maxiter=30, initial_guess=nothing, max_backtracks=12, offset=false) where {N}
    tolerance > 0 || throw(ArgumentError("tolerance must be positive"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    max_backtracks >= 0 || throw(ArgumentError("max_backtracks must be non-negative"))

    T = _coordinate_float_type(eltype(world))
    target = SVector{N,T}(world)
    undefined = (true, SVector{N,T}(ntuple(_ -> T(NaN), N)), zero(SVector{N,T}), 0)
    all(isfinite, target) || return undefined
    if wcs.inverse_split === nothing
        # An analytical inverse may only approximate the forward transform, as a fitted polynomial does, so it only
        # starts Newton on the forward transform itself.
        estimate = initial_guess === nothing ? wcs.backward(target) : SVector{N,T}(initial_guess) .- offset
        all(isfinite, estimate) || wcs.bounding_box === nothing || (estimate = _box_initial_guess(wcs, wcs.bounding_box, T))
        return _newton(wcs.forward, target, _longitude_axes(output_frame(wcs)), estimate, tolerance, maxiter, max_backtracks)
    end
    # NaN from the analytically invertible tail proves that no pixel maps to the target.
    intermediate = wcs.inverse_split.tail_inverse(target)
    any(isnan, intermediate) && return undefined
    prefix = wcs.inverse_split.prefix
    # Like astropy's all_world2pix, first solve only the non-invertible prefix against the target mapped through the
    # invertible tail, which keeps projections and rotations out of the loop. Intermediate coordinates have no known
    # periodicity, so a longitude formed before the tail can sit 360 degrees from its goal; when this solve fails,
    # solve the full forward map in world coordinates, as Python GWCS does.
    if output_dimension(prefix) == N && !(wcs.inverse_split.tail_inverse isa IdentityTransform)
        estimate = _initial_guess(wcs, intermediate, T, initial_guess, prefix, offset)
        result = _newton(prefix, intermediate, ntuple(_ -> false, N), estimate, tolerance, maxiter, max_backtracks)
        first(result) && return result
    end
    estimate = _initial_guess(wcs, intermediate, T, initial_guess, wcs.forward, offset)
    return _newton(wcs.forward, target, _longitude_axes(output_frame(wcs)), estimate, tolerance, maxiter, max_backtracks)
end

# Damped Newton for `solve(pixels) = goal`; returns the root when converged and the best estimate otherwise.
function _newton(solve::AbstractTransform{N,N}, goal::StaticVector{N,T}, longitude_axes, estimate, tolerance, maxiter, max_backtracks) where {N,T}
    # Steps that leave a transform's domain give NaN residuals, so the line search backtracks.
    residual_function = pixels -> _inverse_residual(solve, pixels, goal, longitude_axes)
    best_solution = estimate
    best_residual = residual_function(estimate)
    best_norm = maximum(abs, best_residual)
    accuracy = SVector{N,T}(ntuple(i -> T(Inf), N))
    iterations = 0
    last_step_accepted = false

    for iteration in 1:maxiter
        iterations = iteration
        last_step_accepted = false
        residual = residual_function(estimate)
        residual_norm = maximum(abs, residual)
        if residual_norm < best_norm
            best_solution = estimate
            best_residual = residual
            best_norm = residual_norm
        end
        all(isfinite, residual) || break

        jacobian = _inverse_jacobian(solve, estimate)
        all(isfinite, jacobian) || break
        scale = max(one(T), maximum(abs, jacobian))
        abs(det(jacobian)) > eps(T) * scale^N || break

        step = jacobian \ residual
        accuracy = abs.(step)
        maximum(accuracy) <= tolerance && return (true, estimate - step, accuracy, iterations)

        accepted = false
        candidate = estimate
        factor = one(T)
        for backtrack in 0:max_backtracks
            candidate = estimate - factor * step
            candidate_residual = residual_function(candidate)
            if all(isfinite, candidate_residual) && maximum(abs, candidate_residual) < residual_norm
                accepted = true
                break
            end
            factor /= T(2)
        end
        accepted || break
        estimate = candidate
        last_step_accepted = true
    end

    if last_step_accepted && iterations == maxiter
        residual = residual_function(estimate)
        if all(isfinite, residual)
            jacobian = _inverse_jacobian(solve, estimate)
            scale = max(one(T), maximum(abs, jacobian))
            if all(isfinite, jacobian) && abs(det(jacobian)) > eps(T) * scale^N
                step = jacobian \ residual
                accuracy = abs.(step)
                maximum(accuracy) <= tolerance && return (true, estimate - step, accuracy, iterations)
            end
        end
    end
    return (false, best_solution, accuracy, iterations)
end

function _numerical_result(wcs::WCS{Nin,Nout}, world::StaticVector; kwargs...) where {Nin,Nout}
    throw(ArgumentError("numerical inversion requires equal pixel and world dimensions, got $Nin and $Nout"))
end

function _numerical_inverse(wcs::WCS, world::StaticVector; with_bounding_box=true, fill_value=NaN, strict=false, origin=0, kwargs...)
    offset = _origin_offset(origin)
    return _finish_inverse(wcs, world, _numerical_solve(wcs, world; offset, kwargs...); with_bounding_box, fill_value, strict, offset)
end

"""
    numerical_inverse(wcs, world...; tolerance=1e-5, maxiter=30, initial_guess=nothing, max_backtracks=12, quiet=false,
                      with_bounding_box=true, fill_value=NaN, strict=false, origin=0) -> SVector

Solve `pixel_to_world(wcs, pixel) == world` for the pixel with the damped Newton solver that `world_to_pixel` uses for a
square WCS without an analytical inverse, to `tolerance` pixels. With an analytical inverse, which may only approximate
the forward transform, the solver starts from that inverse's pixel and refines it against the forward transform.
Nonconvergence throws `NoConvergence` unless `quiet=true`; the other keywords act as in `world_to_pixel`.
"""
function numerical_inverse(wcs::WCS{Nin,Nout}, coords::AbstractVector; kwargs...) where {Nin,Nout}
    length(coords) == Nout || throw(DimensionMismatch("world coordinate has length $(length(coords)), expected $Nout"))
    T = _coordinate_float_type(eltype(coords))
    return _numerical_inverse(wcs, SVector{Nout,T}(coords); kwargs...)
end

numerical_inverse(wcs::WCS, coords::Tuple{Vararg{Real}}; kwargs...) = numerical_inverse(wcs, collect(coords); kwargs...)
numerical_inverse(wcs::WCS, coords::Real...; kwargs...) = numerical_inverse(wcs, coords; kwargs...)
