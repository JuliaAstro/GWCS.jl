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
function _numerical_inverse_split(forward::AbstractTransform)
    transforms = forward isa ComposedTransform ? forward.transforms : (forward,)
    first_invertible = length(transforms) + 1
    while first_invertible > 1 && !(_inverse_transform_checked(transforms[first_invertible - 1]) isa NoInverse)
        first_invertible -= 1
    end
    prefix = compose(transforms[1:(first_invertible - 1)]...)
    first_invertible > length(transforms) && return (prefix=prefix, tail_inverse=IdentityTransform{output_dimension(forward)}())
    return (prefix=prefix, tail_inverse=_inverse_transform_checked(compose(transforms[first_invertible:end]...)))
end

function _initial_guess(wcs::WCS{N}, intermediate::StaticVector, ::Type{T}, initial_guess, solve) where {N,T}
    if initial_guess !== nothing
        length(initial_guess) == N || throw(DimensionMismatch("initial guess has length $(length(initial_guess)), expected $N"))
        return SVector{N,T}(initial_guess)
    end
    initialized, used = _inverse_initialize(wcs.inverse_split.prefix, intermediate)
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
    return ntuple(i -> endswith(frame.axis_physical_types[i], ".ra") || endswith(frame.axis_physical_types[i], ".lon"), N)
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
function _numerical_solve(wcs::WCS, world::StaticVector; quiet=false, kwargs...)
    converged, solution, accuracy, iterations = _numerical_result(wcs, world; kwargs...)
    (converged || quiet) && return solution
    message = "numerical inverse failed to converge to the requested accuracy after $iterations iterations"
    throw(NoConvergence(message, solution, accuracy, iterations))
end

# Returns `(converged, solution, accuracy, iterations)`, so callers that only need the outcome, like `in_image`,
# avoid constructing a `NoConvergence` per point. A target without a preimage converges to NaN.
function _numerical_result(wcs::WCS{0,0}, world::StaticVector{0}; tolerance=1e-5, maxiter=30, initial_guess=nothing, max_backtracks=12)
    tolerance > 0 || throw(ArgumentError("tolerance must be positive"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    max_backtracks >= 0 || throw(ArgumentError("max_backtracks must be non-negative"))
    initial_guess === nothing || isempty(initial_guess) || throw(DimensionMismatch("initial guess has length $(length(initial_guess)), expected 0"))
    empty = SVector{0,_coordinate_float_type(eltype(world))}()
    return (true, empty, empty, 0)
end

function _numerical_result(wcs::WCS{N,N}, world::StaticVector{N}; tolerance=1e-5, maxiter=30, initial_guess=nothing, max_backtracks=12) where {N}
    tolerance > 0 || throw(ArgumentError("tolerance must be positive"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    max_backtracks >= 0 || throw(ArgumentError("max_backtracks must be non-negative"))

    T = _coordinate_float_type(eltype(world))
    target = SVector{N,T}(world)
    undefined = (true, SVector{N,T}(ntuple(_ -> T(NaN), N)), zero(SVector{N,T}), 0)
    all(isfinite, target) || return undefined
    # NaN from the analytically invertible tail proves that no pixel maps to the target.
    intermediate = wcs.inverse_split.tail_inverse(target)
    any(isnan, intermediate) && return undefined
    prefix = wcs.inverse_split.prefix
    # Like astropy's all_world2pix, first solve only the non-invertible prefix against the target mapped through the
    # invertible tail, which keeps projections and rotations out of the loop. Intermediate coordinates have no known
    # periodicity, so a longitude formed before the tail can sit 360 degrees from its goal; when this solve fails,
    # solve the full forward map in world coordinates, as Python GWCS does.
    if output_dimension(prefix) == N && !(wcs.inverse_split.tail_inverse isa IdentityTransform)
        estimate = _initial_guess(wcs, intermediate, T, initial_guess, prefix)
        result = _newton(prefix, intermediate, ntuple(_ -> false, N), estimate, tolerance, maxiter, max_backtracks)
        first(result) && return result
    end
    estimate = _initial_guess(wcs, intermediate, T, initial_guess, wcs.forward)
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

function _numerical_inverse(wcs::WCS, world::StaticVector; with_bounding_box=true, fill_value=NaN, strict=false, kwargs...)
    return _finish_inverse(wcs, world, _numerical_solve(wcs, world; kwargs...); with_bounding_box, fill_value, strict)
end

function numerical_inverse(wcs::WCS{Nin,Nout}, coords::AbstractVector; kwargs...) where {Nin,Nout}
    length(coords) == Nout || throw(DimensionMismatch("world coordinate has length $(length(coords)), expected $Nout"))
    T = _coordinate_float_type(eltype(coords))
    return _numerical_inverse(wcs, SVector{Nout,T}(coords); kwargs...)
end

numerical_inverse(wcs::WCS, coords::Tuple{Vararg{Real}}; kwargs...) = numerical_inverse(wcs, collect(coords); kwargs...)
numerical_inverse(wcs::WCS, coords::Real...; kwargs...) = numerical_inverse(wcs, coords; kwargs...)
