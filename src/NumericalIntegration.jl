__precompile__()

module NumericalIntegration

using Compat

export integrate
export Trapezoidal, TrapezoidalEven, TrapezoidalFast, TrapezoidalEvenFast
export SimpsonEven, SimpsonEvenFast

@compat abstract type IntegrationMethod end

immutable Trapezoidal         <: IntegrationMethod end
immutable TrapezoidalEven     <: IntegrationMethod end
immutable TrapezoidalFast     <: IntegrationMethod end
immutable TrapezoidalEvenFast <: IntegrationMethod end
immutable SimpsonEven         <: IntegrationMethod end # https://en.wikipedia.org/wiki/Simpson%27s_rule#Alternative_extended_Simpson.27s_rule
immutable SimpsonEvenFast     <: IntegrationMethod end

function integrate{X<:Real, Y<:Real}(x::AbstractVector{X}, y::AbstractVector{Y}, ::Trapezoidal)
  @assert length(x) == length(y) "x and y vectors must be of the same length!"
  retval = zero(promote(x[1], y[1])[1])
    for i in 1 : length(y)-1
      retval += (x[i+1] - x[i]) * (y[i] + y[i+1])
    end
  return 0.5 * retval
end

function integrate{X<:Real, Y<:Real}(x::AbstractVector{X}, y::AbstractVector{Y}, ::TrapezoidalEven)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    return 0.5 * (x[end] - x[1]) / (length(y) - 1) * (y[1] + y[end] + sum(y[2:end-1]))
end

function integrate{X<:Real, Y<:Real}(x::AbstractVector{X}, y::AbstractVector{Y}, ::TrapezoidalFast)
    retval = zero(promote(x[1], y[1])[1])
    @fastmath @simd for i in 1 : length(y)-1
        @inbounds retval += (x[i+1] - x[i]) * (y[i] + y[i+1])
    end
    return 0.5 * retval
end

function integrate{X<:Real, Y<:Real}(x::AbstractVector{X}, y::AbstractVector{Y}, ::TrapezoidalEvenFast)
    retval = zero(promote(x[1], y[1])[1])
    N = length(y) - 1
    @fastmath @simd for i in 2 : N
        @inbounds retval += y[i]
    end
    @inbounds return (x[end] - x[1]) / N * (retval + 0.5*y[1] + 0.5*y[end])
end

function integrate{X<:Real, Y<:Real}(x::AbstractVector{X}, y::AbstractVector{Y}, ::SimpsonEven)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    retval = (17*y[1] + 59*y[2] + 43*y[3] + 49*y[4] + 49*y[end-3] + 43*y[end-2] + 59*y[end-1] + 17*y[end]) / 48
    for i in 5 : length(y) - 1
        retval += y[i]
    end
    return (x[end] - x[1]) / (length(y) - 1) * retval
end

function integrate{X<:Real, Y<:Real}(x::AbstractVector{X}, y::AbstractVector{Y}, ::SimpsonEvenFast)
    @inbounds retval = 17*y[1] + 59*y[2] + 43*y[3] + 49*y[4]
    @inbounds retval += 49*y[end-3] + 43*y[end-2] + 59*y[end-1] + 17*y[end]
    retval /= 48
    @fastmath @inbounds for i in 5 : length(y)-1
        retval += y[i]
    end
    @inbounds return (x[end] - x[1]) / (length(y) - 1) * retval
end

integrate{X<:Real, Y<:Real}(x::AbstractVector{X}, y::AbstractVector{Y}) = integrate(x, y, TrapezoidalFast())

end
