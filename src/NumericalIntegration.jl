__precompile__()

module NumericalIntegration

export integrate
export Trapezoidal, TrapezoidalEven, TrapezoidalFast, TrapezoidalEvenFast
export SimpsonEven, SimpsonEvenFast
export IntegrationMethod

abstract type IntegrationMethod end

struct Trapezoidal         <: IntegrationMethod end
struct TrapezoidalEven     <: IntegrationMethod end
struct TrapezoidalFast     <: IntegrationMethod end
struct TrapezoidalEvenFast <: IntegrationMethod end
struct SimpsonEven         <: IntegrationMethod end # https://en.wikipedia.org/wiki/Simpson%27s_rule#Alternative_extended_Simpson.27s_rule
struct SimpsonEvenFast     <: IntegrationMethod end

const HALF = 1//2

function _zero(x,y)
    ret = zero(eltype(x)) + zero(eltype(y))
    ret / 2
end

function integrate(x::AbstractVector, y::AbstractVector, ::Trapezoidal)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    retval = _zero(x,y)
    for i in 1 : length(y)-1
        retval += (x[i+1] - x[i]) * (y[i] + y[i+1])
    end
    return HALF * retval
end

function integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEven)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    return (x[2] - x[1]) * (HALF * (y[1] + y[end]) + sum(y[2:end-1]))
end

function integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalFast)
    retval = _zero(x,y)
    @fastmath @simd for i in 1 : length(y)-1
        @inbounds retval += (x[i+1] - x[i]) * (y[i] + y[i+1])
    end
    return HALF * retval
end

function integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEvenFast)
    retval = _zero(x,y)
    N = length(y) - 1
    @fastmath @simd for i in 2 : N
        @inbounds retval += y[i]
    end
    @inbounds return (x[2] - x[1]) * (retval + HALF*y[1] + HALF*y[end])
end

function integrate(x::AbstractVector, y::AbstractVector, ::SimpsonEven)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    retval = (17*y[1] + 59*y[2] + 43*y[3] + 49*y[4] + 49*y[end-3] + 43*y[end-2] + 59*y[end-1] + 17*y[end]) / 48
    for i in 5 : length(y) - 1
        retval += y[i]
    end
    return (x[2] - x[1]) * retval
end

function integrate(x::AbstractVector, y::AbstractVector, ::SimpsonEvenFast)
    @inbounds retval = 17*y[1] + 59*y[2] + 43*y[3] + 49*y[4]
    @inbounds retval += 49*y[end-3] + 43*y[end-2] + 59*y[end-1] + 17*y[end]
    retval /= 48
    @fastmath @inbounds for i in 5 : length(y)-1
        retval += y[i]
    end
    @inbounds return (x[2] - x[1]) * retval
end

integrate(x::AbstractVector, y::AbstractVector) = integrate(x, y, TrapezoidalFast())

end
