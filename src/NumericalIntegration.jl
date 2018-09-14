module NumericalIntegration

using LinearAlgebra
using Logging

export integrate
export Trapezoidal, TrapezoidalEven, TrapezoidalFast, TrapezoidalEvenFast
export SimpsonEven, SimpsonEvenFast
export RombergEven
export IntegrationMethod

abstract type IntegrationMethod end

struct Trapezoidal         <: IntegrationMethod end
struct TrapezoidalEven     <: IntegrationMethod end
struct TrapezoidalFast     <: IntegrationMethod end
struct TrapezoidalEvenFast <: IntegrationMethod end
struct SimpsonEven         <: IntegrationMethod end # https://en.wikipedia.org/wiki/Simpson%27s_rule#Alternative_extended_Simpson.27s_rule
struct SimpsonEvenFast     <: IntegrationMethod end
struct RombergEven{T<:AbstractFloat}      <: IntegrationMethod
    acc::T
end # https://en.wikipedia.org/wiki/Romberg%27s_method
RombergEven() = RombergEven(1e-12)

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

function integrate(x::AbstractVector, y::AbstractVector, m::RombergEven)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    @assert ((length(x) - 1) & (length(x) - 2)) == 0 "Need length of vector to be 2^n + 1"
    maxsteps::Integer = Int(log2(length(x)-1))
    rombaux = zeros(eltype(y), maxsteps, 2)
    prevrow = 1
    currrow = 2
    @inbounds h = x[end] - x[1]
    @inbounds rombaux[prevrow, 1] = (y[1] + y[end])*h*HALF
    @inbounds for i in 1 : (maxsteps-1)
        h *= HALF
        npoints = 1 << (i-1)
        jumpsize = div(length(x)-1, 2*npoints)
        c = 0.0
        for j in 1 : npoints
            c += y[1 + (2*j-1)*jumpsize]
        end
        rombaux[1, currrow] = h*c + HALF*rombaux[1, prevrow]
        for j in 2 : (i+1)
            n_k = 4^(j-1)
            rombaux[j, currrow] = (n_k*rombaux[j-1, currrow] - rombaux[j-1, prevrow])/(n_k - 1)
        end

        if i > maxsteps//3 && norm(rombaux[i, prevrow] - rombaux[i+1, currrow], Inf) < m.acc
            return rombaux[i+1, currrow]
        end

        prevrow, currrow = currrow, prevrow
    end
    finalerr = norm(rombaux[maxsteps-1, prevrow] - rombaux[maxsteps, currrow], Inf)
    @warn "RombergEven :: final step reached, but accuracy not: $finalerr > $(m.acc)"
    @inbounds return rombaux[maxsteps, prevrow]
end

integrate(x::AbstractVector, y::AbstractVector) = integrate(x, y, TrapezoidalFast())

end
