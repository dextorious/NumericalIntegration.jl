module NumericalIntegration

using LinearAlgebra
using Logging

export integrate, cumul_integrate
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
struct RombergEven{T<:AbstractFloat} <: IntegrationMethod
    acc::T
end # https://en.wikipedia.org/wiki/Romberg%27s_method
RombergEven() = RombergEven(1e-12)

const HALF = 1//2


#documentation

"""
    integrate(x,y...)

Compute numerical integral of y(x) from x=x[1] to x=x[end]. Return a scalar of the same type as the input. If not method is supplied, use TrapezdoialFast.
"""
function integrate(x,y...) end


"""
    cumul_integrate(x,y...)

Compute cumulative numerical integral of y(x) from x=x[1] to x=x[end]. Return a vector with elements of the same type as the input. If not method is supplied, use TrapezdoialFast.
"""
function cumul_integrate(x,y...) end

#implementation

function _zero(x,y)
    ret = zero(eltype(x)) + zero(eltype(y))
end

function _zeros(x::AbstractVector,y::AbstractVector)
    ret = zeros(eltype(x),size(x)) + zeros(eltype(y),size(y))
end

"""
    integrate(x::AbstractVector, y::AbstractVector, ::Trapezoidal)

Use Trapezoidal rule.
"""
function integrate(x::AbstractVector, y::AbstractVector, ::Trapezoidal)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    retval = _zero(x,y)
    for i in 1 : length(y)-1
        retval += (x[i+1] - x[i]) * (y[i] + y[i+1])
    end
    return HALF * retval
end

"""
    integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEven)

Use Trapezoidal rule, assuming evenly spaced vector x.
"""
function integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEven)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    return (x[2] - x[1]) * (HALF * (y[1] + y[end]) + sum(y[2:end-1]))
end

"""
    integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalFast)

Use Trapezoidal rule. Unsafe method: no bound checking. This is the default when no method is supplied.
"""
function integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalFast)
    retval = _zero(x,y)
    @fastmath @simd for i in 1 : length(y)-1
        @inbounds retval += (x[i+1] - x[i]) * (y[i] + y[i+1])
    end
    return HALF * retval
end

"""
    integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEvenFast)

Use Trapezoidal rule, assuming evenly spaced vector x. Unsafe method: no bound checking.
"""
function integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEvenFast)
    retval = _zero(x,y)
    N = length(y) - 1
    @fastmath @simd for i in 2 : N
        @inbounds retval += y[i]
    end
    @inbounds return (x[2] - x[1]) * (retval + HALF*y[1] + HALF*y[end])
end

"""
    integrate(x::AbstractVector, y::AbstractVector, ::SimpsonEven)

Use Simpson's rule, assuming evenly spaced vector x.
"""
function integrate(x::AbstractVector, y::AbstractVector, ::SimpsonEven)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    retval = (17*y[1] + 59*y[2] + 43*y[3] + 49*y[4] + 49*y[end-3] + 43*y[end-2] + 59*y[end-1] + 17*y[end]) / 48
    for i in 5 : length(y) - 1
        retval += y[i]
    end
    return (x[2] - x[1]) * retval
end

"""
    integrate(x::AbstractVector, y::AbstractVector, ::SimpsonEven)

Use Simpson's rule, assuming evenly spaced vector x.  Unsafe method: no bound checking.
"""
function integrate(x::AbstractVector, y::AbstractVector, ::SimpsonEvenFast)
    @inbounds retval = 17*y[1] + 59*y[2] + 43*y[3] + 49*y[4]
    @inbounds retval += 49*y[end-3] + 43*y[end-2] + 59*y[end-1] + 17*y[end]
    retval /= 48
    @fastmath @inbounds for i in 5 : length(y)-1
        retval += y[i]
    end
    @inbounds return (x[2] - x[1]) * retval
end

"""
    integrate(x::AbstractVector, y::AbstractMatrix, method; dims=2)

When y is an array, compute integral along dimension specified by dims (default 2: columns).
"""
function integrate(x::AbstractVector, y::AbstractMatrix, M::IntegrationMethod; dims=2)
    out = [integrate(x,selectdim(y,dims,j),M) for j=1:size(y,dims)]
    return out
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



function integrate(X::Tuple{AbstractVector}, Y::AbstractVector{T}, M::IntegrationMethod) :: T where {T}
    return integrate(X[1], Y, M)
end
"""
    integrate(X::NTuple{N,AbstractVector}, Y::AbstractArray{T,N}, method, cache=nothing)

Given an n-dimensional grid of values, compute the total integral along each dim
"""
function integrate(X::NTuple{N,AbstractVector}, Y::AbstractArray{T,N}, M::IntegrationMethod) :: T where {T,N}
    n = last(size(Y))
    cache = Vector{T}(undef, n)
    x = X[1:N-1]
    @inbounds for i in 1:n
        cache[i] = integrate(x, selectdim(Y,N,i), M)
    end
    return integrate(X[end], cache, M)
end



# cumulative integrals

"""
    cumul_integrate(x::AbstractVector, y::AbstractVector, ::Trapezoidal)

Use Trapezoidal rule.
"""
function cumul_integrate(x::AbstractVector, y::AbstractVector, ::Trapezoidal)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    retarr = _zeros(x,y)
    for i in 2 : length(y)
        retarr[i] = retarr[i-1] + (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    return HALF * retarr
end

"""
    cumul_integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEven)

Use Trapezoidal rule, assuming evenly spaced vector x.
"""
function cumul_integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEven)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    retarr = _zeros(x,y)
    for i in 2 : length(y)
        retarr[i] = retarr[i-1] + (y[i-1] + y[i])
    end
    return (x[2] - x[1]) * HALF * retarr
end

"""
    cumul_integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalFast)

Use Trapezoidal rule. Unsafe method: no bound checking. This is the default when no method is supplied.
"""
function cumul_integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalFast)
    retarr = _zeros(x,y)
    @fastmath for i in 2 : length(y) #not sure if @simd can do anything here
        @inbounds retarr[i] = retarr[i-1] + (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    return HALF * retarr
end

"""
    cumul_integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEvenFast)

Use Trapezoidal rule, assuming evenly spaced vector x. Unsafe method: no bound checking.
"""
function cumul_integrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEvenFast)
    retarr = _zeros(x,y)
    @fastmath for i in 2 : length(y)
        @inbounds retarr[i] = retarr[i-1] + (y[i] + y[i-1])
    end
    @inbounds return (x[2] - x[1]) * HALF * retarr
end

"""
    cumul_integrate(x::AbstractVector, y::AbstractMatrix, method; dims=2)

When y is an array, compute integral along each dimension specified by dims (default 2: columns)
"""
function cumul_integrate(x::AbstractVector, y::AbstractMatrix, M::IntegrationMethod; dims=2)
    return hcat([cumul_integrate(x,selectdim(y,dims,j),M) for j=1:size(y,dims)]...)
end

#default behaviour
integrate(x::AbstractVector, y::AbstractVector) = integrate(x, y, Trapezoidal())

integrate(x::AbstractVector, y::AbstractMatrix; dims=2) = integrate(x, y, Trapezoidal(); dims=dims)

integrate(X::NTuple, Y::AbstractArray) = integrate(X, Y, Trapezoidal())

cumul_integrate(x::AbstractVector, y::AbstractVector) = cumul_integrate(x, y, Trapezoidal())

cumul_integrate(x::AbstractVector, y::AbstractMatrix; dims=2) = cumul_integrate(x, y, Trapezoidal(); dims=dims)

end # module
