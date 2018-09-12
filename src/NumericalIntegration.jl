module NumericalIntegration

export integrate, cintegrate
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


#documentation

"""
    integrate(x,y...)

Compute numerical integral of y(x) from x=x[1] to x=x[end]. Return a scalar of the same type as the input. If not method is supplied, use TrapezdoialFast.
"""
function integrate(x,y...) end


"""
    cintegrate(x,y...)

Compute cumulative numerical integral of y(x) from x=x[1] to x=x[end]. Return a vector with elements of the same type as the input. If not method is supplied, use TrapezdoialFast.
"""
function cintegrate(x,y...) end

#implementation

function _zero(x,y)
    ret = zero(eltype(x)) + zero(eltype(y))
    ret / 2
end

function _zeros(x::AbstractVector,y::AbstractVector)
    ret = zeros(eltype(x),size(x)) + zeros(eltype(y),size(y))
    ret / 2
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
    integrate(x::AbstractVector, y::AbstractArray, method)

When y is an array, compute integral for each column.
"""
function integrate(x::AbstractVector, y::AbstractArray, M::IntegrationMethod)
    out = [integrate(x,y[:,j],M) for j=1:size(y,2)]
    return out
end


# cumulative integrals

"""
    cintegrate(x::AbstractVector, y::AbstractVector, ::Trapezoidal)

Use Trapezoidal rule.
"""
function cintegrate(x::AbstractVector, y::AbstractVector, ::Trapezoidal)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    retarr = _zeros(x,y)
    for i in 2 : length(y)
        retarr[i] = retarr[i-1] + (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    return HALF * retarr
end

"""
    cintegrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEven)

Use Trapezoidal rule, assuming evenly spaced vector x.
"""
function cintegrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEven)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    retarr = _zeros(x,y)
    for i in 2 : length(y)
        retarr[i] = retarr[i-1] + (y[i-1] + y[i])
    end
    return (x[2] - x[1]) * HALF * retarr
end

"""
    cintegrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalFast)

Use Trapezoidal rule. Unsafe method: no bound checking. This is the default when no method is supplied.
"""
function cintegrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalFast)
    retarr = _zeros(x,y)
    @fastmath for i in 2 : length(y) #not sure if @simd can do anything here
        @inbounds retarr[i] = retarr[i-1] + (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    return HALF * retarr
end

"""
    cintegrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEvenFast)

Use Trapezoidal rule, assuming evenly spaced vector x. Unsafe method: no bound checking.
"""
function cintegrate(x::AbstractVector, y::AbstractVector, ::TrapezoidalEvenFast)
    retarr = _zeros(x,y)
    @fastmath for i in 2 : length(y)
        @inbounds retarr[i] = retarr[i-1] + (y[i] + y[i-1])
    end
    @inbounds return (x[2] - x[1]) * HALF * retarr
end

"""
    cintegrate(x::AbstractVector, y::AbstractArray, method)

When y is an array, compute integral for each column.
"""
function cintegrate(x::AbstractVector, y::AbstractArray, M::IntegrationMethod)
    return hcat([cintegrate(x,y[:,j],M) for j=1:size(y,2)]...)
end


#default behaviour
integrate(x::AbstractVector, y::AbstractVector) = integrate(x, y, TrapezoidalFast())

integrate(x::AbstractVector, y::AbstractArray) = integrate(x,y,TrapezoidalFast())

cintegrate(x::AbstractVector, y::AbstractVector) = cintegrate(x, y, TrapezoidalFast())

cintegrate(x::AbstractVector, y::AbstractArray) = cintegrate(x,y,TrapezoidalFast())

end
