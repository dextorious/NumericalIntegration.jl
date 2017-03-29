__precompile__()

module NumericalIntegration

abstract IntegrationMethod
immutable Trapezoidal         <: IntegrationMethod end
immutable TrapezoidalEven     <: IntegrationMethod end
immutable TrapezoidalFast     <: IntegrationMethod end
immutable TrapezoidalEvenFast <: IntegrationMethod end
immutable SimpsonEven         <: IntegrationMethod end # https://en.wikipedia.org/wiki/Simpson%27s_rule#Alternative_extended_Simpson.27s_rule
immutable SimpsonEvenFast     <: IntegrationMethod end

integrate{T<:AbstractFloat}(x::Vector{T}, y::Vector{T}, method::IntegrationMethod=Trapezoidal()) = integrate(x, y, method)

function integrate{T<:AbstractFloat}(x::Vector{T}, y::Vector{T}, ::Trapezoidal)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    retval = zero(eltype(x))
    for i in 1 : length(y)-1
        retval += (x[i+1] - x[i]) * (y[i] + y[i+1])
    end
    return 0.5 * retval
end

function integrate{T<:AbstractFloat}(x::Vector{T}, y::Vector{T}, ::TrapezoidalEven)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    return 0.5 * (x[end] - x[1]) / (length(y) - 1) * (y[1] + y[end] + sum(y[2:end-1]))
end

function integrate{T<:AbstractFloat}(x::Vector{T}, y::Vector{T}, ::TrapezoidalFast)
    retval = zero(eltype(x))
    @fastmath @simd for i in 1 : length(y)-1
        @inbounds retval += (x[i+1] - x[i]) * (y[i] + y[i+1])
    end
    return 0.5 * retval
end

function integrate{T<:AbstractFloat}(x::Vector{T}, y::Vector{T}, ::TrapezoidalEvenFast)
    retval = zero(eltype(x))
    N = length(y) - 1
    @fastmath @simd for i in 2 : N
        @inbounds retval += y[i]
    end
    @inbounds return (x[end] - x[1]) / N * (retval + 0.5*y[1] + 0.5*y[end])
end

function integrate{T<:AbstractFloat}(x::Vector{T}, y::Vector{T}, ::SimpsonEven)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    retval = (17*y[1] + 59*y[2] + 43*y[3] + 49*y[4] + 49*y[end-3] + 43*y[end-2] + 59*y[end-1] + 17*y[end]) / 48
    for i in 5 : length(y) - 1
        retval += y[i]
    end
    return (x[end] - x[1]) / (length(y) - 1) * retval
end

function integrate{T<:AbstractFloat}(x::Vector{T}, y::Vector{T}, ::SimpsonEvenFast)
    @inbounds retval = 17*y[1] + 59*y[2] + 43*y[3] + 49*y[4]
    @inbounds retval += 49*y[end-3] + 43*y[end-2] + 59*y[end-1] + 17*y[end]
    retval /= 48
    @fastmath @inbounds for i in 5 : length(y)-1
        retval += y[i]
    end
    @inbounds return (x[end] - x[1]) / (length(y) - 1) * retval
end

export integrate
export Trapezoidal, TrapezoidalEven, TrapezoidalFast, TrapezoidalEvenFast
export SimpsonEven, SimpsonEvenFast

end
