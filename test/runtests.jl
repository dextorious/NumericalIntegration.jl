using NumericalIntegration
using Test
using InteractiveUtils # for subtypes
using StaticArrays
using HCubature # for testing n-dimensional integration

@testset "compare with analytic result" begin
    x = collect(-π : 2*π/2048 : π)
    y = sin.(x)
    p = collect(0 : π/1024 : π)
    q = sin.(p)
    for M in subtypes(IntegrationMethod)
        for T in [Float32, Float64, BigFloat]
            for (xs,ys,val,atol) in [
                                     (x,y,0,1e-4),
                                     (p,q,2,1e-4),
                                  ]
                result = @inferred integrate(T.(xs), T.(ys),M())
                @test isapprox(result, val, atol=atol)
                @test typeof(result) == T
            end
        end
    end
end

@testset "cumulative integration" begin
    x = collect(-π : 2*π/2048 : π)
    y = sin.(x)
    exact = @. -cos(x) + cos(-π)
    for M in subtypes(IntegrationMethod)
        hasmethod(cumul_integrate, Tuple{AbstractVector, AbstractVector, M}) || continue
        for T in [Float32, Float64, BigFloat]
            result = @inferred cumul_integrate(T.(x), T.(y), M())
            @test typeof(result) == Vector{T}
            @test all(isapprox.(result, exact, atol=1e-4))
        end
    end
end

@testset "n-dimensional integration testing" begin
    X = range(0,stop=2π,length=10)
    Y = range(-π,stop=π,length=10)
    Z = range(0,stop=2,length=10)

    A = Array{Float64}(undef,length(X),length(Y),length(Z))
    f(x) = sin(x[1]) + cos(x[2]) + 2x[3]

    for (k,z) in enumerate(Z)
        for (j,y) in enumerate(Y)
            for (i,x) in enumerate(X)
                A[i,j,k] = f([x,y,z])
            end
        end
    end

    @test isapprox(integrate((X,Y,Z),A), hcubature(f,[0,-π,0],[2π,π,2])[1], atol=1e-4)
end

@testset "SVector" begin
    xs = range(0, stop=1, length=9)
    ys1 = randn(9)
    ys2 = randn(9)
    ys = map(SVector, ys1, ys2)
    for M in subtypes(IntegrationMethod)
        m = M()
        res1 = integrate(xs,ys1,m)
        res2 = integrate(xs,ys2,m)
        res = @inferred integrate(xs,ys,m)
        @test res ≈ @SVector [res1, res2]
    end
end

@testset "Raising Warnings" begin
    xs = collect(-1.0 : 0.5 : 1.0)
    ys = xs.^2
    m = RombergEven()
    expwarn = "RombergEven :: final step reached, but accuracy not: 1.3333333333333335 > 1.0e-12"
    @test_logs (:warn, expwarn) integrate(xs, ys, m)
end

using Unitful

@testset "Unitful compatibility" begin
    x = LinRange(1u"s", 10u"s", 9)
    y = rand(9)*u"m/s"
    for M in subtypes(IntegrationMethod)
        @test typeof(integrate(x,y,M())) == typeof(1.0u"m")
        if hasmethod(cumul_integrate, Tuple{AbstractVector, AbstractVector, M})
            @test typeof(cumul_integrate(x,y,M())) == Vector{typeof(1.0u"m")}
        end
    end
end
