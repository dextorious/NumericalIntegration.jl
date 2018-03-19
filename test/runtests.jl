using NumericalIntegration
using Base.Test
using StaticArrays

@testset "compare with analytic result" begin
    x = collect(-π : π/1000 : π)
    y = sin.(x)
    p = collect(0 : π/1000 : π)
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

@testset "SVector" begin
    xs  = linspace(0,1,10)
    ys1 = randn(10)
    ys2 = randn(10)
    ys = map(SVector, ys1, ys2)
    for M in subtypes(IntegrationMethod)
        m = M()
        res1 = integrate(xs,ys1,m)
        res2 = integrate(xs,ys2,m)
        res = @inferred integrate(xs,ys,m)
        @test res ≈ @SVector [res1, res2]
    end
end
