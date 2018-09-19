using NumericalIntegration
using Test
using InteractiveUtils # for subtypes
using StaticArrays

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
    expwarn = "RombergEven :: final step reached, but accuracy not: 1.0 > 1.0e-12"
    @test_logs (:warn, expwarn) integrate(xs, ys, m)
end
