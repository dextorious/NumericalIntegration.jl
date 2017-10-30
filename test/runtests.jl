using NumericalIntegration
using Base.Test

methods = [Trapezoidal(), TrapezoidalEven(), TrapezoidalFast(), TrapezoidalEvenFast(), SimpsonEven(), SimpsonEvenFast()]
x = collect(-π : π/1000 : π)
y = sin.(x)
p = collect(0 : π/1000 : π)
q = sin.(p)
for method in methods
    println(string("Testing method: ", typeof(method)))
    @test abs(integrate(x, y, method)) < 1e-4
    @test abs(integrate(p, q, method)-2.0) < 1e-4
end
