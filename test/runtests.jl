using NumericalIntegration
using Base.Test

methods = [Trapezoidal(), TrapezoidalEven(), TrapezoidalFast(), TrapezoidalEvenFast(), SimpsonEven(), SimpsonEvenFast()]
x = collect(-π : π/1000 : π)
y = sin(x)
for method in methods
    println(string("Testing method: ", typeof(method)))
    @test abs(integrate(x, y, method)) < 1e-4
end
