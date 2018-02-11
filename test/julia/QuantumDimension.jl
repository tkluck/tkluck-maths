using Base.Test

using OrbifoldEquivalence
using PolynomialRings

@testset "QuantumDimension" begin
    @ring! ℚ[x,y]

    Q = [0 x-y; x^2 + x*y + y^2 0]

    @test quantum_dimension(Q,x^3-y^3, [:x], [:y]) == 1
end
