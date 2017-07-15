using Base.Test

@testset "QuantumDimension" begin
    using OrbifoldEquivalence
    using PolynomialRings

    R,(x,y) = polynomial_ring(Rational{Int}, :x, :y)

    Q = [0 x-y; x^2 + x*y + y^2 0]

    @test quantum_dimension(Q,x^3-y^3, [:x], [:y]) == 1
end
