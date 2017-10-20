using Base.Test

using PolynomialRings

@testset "CommutativeAlgebras" begin
    using CommutativeAlgebras
    @testset "Construction and conversion" begin
        R = @ring ℚ[α]

        S = R/Ideal(α^2 - 2)

        @test S(α)^2 - 2 == 0

        Q = NumberField(S)
        @test Q(α)^2 - 2 == 0

        A = @ring Q[x]
        @test (α*x)^2 == 2x^2

        @test_throws ArgumentError @ring Q[α]

        B = @ring ℚ[α, x]
        @test (α*x)^2 != 2x^2
        @test A(α*x)^2 == A(2x^2)

        R = @ring ℚ[α,β]
        S = R/Ideal(α^2 - 2, β^3 - 2)
        Q = NumberField(S)
        @test Q(α^2) == Q(β^3) == 2

        #RR = @ring Q[β]

        #SS = RR/Ideal(β^3 - 2)

        #QQ = NumberField(SS)

        #@test QQ(β + α) == QQ(β) + Q(α)

    end

end
