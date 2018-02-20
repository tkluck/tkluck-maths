using Base.Test

using PolynomialRings

@testset "FlintNumbers" begin
    using FlintNumbers
    @testset "Construction and conversion" begin
        R = @ring! ℚ[α]

        S = R/Ideal(α^2 - 2)

        @test S(α)^2 - 2 == 0

        Q = FlintNumberField(S)
        @test Q(α)^2 - 2 == 0

        @test Q(1+α) // Q(1+α) == 1
        @test 1 // Q(1+α) == -Q(1-α)

        A = @ring! Q[x]
        @test (α*x)^2 == 2x^2

        @test_throws ArgumentError @ring Q[α]

        R = @ring! Q[γ]
        S = R/Ideal(γ^2 - α)
        #Q = FlintNumberField(S)

        #@test Q(γ)^2 == α
        #@test Q(γ)^3 == α*γ
        #@test Q(γ)^4 == 2

        B = @ring! ℚ[α, x]
        @test (α*x)^2 != 2x^2
        @test A(α*x)^2 == A(2x^2)

        R = @ring! ℚ[α,β]
        S = R/Ideal(α^2 - 2, β^3 - 2)
        Q = FlintNumberField(S)
        @test Q(α^2) == Q(β^3) == 2

        #RR = @ring! Q[β]

        #SS = RR/Ideal(β^3 - 2)

        #QQ = NumberField(SS)

        #@test QQ(β + α) == QQ(β) + Q(α)
    end
end
