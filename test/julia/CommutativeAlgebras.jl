using Base.Test

using PolynomialRings

@testset "CommutativeAlgebras" begin
    using CommutativeAlgebras
    @testset "Construction" begin
        R = @ring ℚ[α]

        S = R/Ideal(α^2 - 2)

        @test S(α)^2 - 2 == 0
    end
end
