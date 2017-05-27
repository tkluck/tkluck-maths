using Base.Test
using PolynomialRings: polynomial_ring
using PolynomialRings.Groebner: red

@testset "Groebner" begin

    R,(x,y) = polynomial_ring(Rational{Int}, :x, :y)
    S,(z,) = polynomial_ring(Int, :z)

    @test red(x^2, [x]) == (0, [x]')

    @test red(x + y, [x]) == (y, [1]')

    @test begin
        f, G = (x^23 + y -x*43, [x^3*y^4, x^7])
        f_red, factors = red(f,G)
        f == f_red + (factors * G)[1]
    end

end
