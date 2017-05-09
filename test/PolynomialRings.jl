using Base.Test
using PolynomialRings: polynomial_ring, expansion

@testset "PolynomialRings" begin

    R,(x,y) = polynomial_ring(Rational{Int}, :x, :y)
    S,(z,) = polynomial_ring(Int, :z)


    @test x != 0
    @test 0*x == 0
    @test 1*x != 0
    @test 1*x == x
    @test 2*x == x+x
    @test 2*x*y == y*x*2
    @test x*y*x == x^2*y
    @test x-x == 0

    @test x*z == z*x
    @test x*y*z == x*z*y

    @test (x+z)*(x-z) == x^2 - z^2

    @test expansion(x*y*z + x*z + z^2, :z) == [(z, x*y + x), (z^2, 1)]

end
