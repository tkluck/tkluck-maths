using Base.Test
using Polynomials: Polynomial, Monomial, Exponent

@testset "Polynomials" begin

    x = Polynomial{Int, 2}([Monomial(Exponent{2}((1,0)), 1)])
    y = Polynomial{Int, 2}([Monomial(Exponent{2}((0,1)), 1)])


    @test x != 0
    @test 0*x == 0
    @test 1*x != 0
    @test 1*x == x
    @test 2*x == x+x
    @test 2*x*y == y*x*2
    @test x*y*x == x^2*y
    @test x-x == 0

end
