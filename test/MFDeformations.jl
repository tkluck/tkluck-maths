using PolynomialRings: polynomial_ring
using MFDeformations: diff, H1, deformation

using Base.Test

@testset begin

    A, (x,y,z) = polynomial_ring(Int, :x, :y, :z)

    d1 = [z   y^2 x^3 0 ; y -x*z 0 x^3 ; x 0 -x*z -y^2 ; 0 x -y z  ]
    d0 = [x*z y^2 x^3 0 ; y -z   0 x^3 ; x 0 -z   -y^2 ; 0 x -y x*z]
    zz = zero(d1)
    Q = [ zz d1; d0 zz ]

    dQ, dQ_even, dQ_odd = diff(Q)
    H = H1(Q)

    @test length(H) == 5

    A, (a,b,c,d,e) = polynomial_ring(Int, :a, :b, :c, :d, :e)

    @test deformation(Q, :a, :b, :c, :d, :e)(a=0, b=0, c=0, d=0, e=0) == Q

end
