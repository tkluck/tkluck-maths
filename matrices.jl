using PolynomialRings: polynomial_ring
using MFDeformations: diff, H1, deformation

StrangeDuality() = begin
    A, (x,y,z) = polynomial_ring(Int, :x, :y, :z)

    d1 = [z   y^2 x^3 0 ; y -x*z 0 x^3 ; x 0 -x*z -y^2 ; 0 x -y z  ]
    d0 = [x*z y^2 x^3 0 ; y -z   0 x^3 ; x 0 -z   -y^2 ; 0 x -y x*z]
    zz = zero(d1)

    [ zz d1; d0 zz ]
end

S1() = begin
    A, (x,y) = polynomial_ring(Int, :x, :y)

    q1=[x x*y; -y -x^3]
    q0=[x^3 x*y; -y -x]
    zz = zero(q1)

    [ zz q1; q0 zz ]
end

T2() = begin
    A, (x,y) = polynomial_ring(Int, :x, :y)

    q1 = [x^2 y ; -y -x^2]
    q0 = [x^3 x*y; -x*y -x^3]
    zz = zero(q1)

    [ zz q1; q0 zz ]
end


println("Finding deformations for:")
Q = T2()
display(Q); println()

A, (t,u,v,w) = polynomial_ring(Int, :t, :u, :v, :w)

Qdef =deformation(Q, :t, :u, :v, :w)
println("Deformation is:")
display(Qdef); println()
println("which squares to:")
display(Qdef^2); println()


