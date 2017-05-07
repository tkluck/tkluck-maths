using Polynomials

x,y,z = Polynomials.variables(Polynomials._Polynomial{Int, 3})

println(x*y*z)

d1 = [z   y^2 x^3 0 ; y -x*z 0 x^3 ; x 0 -x*z -y^2 ; 0 x -y z  ]
d0 = [x*z y^2 x^3 0 ; y -z   0 x^3 ; x 0 -z   -y^2 ; 0 x -y x*z]
zz = zeros(typeof(x), 4,4)

Q = [ zz d1; d0 zz ]

display(Q)
println()
display(Q^2)
println()

using Modules

delta(Q, i,j) = begin d=zeros(Q); d[i,j]=one(eltype(Q)); d end
grading(i,j) = (
	(i in 1:4 && j in 1:4) || (i in 5:8 && j in 5:8)
	? 0: 1
)

dQ = [ Q*delta(Q,i,j) - (-1)^grading(i,j)*delta(Q,i,j)*Q for i in indices(Q,1), j in indices(Q,2) ]

Polynomials.iszero(zero(typeof(x)))
Polynomials.iszero.(dQ)

dQ = HomspaceMorphism([ dQ[i,j][k,l] for i=1:8,j=1:8,k=1:8,l=1:8])

K = kernel(dQ)
display( K )
println()

display( Q*K[1] - K[1]*Q ); println()
display( Q*K[end] + K[end]*Q ); println()


import Polynomials: monomials_not_in_ideal
display( monomials_not_in_ideal(typeof(x)[ x^3, y^2 , z^3, x*z]) ); println()
# 1 x x^2 y z z^2 xy yx^2 yz yz^2

import Modules: span
import Polynomials: _reduce, iszero

b,transformation = minimal_groebner_basis(span(dQ))

H1 = filter(K) do k
    (k_red, factors) = _reduce(k, b)
    !iszero(k_red)
end

foreach(H1) do h1
    display(h1); println()
end


Q1 = H1[1]

Q2, obs = lift_and_obstruction(dQ, Q*Q1 + Q1*Q)

display(obs); println()
display(Q2); println()
