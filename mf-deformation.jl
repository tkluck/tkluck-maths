using PolynomialRings: polynomial_ring

Ring, (x,y,z,w,u,v,a,b) = polynomial_ring(Int, :x, :y, :z, :w, :u, :v, :a, :b)

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
dQ_even = [ grading(i,j) == 1 ? (Q*delta(Q,i,j) + delta(Q,i,j)*Q) : zero(Q) for i in indices(Q,1), j in indices(Q,2) ]
dQ_odd =  [ grading(i,j) == 0 ? (Q*delta(Q,i,j) - delta(Q,i,j)*Q) : zero(Q) for i in indices(Q,1), j in indices(Q,2) ]

dQ = HomspaceMorphism([ dQ[i,j][k,l] for i=1:8,j=1:8,k=1:8,l=1:8])
dQ_even = HomspaceMorphism([ dQ_even[i,j][k,l] for i=1:8,j=1:8,k=1:8,l=1:8])
dQ_odd = HomspaceMorphism([ dQ_odd[i,j][k,l] for i=1:8,j=1:8,k=1:8,l=1:8])

#K = kernel(dQ)
#display( K )
#println()

#display( Q*K[1] - K[1]*Q ); println()
#display( Q*K[end] + K[end]*Q ); println()


#import Polynomials: monomials_not_in_ideal
#display( monomials_not_in_ideal(typeof(x)[ x^3, y^2 , z^3, x*z]) ); println()
# 1 x x^2 y z z^2 xy yx^2 yz yz^2

import Modules: span
import Polynomials: _reduce, iszero

b,transformation = minimal_groebner_basis(span(dQ_odd))

H1 = filter(kernel(dQ_even)) do k
    (k_red, factors) = _reduce(k, b)
    !iszero(k_red)
end

foreach(H1) do h1
    display(h1); println()
end


#Q1 = Q + sum(eps*h for (eps,h) in zip([w,u,v,a,b],H1))

import Iterators: groupby
import Polynomials: Polynomial, _Monomial, Exponent
function expansion{R<:Number}(p::Polynomial{R,8})
    res = Tuple{Polynomial{R,3}, Polynomial{R,5}}[]
    for coeffs in groupby(c -> c[1].e[4:8], p.coeffs)
        res_coeffs = _Monomial{R,3}[]
        local exp_wuv
        for coeff in coeffs
            exp_xyz = coeff[1].e[1:3]
            exp_wuv = coeff[1].e[4:8]
            push!(res_coeffs, _Monomial((Exponent(exp_xyz), coeff[2])))
        end
        push!(res, (Polynomial(res_coeffs), Polynomial([_Monomial((Exponent(exp_wuv), one(R)))]) ) )
    end
    return res
end

display(expansion(x*y*z*w + w + x*v)); println()

Q2, obs2 = lift_and_obstruction(dQ_even, -Q1*Q1)
Q3, obs3 = lift_and_obstruction(dQ_even, -(Q2*Q1 + Q1*Q2))

println("obstructions:")
display(obs2 + obs3); println()

println("deformation in Q^2:")
display((Q + w * Q1 + w^2 * Q2 + w^3 * Q3)^2 - Q^2); println()

