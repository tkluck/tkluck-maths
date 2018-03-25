using PolynomialRings
using MatrixFactorizations

R = @ring! ℚ[x,y,z]

A = unit_matrix_factorization(x^3, [:x], [:y])
B = unit_matrix_factorization(y^3, [:y], [:z])

∇W = diff.(y^3, [:y])
Jacobian = @ring(ℚ[y])/Ideal(∇W...)

X = A⨶B
Q = X ⊗ Jacobian

λ = eye(Int,2)⨷diff.(B, :y)//(-3)

# the differential operator taking ``R[var]`` into ``Ω^1_{R[var^exp] / R}``.
function reldiff(f, var, exp)
    res = zero(f)
    varval = typeof(f)(var)
    for ((n,),c) in expansion(f, var)
        d,r = divrem(n, exp)
        r == 0 || throw(ValueError("reldiff needs multiples for $var^$exp"))
        d == 0 && continue
        res += c*d*varval^((d-1)*exp)
    end
    res
end
function matrix_over_t(X, var, exp)
    varval = eltype(X)(var)
    blocks = map(X) do f
        res = zeros(typeof(f), exp, exp)
        for ((n,),c) in expansion(f, var)
            for i = 0:exp-1
                j = mod(i+n, exp)
                res[j+1,i+1] += c*varval^(n + i - j)
            end
        end
        res
    end
    return MatrixFactorizations.flatten_blocks(blocks)
end
∇(v) = reldiff.(v, :y, 2)

XX = matrix_over_t(X, :y, 2)
At = -∇(XX)(y=0)

iszero(XX*At + At*XX)

e = -(λ⊗Jacobian) * At

# e is a morphism
Q*e == e*Q

# e is idempotent up to homotopy
hhh = (matrix_over_t(λ, :y, 2)*∇( matrix_over_t(λ, :y, 2) * ∇(matrix_over_t(X, :y, 2))))(y=0)
e^2 == e + Q*hhh + hhh*Q

"""
Newton's method: define e1 = e + Q*b + b*Q, and solve for b s.t.

e1^2 - e1 = o(b)

This gives

{e - 1//2, b} = hhh

(anti-commutator on the LHS).
"""
function solve_anti_commutator(A,h)
    basis = []
    for i in eachindex(A)
        b = zeros(eltype(A), size(A)...)
        b[i] = one(eltype(A))
        push!(basis, b)
    end
    image_basis = map(b->A*b + b*A, basis)

    gr, T = gröbner_transformation(image_basis)

    factors, r = divrem(h, gr)
    if !iszero(r)
        return nothing
    else
        return (factors*T)*basis
    end
end

b1 = solve_anti_commutator(e - one(e)//2, hhh)
e1 = e + Q*b1 + b1*Q
e1^2  == e1
"""
WOW! We only need one iteration of Newton's method.

e1's image is spanned by the vectors:
"""
Im = hcat(gröbner_basis([e1[:,i] for i in indices(e1,2)])...)
"""
and expressed on this basis, e1 is the identity
"""
e1*Im == Im
"""
The interesting one is Q:
"""
QQ=[
0                  0           -x^2      -3z;
0                  0           -(x+z)//3   1;
x-z            3x*z - 3z^2      0          0;
(x^2-z^2)//3   (x^2*z - x^3)    0          0;
]

Q*Im == Im*QQ

QQ^2
"""
Yay! So QQ is the finite dimensional representative of the homotopy class of
A⨶B.

Next step: find the "unitor", the map identifying A with QQ (since B is the unit
matrix factorization). Note that the dimensions are different, so I'm guessing
that identification up to homotopy is good enough?
"""
block_diagonalization(QQ)
"""
So it's A ⊕ "trivial factorization". Am I okay with that? E.g.
is the trivial factorization contractible (homotopy equiv to zero)?
"""

# Knörrer perdiodicity
S = @ring! ℚ[u,v]
K = [0 u; v 0]
Kdual = [0 v; -u 0]

Y = A ⨶ K
YY = Y ⨶ Kdual

Jac = S/Ideal(u,v)
YYY = YY ⊗ Jac

λ1 = -eye(Y) ⨷ diff(Kdual, :v)
λ2 = -eye(Y) ⨷ diff(Kdual, :u)

M = matrix_over_t(matrix_over_t(YY, :u, 1), :v, 1)
Att = (reldiff.(M, :u, 1)*reldiff.(M, :v, 1) - reldiff.(M, :v, 1)*reldiff.(M, :u, 1))//2

ee = -λ1 * λ2 * Att

ee^2 == ee
iszero(ee*YYY - YYY*ee)
