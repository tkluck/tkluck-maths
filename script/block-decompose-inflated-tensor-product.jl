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

iszero(e^2 - e)

iszero(Q*e - e*Q)

block_diagonalization(Q)

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
