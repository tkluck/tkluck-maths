using PolynomialRings
using MatrixFactorizations

R = @ring! ℚ[x,y,z]

A = unit_matrix_factorization(x^3, [:x], [:y])
B = unit_matrix_factorization(y^3, [:y], [:z])

∇W = diff.(y^3, [:y])
Jacobian = @ring(ℚ[y])/Ideal(∇W...)

X = A⨶B
Q = X ⊗ Jacobian

λ = eye(Int,2)⨷diff.(B, :y)⊗Jacobian

# the differential operator taking ``R[var]`` into ``Ω^1_{R[var^exp] / R}``.
function reldiff(f, var, exp)
    res = zero(f)
    varval = typeof(f)(var)
    for ((n,),c) in expansion(f, var)
        d,r = divrem(n, exp)
        # writing ``t = var^exp``, we now have
        #     term = c * var^r * t^d
        d == 0 && continue
        # so its differential w.r.t. t is equal to
        #     d(term) = c * var^r * d*t^(d-1) * dt
        res += c*d*varval^(r+(d-1)*exp)
    end
    res
end
function matrix_over_t(X, var, exp)
    varval = eltype(X)(var)
    blocks = map(X) do f
        res = zeros(typeof(f), exp, exp)
        for ((n,),c) in expansion(f, var)
            d,r = divrem(n, exp)
            for i = 0:exp-1
                j = mod(i+r, exp)
                res[j+1,i+1] += c*varval^(n + i - j)
            end
        end
        res
    end
    return MatrixFactorizations.flatten_blocks(blocks)
end
∇(v) = reldiff.(v, :y, 2)

XX = matrix_over_t(X, :y, 2)
At = -∇(XX)

iszero(XX*At + At*XX)

e = -λ * At

iszero(e^2 - e)

iszero(Q*e - e*Q)

block_diagonalization(Q)
