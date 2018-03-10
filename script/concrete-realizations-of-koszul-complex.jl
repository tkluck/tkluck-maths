using PolynomialRings
using MatrixFactorizations

R = @ring! ℚ[x,y,z]

A = unit_matrix_factorization(x^3, [:x], [:y])
B = unit_matrix_factorization(y^3, [:y], [:z])

∇W = diff.(y^3, [:y])
Jacobian = @ring(ℚ[y])/Ideal(∇W...)

X = A⨶B

λ = eye(Int,2)⨷diff.(B, :y)//(-3)

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
"""
representing the tensor product with Koszul complex
(single variable only)
"""
function Kosz(X, t)
    d = X⨷eye(eltype(X),2)
    δ = eye(X)⨷[0 t; 0 0]
    d + δ
end

_suspend(e) = e
function _suspend(e::Expr)
    args = map(_suspend, e.args)
    if e.head === :ref && length(args) == 2
        return :( suspend($(args[1]), $(args[2])) )
    else
        res = copy(e)
        res.args = args
        return res
    end
end
macro suspend(expr)
    _suspend(expr)
end

function suspend(X, n)
    iseven(n) && return copy(X)
    isodd(n) && return -X
end
function susp_eye(X)
    @assert size(X,1) == size(X,2)
    T = eltype(X)
    n = size(X,1)÷2
    [zeros(T,n,n) eye(T,n); eye(T,n) zeros(T,n,n)]
end

K = Kosz(X, y^2)

inclusion(X)  = hcat([X[:,i]  for i in indices(X,2) if !iszero(X[:,i])]...)
projection(X) = vcat([X[i,:]' for i in indices(X,1) if !iszero(X[i,:])]...)

"""
    α: X[1] → X
    (represented by identity matrix)
"""
α = eye(X)
@suspend iszero(X*α +α*X[1])

ι = [1 0; 0 0]
dt = [0 0; 1 0]
dt⁻ = [0 1; 0 0]
ϑ′ = inclusion(α⨷dt - (λ*α)⨷ι)
ε = projection(α⨷dt⁻)

# Lemma 4.1
ε * ϑ′ == eye(ε * ϑ′)
@suspend iszero(K*ϑ′ - ϑ′*X[1])
@suspend iszero(ε*K - X[1]*ε)

# projection morphism
π = projection(eye(X)⨷[1 0; 0 0])
iszero(matrix_over_t(π*K - X*π, :y, 2)(y=0))

# composition of π and ϑ′
ϑ = projection(matrix_over_t(π*ϑ′, :y, 2)(y=0))

ϑ == -λ⊗Jacobian
