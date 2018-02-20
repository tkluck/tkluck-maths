module MatrixFactorizations

using PolynomialRings
using PolynomialRings: basering, variablesymbols
using PolynomialRings.QuotientRings: QuotientRing, representation_matrix

import PolynomialRings: ⊗

include("MatrixFactorizations/Library.jl")

"""
    M = flatten_blocks(X)

Construct a matrix from a matrix of matrices by concatenating
them horizontally and vertically.
"""
flatten_blocks(X) = vcat([hcat(X[i,:]...) for i=1:size(X,1)]...)

from_alternating_grades(M::Matrix) = [
    M[1:2:end,1:2:end] M[1:2:end,2:2:end];
    M[2:2:end,1:2:end] M[2:2:end,2:2:end];
]
const to_alternating_grades = from_alternating_grades

"""
    A⨷B

Graded tensor product of ℤ/2 graded block matrices. The result consists
of even/odd blocks just like the input and the operation is associative.

_Graded_ means that we introduce Koszul signs. Writing ``A = a^i_j e_i e^j`` and
``B = b^k_ℓ f_k f^ℓ`` for some (co)basis vectors, we represent the tensor
product as

``A⨷B = (-1)^{|j|(|k|+|ℓ|)} p^i_j q^k_ℓ e_i f_k f^ℓ e^j``

where the sign comes from commuting ``e^j`` with ``f_k f^ℓ``. This ensures that,
upon matrix multiplication between two such tensor products, the contractions
between ``e^j`` and ``e_j``, and between ``f^ℓ`` and ``f_ℓ``, are already
adjacent and do not introduce more signs.

Mapping the multi-indices ``(i,k)`` to rows and ``(ℓ,j)`` to columns requires
some care to ensure that we end up with even/odd blocks. A comment in the code
explains how this is done. We first interlace the rows/columns to get alternating
grades in rows/columns (`to_alternating_grades` and `from_alternating_grades`).
Then we need only to reverse some row/column orders. Finally, we de-interlace
to get even/odd blocks. (I wonder if there is a way to avoid making the detour
through the alternating signs representation, but it seems hard to maintain
associativity.)
"""
function ⨷(A,B)
    n,m = size(A)
    k,l = size(B)
    n%2 == m%2 == k%2 == l%2 == 0 || throw(ValueError("Need ℤ/2 graded matrices, even rank and corank"))

    A,B = to_alternating_grades.((A,B))

    res = [a.*B for a in A]

    for row in indices(res,1), col in indices(res,2)
        inner = res[row,col]
        for inner_row in indices(inner,1), inner_col in indices(inner,2)
            # Koszul signs; see formula in the documentation string
            if col%2 == 0 && (inner_row+inner_col)%2 == 1
                inner[inner_row, inner_col] *= -1
            end
        end
        # make sure gradings alternate in the output matrix: they alternate
        # in both inner and outer, and the result grading is their sum. At a
        # boundary between two 'inner' matrices, both signs change, so the
        # adjacent rows/columns have the same sign. We can fix this by reversing
        # the rows/columns within an inner matrix, depending on the sign of the
        # outer row/column.
        if row%2 == 0
            inner[:,:] = inner[end:-1:1,:]
        end
        if col%2 == 0
            inner[:,:] = inner[:,end:-1:1]
        end
    end
    return from_alternating_grades(flatten_blocks(res))
end

"""
    A⨶B

Tensor product of matrix factorizations.
"""
⨶(A,B) = A⨷eye(B) + eye(A)⨷B

"""
    X⊗quotient_ring

Matrix-representation of the operator obtained from ``X``, a matrix acting on a
module over a ring ``R``, by tensoring (over ``R``) with a quotient of ``R``.

Note that `eltype(X)` is not necessarily equal to ``R``; it may also be
an ``R``-algebra. For example, ``R=k[y]`` and `eltype(X) == @ring(k[x,y])` works.
"""
function ⊗(X::AbstractMatrix{<:Polynomial}, quotient_ring::Type{<:QuotientRing})
    X_inflated = representation_matrix.(quotient_ring, X)
    flatten_blocks(X_inflated)
end


"""
    Q,ϵ = ⨶(A,B,W,vars...)

Finite-rank homotopy representation of ``A⨶B``, where we remove the variables
`vars` from the result.

We return a matrix factorization ``Q`` together with an idempotent ϵ representing
the direct summand of ``Q`` that represents ``A⨶B``.

See the pushforward paper by Dykerhoff&Murfet.
"""
function ⨶(A,B,W,vars...)
    R,_ = polynomial_ring(vars...;basering=basering(W))
    ∇W = diff.(W, vars)
    Jacobian = R/Ideal(∇W...)

    Q = (A⨶B) ⊗ Jacobian
end

"""
    unit_matrix_factorization(f, left_vars, right_vars)

A ℤ/2-graded matrix that squares to `f - f(;left_vars => right_vars)` times
the identity matrix.

The source for this formulation is

> Adjunctions and defects in Landau-Ginzburg models, Nils Carqueville and Daniel Murfet
"""
function unit_matrix_factorization(f, left_vars, right_vars)
    R = typeof(f)
    function ∂(f, n)
        for i in 1:n-1
            f = f(; left_vars[i] => R(right_vars[i]))
        end
        factors = div(R(f - f(; left_vars[n] => R(right_vars[n]))), [R(left_vars[n]) - R(right_vars[n])])
        factors[1]
    end

    # x represents, through its bit-representation, a basis element of the exterior
    # algebra. To be precise, x represents the element theta_i_1 \wedge ... \wedge theta_i_n
    # where i_1 ... i_n are the bits set in x.
    #
    # The use of 'gray code' (see wikipedia) ensures that subsequent elements differ by
    # exactly one bit. This way, rows/columns of our result matrix have _alternating_ signs.
    N = length(left_vars)
    gray_code(x) = xor(x, x>>1)
    permutation = map(n->gray_code(n)+1, 0:2^N-1)
    inv_perm = invperm(permutation)
    to_index(x) = inv_perm[x+1]

    function wedge_product_matrix(T, i)
        result = zeros(T, 2^N,2^N)
        for j in 0:2^N-1
            j&(1<<(i-1)) != 0 && continue

            k = j | (1 << (i-1))
            sign = (-1)^count_ones(j & (1<<i - 1))
            result[to_index(j), to_index(k)] = T(sign)
        end
        return result
    end

    function lift_matrix(T, i)
        result = zeros(T, 2^N, 2^N)
        for j in 0:2^N-1
            j&(1<<(i-1)) == 0 && continue

            k = j & ~(1<<(i-1))
            sign = (-1)^count_ones(j & (1<<i - 1))
            result[to_index(j), to_index(k)] = T(sign)
        end
        return result
    end

    delta_plus = sum( ∂(f, i) * wedge_product_matrix(R, i) for i=1:N )
    delta_minus = sum( (R(left_vars[i]) - R(right_vars[i])) * lift_matrix(R, i) for i = 1:N )
    return from_alternating_grades(delta_plus + delta_minus)
end

"""
    D, A = block_diagonalization(X)

Decompose the matrix factorization X into a direct sum of irreducible
matrix factorizations, and represent this direct sum as a block-diagonal
matrix factorization `D` such that ``A^{-1} D A = X``.

NOTE: for now, this function only returns D, and not yet A!
"""
function block_diagonalization(X)
    D = copy(X)
    A = eye(X)
    top_right = @view D[1:end÷2, end÷2+1:end]

    # the following functions take parameters indexing `top_right`, but they
    # operate simultaneously on the bottom left of `D` as well. This happens
    # in such a way that X^2 remains the same.
    function rowop(i, factor, j)
        D[i,:] += factor * D[j,:]
        D[:,j] -= factor * D[:,i]
    end
    function colop(i, factor, j)
        D[:,end÷2 + i] += factor * D[:,end÷2 + j]
        D[end÷2 + j,:] -= factor * D[end÷2 + i,:]
    end

    for _=1:100
        for row in indices(top_right,1), col in indices(top_right,2)
            for row2 in indices(top_right,1)
                row2 == row && continue
                iszero(top_right[row2,col]) && continue
                (d,),r = divrem(top_right[row2,col], [top_right[row,col]])
                if iszero(r)
                    rowop(row2, -d, row)
                    @goto more_loops
                end
            end
            for col2 in indices(top_right,2)
                col2 == col && continue
                iszero(top_right[row,col2]) && continue
                (d,),r = divrem(top_right[row,col2], [top_right[row,col]])
                if iszero(r)
                    colop(col2, -d, col)
                    @goto more_loops
                end
            end
        end
        break
        @label more_loops
    end

    # most inefficient algorithm I can think of
    blocks = []
    for row in indices(top_right,1), col in indices(top_right,2)
        if !iszero(top_right[row, col])
            push!(blocks, (Set([row]), Set([col])))
        end
    end
    for i in 1:length(blocks)
        j = i+1
        while j <= length(blocks)
            if !isempty(blocks[i][1] ∩ blocks[j][1]) || !isempty(blocks[i][2] ∩ blocks[j][2])
                blocks[i] = (blocks[i][1] ∪ blocks[j][1], blocks[i][2] ∪ blocks[j][2])
                deleteat!(blocks, j)
                j = i+1
            else
                j += 1
            end
        end
    end

    rowperm = vcat(sort.(collect.(getindex.(blocks, 1)))...)
    colperm = vcat(sort.(collect.(getindex.(blocks, 2)))...)
    D[1:end÷2,:] = D[rowperm,:]
    D[:,end÷2+1:end] = D[:,colperm .+ end÷2]

    D[end÷2+1:end,:] = D[colperm .+ end÷2,:]
    D[:,1:end÷2] = D[:,rowperm]

    #D,A
    D
end

export ⨷, ⨶
export unit_matrix_factorization
export block_diagonalization

end
