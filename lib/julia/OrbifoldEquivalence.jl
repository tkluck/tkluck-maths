module OrbifoldEquivalence

using Base.Iterators
using Nulls

using PolynomialRings
using QuasiHomogeneous

function supertrace(Q::Matrix)
    n,m = size(Q)
    n == m && n%2 ==0 || throw(ArgumentError("Cannot compute supertrace of $n x $m matrix"))

    k = div(n,2)

    return sum(Q[i,i] for i=1:k) - sum(Q[i,i] for i=(k+1):2k)
end

function multivariate_residue(g, f, vars...)
    R = eltype(f)
    G,tr = groebner_basis(f)

    # TODO: compute that R/G is finite dimensional; otherwise, this computation
    # does not terminate
    M = spzeros(R, length(vars), length(f))
    for (row,v) in enumerate(vars)
        x = convert(R,v)
        x_red, factors = red(x,G)
        while !iszero(x_red)
            x = x^2
            x_red, factors = red(x,G)
        end
        M[row,:] = factors*tr
    end

    # workaround: convert M and f to dense matrix/array for the multiplication;
    # for some reason, we get a method error otherwise.
    f_transformed = [m_i for m_i in M] * [f_i for f_i in f]
    g_transformed = g * det(M)

    term, r = divrem(prod(f_transformed), prod(convert(R,v) for v in vars))
    assert(iszero(r))

    return coefficient(g_transformed, term, vars...)
end

function quantum_dimension(Q::Matrix, W, left_vars, right_vars)
    g = supertrace(prod( diff(Q, v) for v in left_vars) * prod( diff(Q, v) for v in right_vars) )
    f = [diff(W,v) for v in left_vars]

    return multivariate_residue(g, f, left_vars...)
end

function equivalence_exists(R, W, Wvars, V, Vvars, rank, a, b)

    R,allvars = polynomial_ring(Rational{BigInt}, Wvars..., Vvars...)

    total_grading, vgr = QuasiHomogeneous.find_quasihomogeneous_degrees(W - V, Wvars..., Vvars...)

    for (next_coeff,Q) in QuasiHomogeneous.generic_matrix_factorizations(rank, a, b, total_grading, vgr, R, :c)

        c1 = take!(next_coeff)
        C = [coefficient(t) for entry in (Q^2 - c1^2*(V-W)*eye(Int,size(Q)...)) for t in entry.p.terms]

        if 1 in C
            continue
        end

        qdim1 = constant_coefficient(quantum_dimension(Q,W,Wvars,Vvars), Vvars...)
        qdim2 = constant_coefficient(quantum_dimension(Q,V,Vvars,Wvars), Wvars...)

        if iszero(qdim1) || iszero(qdim2)
            continue
        end

        # to dense monomials
        converted = to_dense_monomials([qdim1; qdim2; C])
        qdim1 = converted[1]
        qdim2 = converted[2]
        C = converted[3:end]

        info("Found a potentially interesting matrix: doing the full computation.")

        CC = groebner_basis(C, Val{false}, max_degree=max(deg(qdim1),deg(qdim2)))

        qdim1_red,_ = red(qdim1, CC)
        qdim2_red,_ = red(qdim2, CC)

        if !iszero(qdim1_red) && !iszero(qdim2_red)
            return true
        end
    end
    return false
end

function is_orbifold_equivalent(R, W, Wvars, V, Vvars, N=Inf)

    for S = Base.Iterators.countfrom(1)
        for rank = 1:S
            for a = 0:(S-rank)
                for b = 0:(S-rank-a)
                    info("Trying rank=$rank, graded module isomorphism class = ($a, $b)")
                    if equivalence_exists(R, W, Wvars, V, Vvars, rank, a, b)
                        return true
                    end
                    if (N-=1)<=0
                        return null
                    end
                end
            end
        end
    end
end

export supertrace, multivariate_residue, quantum_dimension, equivalence_exists, is_orbifold_equivalent


end
