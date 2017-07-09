module QuantumDimension

using PolynomialRings

function supertrace(Q::Matrix)
    n,m = size(Q)
    n == m && n%2 ==0 || throw(ArgumentError("Cannot compute supertrace of $n x $m matrix"))

    k = div(n,2)

    return sum(Q[i,i] for i=1:k) - sum(Q[i,i] for i=(k+1):2k)
end

function multivariate_residue(g, f, vars...)
    R = typeof(g)
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

export supertrace, multivariate_residue, quantum_dimension


end
