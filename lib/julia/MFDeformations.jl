module MFDeformations

using PolynomialRings: polynomial_ring, expansion
using PolynomialRings.Polynomials: Polynomial, basering
using PolynomialRings.Modules: HomspaceMorphism, kernel, span, lift_and_obstruction
using PolynomialRings.Groebner: minimal_groebner_basis, reduce


function diff{P <: Polynomial}(Q::Matrix{P})
    two_n,two_m = size(Q)
    if two_n != two_m || two_n % 2 != 0
        throw(ArgumentError("The matrix Q needs to be a 2n x 2n block matrix"))
    end
    n = two_n >> 1

    delta(Q, i,j) = ( d=zeros(Q); d[i,j]=one(eltype(Q)); d )
    grading(i,j) = (
        (i in 1:n && j in 1:n) || (i in (n+1):(2*n) && j in (n+1):(2*n))
        ? 0: 1
    )
    dQ      = [ Q*delta(Q,i,j) - (-1)^grading(i,j)*delta(Q,i,j)*Q for i in indices(Q,1), j in indices(Q,2) ]
    dQ_even = [ grading(i,j) == 1 ? (Q*delta(Q,i,j) + delta(Q,i,j)*Q) : zero(Q) for i in indices(Q,1), j in indices(Q,2) ]
    dQ_odd  = [ grading(i,j) == 0 ? (Q*delta(Q,i,j) - delta(Q,i,j)*Q) : zero(Q) for i in indices(Q,1), j in indices(Q,2) ]

    dQ      = HomspaceMorphism([ dQ[i,j][k,l] for i      = 1:two_n,j = 1:two_n,k = 1:two_n,l = 1:two_n])
    dQ_even = HomspaceMorphism([ dQ_even[i,j][k,l] for i = 1:two_n,j = 1:two_n,k = 1:two_n,l = 1:two_n])
    dQ_odd  = HomspaceMorphism([ dQ_odd[i,j][k,l] for i  = 1:two_n,j = 1:two_n,k = 1:two_n,l = 1:two_n])

    return dQ, dQ_even, dQ_odd
end

function colspan{R <:Number}(M::Matrix{R})
    M = copy(M)
    for j in indices(M, 2)
        for i in indices(M, 1)
            if M[i,j] != 0
                pivot = findfirst(M[i,(j+1):end])
                if pivot == 0
                    break
                end
                M[:,j] = M[i,pivot] * M[:,j] - M[i,j] * M[:,pivot]
            end
        end
    end

    M = M[:, [ j for j in indices(M, 2) if any(m != 0 for m in M[:,j])]]

    return M
end

import PolynomialRings.Polynomials: Term, Monomial
function H1{P <: Polynomial}(Q::Matrix{P})
    dQ, dQ_even, dQ_odd = diff(Q)

    groeb,transformation = minimal_groebner_basis(span(dQ_odd))

    H1 = filter(kernel(dQ_even)) do k
        (k_red, factors) = reduce(k, groeb)
        any(c != 0 for c in k_red)
    end

    finite_basis_set = Set(
        (i, term[1].e)
        for h in H1
        for (i, c) in enumerate(h)
        for term in c.terms
    )
    index = Dict( m => i for (i,m) in enumerate(finite_basis_set) )
    reverse_index = Dict( i => m for (i,m) in enumerate(finite_basis_set) )

    M = zeros(basering(P), (length(finite_basis_set), length(H1)))
    for (k,h) in enumerate(H1)
        for (j, c) in enumerate(h)
            for term in c.terms
                M[index[j, term[1].e], k] = term[2]
            end
        end
    end

    N = colspan(M)
    return map(1:size(N, 2)) do j
        h = zeros(P, size(Q))
        for (i,v) in enumerate(N[:, j])
            j,e  = reverse_index[i]
            h[j] += P([ Term(Monomial(e), v) ])
        end
        h
    end

end


function deformation(Q, var_symbols...; max_order=20)
    dQ, dQ_even, dQ_odd = diff(Q)
    H = H1(Q)

    if length(var_symbols) == 0
        var_symbols = [ gensym() for _ in H ]
    end
    _, vars = polynomial_ring(Int, var_symbols...)
    if length(vars) != length(H)
        throw(ArgumentError("Need $( length(H) ) variables for this deformation"))
    end

    Q1 = sum(w * h for (w,h) in zip(vars, H))
    Qdef = Q1
    sumobs = zero(Q)

    Qs = [Q1]
    obs = []

    to_vanish = (Q + Qdef)^2 - Q^2

    for ord = 2:max_order
        tic()
        # it already vanishes up to (not including) ord; if it vanishes
        # completely, we are done
        if to_vanish == zero(to_vanish)
            break
        end

        N = length(Q)
        MC = sum(Q_i * Q_j for (Q_i, Q_j) in zip(Qs, reverse(Qs)))
        Q_next, obs_next = mapreduce(
            (a,b)->(a[1]+b[1],a[2]+b[2]),
            (zero(Q), zero(Q)),
            expansion(MC, var_symbols...)) do x
            w, MC_w = x
            Q_w, obs_w = lift_and_obstruction(dQ_even, -MC_w)
            (w*Q_w, w*obs_w)
        end
        to_vanish += Q_next^2 + Q_next*(Q + Qdef) + (Q + Qdef)*Q_next + obs_next
        push!(Qs, Q_next)
        push!(obs, obs_next)
        Qdef += Q_next
        sumobs += obs_next

        println(STDERR, "Step $(ord): $( toq() )")
    end

    return Q + Qdef

end

end
