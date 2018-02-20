module MFDeformations

using PolynomialRings
using HomspaceMorphisms: HomspaceMorphism, kernel, span, lift_and_obstruction

using QuasiHomogeneous

function diff(Q::Matrix{P}) where P <: Polynomial
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

using PolynomialRings: monomialtype, variablesymbols

function finite_subspace_conversion(arrays::AbstractArray{<:AbstractArray{<:Polynomial}}, vars...)
    P = eltype(eltype(arrays))
    MonomialType, CoeffType = PolynomialRings.Expansions._expansion_types(P, vars...)
    finite_basis_set = Set(
                           (i, w)
                           for array in arrays
                           for (i, c) in enumerate(array)
                           for (w, p) in expansion(c, vars...)
                          )
    index = Dict( m => i for (i,m) in enumerate(finite_basis_set) )
    reverse_index = Dict( i => m for (i,m) in enumerate(finite_basis_set) )

    to_vector(array) = begin
        res = zeros(CoeffType, length(finite_basis_set))
        for (i, c) in enumerate(array)
            for (w, p) in expansion(c, vars...)
                res[index[i, w]] = p
            end
        end
        res
    end
    unexpand(w) = prod(convert(P, v)^i for (v,i) in zip(vars, w))
    to_polynomial_array(vec) = begin
        res = zeros(P, size(arrays[1]))
        for (j, p) in enumerate(vec)
            i, w = reverse_index[j]
            res[i] += p * unexpand(w)
        end
        return res
    end

    to_vector, to_polynomial_array
end

function H1(Q, dQ_even::HomspaceMorphism{P}, dQ_odd::HomspaceMorphism{P}) where P <: Polynomial
    groeb = gröbner_basis(dQ_odd)

    H1 = map(k->rem(k, groeb), kernel(dQ_even))

    to_vector, to_polynomial_array = finite_subspace_conversion(H1, variablesymbols(P)...)

    M = hcat(map(to_vector, H1)...)
    N = PolynomialRings.Util.LinAlgUtil.colspan(M)

    return [to_polynomial_array(N[:,j]) for j=1:size(N,2)]
end

function H1(Q::Matrix{P}) where P <: Polynomial
    dQ, dQ_even, dQ_odd = diff(Q)
    return H1(Q, dQ_even, dQ_odd)
end


function graded_implicit_tangent_space(f, Q, vars::Gradings)
    info("Computing gradings")
    gr = map(q_i->quasidegree(q_i, vars), Q)

    info("Creating deformation vector")
    ch = formal_coefficients(eltype(Q), :c)
    N = generic_quasihomogeneous_map(gr, vars, ch)

    info("Applying function")
    CC = flat_coefficients(f(Q+N) - f(Q), symbols(vars)...)

    info("Getting linear coefficients")
    coeffs = [@linear_coefficients(cc, c[]) for cc in CC]
    rows = length(coeffs)
    cols = maximum(length, coeffs)
    info("Computing matrix")
    M = zeros(eltype(eltype(coeffs)), rows, cols)
    for (i,c) in enumerate(coeffs)
        for (j,cc) in enumerate(c)
            M[i,j] = cc
        end
    end

    info("Computing kernel")
    K = PolynomialRings.Util.LinAlgUtil.kernel(M)

    info("Substituting kernel back into deformation vector")
    result = map(indices(K,2)) do j
        N(c = i->K[i,j])
    end

    # workaround implementation detail: generic_quasihomogeneous map sometimes
    # skips variables
    filter!(!iszero, result)

    info("Done")

    return result
end

function first_order_deformation(Q, vars::Gradings, ε::Symbol)
    T = graded_implicit_tangent_space(Q->Q^2,Q,vars)
    Q1 = sum(prod, zip(formal_coefficients(eltype(Q),ε), T))
    return Q1
end

function deformation(Q, vars::Gradings; max_order=20)
    Q1 = first_order_deformation(Q, vars, :ε)
    Qdef = Q1
    sumobs = zero(Q)

    Qs = [Q1]
    obs = []

    to_vanish = (Q + Qdef)^2 - Q^2

    for ord = 2:max_order
        tic()
        # it already vanishes up to (not including) ord; if it vanishes
        # completely, we are done
        #if to_vanish == zero(to_vanish)
        #    break
        #end

        MC = sum(Q_i * Q_j for (Q_i, Q_j) in zip(Qs, reverse(Qs)))
        Q_next, obs_next = mapreduce(
            (a,b)->(a[1]+b[1],a[2]+b[2]),
            (zero(Q), zero(Q)),
            @expansion(MC, ε[])) do x
            exps, MC_w = x
            w = prod(v^e for (v,e) in zip(vars,exps))
            Q_w, obs_w = lift_and_obstruction(dQ_even, -MC_w)
            (w*Q_w, w*obs_w)
        end
        #to_vanish += Q_next^2 + Q_next*(Q + Qdef) + (Q + Qdef)*Q_next + obs_next
        push!(Qs, Q_next)
        push!(obs, obs_next)
        Qdef += Q_next
        sumobs += obs_next

        info("Step $(ord): $( toq() )")
    end

    return Q + Qdef

end

end
