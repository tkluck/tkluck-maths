module Groebner

using PolynomialRings.Polynomials: Polynomial, AbstractModuleElement
using PolynomialRings.Polynomials: _div_with_remainder, leading_term, iszero, _maybe_lcm_multipliers, _ModuleElement, modulebasering

function reduce{M <: AbstractModuleElement}(f::M, G::AbstractArray{M})
    factors = transpose(spzeros(modulebasering(M), length(G)))
    frst = true
    more_loops = false
    f_red = f
    while frst || more_loops
        frst = false
        more_loops = false
        for (i, g) in enumerate(G)
            q, f_red = _div_with_remainder(f_red, g)
            if !isnull(q)
                factors[1, i] += get(q)
                more_loops = true
            end
            if iszero(f_red)
                return f_red, factors
            end
        end
    end

    return f_red, factors
end

function groebner_basis{M <: AbstractModuleElement}(polynomials::AbstractVector{M})

    P = modulebasering(M)
    nonzero_indices = find(p->!iszero(p), polynomials)
    result = polynomials[nonzero_indices]
    transformation =[ P[ i==nonzero_indices[j] ? 1 : 0 for i in eachindex(polynomials)] for j in eachindex(result)]
    pairs_to_consider = [
        (i,j) for i in eachindex(result) for j in eachindex(result) if i < j
    ]

    while length(pairs_to_consider) > 0
        (i,j) = pop!(pairs_to_consider)
        a = result[i]
        b = result[j]

        lt_a = leading_term(a)
        lt_b = leading_term(b)

        maybe_multipliers = _maybe_lcm_multipliers(lt_a, lt_b)
        if !isnull(maybe_multipliers)
            m_a, m_b = get(maybe_multipliers)
            S = m_a * a - m_b * b

            # potential speedup: wikipedia says that in all but the 'last steps'
            # (whichever those may be), we can get away with a version of reduce
            # that only does lead division
            (S_red, factors) = reduce(S, result)

            factors[1, i] -= m_a
            factors[1, j] += m_b

            if !iszero(S_red)
                new_j = length(result)+1
                append!(pairs_to_consider, [(new_i, new_j) for new_i in eachindex(result)])
                push!(result, S_red)

                nonzero_factors = find(factors)
                tr = [ sum(-factors[x] * transformation[x][y] for x in nonzero_factors) for y in eachindex(polynomials) ]
                push!(transformation, tr)
            end
        end
    end

    flat_tr = sparse([ transformation[x][y] for x=eachindex(result), y=eachindex(polynomials) ])

    return result, flat_tr

end

function syzygies{M <: AbstractModuleElement}(polynomials::AbstractVector{M})
    pairs_to_consider = [
        (i,j) for i in eachindex(polynomials) for j in eachindex(polynomials) if i < j
    ]

    result = Vector{Vector{modulebasering(M)}}()

    for (i,j) in pairs_to_consider
        a = polynomials[i]
        b = polynomials[j]
        lt_a = leading_term(a)
        lt_b = leading_term(b)

        maybe_multipliers = _maybe_lcm_multipliers(lt_a, lt_b)
        if !isnull(maybe_multipliers)
            m_a, m_b = get(maybe_multipliers)
            S = m_a * a - m_b * b

            (S_red, syzygy) = reduce(S, polynomials)
            if !iszero(S_red)
                throw(ArgumentError("syzygies(...) expects a Groebner basis, so S_red = $( S_red ) should be zero"))
            end
            syz_vector= vec(syzygy)
            syz_vector[i] -= m_a
            syz_vector[j] += m_b

            (syz_red, _) = reduce(syz_vector, result)
            if !iszero(syz_red)
                push!(result, syz_red)
            end
        end
    end

    flat_result = [ result[x][y] for x=eachindex(result), y=eachindex(polynomials) ]

    return flat_result
end


end
