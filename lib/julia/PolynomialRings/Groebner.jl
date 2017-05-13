module Groebner

using PolynomialRings.Polynomials: Polynomial, _AbstractModuleElement, _AbstractModuleElementVector
using PolynomialRings.Polynomials: _div_with_remainder, leading_term, iszero, _maybe_lcm_multipliers, _ModuleElement

function reduce{P <: Polynomial}(f::_AbstractModuleElement{P}, G::_AbstractModuleElementVector{P})
    factors = transpose(zeros(P, length(G)))
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

function groebner_basis{P <: Polynomial}(polynomials::_AbstractModuleElementVector{P})

    result = copy(polynomials)
    transformation =[ P[ i==j ? 1 : 0 for i in eachindex(polynomials)] for j in eachindex(polynomials)]
    pairs_to_consider = [
        (i,j) for i in eachindex(result) for j in eachindex(result) if i < j
    ]

    while length(pairs_to_consider) > 0
        (i,j) = pop!(pairs_to_consider)
        a = result[i]
        b = result[j]
        if iszero(a) || iszero(b)
            continue
        end
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

                tr = [ sum(-f * transformation[x][y] for (x,f) in enumerate(factors)) for y in eachindex(polynomials) ]
                push!(transformation, tr)
            end
        end
    end

    flat_tr = [ transformation[x][y] for x=eachindex(result), y=eachindex(polynomials) ]

    return result, flat_tr
end

function minimal_groebner_basis{P <: Polynomial}(polynomials::_AbstractModuleElementVector{P})

    (basis, transformation) = groebner_basis(polynomials)

    redundant = Set{Int}()
    for i in eachindex(basis)
        #(p_red, factors) = reduce(basis[i], [b for (j,b) in enumerate(basis) if j!=i && !(j in redundant)])
        if iszero(basis[i])
            push!(redundant, i)
        end
    end

    necessary = [ i for i in eachindex(basis) if !(i in redundant) ]
    return basis[ necessary ], transformation[ necessary, : ]

end

function syzygies{P <: Polynomial}(polynomials::_AbstractModuleElementVector{P})
    pairs_to_consider = [
        (i,j) for i in eachindex(polynomials) for j in eachindex(polynomials) if i < j
    ]

    result = Vector{_ModuleElement{P}}()

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

    (result, _) = minimal_groebner_basis(result)
    flat_result = [ result[x][y] for x=eachindex(result), y=eachindex(polynomials) ]

    return flat_result
end


end
