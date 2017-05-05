module Polynomials

export groebner_basis


import Base: +,==,*,//,-,convert,promote_rule,show,cmp,isless

type Exponent{NumVars}
    e::NTuple{NumVars, Int}
end

function =={E <: Exponent}(a::E, b::E)
    return a.e == b.e
end

function +{NumVars}(a::Exponent{NumVars}, b::Exponent{NumVars})
    return Exponent{NumVars}(
        ntuple(Val{NumVars}) do i
            return a.e[i] + b.e[i]
        end
    )
end

function cmp{E <: Exponent}(a::E, b::E)
    # degrevlex
    if(sum(a.e) == sum(b.e))
        return -cmp(a.e, b.e)
    else
        return cmp(sum(a.e), sum(b.e))
    end
end

isless{E <: Exponent}(a::E, b::E) = cmp(a, b)<0

typealias _Monomial{R<:Number, NumVars} Tuple{Exponent{NumVars}, R}

_Monomial{R<:Number,NumVars}(a::Exponent{NumVars}, b::R)= _Monomial((a,b))

coefficient{R,NumVars}(a::_Monomial{R, NumVars})::R = a[2]

function //{R,NumVars,S}(a::_Monomial{R, NumVars}, b::S)
    exponent, coeff = a
    T = promote_type(R,S)
    return _Monomial(exponent, coeff // b)
end

immutable _Polynomial{R<:Number, NumVars} <: Number
    coeffs::Vector{ _Monomial{R, NumVars} }
end


convert{R<:Number, NumVars}(::Type{_Polynomial{R, NumVars}}, c0::R) = (
    c0 != 0
       ? _Polynomial{R, NumVars}([_Monomial(Exponent(ntuple(i->0, Val{NumVars})), c0)])
       : _Polynomial{R, NumVars}([])
)
promote_rule{R<:Number,NumVars,S<:Number}(::Type{_Polynomial{R, NumVars}}, ::Type{S}) =
    _Polynomial{promote_type(R,S), NumVars}

convert{R<:Number, NumVars}(::Type{_Polynomial{R, NumVars}}, c::_Monomial{R,NumVars}) =
    _Polynomial{R, NumVars}([c])
promote_rule{R<:Number,NumVars}(::Type{_Polynomial{R, NumVars}}, ::Type{_Monomial{R, NumVars}}) =
    _Polynomial{R, NumVars}

+{R,NumVars}(a::_Monomial{R,NumVars},b::_Polynomial{R,NumVars})=+(promote(a,b)...)
*{R,NumVars}(a::_Monomial{R,NumVars},b::_Polynomial{R,NumVars})=*(promote(a,b)...)
+{R,NumVars}(a::_Polynomial{R,NumVars},b::_Monomial{R,NumVars})=+(promote(a,b)...)
*{R,NumVars}(a::_Polynomial{R,NumVars},b::_Monomial{R,NumVars})=*(promote(a,b)...)

function +{R, T, NumVars}(a::_Polynomial{R, NumVars}, b::_Polynomial{T, NumVars})
    S = promote_type(R, T)
    res = Vector{ _Monomial{R, NumVars} }()

    state_a = start(a.coeffs)
    state_b = start(b.coeffs)
    while !done(a.coeffs, state_a) && !done(b.coeffs, state_b)
        ((exponent_a, coefficient_a), next_state_a) = next(a.coeffs, state_a)
        ((exponent_b, coefficient_b), next_state_b) = next(b.coeffs, state_b)

        if exponent_a < exponent_b
            append!(res, [_Monomial(exponent_a, coefficient_a)])
            state_a = next_state_a
        elseif exponent_b < exponent_a
            append!(res, [_Monomial(exponent_b, coefficient_b)])
            state_b = next_state_b
        else
            coeff = coefficient_a + coefficient_b
            if coeff != 0
                append!(res, [_Monomial(exponent_a, coeff)])
            end
            state_b = next_state_b
            state_a = next_state_a
        end
    end

    append!(res, collect(rest(a.coeffs, state_a)))
    append!(res, collect(rest(b.coeffs, state_b)))

    return _Polynomial{S, NumVars}(res)
end

function  *{R,T,NumVars}(a::_Polynomial{R,NumVars}, b::_Polynomial{T,NumVars})
    S = promote_type(R,T)
    res = Vector{ _Monomial{R, NumVars} }()

    # the following seems to be implemented through a very naive version
    # of append! that does a reallocation at every step. So implement
    # it manually below
    #summands = [
    #    _Monomial(exp_a + exp_b, coeff_a * coeff_b)
    #    for (exp_a, coeff_a) in a.coeffs for (exp_b, coeff_b) in b.coeffs
    #]
    summands = Vector{_Monomial{S,NumVars}}(length(a.coeffs) * length(b.coeffs))
    ix = 1
    for (exp_a, coeff_a) in a.coeffs
        for (exp_b, coeff_b) in b.coeffs
            summands[ix] = _Monomial(exp_a + exp_b, coeff_a * coeff_b)
            ix += 1
        end
    end
    assert( ix == length(summands)+1)
    sort!(summands)

    last_exp = Union{}
    for (exponent, coef) in summands
        if exponent == last_exp
            _, cur_coef = res[end]
            res[end] = exponent => cur_coef + coef
        else
            append!(res, [ _Monomial(exponent, coef) ])
            last_exp = exponent
        end
    end

    return _Polynomial{S, NumVars}(res)
end

function -{P <: _Polynomial}(f::P)
    return P([_Monomial(exponent, -coeff) for (exponent, coeff) in f.coeffs])
end

function -{P <: _Polynomial}(a::P, b::P)
    return a + -b
end

function =={P <: _Polynomial}(a::P, b::P)
    return a.coeffs == b.coeffs
end

function leading_term{P <: _Polynomial}(p::P)
    if length(p.coeffs) > 0
        return p.coeffs[end]
    else
        throw(ArgumentError("The zero polynomial $( p ) does not have a leading term"))
    end
end

typealias _ModuleElement{P <: _Polynomial} Union{P, Vector{P}}
typealias _ModuleElementVector{P <: _Polynomial} Union{AbstractVector{P}, AbstractVector{Vector{P}}}

function _lcm_multipliers{NumVars}(a::Exponent{NumVars}, b::Exponent{NumVars})
     _lcm = Exponent{NumVars}(
        ntuple(Val{NumVars}) do i
            return max(a.e[i], b.e[i])
        end
    )

    multiplier_a = Exponent{NumVars}(
        ntuple(Val{NumVars}) do i
            return _lcm.e[i] - a.e[i]
        end
    )
    multiplier_b = Exponent{NumVars}(
        ntuple(Val{NumVars}) do i
            return _lcm.e[i] - b.e[i]
        end
    )

    return multiplier_a, multiplier_b
end

function _is_constant{M <: _Monomial}(a::M)
    exponent, coeff = a
    return sum(exponent.e) == 0
end

function _lcm_multipliers{M <: _Monomial}(a::M, b::M)
    exp_a, coeff_a = a
    exp_b, coeff_b = b

    m_exp_a, m_exp_b = _lcm_multipliers(exp_a, exp_b)

    return _Monomial(m_exp_a, coeff_b), _Monomial(m_exp_b, coeff_a)
end

function _monomial_div{M<: _Monomial}(a::M, b::M)::Nullable{M}
    mul_a, mul_b = _lcm_multipliers(a, b)

    if _is_constant(mul_a)
        return mul_b // coefficient(mul_a)
    else
        return nothing
    end
end

function _lead_div_with_remainder{R,NumVars}(f::_Polynomial{R,NumVars}, g::_Polynomial{R,NumVars})::Tuple{Nullable{_Polynomial{R,NumVars}}, _Polynomial{R,NumVars}}
    maybe_factor = _monomial_div(leading_term(f), leading_term(g))

    if isnull(maybe_factor)
        return nothing, f
    else
        factor = get(maybe_factor)
        return factor, f - g*factor
    end
end

function _reduce{P <: _Polynomial}(f::_ModuleElement{P}, G::_ModuleElementVector{P})
    factors = zero(G)

    frst = true
    more_loops = false
    while frst || more_loops
        frst = false
        more_loops = false
        for (i, g) in enumerate(G)
            q, f = _lead_div_with_remainder(f, g)
            if !isnull(q)
                factors[i] += get(q)
                more_loops = true
            end
            if f == 0
                return f, factors
            end
        end
    end

    return f, factors
end



function groebner_basis{P <: _Polynomial}(polynomials::_ModuleElementVector{P})

    result = copy(polynomials)
    transformation =[ P[ i==j ? 1 : 0 for i in eachindex(polynomials)] for j in eachindex(polynomials)]


    pairs_to_consider = [
        (i,j) for i in eachindex(result) for j in eachindex(result) if i < j
    ]

    while length(pairs_to_consider) > 0
        (i,j) = pop!(pairs_to_consider)
        a = result[i]
        b = result[j]
        lt_a = leading_term(a)
        lt_b = leading_term(b)

        m_a, m_b = _lcm_multipliers(lt_a, lt_b)
        S = m_a * a - m_b * b

        if S != 0
            (S_red, factors) = _reduce(S, result)

            factors[i] += m_a
            factors[j] += m_b

            if S_red != 0
                new_j = length(result)+1
                append!(pairs_to_consider, [(new_i, new_j) for new_i in eachindex(result)])
                append!(result, [S_red])

                tr = [ sum(-f * transformation[x][y] for (x,f) in enumerate(factors)) for y in eachindex(polynomials) ]
                append!(transformation, [tr])
            end
        end
    end

    flat_tr = [ transformation[x][y] for x=eachindex(result), y=eachindex(polynomials) ]

    return result, flat_tr
end

function show{R}(io::IO, p::_Polynomial{R, 1})
    frst = true
    for (e, c) in p.coeffs
        (i,) = e.e
        if !frst
            print(io, " + ")
        else
            frst = false
        end
        print(io, c)
        if i == 1
            print(io, " x")
        elseif i > 1
            print(io, " x^$(i)")
        end
    end
end

function show{R}(io::IO, p::_Polynomial{R, 2})
    frst = true
    if length(p.coeffs) == 0
        show(io, zero(R))
    end
    for (e, c) in p.coeffs
        (i,j) = e.e
        if !frst
            print(io, " + ")
        else
            frst = false
        end
        if c != 1 || i == 0
            print(io, c)
        end
        if i == 1
            print(io, " x")
        elseif i > 1
            print(io, " x^$(i)")
        end
        if j == 1
            print(io, " y")
        elseif j > 1
            print(io, " y^$(j)")
        end
    end
end

x = _Polynomial{Int, 2}([_Monomial(Exponent{2}((1,0)), 1)])
y = _Polynomial{Int, 2}([_Monomial(Exponent{2}((0,1)), 1)])

function random_polynomial()
    res = sum([ rand(0:100) * x^i for i = 0:10 ])
end

function random_matrix()
    return Matrix([ random_polynomial() for i = 1:10, j = 1:10 ])
end

end
