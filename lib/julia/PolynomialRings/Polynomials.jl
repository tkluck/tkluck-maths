module Polynomials

import Base: +,==,!=,*,//,-,convert,promote_rule,show,cmp,isless,zero,eltype

type Monomial{NumVars}
    e::NTuple{NumVars, Int}
end

function =={M <: Monomial}(a::M, b::M)
    return a.e == b.e
end

function +{NumVars}(a::Monomial{NumVars}, b::Monomial{NumVars})
    return Monomial{NumVars}(
        ntuple(Val{NumVars}) do i
            return a.e[i] + b.e[i]
        end
    )
end

function cmp{M <: Monomial}(a::M, b::M)
    # degrevlex
    if(sum(a.e) == sum(b.e))
        return -cmp(a.e, b.e)
    else
        return cmp(sum(a.e), sum(b.e))
    end
end

isless{M <: Monomial}(a::M, b::M) = cmp(a, b)<0

typealias Term{R<:Number, NumVars} Tuple{Monomial{NumVars},R}
Term{R<:Number,NumVars}(m::Monomial{NumVars},r::R) = Term((m,r))
coefficient{R,NumVars}(a::Term{R, NumVars})::R = a[2]

function //{R,NumVars,S}(a::Term{R, NumVars}, b::S)
    exponent, coeff = a
    T = promote_type(R,S)
    return Term(exponent, coeff // b)
end

immutable Polynomial{R <: Number, NumVars, T <: Tuple} <: Number
    terms::Vector{ Term{R, NumVars} }
end

basering{R <: Number, NumVars, T <: Tuple}(::Type{Polynomial{R, NumVars, T}}) = R
basering{R <: Number, NumVars, T <: Tuple}(::Polynomial{R, NumVars, T}) = R

num_variables{R<:Number, NumVars, T <: Tuple}(::Type{Polynomial{R,NumVars,T}}) = NumVars
variables{R<:Number, NumVars, T<:Tuple}(::Type{Polynomial{R,NumVars,T}}) = ntuple(Val{NumVars}) do i
    exponent = ntuple(Val{NumVars}) do j
        i == j ? 1 : 0
    end
    Polynomial{R,NumVars,T}([Term(Monomial(exponent), one(R))])
end

num_variables{P <: Polynomial}(p::P) = num_variables(typeof(p))
variables{P <: Polynomial}(p::P) = variables(typeof(p))

convert{R<:Number, NumVars, T<:Tuple, S<:Number}(::Type{Polynomial{R, NumVars, T}}, c0::S) = (
    c0 != 0
       ? Polynomial{promote_type(R,S), NumVars, T}([(Monomial(ntuple(i->0, Val{NumVars})), promote_type(R,S)(c0))])
       : Polynomial{promote_type(R,S), NumVars, T}([])
)

promote_rule{R<:Number,NumVars,T<:Tuple,S<:Number}(::Type{Polynomial{R, NumVars, T}}, ::Type{S}) =
    Polynomial{promote_type(R,S), NumVars, T}

convert{R<:Number, NumVars,T<:Tuple}(::Type{Polynomial{R, NumVars,T}}, c::Term{R,NumVars}) = (
    coefficient(c) != 0
       ?  Polynomial{R, NumVars,T}([c])
       : Polynomial{R, NumVars,T}([])
)

promote_rule{R<:Number,NumVars,T}(::Type{Polynomial{R, NumVars,T}}, ::Type{Term{R, NumVars}}) =
    Polynomial{R, NumVars,T}

promote_rule{S<:Polynomial}(::Type{S}, ::Type{S}) = S
convert{S<:Polynomial}(::Type{S}, p::S) = p

+{R,NumVars,T}(a::Term{R,NumVars},b::Polynomial{R,NumVars,T})=+(promote(a,b)...)
*{R,NumVars,T}(a::Term{R,NumVars},b::Polynomial{R,NumVars,T})=*(promote(a,b)...)
-{R,NumVars,T}(a::Term{R,NumVars},b::Polynomial{R,NumVars,T})=-(promote(a,b)...)
+{R,NumVars,T}(a::Polynomial{R,NumVars,T},b::Term{R,NumVars})=+(promote(a,b)...)
*{R,NumVars,T}(a::Polynomial{R,NumVars,T},b::Term{R,NumVars})=*(promote(a,b)...)
-{R,NumVars,T}(a::Polynomial{R,NumVars,T},b::Term{R,NumVars})=-(promote(a,b)...)

iszero{P <: Polynomial}(p::P)= length(p.terms) == 0
iszero{P <: Polynomial}(p::Vector{P}) = all(iszero(p_i) for p_i in p)
iszero{P <: Polynomial}(p::Matrix{P}) = all(iszero(p_i) for p_i in p)

function +{R1, R2, NumVars, T<:Tuple}(a::Polynomial{R1, NumVars, T}, b::Polynomial{R2, NumVars, T})
    S = promote_type(R1, R2)
    res = Vector{ Term{S, NumVars} }()

    state_a = start(a.terms)
    state_b = start(b.terms)
    while !done(a.terms, state_a) && !done(b.terms, state_b)
        ((exponent_a, coefficient_a), next_state_a) = next(a.terms, state_a)
        ((exponent_b, coefficient_b), next_state_b) = next(b.terms, state_b)

        if exponent_a < exponent_b
            push!(res, Term(exponent_a, coefficient_a))
            state_a = next_state_a
        elseif exponent_b < exponent_a
            push!(res, Term(exponent_b, coefficient_b))
            state_b = next_state_b
        else
            coeff = coefficient_a + coefficient_b
            if coeff != 0
                push!(res, Term(exponent_a, coeff))
            end
            state_b = next_state_b
            state_a = next_state_a
        end
    end

    append!(res, collect(rest(a.terms, state_a)))
    append!(res, collect(rest(b.terms, state_b)))

    return Polynomial{S, NumVars,T}(res)
end

function  *{R1,R2,NumVars,T<:Tuple}(a::Polynomial{R1,NumVars,T}, b::Polynomial{R2,NumVars,T})
    S = promote_type(R1,R2)
    res = Vector{ Term{S, NumVars} }()

    # the following seems to be implemented through a very naive version
    # of push! that does a reallocation at every step. So implement
    # it manually below
    #summands = [
    #    Term(exp_a + exp_b, coeff_a * coeff_b)
    #    for (exp_a, coeff_a) in a.terms for (exp_b, coeff_b) in b.terms
    #]
    summands = Vector{Term{S,NumVars}}(length(a.terms) * length(b.terms))
    ix = 1
    for (exp_a, coeff_a) in a.terms
        for (exp_b, coeff_b) in b.terms
            summands[ix] = Term(exp_a + exp_b, coeff_a * coeff_b)
            ix += 1
        end
    end
    assert( ix == length(summands)+1)
    sort!(summands)

    last_exp = Union{}
    for (exponent, coef) in summands
        if exponent == last_exp
            _, cur_coef = res[end]
            res[end] = Term(exponent, cur_coef + coef)
        else
            push!(res,  Term(exponent, coef))
            last_exp = exponent
        end
    end
    filter!(m -> coefficient(m) != 0, res)
    return Polynomial{S, NumVars,T}(res)
end

-{P <: Polynomial}(f::P) = P([Term(exponent, -coeff) for (exponent, coeff) in f.terms])
-{P <: Polynomial}(a::P, b::P) = a + -b
=={P <: Polynomial}(a::P, b::P) = a.terms == b.terms
!={P <: Polynomial}(a::P, b::P) = !(a == b)

function leading_term{P <: Polynomial}(p::P)
    if length(p.terms) > 0
        return p.terms[end]
    else
        throw(ArgumentError("The zero polynomial $( p ) does not have a leading term"))
    end
end

typealias _ModuleElement{P <: Polynomial} Vector{P}
typealias _HomModuleElement{P <: Polynomial} Matrix{P}
typealias _AbstractModuleElement{P <: Polynomial} Union{P, _ModuleElement{P}, _HomModuleElement{P}}
typealias _AbstractModuleElementVector{P <: Polynomial} Union{AbstractVector{P}, AbstractVector{_ModuleElement{P}}, AbstractVector{_HomModuleElement{P}}}

zero{P <: Polynomial}(a::AbstractVector{Vector{P}}) = [[0 for _ in a_i] for a_i in a]

type _ModuleTerm{T <: Term}
    t::T
    pos::Int
end

*{P<:Polynomial}(a::P, x::_ModuleElement{P})= P[ a*x_i for x_i in x ]
*{P<:Polynomial}(x::_ModuleElement{P}, ::P)= P[ a*x_i for x_i in x ]
*{R<:Number, NumVars,T<:Tuple}(x::_ModuleElement{Polynomial{R, NumVars,T}}, a::Term{R, NumVars})= convert(Polynomial{R,NumVars,T}, a) * x
*{R<:Number, NumVars,T<:Tuple}(a::Term{R, NumVars}, x::_ModuleElement{Polynomial{R, NumVars,T}})= convert(Polynomial{R,NumVars,T}, a) * x

*{P<:Polynomial}(a::P, x::_HomModuleElement{P})= P[ a*x_i for x_i in x ]
*{R<:Number, NumVars,T<:Tuple}(x::_HomModuleElement{Polynomial{R, NumVars,T}}, a::Term{R, NumVars})= convert(Polynomial{R,NumVars,T}, a) * x
*{R<:Number, NumVars,T<:Tuple}(a::Term{R, NumVars}, x::_HomModuleElement{Polynomial{R, NumVars,T}})= convert(Polynomial{R,NumVars,T}, a) * x

function leading_term{P<:Polynomial}(a::_ModuleElement{P})
    for (i, f_i) in enumerate(a)
        if !iszero(f_i)
            return _ModuleTerm(leading_term(f_i), i)
        end
    end
    throw(ArgumentError("The zero element $( a ) does not have a leading term"))
end

function leading_term{P<:Polynomial}(a::_HomModuleElement{P})
    for (i, f_i) in enumerate(a)
        if !iszero(f_i)
            return _ModuleTerm(leading_term(f_i), i)
        end
    end
    throw(ArgumentError("The zero element $( a ) does not have a leading term"))
end

function _lcm_multipliers{NumVars}(a::Monomial{NumVars}, b::Monomial{NumVars})
     _lcm = Monomial{NumVars}(
        ntuple(Val{NumVars}) do i
            return max(a.e[i], b.e[i])
        end
    )

    multiplier_a = Monomial{NumVars}(
        ntuple(Val{NumVars}) do i
            return _lcm.e[i] - a.e[i]
        end
    )
    multiplier_b = Monomial{NumVars}(
        ntuple(Val{NumVars}) do i
            return _lcm.e[i] - b.e[i]
        end
    )

    return multiplier_a, multiplier_b
end

function _is_constant{M <: Term}(a::M)
    exponent, coeff = a
    return sum(exponent.e) == 0
end

function _lcm_multipliers{M <: Term}(a::M, b::M)
    exp_a, coeff_a = a
    exp_b, coeff_b = b

    m_exp_a, m_exp_b = _lcm_multipliers(exp_a, exp_b)

    return Term(m_exp_a, coeff_b), Term(m_exp_b, coeff_a)
end

function _monomial_div{M<: Term}(a::M, b::M)::Nullable{M}
    mul_a, mul_b = _lcm_multipliers(a, b)

    if _is_constant(mul_a)
        return mul_b // coefficient(mul_a)
    else
        return nothing
    end
end

function _monomial_div{M<: Term}(a::_ModuleTerm{M}, b::_ModuleTerm{M})::Nullable{M}
    if a.pos != b.pos
        return nothing
    else
        return _monomial_div(a.t, b.t)
    end
end

_maybe_lcm_multipliers{M <: Term}(a::M, b::M)::Nullable{Tuple{M,M}} = _lcm_multipliers(a,b)
function _maybe_lcm_multipliers{M <: Term}(a::_ModuleTerm{M}, b::_ModuleTerm{M})::Nullable{Tuple{M,M}}
    if a.pos != b.pos
        return nothing
    else
        return _lcm_multipliers(a.t, b.t)
    end
end

function _lead_div_with_remainder{R,NumVars}(f::Polynomial{R,NumVars}, g::Polynomial{R,NumVars})::Tuple{Nullable{Polynomial{R,NumVars}}, Polynomial{R,NumVars}}
    maybe_factor = _monomial_div(leading_term(f), leading_term(g))

    if isnull(maybe_factor)
        return nothing, f
    else
        factor = get(maybe_factor)
        return factor, f - g*factor
    end
end

_monomials{R, NumVars, T<:Tuple}(f::Polynomial{R, NumVars,T}) = reverse(f.terms)

function _monomials{R<:Number, NumVars,T<:Tuple}(f::_ModuleElement{Polynomial{R, NumVars,T}})
    return [
        _ModuleTerm{Term{R, NumVars}}(m, i)
        for (i, f_i) in enumerate(f)
        for m in _monomials(f_i)
    ]
end

function _monomials{R<:Number, NumVars, T<:Tuple}(f::_HomModuleElement{Polynomial{R, NumVars,T}})
    return [
        _ModuleTerm{Term{R, NumVars}}(m, i)
        for (i, f_i) in enumerate(f)
        for m in _monomials(f_i)
    ]
end

function _div_with_remainder{P <: Polynomial}(f::_AbstractModuleElement{P}, g::_AbstractModuleElement{P})::Tuple{Nullable{P}, _AbstractModuleElement{P}}
    if iszero(f)
        return zero(P), f
    elseif iszero(g)
        return nothing, f
    else
        for monomial in _monomials(f)
            maybe_factor = _monomial_div(monomial, leading_term(g))
            if !isnull(maybe_factor)
                factor = get(maybe_factor)
                return factor, f - g*factor
            end
        end
        return nothing, f
    end
end


import Iterators: product

function monomials_not_in_ideal{R <: Number, NumVars}(monomials::Vector{Polynomial{R, NumVars}})
    vars = variables(Polynomial{R, NumVars})

    degree_bounds = [0 for _ in vars]
    for monom in monomials
        exponent_tuple = monom.terms[1][1].e
        variables_appearing = [(i,e) for (i,e) in enumerate(exponent_tuple) if e != 0]
        if length(variables_appearing) == 1
            (variable_index,variable_exp), = variables_appearing
            if degree_bounds[ variable_index ] == 0
                degree_bounds[ variable_index ] = variable_exp
            elseif variable_exp < degree_bounds[ variable_index ]
                degree_bounds[ variable_index ] = variable_exp
            end
        end
    end
    if any(b == 0 for b in degree_bounds)
        throw(ArgumentError("monomials_not_in_ideal: result is not a finite set"))
    end

    result = eltype(monomials)[]
    for p in product([0:(b-1) for b in degree_bounds]...)
        m = prod(v^p[i] for (i,v) in enumerate(vars))
        if !any(!isnull(_monomial_div(m.terms[1], d.terms[1])) for d in monomials)
            push!(result, m)
        end
    end
    return result

end


function random_polynomial()
    res = sum([ rand(0:100) * x^i for i = 0:10 ])
end

function random_matrix()
    return Matrix([ random_polynomial() for i = 1:10, j = 1:10 ])
end

end
