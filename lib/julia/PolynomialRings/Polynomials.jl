module Polynomials

import Base: +,==,!=,*,//,-,convert,promote_rule,show,cmp,isless,zero,eltype

immutable Monomial{NumVars}
    e::NTuple{NumVars, Int}
    t::Int
end
Monomial{NumVars}(e::NTuple{NumVars, Int}) = Monomial(e, sum(e))

function =={M <: Monomial}(a::M, b::M)
    return a.t == b.t && a.e == b.e
end

@generated function +{NumVars}(a::Monomial{NumVars}, b::Monomial{NumVars})
    expr = :( tuple() )
    for i=1:NumVars
        push!(expr.args, :( a.e[$i] + b.e[$i] ))
    end
    quote
        Monomial{$NumVars}($expr, a.t + b.t)
    end
end

function cmp{M <: Monomial}(a::M, b::M)
    # degrevlex
    degcmp = cmp(a.t, b.t)
    if degcmp == 0
        i = length(a.e)
        @inbounds while i >= 1
            varcmp = cmp(a.e[i], b.e[i])
            if varcmp != 0
                return -varcmp
            end
            i -= 1
        end
        return 0
    else
        return degcmp
    end
end

isless{M <: Monomial}(a::M, b::M) = cmp(a, b)<0

abstract AbstractBaseRing
typealias BaseRing{N <: Number}  Union{N, AbstractBaseRing}
typealias Term{R<:BaseRing, NumVars} Tuple{Monomial{NumVars},R}
Term{R<:BaseRing,NumVars}(m::Monomial{NumVars},r::R) = Term((m,r))
coefficient{R,NumVars}(a::Term{R, NumVars})::R = a[2]

isless{T <: Term}(a::T, b::T) = cmp(a[1], b[1])<0

fieldtypes{T <: Tuple}(t::Type{T}) = Symbol[fieldtype(T, i) for i in 1:nfields(T)]
import PolynomialRings
function termmul{T1 <: Term, T2 <: Term, Vars1 <: Tuple, Vars2 <: Tuple}(a::T1, b::T2, ::Type{Vars1}, ::Type{Vars2})

    all_names = Set()
    union!(all_names, fieldtypes(Vars1))
    union!(all_names, fieldtypes(Vars2))
    names = sort(collect(Symbol(s) for s in all_names))
    Vars = Tuple{names...}
    NumVars = length(all_names)

    f = PolynomialRings._converter(Vars, Vars1, true)
    g = PolynomialRings._converter(Vars, Vars2, true)

    exp_a, coef_a = a
    exp_b, coef_b = b
    return Term(Monomial(f(exp_a.e)) + Monomial(g(exp_b.e)), coef_a*coef_b)

end

function termmul{T <: Term, Vars <: Tuple}(a::T, b::T, ::Type{Vars}, ::Type{Vars})
    exp_a, coef_a = a
    exp_b, coef_b = b
    Term(exp_a + exp_b, coef_a * coef_b)
end


function //{R,NumVars,S}(a::Term{R, NumVars}, b::S)
    exponent, coeff = a
    T = promote_type(R,S)
    return Term(exponent, coeff // b)
end

immutable Polynomial{R <: BaseRing, NumVars, T <: Tuple} <: Number
    terms::Vector{ Term{R, NumVars} }
end

basering{R <: BaseRing, NumVars, T <: Tuple}(::Type{Polynomial{R, NumVars, T}}) = R
basering{R <: BaseRing, NumVars, T <: Tuple}(::     Polynomial{R, NumVars, T} ) = R
termtype{R <: BaseRing, NumVars, T <: Tuple}(::Type{Polynomial{R, NumVars, T}}) = Term{R, NumVars}
termtype{R <: BaseRing, NumVars, T <: Tuple}(::     Polynomial{R, NumVars, T} ) = Term{R, NumVars}
monomialtype{R <: BaseRing, NumVars, T <: Tuple}(::Type{Polynomial{R, NumVars, T}}) = Monomial{NumVars}
monomialtype{R <: BaseRing, NumVars, T <: Tuple}(::     Polynomial{R, NumVars, T} ) = Monomial{NumVars}
_varsymbols{R <: BaseRing, NumVars, T <: Tuple}(::Type{Polynomial{R, NumVars, T}}) = T
_varsymbols{R <: BaseRing, NumVars, T <: Tuple}(::     Polynomial{R, NumVars, T} ) = T

num_variables{R<:BaseRing, NumVars, T <: Tuple}(::Type{Polynomial{R,NumVars,T}}) = NumVars
variables{R<:BaseRing, NumVars, T<:Tuple}(::Type{Polynomial{R,NumVars,T}}) = map(1:NumVars) do i
    exponent = ntuple(Val{NumVars}) do j
        i == j ? 1 : 0
    end
    Polynomial{R,NumVars,T}([Term(Monomial(exponent), one(R))])
end

num_variables{P <: Polynomial}(p::P) = num_variables(typeof(p))
variables{P <: Polynomial}(p::P) = variables(typeof(p))

convert{R<:BaseRing, NumVars, T<:Tuple, S<:BaseRing}(::Type{Polynomial{R, NumVars, T}}, c0::S) = (
    c0 != 0
       ? Polynomial{promote_type(R,S), NumVars, T}([(Monomial(ntuple(i->0, Val{NumVars})), promote_type(R,S)(c0))])
       : Polynomial{promote_type(R,S), NumVars, T}([])
)

promote_rule{R<:Number,NumVars,T<:Tuple,S<:Real}(::Type{Polynomial{R, NumVars, T}}, ::Type{S}) =
    Polynomial{promote_type(R,S), NumVars, T}

convert{R<:BaseRing, NumVars,T<:Tuple}(::Type{Polynomial{R, NumVars,T}}, c::Term{R,NumVars}) = (
    coefficient(c) != 0
       ? Polynomial{R, NumVars,T}([c])
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
iszero{P <: Polynomial}(p::AbstractArray{P}) = all(iszero(p_i) for p_i in p)

import Base: zero, one
zero{R, NumVars, T}(p::Polynomial{R, NumVars, T}) = Polynomial{R, NumVars, T}([])
one{R, NumVars, T}(p::Polynomial{R, NumVars, T}) = (
    Polynomial{R, NumVars, T}([(Monomial(ntuple(i->0, Val{NumVars})), one(R))])
)

function +{R1, R2, NumVars, T<:Tuple}(a::Polynomial{R1, NumVars, T}, b::Polynomial{R2, NumVars, T})
    S = promote_type(R1, R2)
    res = Vector{ Term{S, NumVars} }(length(a.terms) + length(b.terms))
    n = 0

    state_a = start(a.terms)
    state_b = start(b.terms)
    while !done(a.terms, state_a) && !done(b.terms, state_b)
        ((exponent_a, coefficient_a), next_state_a) = next(a.terms, state_a)
        ((exponent_b, coefficient_b), next_state_b) = next(b.terms, state_b)

        if exponent_a < exponent_b
            @inbounds res[n+=1] = Term(exponent_a, coefficient_a)
            state_a = next_state_a
        elseif exponent_b < exponent_a
            @inbounds res[n+=1] = Term(exponent_b, coefficient_b)
            state_b = next_state_b
        else
            coeff = coefficient_a + coefficient_b
            if coeff != 0
                @inbounds res[n+=1] = Term(exponent_a, coeff)
            end
            state_b = next_state_b
            state_a = next_state_a
        end
    end

    for t in rest(a.terms, state_a)
        @inbounds res[n+=1] = t
    end
    for t in rest(b.terms, state_b)
        @inbounds res[n+=1] = t
    end

    resize!(res, n)
    return Polynomial{S, NumVars,T}(res)
end

function -{R1, R2, NumVars, T<:Tuple}(a::Polynomial{R1, NumVars, T}, b::Polynomial{R2, NumVars, T})
    S = promote_type(R1, R2)
    res = Vector{ Term{S, NumVars} }(length(a.terms) + length(b.terms))
    n = 0

    state_a = start(a.terms)
    state_b = start(b.terms)
    while !done(a.terms, state_a) && !done(b.terms, state_b)
        ((exponent_a, coefficient_a), next_state_a) = next(a.terms, state_a)
        ((exponent_b, coefficient_b), next_state_b) = next(b.terms, state_b)

        if exponent_a < exponent_b
            @inbounds res[n+=1] = Term(exponent_a, coefficient_a)
            state_a = next_state_a
        elseif exponent_b < exponent_a
            @inbounds res[n+=1] = Term(exponent_b, -coefficient_b)
            state_b = next_state_b
        else
            coeff = coefficient_a - coefficient_b
            if coeff != 0
                @inbounds res[n+=1] = Term(exponent_a, coeff)
            end
            state_b = next_state_b
            state_a = next_state_a
        end
    end

    for t in rest(a.terms, state_a)
        @inbounds res[n+=1] = t
    end
    for t in rest(b.terms, state_b)
        (exp, c) = t
        @inbounds res[n+=1] = Term(exp, -c)
    end

    resize!(res, n)
    return Polynomial{S, NumVars,T}(res)
end

macro _enqueue_term(i,j)
    quote
        @inbounds t = termmul(a.terms[$i], b.terms[$j], _varsymbols(P1), _varsymbols(P2))
        enqueue!(minimal_corners, ($i,$j), t)
    end
end
import Util: BoundedPriorityQueue, enqueue!, dequeue!, peek
function *{P1 <: Polynomial, P2 <: Polynomial}(a::P1, b::P2)
    PP = promote_type(P1, P2)

    if iszero(a) || iszero(b)
        return zero(PP)
    end

    summands = Vector{termtype(PP)}(length(a.terms) * length(b.terms))
    k = 0

    row_indices= zeros(Int, length(a.terms))
    col_indices= zeros(Int, length(b.terms))

    # using a bounded queue not to drop items when it gets too big, but to allocate it
    # once to its maximal theoretical size and never reallocate.
    minimal_corners = BoundedPriorityQueue{Tuple{Int, Int}, termtype(PP)}(min(length(a.terms), length(b.terms)))
    @_enqueue_term(1,1)
    @inbounds while length(minimal_corners)>0
        # I don't understand the type inference breakage here, but making
        # it explicit speeds things up
        (row, col), t = peek(minimal_corners)::Pair{Tuple{Int,Int}, termtype(PP)}
        dequeue!(minimal_corners)
        summands[k+=1] = t
        row_indices[row] = col
        col_indices[col] = row
        if row < length(a.terms) && row_indices[row+1] == col - 1
            @_enqueue_term(row+1, col)
        end
        if col < length(b.terms) && col_indices[col+1] == row - 1
            @_enqueue_term(row, col+1)
        end
    end

    assert(k == length(summands))
    #assert(issorted(summands, by=t -> t[1]))

    if length(summands) > 0
        last_exp, _ = summands[1]
        n = 1
        for j in 2:length(summands)
            exponent, coef = summands[j]
            if exponent == last_exp
                _, cur_coef = summands[n]
                @inbounds summands[n] = Term(exponent, cur_coef + coef)
            else
                @inbounds summands[n+=1] = summands[j]
                last_exp = exponent
            end
        end
        resize!(summands, n)
        filter!(m -> coefficient(m) != 0, summands)
    end
    return PP(summands)
end

-{P <: Polynomial}(f::P) = P([Term(exponent, -coeff) for (exponent, coeff) in f.terms])
=={P <: Polynomial}(a::P, b::P) = a.terms == b.terms
!={P <: Polynomial}(a::P, b::P) = !(a == b)

function leading_term{P <: Polynomial}(p::P)
    if length(p.terms) > 0
        return p.terms[end]
    else
        throw(ArgumentError("The zero polynomial $( p ) does not have a leading term"))
    end
end

typealias AbstractModuleElement{P <: Polynomial} Union{P, AbstractArray{P}}
modulebasering{M <: Polynomial}(::Type{M}) = M
modulebasering{M <: AbstractArray}(::Type{M}) = eltype(M)

*{R<:BaseRing, NumVars,T<:Tuple}(x::AbstractArray{Polynomial{R, NumVars,T}}, a::Term{R, NumVars})= eltype(x)[ x_i*a for x_i in x]
*{R<:BaseRing, NumVars,T<:Tuple}(a::Term{R, NumVars}, x::AbstractArray{Polynomial{R, NumVars,T}})= eltype(x)[ a*x_i for x_i in x]


immutable _ModuleTerm{T <: Term}
    t::T
    pos::Int
end

isless{M <: _ModuleTerm}(a::M, b::M) = a.pos > b.pos || a.t < b.t

function leading_term{P <: Polynomial}(a::AbstractArray{P})
    i = findfirst(x->!iszero(x), a)
    if i>0
        return _ModuleTerm(leading_term(a[i]), i)
    else
        throw(ArgumentError("The zero element $( a ) does not have a leading term"))
    end
end

function _lcm_multipliers{NumVars}(a::Monomial{NumVars}, b::Monomial{NumVars})
     _lcm = Monomial(
        ntuple(Val{NumVars}) do i
            return max(a.e[i], b.e[i])
        end
    )

    multiplier_a = Monomial(
        ntuple(Val{NumVars}) do i
            return _lcm.e[i] - a.e[i]
        end
    )
    multiplier_b = Monomial(
        ntuple(Val{NumVars}) do i
            return _lcm.e[i] - b.e[i]
        end
    )

    return multiplier_a, multiplier_b
end

function _is_constant{T <: Term}(a::T)
    exponent, coeff = a
    return sum(exponent.e) == 0
end

function _lcm_multipliers{T <: Term}(a::T, b::T)
    exp_a, coeff_a = a
    exp_b, coeff_b = b

    m_exp_a, m_exp_b = _lcm_multipliers(exp_a, exp_b)

    return Term(m_exp_a, coeff_b), Term(m_exp_b, coeff_a)
end

function _monomial_div{T<: Term}(a::T, b::T)::Nullable{T}
    mul_a, mul_b = _lcm_multipliers(a, b)

    if _is_constant(mul_a)
        return mul_b // coefficient(mul_a)
    else
        return nothing
    end
end

function _monomial_div{T<: Term}(a::_ModuleTerm{T}, b::_ModuleTerm{T})::Nullable{T}
    if a.pos != b.pos
        return nothing
    else
        return _monomial_div(a.t, b.t)
    end
end

_maybe_lcm_multipliers{T <: Term}(a::T, b::T)::Nullable{Tuple{T,T}} = _lcm_multipliers(a,b)
function _maybe_lcm_multipliers{T <: Term}(a::_ModuleTerm{T}, b::_ModuleTerm{T})::Nullable{Tuple{T,T}}
    if a.pos != b.pos
        return nothing
    else
        return _lcm_multipliers(a.t, b.t)
    end
end

_monomials{P <: Polynomial}(f::P) = reverse(f.terms)
_monomials{P <: Polynomial}(x::AbstractArray{P}) = _MonomialsIter{P, typeof(x)}(x)
immutable _MonomialsIter{P <: Polynomial, M <: AbstractArray}
    f::M
    _MonomialsIter(f::AbstractArray{P}) = new(f)
end
import Base: start, next, done, eltype, length
start{M <: _MonomialsIter}(::M) = (1,0)
function done{M <: _MonomialsIter}(x::M, state::Tuple{Int,Int})
    row, term = state
    if row > length(x.f)
        return true
    end
    if term < length(x.f[row].terms)
        return false
    end
    for i in (row+1):length(x.f)
        if length(x.f[i].terms) > 0
            return false
        end
    end
    return true
end
function next{M <: _MonomialsIter}(x::M, state::Tuple{Int,Int})
    row, term = state
    if term >= length(x.f[row].terms)
        term = 0
        row += 1
        while length(x.f[row].terms) == 0
            row += 1
        end
    end
    return _ModuleTerm(x.f[row].terms[end-term], row), (row, term+1)
end

function _div_with_remainder{M <: AbstractModuleElement}(f::M, g::M)::Tuple{Nullable{modulebasering(M)}, M}
    if iszero(f)
        return zero(modulebasering(M)), f
    else
        lt_g = try
            leading_term(g)
        catch ArgumentError # g is zero
            return nothing, f
        end
        for monomial in _monomials(f)
            maybe_factor = _monomial_div(monomial, lt_g)
            if !isnull(maybe_factor)
                factor = get(maybe_factor)
                return factor, f - g*factor
            end
        end
        return nothing, f
    end
end

function _lead_div_with_remainder{M <: AbstractModuleElement}(f::M, g::M)::Tuple{Nullable{modulebasering(M)}, M}
    lt_f = try
        leading_term(f)
    catch ArgumentError
        return zero(modulebasering(M)), f
    end
    lt_g = try
        leading_term(g)
    catch ArgumentError # g is zero
        return nothing, f
    end
    maybe_factor = _monomial_div(lt_f, lt_g)
    if !isnull(maybe_factor)
        factor = get(maybe_factor)
        return factor, f - g*factor
    end
    return nothing, f
end


end
