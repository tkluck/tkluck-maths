module Polynomials

import Base: +,==,!=,*,//,-,convert,promote_rule,show,cmp,isless,zero,eltype

immutable Monomial{NumVars}
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

macro myplus(numvars)
    expr = :( tuple() )
    for i=1:numvars
        push!(expr.args, :( a.e[$i] + b.e[$i] ))
    end
    return quote
        Monomial{$numvars}($expr)
    end
end

+(a::Monomial{1}, b::Monomial{1}) = @myplus 1
+(a::Monomial{2}, b::Monomial{2}) = @myplus 2
+(a::Monomial{3}, b::Monomial{3}) = @myplus 3
+(a::Monomial{4}, b::Monomial{4}) = @myplus 4
+(a::Monomial{5}, b::Monomial{5}) = @myplus 5
+(a::Monomial{6}, b::Monomial{6}) = @myplus 6
+(a::Monomial{7}, b::Monomial{7}) = @myplus 7
+(a::Monomial{8}, b::Monomial{8}) = @myplus 8
+(a::Monomial{9}, b::Monomial{9}) = @myplus 9
+(a::Monomial{10}, b::Monomial{10}) = @myplus 10
+(a::Monomial{11}, b::Monomial{11}) = @myplus 11
+(a::Monomial{12}, b::Monomial{12}) = @myplus 12
+(a::Monomial{13}, b::Monomial{13}) = @myplus 13
+(a::Monomial{14}, b::Monomial{14}) = @myplus 14
+(a::Monomial{15}, b::Monomial{15}) = @myplus 15
+(a::Monomial{16}, b::Monomial{16}) = @myplus 16
+(a::Monomial{17}, b::Monomial{17}) = @myplus 17
+(a::Monomial{18}, b::Monomial{18}) = @myplus 18
+(a::Monomial{19}, b::Monomial{19}) = @myplus 19
+(a::Monomial{20}, b::Monomial{20}) = @myplus 20
+(a::Monomial{21}, b::Monomial{21}) = @myplus 21
+(a::Monomial{22}, b::Monomial{22}) = @myplus 22
+(a::Monomial{23}, b::Monomial{23}) = @myplus 23
+(a::Monomial{24}, b::Monomial{24}) = @myplus 24
+(a::Monomial{25}, b::Monomial{25}) = @myplus 25
+(a::Monomial{26}, b::Monomial{26}) = @myplus 26
+(a::Monomial{27}, b::Monomial{27}) = @myplus 27
+(a::Monomial{28}, b::Monomial{28}) = @myplus 28
+(a::Monomial{29}, b::Monomial{29}) = @myplus 29
+(a::Monomial{30}, b::Monomial{30}) = @myplus 30
+(a::Monomial{31}, b::Monomial{31}) = @myplus 31
+(a::Monomial{32}, b::Monomial{32}) = @myplus 32
+(a::Monomial{33}, b::Monomial{33}) = @myplus 33
+(a::Monomial{34}, b::Monomial{34}) = @myplus 34
+(a::Monomial{35}, b::Monomial{35}) = @myplus 35
+(a::Monomial{36}, b::Monomial{36}) = @myplus 36
+(a::Monomial{37}, b::Monomial{37}) = @myplus 37
+(a::Monomial{38}, b::Monomial{38}) = @myplus 38
+(a::Monomial{39}, b::Monomial{39}) = @myplus 39
+(a::Monomial{40}, b::Monomial{40}) = @myplus 40


function cmp{M <: Monomial}(a::M, b::M)
    # degrevlex
    degcmp = cmp(sum(a.e), sum(b.e))
    if degcmp == 0
        return -cmp(a.e, b.e)
    else
        return degcmp
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

# ------------------------------------------------------
# utility iteration for iterating over a 2-dimensional cartesian product
# in the following pattern:
#
#  1  3  6 10
#  2  5  9 13
#  4  8 12 15
#  7 11 14 16
#
# i.e. in diagonals. The usefulness for polynomial multiplication is that the
# terms on both axes are sorted by total degree, so on the diagonal, we can hope
# to have similar total degrees of the products of the terms. That means that the
# sort step that immediately follows has less sorting to do.
#
# The result is a very big win, as it turns out that the sorting is a major bottleneck
# for the kind of workload I currently have.
immutable DiagonalIterState
    start_row::Int
    diagonal_index::Int
end

immutable DiagonalIter{A <: AbstractArray, B <: AbstractArray}
    rows::A
    cols::B
end

import Base: start, next, done, eltype, length
start{D <: DiagonalIter}(x::D)= DiagonalIterState(0,0)
function next{A <: AbstractArray, B <: AbstractArray}(x::DiagonalIter{A,B}, state::DiagonalIterState)::Tuple{ Tuple{eltype(A),eltype(B)}, DiagonalIterState }
    next_diag_index = state.diagonal_index + 1
    next_start_row = state.start_row

    tentative_next_row = next_start_row - next_diag_index
    tentative_next_col = 1 + next_diag_index

    if tentative_next_row < 1 || tentative_next_col > length(x.cols)
        next_start_row = state.start_row + 1

        if next_start_row > length(x.rows)
            next_diag_index = next_start_row - length(x.rows)
        else
            next_diag_index = 0
        end

    end

    a = x.rows[next_start_row - next_diag_index]
    b = x.cols[1 + next_diag_index]
    item = (a,b)
    next_state = DiagonalIterState(next_start_row, next_diag_index)
    return item, next_state
end
done{D <: DiagonalIter}(x::D, state::DiagonalIterState)= length(x.rows) == 0 || length(x.cols) == 0 || state.start_row >= length(x.rows) + length(x.cols) - 1
eltype{A <: AbstractArray, B <: AbstractArray}(::Type{DiagonalIter{A,B}}) = Tuple{eltype(A), eltype(B)}
length{D <: DiagonalIter}(x::D)= length(x.rows) * length(x.cols)
# end of utility iterator
# ------------------------------------------------------

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
    for ((exp_a,coeff_a),(exp_b,coeff_b)) in DiagonalIter(a.terms, b.terms)
        summands[ix] = Term(exp_a + exp_b, coeff_a * coeff_b)
        ix += 1
    end
    assert( ix == length(summands)+1)
    sort!(summands, by=t -> t[1], alg=QuickSort)

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

function leading_term{M <: Union{_HomModuleElement,_ModuleElement}}(a::M)
    i = findfirst(x->!iszero(x), a)
    if i>0
        return _ModuleTerm(leading_term(a[i]), i)
    else
        throw(ArgumentError("The zero element $( a ) does not have a leading term"))
    end
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

_monomials{P <: Polynomial}(f::P) = reverse(f.terms)
_monomials{M <: Union{_HomModuleElement,_ModuleElement}}(x::M) = _MonomialsIter(x)

immutable _MonomialsIter{M <: Union{_HomModuleElement,_ModuleElement}}
    f::M
end
start{M <: _MonomialsIter}(::M) = (1,0)
function done{M <: _MonomialsIter}(x::M, state::Tuple{Int,Int})
    row, term = state
    return row > length(x.f) || (
        term == length(x.f[row].terms) && all(length(x.f[i].terms) == 0 for i in (row+1):length(x.f))
    )
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

function random_polynomial()
    res = sum([ rand(0:100) * x^i for i = 0:10 ])
end

function random_matrix()
    return Matrix([ random_polynomial() for i = 1:10, j = 1:10 ])
end

end
