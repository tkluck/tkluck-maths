module CommutativeAlgebras

using Nulls

using PolynomialRings.Polynomials: Polynomial
using PolynomialRings.Groebner: groebner_basis
using PolynomialRings.Util: lazymap

abstract type _AbstractCommutativeAlgebra{P<:Polynomial} end
AbstractCommutativeAlgebra{P<:Polynomial} = Union{P, _AbstractCommutativeAlgebra{P}}

# -----------------------------------------------------------------------------
#
# Imports for overloading
#
# -----------------------------------------------------------------------------
import Base: promote_rule, convert
import Base: zero, one, in, rem, issubset
import Base: +,-,*,/,//,==,!=
import PolynomialRings: generators, expansion
import PolynomialRings: allvariablesymbols, fraction_field
import PolynomialRings.Expansions: _expansion

# -----------------------------------------------------------------------------
#
# Ideals
#
# -----------------------------------------------------------------------------

mutable struct Ideal{P<:Polynomial}
    generators::AbstractVector{P}
    _grb::Union{Null, AbstractVector{P}}
    _trns::Union{Null, AbstractMatrix{P}}
end
Ideal(generators::AbstractVector{<:Polynomial}) = Ideal(generators, null, null)
Ideal(generators::Polynomial...) = Ideal(collect(generators), null, null)
generators(I::Ideal) = I.generators
function _grb(I::Ideal)
    if isnull(I._grb)
        I._grb, I._trns = groebner_basis(I.generators)
    end
    I._grb
end
function _trns(I::Ideal)
    if isnull(I._grb)
        I._grb, I._trns = groebner_basis(I.generators)
    end
    I._trns
end
ring(I::Ideal{P}) where P<:Polynomial = P

rem(f, I::Ideal) = rem(ring(I)(f), _grb(I))
in(f, I::Ideal) = iszero(rem(f, I))

issubset(I::Ideal{P}, J::Ideal{P}) where P<:Polynomial = all(g in J for g in generators(I))
==(I::Ideal{P}, J::Ideal{P}) where P<:Polynomial = I⊆J && J⊆I

# -----------------------------------------------------------------------------
#
# Quotient rings
#
# -----------------------------------------------------------------------------
struct QuotientRing{P<:Polynomial, ID} <: _AbstractCommutativeAlgebra{P}
    f::P
    QuotientRing{P,ID}(f::P) where {P,ID} = new(rem(f, _ideal(QuotientRing{P,ID})))
end
ring(::Type{Q}) where Q<:QuotientRing{P} where P = P

function _ideal end
_id_counter = 0
function /(::Type{P}, I::Ideal{P}) where P<:Polynomial
    global _id_counter
    ID = (_id_counter+=1)

    R = QuotientRing{P, ID}
    @eval _ideal(::Type{$R}) = $I

    return R
end


# -----------------------------------------------------------------------------
#
# Number fields
#
# -----------------------------------------------------------------------------
using PolynomialRings: Term,AbstractMonomial, leading_term, exptype, basering, variablesymbols
using PolynomialRings.Terms: monomial, coefficient
using ExactLinearAlgebra: kernel
TermOver{C} = Term{<:AbstractMonomial,C}
PolynomialOver{C} = Polynomial{<:AbstractVector{<:TermOver{C}}}
struct NumberField{P<:Polynomial, ID} <: Number #_AbstractCommutativeAlgebra{P}
    #coeffs::NTuple{N, C}
    f::P
    NumberField{P, ID}(f::P) where P<:Polynomial where ID = new(rem(f, _ideal(NumberField{P, ID})))
end
ring(::Type{F}) where F<:NumberField{P} where P = P

function _monomial_basis end
function _monomial_coeffs end
function _primitive_matrix end
function _primitive_coeffs end

function NumberField(::Type{Q}) where Q<:QuotientRing
    P = ring(Q)
    I = exptype(P)
    MAX = typemax(I)

    leading_monomials = [monomial(leading_term(f)).e for f in _grb(_ideal(Q))]

    # lattice computation: find the monomials that are not divisible
    # by any of the leading terms
    nvars = nfields(eltype(leading_monomials))
    rectangular_bounds = fill(MAX, nvars)
    for m in leading_monomials
        if count(!iszero, m) == 1
            i = findfirst(m)
            rectangular_bounds[i] = min( m[i], rectangular_bounds[i] )
        end
    end

    if !all(b->b<MAX, rectangular_bounds)
        throw("$Q is not a number field; it is infinite dimensional")
    end

    divisible = BitArray(rectangular_bounds...)
    for m in leading_monomials
        block = [(m[i]+1):b for (i,b) in enumerate(rectangular_bounds)]
        setindex!(divisible, true, block...)
    end

    monomials = [tuple([i-1 for i in ind2sub(divisible,i)]...) for i in eachindex(divisible) if !divisible[i]]
    coeffs(f) = (ff = rem(f,_ideal(Q)); [coefficient(ff,m,variablesymbols(P)...) for m in monomials])
    N = length(monomials)

    global _id_counter
    ID = (_id_counter+=1)
    F = NumberField{P, ID}

    # let's hope this is a primitive element
    α = sum(generators(P))
    M = hcat((coeffs(α^n) for n=0:N)...)

    K = kernel(M)
    if size(K,2) != 1
        throw(AssertionError("OOPS! My naive guess for a primitive element doesn't work. Maybe this is not a number field?"))
    end

    @eval _ideal(::Type{$F}) = _ideal($Q)

    return F
end

#function convert(::Type{F}, f::P) where F<:NumberField{P, ID} where {P<:Polynomial, ID}
#    return F(f)
#end

# -----------------------------------------------------------------------------
#
# Operations for QuotientRings (mostly pass-through)
#
# -----------------------------------------------------------------------------

zero(::Type{Q}) where Q<:QuotientRing{P} where P<:Polynomial = Q(zero(P))
one(::Type{Q})  where Q<:QuotientRing{P} where P<:Polynomial = Q(one(P))
+(a::QuotientRing) = a
-(a::Q) where Q<:QuotientRing = Q(-a.f)
+(a::Q, b::Q)  where Q<:QuotientRing = Q(a.f+b.f)
-(a::Q, b::Q)  where Q<:QuotientRing = Q(a.f-b.f)
*(a::Q, b::Q)  where Q<:QuotientRing = Q(a.f*b.f)
//(a::Q, b::Q) where Q<:QuotientRing = Q(a.f//b.f)
==(a::Q, b::Q) where Q<:QuotientRing = a.f == b.f
allvariablesymbols(::Type{Q}) where Q<:QuotientRing = allvariablesymbols(ring(Q))

# -----------------------------------------------------------------------------
#
# Operations for NumberFields (mostly pass-through)
#
# -----------------------------------------------------------------------------
zero(::Type{Q}) where Q<:NumberField{P} where P<:Polynomial = Q(zero(P))
one(::Type{Q})  where Q<:NumberField{P} where P<:Polynomial = Q(one(P))
+(a::NumberField) = a
-(a::Q) where Q<:NumberField = Q(-a.f)
+(a::Q, b::Q)  where Q<:NumberField = Q(a.f+b.f)
-(a::Q, b::Q)  where Q<:NumberField = Q(a.f-b.f)
*(a::Q, b::Q)  where Q<:NumberField = Q(a.f*b.f)
//(a::Q, b::Q) where Q<:NumberField = Q(a.f//b.f)
==(a::Q, b::Q) where Q<:NumberField = a.f == b.f
allvariablesymbols(::Type{Q}) where Q<:NumberField = allvariablesymbols(ring(Q))

fraction_field(::Type{N}) where N<:NumberField = N

# -----------------------------------------------------------------------------
#
# Operations through promotion
#
# -----------------------------------------------------------------------------
_Q{P<:Polynomial} = Union{QuotientRing{P}, NumberField{P}}
function promote_rule(::Type{Q}, ::Type{C}) where Q<:_Q{P} where {P<:Polynomial,C}
    rule_for_P = typejoin( promote_rule(P,C), promote_rule(C,P) )
    if rule_for_P === P
        return Q
    else
        return Union{}
    end
end

function convert(::Type{Q}, c::C) where Q<:_Q{P} where {P<:Polynomial,C}
    Q(convert(P, c))
end

convert(::Type{Q}, q::Q) where Q<:_Q{P} where P<:Polynomial = q

# -----------------------------------------------------------------------------
#
# Operations through promotion
#
# -----------------------------------------------------------------------------
+(a::QuotientRing, b) = +(promote(a,b)...)
+(a, b::QuotientRing) = +(promote(a,b)...)
-(a::QuotientRing, b) = -(promote(a,b)...)
-(a, b::QuotientRing) = -(promote(a,b)...)
*(a::QuotientRing, b) = *(promote(a,b)...)
*(a, b::QuotientRing) = *(promote(a,b)...)
//(a::QuotientRing, b) = //(promote(a,b)...)
//(a, b::QuotientRing) = //(promote(a,b)...)
==(a::QuotientRing, b) = ==(promote(a,b)...)
==(a, b::QuotientRing) = ==(promote(a,b)...)
!=(a::QuotientRing, b) = !=(promote(a,b)...)
!=(a, b::QuotientRing) = !=(promote(a,b)...)

+(a::NumberField, b::Number) = +(promote(a,b)...)
+(a::Number, b::NumberField) = +(promote(a,b)...)
-(a::NumberField, b::Number) = -(promote(a,b)...)
-(a::Number, b::NumberField) = -(promote(a,b)...)
*(a::NumberField, b::Number) = *(promote(a,b)...)
*(a::Number, b::NumberField) = *(promote(a,b)...)
//(a::NumberField, b::Number) = //(promote(a,b)...)
//(a::Number, b::NumberField) = //(promote(a,b)...)
==(a::NumberField, b::Number) = ==(promote(a,b)...)
==(a::Number, b::NumberField) = ==(promote(a,b)...)
!=(a::NumberField, b::Number) = !=(promote(a,b)...)
!=(a::Number, b::NumberField) = !=(promote(a,b)...)


# -----------------------------------------------------------------------------
#
# Exports
#
# -----------------------------------------------------------------------------
export Ideal, NumberField


end
