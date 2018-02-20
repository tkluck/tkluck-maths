module FlintNumbers

using PolynomialRings
using Nemo

using PolynomialRings: exptype, leading_term
using PolynomialRings.Terms: monomial
import PolynomialRings.QuotientRings: QuotientRing, ring, _grb, _ideal, monomial_basis
import PolynomialRings.Util.LinAlgUtil: nullspace, AbstractExactNumber

struct FlintNumber{P,Q} <: AbstractExactNumber
    x::Nemo.nf_elem
    FlintNumber{P,Q}(x::Nemo.nf_elem) where {P,Q} = new(x)
end

function _nemo_number_field end
function _named_values end
function _primitive_element end

import Base: copy
copy(f::Nemo.fmpq) = f # this needs to be defined for power_by_squaring to work
copy(f::Nemo.fmpq_poly) = f # this needs to be defined for power_by_squaring to work

function FlintNumberField(::Type{Q}) where Q<:QuotientRing
    P = ring(Q)

    monomials = monomial_basis(Q)
    coeffs(f) = (ff = rem(f,_ideal(Q)); [coefficient(ff,m,variablesymbols(P)...) for m in monomials])
    N = length(monomials)

    F = FlintNumber{P, Q}

    # let's hope this is a primitive element
    α = sum(generators(P))
    M = hcat((coeffs(α^n) for n=0:N)...)

    K = nullspace(M)
    if size(K,2) != 1 || iszero(K[end])
        throw(AssertionError("OOPS! My naive guess for a primitive element doesn't work. Maybe this is not a number field?"))
    end

    R,x = Nemo.PolynomialRing(Nemo.QQ, "x")
    minpoly = sum(Nemo.QQ(K[i+1])*x^i for i=0:N)
    NNF,x = Nemo.NumberField(minpoly, "x")

    named_values = Dict{Symbol, F}()
    for β in variablesymbols(P)
        MM = copy(M)
        MM[:,end] = coeffs(P(β))
        KK = nullspace(MM)
        named_values[β] = F( -sum(KK[i+1]*x^i for i=0:(N-1))//KK[end] )
    end

    @eval _nemo_number_field(::Type{$F}) = $NNF
    @eval _named_values(::Type{$F}) = $named_values
    @eval _primitive_element(::Type{$F}) = $α
    @eval _ideal(::Type{$F}) = _ideal($Q)

    return F
end

import Base: zero, one, +,-,*,//, convert, promote_rule, ==, !=
import PolynomialRings: allvariablesymbols, fraction_field

zero(::Type{F}) where F<:FlintNumber = F(zero(_nemo_number_field(F)))
one(::Type{F})  where F<:FlintNumber = F(one(_nemo_number_field(F)))
+(a::FlintNumber) = a
-(a::F) where F<:FlintNumber = F(-a.x)
+(a::F, b::F) where F<:FlintNumber = F(a.x+b.x)
-(a::F, b::F) where F<:FlintNumber = F(a.x-b.x)
*(a::F, b::F) where F<:FlintNumber = F(a.x*b.x)
//(a::F, b::F) where F<:FlintNumber = F(a.x//b.x)
==(a::F, b::F) where F<:FlintNumber = a.x == b.x
!=(a::F, b::F) where F<:FlintNumber = a.x != b.x
ring(::Type{F}) where F<:FlintNumber{P} where P = P
allvariablesymbols(::Type{F}) where F<:FlintNumber = allvariablesymbols(ring(F))
copy(a::FlintNumber) = a
fraction_field(::Type{F}) where F<:FlintNumber = F

promote_rule(::Type{F}, ::Type{N}) where {F<:FlintNumber,N<:Number} = F
promote_rule(::Type{F}, ::Type{C}) where {F<:FlintNumber,C<:PolynomialRings.Constants.Constant} = F
convert(::Type{F}, a::Number) where F<:FlintNumber = F(_nemo_number_field(F)(a))
convert(::Type{F}, a::Rational) where F<:FlintNumber = F(numerator(a)) // F(denominator(a))
convert(::Type{F}, a::F) where F<:FlintNumber = a
convert(::Type{F}, f::P) where F<:FlintNumber{P} where P<:Polynomial  = f(;_named_values(F)...)

+(a::Number,b::FlintNumber) = +(promote(a,b)...)
+(a::FlintNumber,b::Number) = +(promote(a,b)...)
-(a::Number,b::FlintNumber) = -(promote(a,b)...)
-(a::FlintNumber,b::Number) = -(promote(a,b)...)
*(a::Number,b::FlintNumber) = *(promote(a,b)...)
*(a::FlintNumber,b::Number) = *(promote(a,b)...)
//(a::Number,b::FlintNumber) = //(promote(a,b)...)
//(a::FlintNumber,b::Number) = //(promote(a,b)...)
==(a::Number,b::FlintNumber) = ==(promote(a,b)...)
==(a::FlintNumber,b::Number) = ==(promote(a,b)...)
!=(a::Number,b::FlintNumber) = !=(promote(a,b)...)
!=(a::FlintNumber,b::Number) = !=(promote(a,b)...)


export FlintNumberField

function nullspace(M::Matrix{F}) where F <: FlintNumber
    MS = Nemo.MatrixSpace(_nemo_number_field(F), size(M)...)

    MM = MS(map(f->f.x, M))

    nullity, K = nullspace(MM)

    map(F, K.entries)
end

import Base:show

function show(io::IO, x::FlintNumber)
    N = Nemo.degree(_nemo_number_field(typeof(x)))
    α = _primitive_element(typeof(x))
    f = rem(
            sum(0:(N-1)) do n
                c = Nemo.coeff(x.x, n)
                if !iszero(c)
                    Rational{BigInt}(c) * α^n
                else
                    zero(Rational{BigInt})
                end
            end,
            _ideal(typeof(x)),
           )
    print(io, f)
end


end
