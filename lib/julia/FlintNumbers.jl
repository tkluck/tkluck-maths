module FlintNumbers

using PolynomialRings
using CommutativeAlgebras
using Nemo

using PolynomialRings: exptype, leading_term
using PolynomialRings.Terms: monomial
import CommutativeAlgebras: _grb, _ideal
using ExactLinearAlgebra: kernel

_id_counter = 0

struct FlintNumber{P,N} <: Number
    x::Nemo.nf_elem
    FlintNumber{P,N}(x::Nemo.nf_elem) where {P,N} = new(x)
end

function _nemo_number_field end
function _named_values end
function _primitive_element end

import Base: copy
copy(f::Nemo.fmpq) = f # this needs to be defined for power_by_squaring to work
copy(f::Nemo.fmpq_poly) = f # this needs to be defined for power_by_squaring to work

function FlintNumberField(::Type{Q}) where Q<:QuotientRing
    P = ring(Q)
    MAX = typemax(exptype(P))

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
    F = FlintNumber{P, ID}

    # let's hope this is a primitive element
    α = sum(generators(P))
    M = hcat((coeffs(α^n) for n=0:N)...)

    K = kernel(M)
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
        KK = kernel(MM)
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
import CommutativeAlgebras: ring

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

import ExactLinearAlgebra: kernel
function kernel(M::AbstractMatrix{F}) where F <: FlintNumber
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
