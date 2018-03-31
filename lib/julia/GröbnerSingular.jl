module GröbnerSingular

using PolynomialRings.Monomials: AbstractMonomial, num_variables, enumeratenz
import PolynomialRings.MonomialOrderings: MonomialOrder
using PolynomialRings.Terms: Term, coefficient, monomial
using PolynomialRings.Polynomials: Polynomial, PolynomialOver, monomialtype, terms
using PolynomialRings: generators

import Expect: ExpectProc, expect!

import Base: write, print
import PolynomialRings: gröbner_basis
import PolynomialRings.Backends.Gröbner: Backend

struct SingularExpect <: Backend end

struct SingularProc{Stream} <: IO
    s::Stream
end

print(sp::SingularProc, x::String) = print(sp.s, x)
print(sp::SingularProc, x::Char) = print(sp.s, x)
print(ep::ExpectProc, x::Char) = print(ep.in_stream, x)
print(sp::SingularProc, x) = print(sp.s, x)
expect!(sp::SingularProc{<:ExpectProc}, x) = expect!(sp.s, x)
expect!(sp::SingularProc, x) = ""

function set_ring!(singular::SingularProc, ::Type{P}) where P<:Polynomial
    println(singular, "ring R = 0, x(1..$(num_variables(monomialtype(P)))), dp;")
end

print(singular::SingularProc, m::AbstractMonomial) = join(singular, ["x($ix)^$e" for (ix, e) in enumeratenz(m)], "*")

function print(singular::SingularProc, t::Term)
    print(singular, coefficient(t))
    print(singular, "*")
    print(singular, monomial(t))
end

print(singular::SingularProc, f::Polynomial) = join(singular, terms(f), " + ")

function parse_array(::Type{P}, a::AbstractString) where P<:Polynomial
    map(split(rstrip(a), "\n")) do line
        (_, poly) = split(line, "=", limit=2)
        mapreduce(+, zero(P), split(poly, "+")) do positive_term
            mapfoldl(-, split(positive_term, "-")) do term
                mapreduce(*, one(P), split(term, "*")) do var
                    m = match(r"
                        x\(  (?<varnum>[0-9]+)  \)  (?:\^  (?<exp>[0-9]+)  )?
                        |
                        (?<num>[0-9]+)
                    "x, var)
                    if m !== nothing && m[:exp] !== nothing
                        varnum = parse(Int, m[:varnum])
                        exp = parse(Int, m[:exp])
                        generators(P)[varnum]^exp
                    elseif m !== nothing && m[:varnum] !== nothing
                        varnum = parse(Int, m[:varnum])
                        generators(P)[varnum]
                    elseif m !== nothing && m[:num] !== nothing
                        num = parse(BigInt, m[:num])
                        P(num)
                    else
                        throw(ErrorException("Can't parse Singular's output at $var (in $poly)"))
                    end
                end
            end
        end
    end
end

function singular_std(G::AbstractArray{P}) where P<:Polynomial
    singular = SingularProc( ExpectProc(`Singular -q`, 86400) )

    expect!(singular, "> ")
    set_ring!(singular, P)
    expect!(singular, "> ")

    print(singular, "ideal I = ")
    join(singular, G, ", ")
    println(singular, ";")

    expect!(singular, "> ")

    println(singular, "std(I);")

    result = expect!(singular, "> ")

    return parse_array(P, result)
end

const ApplicableBaserings = Union{BigInt}
const ApplicablePolynomial = PolynomialOver{<:ApplicableBaserings}
function gröbner_basis(::SingularExpect, ::MonomialOrder{:degrevlex}, polynomials::AbstractArray{<:ApplicablePolynomial}; kwds...)
    return singular_std(polynomials)
end

import PolynomialRings
PolynomialRings.Backends.Gröbner.set_default(SingularExpect())

end
