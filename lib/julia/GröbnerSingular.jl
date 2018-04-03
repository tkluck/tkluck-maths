module GröbnerSingular

using PolynomialRings.Monomials: AbstractMonomial, num_variables, enumeratenz
import PolynomialRings.MonomialOrderings: MonomialOrder
using PolynomialRings.Terms: Term, coefficient, monomial
using PolynomialRings.Polynomials: Polynomial, PolynomialOver, monomialtype, terms
using PolynomialRings: generators, integral_fraction, base_extend

import Expect: ExpectProc, expect!

import Base: write, print
import PolynomialRings: gröbner_basis, gröbner_transformation
import PolynomialRings.Backends.Gröbner: Backend

struct SingularExpect <: Backend end

struct SingularProc{Stream} <: IO
    s::Stream
end

function singularproc(f)
    # -q suppresses banners on stderr
    # -t avoids a segfault when inputting large input
    # pty=true is default, but it's good to be explicit
    # because we really need the "> " prompts
    stream = ExpectProc(`Singular -q -t`, 86400, pty=true)
    singular = SingularProc(stream)
    try
        return f(singular)
    finally
        close(stream)
    end
end

print(sp::SingularProc, x::String) = print(sp.s, x)
print(sp::SingularProc, x::Char) = print(sp.s, x)
print(sp::SingularProc, x) = print(sp.s, x)

expect!(sp::SingularProc{<:ExpectProc}, x) = expect!(sp.s, x)
expect!(sp::SingularProc, x) = ""

print(ep::ExpectProc, x::Char) = print(ep.in_stream, x)

function set_ring!(singular::SingularProc, ::Type{P}) where P<:Polynomial
    println(singular, "ring R = 0, x(1..$(num_variables(monomialtype(P)))), (c,dp);")
end

print(singular::SingularProc, m::AbstractMonomial) = join(singular, ["x($ix)^$e" for (ix, e) in enumeratenz(m)], "*")

function print(singular::SingularProc, t::Term)
    print(singular, coefficient(t))
    print(singular, "*")
    print(singular, monomial(t))
end

print(singular::SingularProc, f::Polynomial) = iszero(f) ? print(singular, "0") : join(singular, terms(f), " + ")

function print(singular::SingularProc, a::AbstractArray{<:Polynomial})
    print(singular, "[")
    join(singular, a, ",")
    print(singular, "]")
end

function parse_polynomial(::Type{P}, a::AbstractString) where P<:Polynomial
    if startswith(a, "-")
        a = "0$a"
    end
    mapreduce(+, zero(P), split(a, "+")) do positive_term
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
                    throw(ErrorException("Can't parse Singular's output at $var (in $a)"))
                end
            end
        end
    end
end

function parse_vector(::Type{P}, a::AbstractString, module_dims) where P<:Polynomial
    if startswith(a, "[") && endswith(a, "]")
        entries = map(s->parse_polynomial(P, s), split(a[2:end-1], ","))
        res = zeros(P, module_dims)
        res[1:length(entries)] = entries
        res
    else
        throw(ErrorException("Can't parse Singular's output at $a"))
    end
end

function parse_array(::Type{P}, a::AbstractString, module_dims=nothing) where P<:Polynomial
    map(split(rstrip(a), "\n")) do line
        (_, poly) = split(line, "=", limit=2)
        module_dims === nothing ? parse_polynomial(P, poly) : parse_vector(P, poly, module_dims)
    end
end

function parse_matrix(::Type{P}, a::AbstractString) where P<:Polynomial
    entries = map(split(rstrip(a), "\n")) do line
        (lhs, poly) = split(line, "=", limit=2)
        m = match(r".*\[(?<row>[0-9]+),(?<col>[0-9]+)\]", lhs)
        if m === nothing
            throw(ErrorException("Can't parse Singular's output at $lhs (in $a)"))
        end
        (parse(Int, m[:row]), parse(Int, m[:col]), parse_polynomial(P, poly))
    end
    rows = maximum(e->e[1], entries)
    cols = maximum(e->e[2], entries)
    res = zeros(P, rows, cols)
    for (row, col, p) in entries
        res[row, col] = p
    end
    res
end

function singular_std(G::AbstractArray{P}) where P<:Polynomial
    singularproc() do singular
        expect!(singular, "> ")
        set_ring!(singular, P)
        expect!(singular, "> ")

        print(singular, "ideal I = ")
        join(singular, G, ", ")
        println(singular, ";")
        expect!(singular, "> ")

        println(singular, "std(I);")
        result = expect!(singular, "> ")

        parse_array(P, result)
    end
end

function singular_std(G::AbstractArray{<:AbstractArray{P}}) where P<:Polynomial
    isempty(G) && return copy(G)
    module_dims = size(G[1])

    singularproc() do singular
        expect!(singular, "> ")
        set_ring!(singular, P)
        expect!(singular, "> ")

        print(singular, "module M = ")
        join(singular, G, ", ")
        println(singular, ";")
        expect!(singular, "> ")

        println(singular, "std(M);")
        result = expect!(singular, "> ")

        sparse_gr = parse_array(P, result, module_dims)
        gr = issparse(G[1]) ? sparse_gr : collect.(sparse_gr)

        gr
    end
end

function singular_liftstd(G::AbstractArray{P}) where P<:Polynomial
    singularproc() do singular
        expect!(singular, "> ")
        set_ring!(singular, P)
        expect!(singular, "> ")

        print(singular, "ideal I = ")
        join(singular, G, ", ")
        println(singular, ";")
        expect!(singular, "> ")

        println(singular, "matrix T;")
        expect!(singular, "> ")
        println(singular, "liftstd(I, T);")
        result = expect!(singular, "> ")

        println(singular, "T;")
        result_matrix = expect!(singular, "> ")

        gr = parse_array(P, result)
        tr = parse_matrix(P, result_matrix)

        gr, tr
    end
end

function singular_liftstd(G::AbstractArray{<:AbstractArray{P}}) where P<:Polynomial
    isempty(G) && return copy(G)
    module_dims = size(G[1])

    singularproc() do singular
        expect!(singular, "> ")
        set_ring!(singular, P)
        expect!(singular, "> ")

        print(singular, "module M = ")
        join(singular, G, ", ")
        println(singular, ";")
        expect!(singular, "> ")

        println(singular, "matrix T;")
        expect!(singular, "> ")
        println(singular, "liftstd(M, T);")
        result = expect!(singular, "> ")

        println(singular, "T;")
        result_matrix = expect!(singular, "> ")

        sparse_gr = parse_array(P, result, module_dims)
        tr = parse_matrix(P, result_matrix)

        gr = issparse(G[1]) ? sparse_gr : collect.(sparse_gr)

        gr, tr
    end
end

const ApplicableBaserings = Union{BigInt}
const ApplicablePolynomial = PolynomialOver{<:ApplicableBaserings}
const ApplicableModuleElement = Union{ApplicablePolynomial, AbstractArray{<:ApplicablePolynomial}}
function gröbner_basis(::SingularExpect, ::MonomialOrder{:degrevlex}, polynomials::AbstractArray{<:ApplicableModuleElement}; kwds...)
    return singular_std(polynomials)
end
function gröbner_transformation(::SingularExpect, ::MonomialOrder{:degrevlex}, polynomials::AbstractArray{<:ApplicableModuleElement}; kwds...)
    gr, tr = singular_liftstd(polynomials)
    return gr, sparse(transpose(tr)) # opposite convention for matrix multiplication in Singular compared to us
end

const RationalPolynomial = PolynomialOver{Rational{BigInt}}
const RationalModuleElement = Union{RationalPolynomial, AbstractArray{<:RationalPolynomial}}
function gröbner_basis(::SingularExpect, ::MonomialOrder{:degrevlex}, polynomials::AbstractArray{<:RationalModuleElement}; kwds...)
    integral_polynomials = [p for (p, _) in integral_fraction.(polynomials)]
    return base_extend.(gröbner_basis(SingularExpect(), MonomialOrder{:degrevlex}(), integral_polynomials))
end
function gröbner_transformation(::SingularExpect, ::MonomialOrder{:degrevlex}, polynomials::AbstractArray{<:RationalModuleElement}; kwds...)
    integral_polynomials = [p for (p, _) in integral_fraction.(polynomials)]
    multipliers          = [D for (_, D) in integral_fraction.(polynomials)]
    gr, tr = gröbner_transformation(SingularExpect(), MonomialOrder{:degrevlex}(), integral_polynomials)
    for (i,D) in enumerate(multipliers)
        tr[:,i] *= D
    end
    return base_extend.(gr), base_extend.(tr)
end


import PolynomialRings
PolynomialRings.Backends.Gröbner.set_default(SingularExpect())

end
