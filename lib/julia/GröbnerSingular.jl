module GröbnerSingular

using SparseArrays: AbstractSparseArray, issparse, sparse

using PolynomialRings.Monomials: AbstractMonomial, num_variables, enumeratenz
import PolynomialRings.MonomialOrderings: MonomialOrder
using PolynomialRings.Terms: Term, coefficient, monomial
using PolynomialRings.Polynomials: Polynomial, PolynomialOver, monomialtype, terms
using PolynomialRings: generators, integral_fraction, base_extend, basering
using GaloisFields: AbstractGaloisField, char

import Expect: ExpectProc, expect!

import Base: write, print
import PolynomialRings: gröbner_basis, gröbner_transformation, lift
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
    stream = ExpectProc(`Singular -q -t --cntrlc=q`, 86400, pty=true)
    singular = SingularProc(stream)
    try
        return f(singular)
    finally
        # be aggressive, as it doesn't seem to listen to friendly
        # requests to shut down
        kill(stream, 9)
    end
end

print(sp::SingularProc, x::String) = print(sp.s, x)
print(sp::SingularProc, x::Char)   = print(sp.s, x)
print(sp::SingularProc, x)         = print(sp.s, x)

struct SingularError
    msg::String
end

function expect_output!(sp::SingularProc{<:ExpectProc})
    output = expect!(sp.s, "> ")
    if "// **" ⊆ output
        msg = mapreduce(m->m[:msg], *, eachmatch(r"// \*\* (?<msg>.*)$"m, output))
        throw(SingularError("Error from Singular: $msg"))
    elseif "   ?" ⊆ output
        msg = mapreduce(m->m[:msg], *, eachmatch(r"   \? (?<msg>.*)$"m, output))
        throw(SingularError("Error from Singular: $msg"))
    else
       return output
    end
end

expect!(sp::SingularProc, x) = ""

print(ep::ExpectProc, x::Char) = print(ep.in_stream, x)

function set_ring!(singular::SingularProc, ::Type{P}) where P<:Polynomial
    q = char(basering(P))
    println(singular, "ring R = $q, x(1..$(num_variables(monomialtype(P)))), (c,dp);")
end

function print(singular::SingularProc, x::Rational)
    print(singular, numerator(x))
    print(singular, "/")
    print(singular, denominator(x))
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

function print(singular::SingularProc, a::AbstractSparseArray{<:Polynomial})
    print(singular,"0")
    for i in LinearIndices(a)[findall(!iszero, a)]
        print(singular, "+")
        print(singular, "gen($i)*(")
        print(singular, a[i])
        print(singular, ")")
    end
end

function parse_polynomial(::Type{P}, a::AbstractString) where P<:Polynomial
    if startswith(a, "-")
        a = "0$a"
    end
    mapreduce(+, split(a, "+"), init=zero(P)) do positive_term
        mapfoldl(-, split(positive_term, "-")) do term
            mapreduce(*, split(term, "*"), init=one(P)) do var
                m = match(r"
                    x\(  (?<varnum>[0-9]+)  \)  (?:\^  (?<exp>[0-9]+)  )?
                    |
                    ( (?<num>[0-9]+) (/(?<denom>[0-9]+))? )
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
                    if m[:denom] !== nothing
                        denom = parse(BigInt, m[:denom])
                        P(num//denom)
                    else
                        P(num)
                    end
                else
                    throw(SngularError("Can't parse Singular's output at $var (in $a)"))
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
        throw(SngularError("Can't parse Singular's output at $a"))
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
            throw(SngularError("Can't parse Singular's output at $lhs (in $a)"))
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
        expect_output!(singular)
        set_ring!(singular, P)
        expect_output!(singular)

        print(singular, "ideal I = ")
        join(singular, G, ", ")
        println(singular, ";")
        expect_output!(singular)

        println(singular, "std(I);")
        result = expect_output!(singular)

        parse_array(P, result)
    end
end

function singular_std(G::AbstractArray{<:AbstractArray{P}}) where P<:Polynomial
    module_dims = size(G[1])

    singularproc() do singular
        expect_output!(singular)
        set_ring!(singular, P)
        expect_output!(singular)

        print(singular, "module M = ")
        join(singular, G, ", ")
        println(singular, ";")
        expect_output!(singular)

        println(singular, "std(M);")
        result = expect_output!(singular)

        sparse_gr = parse_array(P, result, module_dims)
        gr = issparse(G[1]) ? sparse_gr : collect.(sparse_gr)

        gr
    end
end

function singular_liftstd(G::AbstractArray{P}) where P<:Polynomial
    singularproc() do singular
        expect_output!(singular)
        set_ring!(singular, P)
        expect_output!(singular)

        print(singular, "ideal I = ")
        join(singular, G, ", ")
        println(singular, ";")
        expect_output!(singular)

        println(singular, "matrix T;")
        expect_output!(singular)
        println(singular, "liftstd(I, T);")
        result = expect_output!(singular)

        println(singular, "T;")
        result_matrix = expect_output!(singular)

        gr = parse_array(P, result)
        tr = parse_matrix(P, result_matrix)

        gr, tr
    end
end

function singular_liftstd(G::AbstractArray{<:AbstractArray{P}}) where P<:Polynomial
    module_dims = size(G[1])

    singularproc() do singular
        expect_output!(singular)
        set_ring!(singular, P)
        expect_output!(singular)

        print(singular, "module M = ")
        join(singular, G, ", ")
        println(singular, ";")
        expect_output!(singular)

        println(singular, "matrix T;")
        expect_output!(singular)
        println(singular, "liftstd(M, T);")
        result = expect_output!(singular)

        println(singular, "T;")
        result_matrix = expect_output!(singular)

        sparse_gr = parse_array(P, result, module_dims)
        tr = parse_matrix(P, result_matrix)

        gr = issparse(G[1]) ? sparse_gr : collect.(sparse_gr)

        gr, tr
    end
end

function singular_lift(G::AbstractArray{P}, y::P) where P<:Polynomial
    singularproc() do singular
        expect_output!(singular)
        set_ring!(singular, P)
        expect_output!(singular)

        print(singular, "ideal I = ")
        join(singular, G, ", ")
        println(singular, ";")
        expect_output!(singular)

        print(singular, "poly y = ")
        print(singular, y)
        println(singular, ";")
        expect_output!(singular)

        try
            println(singular, "lift(I, y);")
            result = expect_output!(singular)
            return parse_matrix(P, result)
        catch e
            if isa(e, SingularError) && (
                "not a proper submodule" ⊆ e.msg ||
                "2nd module does not lie in the first" ⊆ e.msg
            )
                return nothing
            else
                rethrow(e)
            end
        end
    end
end

function singular_lift(G::AbstractArray{<:A}, y::A) where A<:AbstractArray{P} where P<:Polynomial
    module_dims = size(G[1])

    singularproc() do singular
        expect_output!(singular)
        set_ring!(singular, P)
        expect_output!(singular)

        print(singular, "module M = ")
        join(singular, G, ", ")
        println(singular, ";")
        expect_output!(singular)

        print(singular, "vector y = ")
        print(singular, y)
        println(singular, ";")
        expect_output!(singular)

        try
            println(singular, "lift(M, y);")
            result = expect_output!(singular)
            return parse_matrix(P, result)
        catch e
            if isa(e, SingularError) && (
                "not a proper submodule" ⊆ e.msg ||
                "2nd module does not lie in the first" ⊆ e.msg
            )
                return nothing
            else
                rethrow(e)
            end
        end
    end
end

const ApplicableBaserings = Union{BigInt,Rational{BigInt},AbstractGaloisField}
const ApplicablePolynomial = PolynomialOver{<:ApplicableBaserings}
const ApplicableModuleElement{P<:ApplicablePolynomial} = Union{P, AbstractArray{<:P}}
function gröbner_basis(::SingularExpect, ::MonomialOrder{:degrevlex}, polynomials::AbstractArray{<:ApplicableModuleElement}; kwds...)
    isempty(polynomials) && return copy(polynomials)
    return singular_std(polynomials)
end
function gröbner_transformation(::SingularExpect, ::MonomialOrder{:degrevlex}, polynomials::AbstractArray{<:ApplicableModuleElement}; kwds...)
    isempty(polynomials) && return copy(polynomials), eye(P, 0)
    gr, tr = singular_liftstd(polynomials)
    return gr, sparse(transpose(tr)) # opposite convention for matrix multiplication in Singular compared to us
end

function lift(::SingularExpect, G::AbstractArray{<:ApplicableModuleElement}, y::ApplicableModuleElement; kwds...)
    isempty(G) && return iszero(y) ? transpose(P[]) : nothing
    res = singular_lift(G, y)
    if res !== nothing
        res = transpose(res) # opposite convention for matrix multiplication in Singular compared to us
    end
    return res
end


import PolynomialRings
function enable()
    PolynomialRings.Backends.Gröbner.set_default(SingularExpect())
end

end
