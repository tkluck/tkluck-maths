module PolynomialRings

include("PolynomialRings/Polynomials.jl")
include("PolynomialRings/Groebner.jl")
include("PolynomialRings/Modules.jl")

import PolynomialRings.Polynomials: Polynomial, Term, Monomial, variables, _varsymbols, basering

import Iterators: groupby

immutable PolynomialRing
    basering :: DataType
    variable_names :: Any
    datatype :: DataType
end


function polynomial_ring{R <: Number}(::Type{R}, variable_names::Symbol...)
    T = Tuple{variable_names...}
    NumVars = nfields(T)
    datatype = Polynomial{R, NumVars, T}
    return PolynomialRing(R, T, datatype), variables(datatype)
end

import Base: show
function show{R, NumVars, T <: Tuple}(io::IO, p::Polynomial{R, NumVars, T})
    frst = true
    if length(p.terms) == 0
        print(io, zero(R))
    end
    for (e, c) in p.terms
        if !frst
            print(io, " + ")
        else
            frst = false
        end
        print(io, c)
        for (ix, i) in enumerate(e.e)
            varname = repr(fieldtype(T, ix))[2:end]
            if i == 1
                print(io, " $varname")
            elseif i > 1
                print(io, " $varname^$i")
            end
        end
    end
end

import Base: promote_rule, promote_type, convert

fieldtypes{T <: Tuple}(t::Type{T}) = Symbol[fieldtype(T, i) for i in 1:nfields(T)]
import Base: convert, promote_rule
function promote_rule{R <: Number, S <: Number, NumVars1, NumVars2, T <: Tuple, U <: Tuple}(
        ::Type{Polynomial{R, NumVars1, T}},
        ::Type{Polynomial{S, NumVars2, U}},
    )
    RS = promote_type(R,S)
    all_names = Set()
    union!(all_names, fieldtypes(U))
    union!(all_names, fieldtypes(T))
    names = sort(collect(Symbol(s) for s in all_names))
    TU = Tuple{names...}
    NumVars = length(all_names)
    return Polynomial{RS, NumVars, TU}
end
convert{R <: Number, NumVars, T <: Tuple}(
        ::Type{Polynomial{R, NumVars, T}},
        x::Polynomial{R, NumVars, T},
    ) = x


# some black magic: create a function that does the re-indexing from
# the order in U to the order in T. It is actually about as easy to
# do it in this way as it is to write the appropriate data structures.
#
# By memoizing the result, we ensure that we only need to compile the
# function once.
_converter_cache = Dict{Tuple{DataType, DataType}, Function}()
function _converter{T <: Tuple, U <: Tuple}(::Type{T}, ::Type{U}, safe::Bool)
    if (T,U) in keys(_converter_cache)
        return _converter_cache[T,U]
    end
    # ensure that we do not throw away data
    if safe
        for j in 1:nfields(U)
            if !any(fieldtype(T,i) == fieldtype(U,j) for i in 1:nfields(T))
                throw(ArgumentError("Cannot convert variables $U to variables $T"))
            end
        end
    end
    # create an expression that calls the tuple constructor. No arguments -- so far
    converter = :( tuple() )
    for i in 1:nfields(T)
        # for every result field, add the constant 0 as an argument
        push!(converter.args, 0)
        for j in 1:nfields(U)
            if fieldtype(T, i) == fieldtype(U,j)
                # HOWEVER, if it actually also exists in U, then replace the 0
                # by reading from exponent_tuple
                converter.args[end]= :( exponent_tuple[$j] )
                break
            end
        end
    end
    # now specify that exponent_tuple is an argument
    converter = :( exponent_tuple -> $converter )
    # and make the result into a callable Function
    f = eval(converter)
    _converter_cache[T,U] = f
    return f
end

function convert{R <: Number, S <: Number, NumVars1, NumVars2, T <: Tuple, U <: Tuple}(
        ::Type{Polynomial{R, NumVars1, T}},
        x::Polynomial{S, NumVars2, U},
    )

    f = _converter(T, U, true)
    new_terms = map(x.terms) do term
        exponent, c = term
        new_exponent = f(exponent.e)
        Term(Monomial(new_exponent), c)
    end

    return Polynomial{R, NumVars1, T}(new_terms)
end

function expansion{P <: Polynomial}(x::P, vars::Symbol...)
    if length(vars) == 0
        throw(ArgumentError("Need to pass at least one variable for expansion"))
    end
    R = basering(P)
    T = _varsymbols(P)
    other_vars = Symbol[fieldtype(T, i) for i in 1:nfields(T) if !(fieldtype(T,i) in vars)]
    if length(other_vars) == 0
        return [ (Polynomial([Term(Monomial(exp.e), one(R))]), coef) for (exp, coef) in x.terms ]
    end

    f = _converter(Tuple{vars...},       T, false)
    g = _converter(Tuple{other_vars...}, T, false)

    res = []
    separated_terms = [(f(exp.e), g(exp.e), coeff) for (exp, coeff) in x.terms]
    sort!(separated_terms, by=x->x[1])
    for term_group in groupby(x->x[1], separated_terms)
        w_term = Term(Monomial(term_group[1][1]), one(R))
        p_terms = [Term(Monomial(t[2]), t[3]) for t in term_group]
        sort!(p_terms, by=t->t[1])
        w = Polynomial{R, length(vars),       Tuple{vars...}      }([ w_term ])
        p = Polynomial{R, length(other_vars), Tuple{other_vars...}}(p_terms)

        push!(res, (w, p))
    end

    return res

end

function expansion{P <: Polynomial}(a::AbstractArray{P}, vars::Symbol...)
    array_of_expansions = [ (w,p,i) for (i, a_i) in enumerate(a) for (w,p) in expansion(a_i, vars...)]
    sort!(array_of_expansions, by=a->a[1].terms[1][1])

    res = []
    for group in groupby(x -> x[1], array_of_expansions)
        r = zeros(typeof(group[1][2]), size(a))
        for (w,p,i) in group
            r[i] = p
        end
        push!(res, (group[1][1], r))
    end

    return res
end

function (p::Polynomial{R, NumVars, T}){R <: Number, NumVars, T <: Tuple}(; kwargs...)
    vars = [k for (k, v) in kwargs]
    values = [v for (k, v) in kwargs]

    remaining_vars = Symbol[fieldtype(T, i) for i in 1:nfields(T) if !(fieldtype(T,i) in vars)]

    f = _converter(Tuple{vars...}, T, false)

    if length(remaining_vars) == 0
        res = zero(R)
        for term in p.terms
            m, c = term
            res += c * prod(v^e_i for (v, e_i) in zip(values, f(m.e)))
        end
        return res
    else
        P = Polynomial{R, length(remaining_vars), Tuple{remaining_vars...}}
        res = zero(P)
        g = _converter(Tuple{remaining_vars...}, T, false)
        for term in p.terms
            m, c = term
            res += prod(v^e_i for (v, e_i) in zip(values, f(m.e))) * P([ Term(Monomial(g(m.e)), c) ])
        end
        return res
    end

end

function (p::Array{Polynomial{R, NumVars, T}}){R <: Number, NumVars, T <: Tuple}(; kwargs...)

    return [p_i(;kwargs...) for p_i in p]

end


end
