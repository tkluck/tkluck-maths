module PolynomialRings

include("PolynomialRings/Polynomials.jl")
include("PolynomialRings/Groebner.jl")
include("PolynomialRings/Modules.jl")

import PolynomialRings.Polynomials: BaseRing, Polynomial, Term, Monomial, variables, _varsymbols, basering

import Iterators: groupby

immutable PolynomialRing
    basering :: DataType
    variable_names :: Any
    datatype :: DataType
end



function polynomial_ring{R <: BaseRing}(::Type{R}, variable_names::Symbol...)
    T = Tuple{variable_names...}
    NumVars = nfields(T)
    datatype = Polynomial{R, NumVars, T}
    return PolynomialRing(R, T, datatype), variables(datatype)
end

function polynomial_ring{A <: AbstractArray}(dummy_element::A, variable_names::Symbol...)
    T = Tuple{variable_names...}
    NumVars = nfields(T)
    R = typeof(Coefficients.ArrayCoefficient(dummy_element))
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
function promote_rule{R <: BaseRing, S <: BaseRing, NumVars1, NumVars2, T <: Tuple, U <: Tuple}(
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
convert{R <: BaseRing, NumVars, T <: Tuple}(
        ::Type{Polynomial{R, NumVars, T}},
        x::Polynomial{R, NumVars, T},
    ) = x

@generated function _convert(::Type{T}, ::Type{U}, exponent_tuple::Tuple) where T <: Tuple where U <: Tuple
    for j in 1:nfields(U)
        if !any(fieldtype(T,i) == fieldtype(U,j) for i in 1:nfields(T))
            throw(ArgumentError("Cannot convert variables $U to variables $T"))
        end
    end
    :( _lossy_convert(T, U, exponent_tuple) )
end

@generated function _lossy_convert(::Type{T}, ::Type{U}, exponent_tuple::Tuple) where T <: Tuple where U <: Tuple
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
    return converter
end

function convert{R <: BaseRing, S <: BaseRing, NumVars1, NumVars2, T <: Tuple, U <: Tuple}(
        ::Type{Polynomial{R, NumVars1, T}},
        x::Polynomial{S, NumVars2, U},
    )

    new_terms = map(x.terms) do term
        exponent, c = term
        new_exponent = _convert(T, U, exponent.e)
        Term(Monomial(new_exponent), c)
    end

    return Polynomial{R, NumVars1, T}(new_terms)
end

expansion{P <: Polynomial}(x::P, vars::Symbol...) = _expansion_impl(x, vars...)

function _expansion_impl{P <: Polynomial}(x::P, vars::Symbol...)
    if length(vars) == 0
        throw(ArgumentError("Need to pass at least one variable for expansion"))
    end
    R = basering(P)
    T = _varsymbols(P)
    other_vars = Symbol[fieldtype(T, i) for i in 1:nfields(T) if !(fieldtype(T,i) in vars)]
    if length(other_vars) == 0
        return [ (P([Term(Monomial(exp.e), one(R))]), coef) for (exp, coef) in x.terms ]
    end

    f = exponent_tuple -> _lossy_convert(Tuple{vars...},       T, exponent_tuple)
    g = exponent_tuple -> _lossy_convert(Tuple{other_vars...}, T, exponent_tuple)

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
    array_of_expansions = [ (w,p,i) for (i, a_i) in enumerate(a) for (w,p) in _expansion_impl(a_i, vars...)]
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

function (p::Polynomial{R, NumVars, T}){R <: BaseRing, NumVars, T <: Tuple}(; kwargs...)
    vars = [k for (k, v) in kwargs]
    values = [v for (k, v) in kwargs]

    remaining_vars = Symbol[fieldtype(T, i) for i in 1:nfields(T) if !(fieldtype(T,i) in vars)]

    f = exponent_tuple -> _lossy_convert(Tuple{vars...}, T, exponent_tuple)

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
        g = exponent_tuple -> _lossy_convert(Tuple{remaining_vars...}, T, exponent_tuple)
        for term in p.terms
            m, c = term
            res += prod(v^e_i for (v, e_i) in zip(values, f(m.e))) * P([ Term(Monomial(g(m.e)), c) ])
        end
        return res
    end

end

function (p::Array{Polynomial{R, NumVars, T}}){R <: BaseRing, NumVars, T <: Tuple}(; kwargs...)

    return [p_i(;kwargs...) for p_i in p]

end

function random_term(total_degree::Int)
    A, (x,y,z) = polynomial_ring(Int, :x, :y, :z)

    a = rand(0:total_degree)
    b = rand(0:(total_degree-a))
    c = rand(0:(total_degree-a-b))

    return rand(0:100) * x^a*y^b*z^c

end
function random_polynomial()
    res = sum([ random_term(i) for i = 0:10 ])
end

function random_matrix()
    return Matrix([ random_polynomial() for i = 1:10, j = 1:10 ])
end

include("PolynomialRings/Coefficients.jl")


end
