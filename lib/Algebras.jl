module Algebras

import Base: show

import Polynomials: _Polynomial, _Monomial, Exponent, variables

immutable AlgebraElement{R <: Number, NumVars, T <: Tuple} <: Number
    p::_Polynomial{R,NumVars}
end

immutable Algebra
    basering :: DataType
    variable_names :: Any
    datatype :: DataType
end


function algebra{R <: Number}(::Type{R}, variable_names...)
    T = Tuple{variable_names...}
    numvars = nfields(T)

    vars = variables(_Polynomial{R, numvars})

    datatype = AlgebraElement{R, numvars, T}
    return Algebra(R, T, datatype), [datatype(v) for v in vars]
end

function show{R, NumVars, T <: Tuple}(io::IO, p::AlgebraElement{R, NumVars, T})
    frst = true
    if length(p.p.coeffs) == 0
        print(io, zero(R))
    end
    for (e, c) in p.p.coeffs
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

import Base: +,*,-,==,!=
+{A <: AlgebraElement}(a::A, b::A)  = A(a.p * b.p)
*{A <: AlgebraElement}(a::A, b::A)  = A(a.p * b.p)
-{A <: AlgebraElement}(a::A, b::A)  = A(a.p * b.p)
-{A <: AlgebraElement}(a::A)        = A(-a.p)
=={A <: AlgebraElement}(a::A, b::A) = a.p == b.p
!={A <: AlgebraElement}(a::A, b::A) = a.p != b.p

_symname(s::Symbol)=repr(s)[2:end]
fieldtypes{T <: Tuple}(t::Type{T}) = [_symname(fieldtype(T, i)) for i in 1:nfields(T)]
import Base: convert, promote_rule
function promote_rule{R <: Number, S <: Number, NumVars1, NumVars2, T <: Tuple, U <: Tuple}(
        ::Type{AlgebraElement{R, NumVars1, T}},
        ::Type{AlgebraElement{S, NumVars2, U}},
    )
    RS = promote_type(R,S)
    all_names = Set()
    union!(all_names, fieldtypes(U))
    union!(all_names, fieldtypes(T))
    names = sort(collect(Symbol(s) for s in all_names))
    TU = Tuple{names...}
    NumVars = length(all_names)
    return AlgebraElement{RS, NumVars, TU}
end
convert{R <: Number, NumVars, T <: Tuple}(
        ::Type{AlgebraElement{R, NumVars, T}},
        x::AlgebraElement{R, NumVars, T},
    ) = x


# some black magic: create a function that does the re-indexing from
# the order in U to the order in T. It is actually about as easy to
# do it in this way as it is to write the appropriate data structures.
#
# By memoizing the result, we ensure that we only need to compile the
# function once.
_converter_cache = Dict{Tuple{DataType, DataType}, Function}()
function _converter{T <: Tuple, U <: Tuple}(::Type{T}, ::Type{U})
    if (T,U) in keys(_converter_cache)
        return _converter_cache[T,U]
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
        ::Type{AlgebraElement{R, NumVars1, T}},
        x::AlgebraElement{S, NumVars2, U},
    )

    f = _converter(T, U)
    new_terms = map(x.p.coeffs) do term
        exponent, c = term
        new_exponent = f(exponent.e)
        _Monomial{R, NumVars1}((Exponent(new_exponent), c))
    end

    return AlgebraElement{R, NumVars1, T}(_Polynomial(new_terms))
end


A,(x,y) = algebra(Int64, :x, :y)
B,(z,) = algebra(Int64, :z)

println(typeof(x))
println(typeof(z))

println(x*z)


end
