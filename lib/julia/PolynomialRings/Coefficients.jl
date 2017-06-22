module Coefficients

import PolynomialRings.Polynomials: AbstractBaseRing, Polynomial, Term, Monomial, num_variables

immutable ArrayCoefficient{ArrayType <: AbstractArray, Dims <: Val} <: AbstractBaseRing
    c::ArrayType
    ArrayCoefficient{ArrayType,Dims}(c::ArrayType) where {ArrayType,Dims} = new(c)
end
ArrayCoefficient{ArrayType <: AbstractArray}(c::ArrayType) = ArrayCoefficient{ArrayType, Val{size(c)}}(c)

import Base: size, eltype
size{ArrayType <: AbstractArray, dims}(::Type{ArrayCoefficient{ArrayType, Val{dims}}}) = dims
eltype{ArrayType <: AbstractArray, Dims <: Val}(::Type{ArrayCoefficient{ArrayType, Dims}}) = eltype(ArrayType)
matrixorder{ArrayType <: AbstractArray, Dims}(::Type{ArrayCoefficient{ArrayType, Dims}}) = matrixorder(Dims)
matrixorder{dims}(::Type{Val{dims}}) = (assert(length(dims) == 2 && dims[1] == dims[2]); dims[1])


import Base: zero, one, +, -, *, ==, !=
zero{A <: ArrayCoefficient}(::Type{A}) = A(zeros(eltype(A), size(A)))
one{A <: ArrayCoefficient}(::Type{A}) = A(eye(eltype(A), matrixorder(A)))
+{A <: ArrayCoefficient}(a::A, b::A) = A(a.c + b.c)
-{A <: ArrayCoefficient}(a::A, b::A) = A(a.c - b.c)
-{A <: ArrayCoefficient}(a::A) = A(-a.c)
*{A <: ArrayCoefficient}(a::A, b::A) = A(a.c * b.c)
=={A <: ArrayCoefficient}(a::A, b::A) = a.c == b.c
!={A <: ArrayCoefficient}(a::A, b::A) = a.c != b.c

*{A <: ArrayCoefficient}(a::A, b::eltype(A))= A(a.c * b)
*{A <: ArrayCoefficient}(a::eltype(A), b::A)= A(a * b.c)

import Base: promote_rule
_promote{T,N,R}(::Type{Array{T,N}}, ::Type{R}) = Array{promote_type(T,R),N}
_promote{T,R}(::Type{SparseVector{T}}, ::Type{R}) = SparseVector{promote_type(T,R)}
_promote{T,R}(::Type{SparseMatrixCSC{T}}, ::Type{R}) = SparseMatrixCSC{promote_type(T,R)}
promote_rule{R, ArrayType <: AbstractArray, Dims <: Val}(::Type{R}, ::Type{ArrayCoefficient{ArrayType, Dims}}) = ArrayCoefficient{_promote(ArrayType, R), Dims}

function +{A<:ArrayCoefficient}(a::A, b::Polynomial{A})
    Term(Monomial(ntuple(i->0, Val{NumVars})), a) + b
end
function +{A<:ArrayCoefficient}(b::Polynomial{A}, a::A)
    b + Term(Monomial(ntuple(i->0, Val{NumVars})), a)
end
*{A<:ArrayCoefficient}(a::A, b::Polynomial{A}) = (
    Polynomial([ Term(exp, a*coeff) for (exp,coeff) in b.terms ])
)
*{A<:ArrayCoefficient}(a::Polynomial{A}, b::A) = (
    Polynomial{R,NumVars,T}([ Term(exp, coeff*b) for (exp,coeff) in a.terms ])
)
function -{A<:ArrayCoefficient}(a::A, b::Polynomial{A})
    Term(Monomial(ntuple(i->0, Val{NumVars})), a) + -b
end
function -{A<:ArrayCoefficient}(b::Polynomial{A}, a::A)
    b + Term(Monomial(ntuple(i->0, Val{NumVars})), -a)
end


function *{P <: Polynomial, N, Dims}(a::Polynomial{ArrayCoefficient{Array{P, N}, Dims}}, b::P)
    typeof(a)([Term(exponent, coeff*b) for (exponent, coeff) in a.terms])
end
function *{P <: Polynomial, N, Dims}(a::P, b::Polynomial{ArrayCoefficient{Array{P, N}, Dims}})
    typeof(b)([Term(exponent, a*coeff) for (exponent, coeff) in b.terms])
end

function +{P <: Polynomial, N, Dims}(a::Polynomial{ArrayCoefficient{Array{P, N}, Dims}}, b::Array{P,N})
    P2 = typeof(a)
    return a+P2([Term(Monomial(ntuple(n->0,Val{num_variables(a)})), ArrayCoefficient(b))])
end
function +{P <: Polynomial, N, Dims}(a::Array{P,N}, b::Polynomial{ArrayCoefficient{Array{P, N}, Dims}})
    P2 = typeof(b)
    return P2([Term(Monomial(ntuple(n->0,Val{num_variables(b)})), ArrayCoefficient(a))])+b
end

function *{P <: Polynomial, Dims}(a::Polynomial{ArrayCoefficient{Matrix{P}, Dims}}, b::Matrix{P})
    P2 = typeof(a)
    P2([Term(exponent, coeff*b) for (exponent, coeff) in a.terms])
end
function *{P <: Polynomial, Dims}(a::Matrix{P}, b::Polynomial{ArrayCoefficient{Matrix{P}, Dims}})
    P2 = typeof(b)
    P2([Term(exponent, a*coeff) for (exponent, coeff) in b.terms])
end

function -{P <: Polynomial, N, Dims}(a::Polynomial{ArrayCoefficient{Array{P, N}, Dims}}, b::Array{P,N})
    return a+ -b
end
function -{P <: Polynomial, N, Dims}(a::Array{P,N}, b::Polynomial{ArrayCoefficient{Array{P, N}, Dims}})
    return a+ -b
end



import PolynomialRings: expansion, _expansion_impl
function expansion{A <: ArrayCoefficient}(x::Polynomial{A}, vars::Symbol...)
    return [(w, p.c) for (w, p) in _expansion_impl(x, vars...)]
end

end
