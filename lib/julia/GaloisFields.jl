module GaloisFields

import Base: zero, one, +,-,*,//,inv
import Base: show
import Base: convert, promote_rule, promote_type

import PolynomialRings: fraction_field

struct Reduced end

struct GaloisField{I,p} <: Number
    n::I
    GaloisField{I,p}(n::I) where {I,p} = new(rem(n,p))
    GaloisField{I,p}(::Reduced, n::I) where {I,p} = new(n)
end

const GF{I,p} = GaloisField{I,p}

GaloisField(p::Integer) = GaloisField{typeof(p), p}



char(x) = char(typeof(x))
char(::Type{GF{I,p}}) where {I,p} = p


zero(x::Type{<:GF}) = x(0)
one(x::Type{<:GF}) = x(1)
+(a::F,b::F) where F<:GF = F(a.n+b.n)
-(a::F,b::F) where F<:GF = F(a.n-b.n)
+(a::F) where F<:GF = a
-(a::F) where F<:GF = F(Reduced(), char(F)-a.n)
*(a::F,b::F) where F<:GF = F(a.n*b.n)
inv(a::F) where F<:GF = F(Reduced(), invmod(a.n,char(F)))
//(a::F,b::F) where F<:GF = a*inv(b)

show(io::IO, a::GF) = show(io, a.n)
show(io::IO, ::Type{GF{I,p}}) where {I,p} = write(io, "𝔽$p")


promote_rule(F::Type{<:GF}, ::Type{<:Integer}) = F
convert(F::Type{GF{I,p}}, i::Integer) where {I,p} = F(Reduced(), rem(i,p))

convert(::Type{F}, n::F) where F<:GF = n
(::Type{F})(n::F) where F<:GF = n

fraction_field(F::Type{<:GF}) = F


export GaloisField

end
