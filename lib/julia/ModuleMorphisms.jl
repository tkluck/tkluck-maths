module ModuleMorphisms

import PolynomialRings


import Base: zero
import PolynomialRings.Modules: modulebasering

struct FreeModule{ArrayType,Size}
end

function FreeModule(P::Type{<:Polynomial}, size...; sparse=false)
    if sparse
        return FreeModule{Array{P},size}
    elseif length(size) == 1
        return FreeModule{SparseVector{P,Int}, size}
    elseif length(size) == 2
        return FreeModule{SparseMatrixCSC{P,Int}, size}
    else
        throw(ValueError("Cannot create sparse free module with size $size"))
    end
end

arraytype(::Type{<:FreeModule{ArrayType}}) where ArrayType = ArrayType
modulebasering(f::Type{<:FreeModule}) = eltype(arraytype(f))
rank(::Type{<:FreeModule{ArrayType,Size}}) = prod(Size)

zero(f::Type{<:FreeModule{ArrayType,Size}}) = (res = ArrayType(Size...); res .= zero(modulebasering(f)); res)
function basis(f::Type{<:FreeModule{ArrayType,Size}})
    z0 = zero(f)
    res = ArrayType[]
    for i in eachindex(z)
        b = zero(f)
        b[i] = one(modulebasering(f))
        push!(res, b)
    end
    return res
end

export FreeModule, rank, basis

mutable struct Submodule{F,A,P}
    generators::Vector{A}
    _groebner_basis::Nullable{Vector{A}}
    _groebner_transformation::Nullable{Matrix{P}}
end

mutable struct ModuleMorphism{Domain,Codomain,A}
    f::Function
    image::Submodule{Codomain}
    _kernel::Nullable{Submodule{Domain}}
end


end
