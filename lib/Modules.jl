module Modules

using Polynomials: _Polynomial, _reduce, groebner_basis

export Morphism

type Morphism{P <: _Polynomial, DimDomain, DimCodomain}
    m::Matrix{P}
end


function lift{P <: _Polynomial, DimDomain, DimCodomain}(F::Morphism{P, DimDomain, DimCodomain}, x::Vector{P})::Nullable{Matrix{P}}

    assert(length(x) == DimCodomain)

    F_vectors = [
        getindex(F.m, :, i)
        for i in 1:DimDomain
    ]
    (basis, transformation) = groebner_basis(F_vectors)
    (x_red, factors) = _reduce(x, basis)

    if any(x_i != 0 for x_i in x_red)
        return nothing
    else
        return factors * transformation
    end



end

function lift{P <: _Polynomial}(F::Vector{P}, x::P)::Nullable{Matrix{P}}

    (basis, transformation) = groebner_basis(F)
    (x_red, factors) = _reduce(x, basis)

    if x_red != 0
        return nothing
    else
        return factors * transformation
    end

end


end
