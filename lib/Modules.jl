module Modules

using Polynomials: _Polynomial, _reduce, groebner_basis, syzygies

export ModuleMorphism, HomspaceMorphism, lift, kernel

abstract Morphism{P <: _Polynomial}

type ModuleMorphism{P <: _Polynomial} <: Morphism{P}
    m::Matrix{P}
end

type HomspaceMorphism{P <: _Polynomial} <: Morphism{P}
    m::Array{P, 4}
end


function lift{P <: _Polynomial}(F::Morphism{P}, x::Vector{P})::Nullable{Matrix{P}}

    F_vectors = [
        getindex(F.m, :, i)
        for i in indices(F.m, 2)
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

function kernel{P <: _Polynomial}(F::Morphism{P})
    F_vectors = [
        getindex(F.m, :, i)
        for i in indices(F.m, 2)
    ]
    (basis, transformation) = groebner_basis(F_vectors)
    S = syzygies(basis)

    return transpose(S * transformation)

end

end
