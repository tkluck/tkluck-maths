module Modules

using Polynomials: _Polynomial, _reduce, minimal_groebner_basis, syzygies, groebner_basis

export ModuleMorphism, HomspaceMorphism, lift, kernel

abstract Morphism{P <: _Polynomial}

type ModuleMorphism{P <: _Polynomial} <: Morphism{P}
    m::Matrix{P}
end

type HomspaceMorphism{P <: _Polynomial} <: Morphism{P}
    m::Array{P, 4}
end

span(F::ModuleMorphism) = [
    getindex(F.m, :, i)
    for i in indices(F.m, 2)
]

span(F::HomspaceMorphism) = [
    getindex(F.m, i, j, :, :)
    # this 'flattens' the matrix basis into a single index.
    # the corresponding 'unflatten' operation happens elsewhere
    # in this file; search for 'unflatten'.
    for j in indices(F.m, 2)
    for i in indices(F.m, 1)
]

function lift{P <: _Polynomial}(F::Morphism{P}, x::Vector{P})::Nullable{Matrix{P}}

    (basis, transformation) = minimal_groebner_basis(span(F))
    (x_red, factors) = _reduce(x, basis)

    if any(x_i != 0 for x_i in x_red)
        return nothing
    else
        return factors * transformation
    end



end

function lift{P <: _Polynomial}(F::Vector{P}, x::P)::Nullable{Matrix{P}}

    (basis, transformation) = minimal_groebner_basis(F)
    (x_red, factors) = _reduce(x, basis)

    if x_red != 0
        return nothing
    else
        return factors * transformation
    end

end

function kernel{P <: _Polynomial}(F::ModuleMorphism{P})
    (basis, transformation) = minimal_groebner_basis(span(F))
    S = syzygies(basis)

    coefficients = S * transformation

    return transpose(coefficients)

end

function kernel{P <: _Polynomial}(F::HomspaceMorphism{P})
    (basis, transformation) = minimal_groebner_basis(span(F))
    S = syzygies(basis)

    coefficients = S * transformation
    return [
        begin
            k = zeros(P, size(F.m, 1), size(F.m, 2))
            # this is where we 'unflatten' the basis matrices
            for j in indices(coefficients,2)
                k[j] = coefficients[i,j]
            end
            k
        end
        for i in indices(coefficients, 1)
    ]

end


end
