module Modules

using Polynomials: _Polynomial, _reduce, minimal_groebner_basis, syzygies, groebner_basis

export ModuleMorphism, HomspaceMorphism, lift, kernel, lift_and_obstruction

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
function unflatten(F::HomspaceMorphism, coefficients)
    k = zeros(eltype(F.m), size(F.m, 1), size(F.m, 2))
    # this is where we 'unflatten' the basis matrices
    for (j,c) in enumerate(coefficients)
        k[j] = c
    end
    return k
end
function lift{P <: _Polynomial}(F::ModuleMorphism{P}, x::Vector{P})::Nullable{Vector{P}}

    (basis, transformation) = minimal_groebner_basis(span(F))
    (x_red, factors) = _reduce(x, basis)

    if any(x_i != 0 for x_i in x_red)
        return nothing
    else
        return vec(factors * transformation)
    end
end

function lift_and_obstruction{P <: _Polynomial}(F::HomspaceMorphism{P}, x::Matrix{P})

    (basis, transformation) = minimal_groebner_basis(span(F))
    (x_red, factors) = _reduce(x, basis)

    return unflatten(F, factors * transformation), x_red
end

function lift{P <: _Polynomial}(F::HomspaceMorphism{P}, x::Matrix{P})::Nullable{Matrix{P}}

    lift, obstruction = lift_and_obstruction(F, x)

    if any(x_i != 0 for x_i in obstruction)
        return nothing
    else
        return lift
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
        unflatten(F, coefficients[i, :])
        for i in indices(coefficients, 1)
    ]

end


end
