module HomspaceMorphisms

using PolynomialRings

import PolynomialRings: gröbner_basis, gröbner_transformation

export ModuleMorphism, HomspaceMorphism, lift, kernel, lift_and_obstruction

abstract type Morphism{P <: Polynomial} end

type ModuleMorphism{P <: Polynomial} <: Morphism{P}
    m::Matrix{P}
end

type HomspaceMorphism{P <: Polynomial} <: Morphism{P}
    m::Array{P, 4}
    _image_basis::Nullable{AbstractVector}
    _image_basis_transformation::Nullable{AbstractMatrix}
end
HomspaceMorphism(m::Array{P,4}) where P <: Polynomial = HomspaceMorphism(m, Nullable{Vector{Array{P,2}}}(nothing), Nullable{AbstractMatrix{P}}(nothing))

(F::ModuleMorphism{P})(x::Vector{P}) where P <: Polynomial= F.m * x

span(F::ModuleMorphism) = [
    getindex(F.m, :, i)
    for i in indices(F.m, 2)
]

span(F::HomspaceMorphism) = [
    getindex(F.m, i, j, :, :)
    # this 'flattens' the matrix basis into a single index.
    # the corresponding 'unflatten' operation happens elsewhere
    # in this file; search for 'unflatten'.
    for j in indices(F.m, 2) for i in indices(F.m, 1)
]
function unflatten(F::HomspaceMorphism, coefficients)
    k = zeros(eltype(F.m), size(F.m, 1), size(F.m, 2))
    # this is where we 'unflatten' the basis matrices
    for (j,c) in enumerate(coefficients)
        k[j] = c
    end
    return k
end
function lift(F::ModuleMorphism{P}, x::Vector{P})::Nullable{Vector{P}} where P <: Polynomial

    (basis, transformation) = gröbner_transformation(span(F))
    factors, x_red = divrem(x, basis)

    if any(x_i != 0 for x_i in x_red)
        return nothing
    else
        return vec(factors * transformation)
    end
end

function lift_and_obstruction(F::HomspaceMorphism{P}, x::Matrix{P}) where P <: Polynomial

    (basis, transformation) = gröbner_transformation(span(F))
    factors, x_red = divrem(x, basis)

    return unflatten(F, factors * transformation), x_red
end

function lift(F::HomspaceMorphism{P}, x::Matrix{P})::Nullable{Matrix{P}} where P <: Polynomial

    lift, obstruction = lift_and_obstruction(F, x)

    if any(x_i != 0 for x_i in obstruction)
        return nothing
    else
        return lift
    end
end

function lift(F::Vector{P}, x::P)::Nullable{Matrix{P}} where P <: Polynomial

    (basis, transformation) = gröbner_transformation(span(F))
    factors, x_red = divrem(x, basis)

    if x_red != 0
        return nothing
    else
        return factors * transformation
    end

end

function kernel(F::ModuleMorphism{P}) where P <: Polynomial
    (basis, transformation) = gröbner_transformation(span(F))
    S = syzygies(basis)

    coefficients = S * transformation

    return transpose(coefficients)

end

function gröbner_transformation(F::HomspaceMorphism)
    if isnull(F._image_basis)
        basis, transformation = gröbner_transformation(span(F))
        F._image_basis = basis
        F._image_basis_transformation = transformation
    end
    return (get(F._image_basis), get(F._image_basis_transformation))
end

gröbner_basis(F::HomspaceMorphism) = gröbner_transformation(F)[1]

function kernel(F::HomspaceMorphism{P}) where P <: Polynomial
    basis, transformation = gröbner_transformation(F)
    S = syzygies(basis)

    coefficients = S * transformation
    return [
        unflatten(F, coefficients[i, :])
        for i in indices(coefficients, 1)
    ]

end

function cokernel_basis(F::Morphism{P}) where P <: Polynomial
    (basis, transformation) = gröbner_transformation(span(F))
    leading_terms = [leading_term(f) for f in basis]

end


end
