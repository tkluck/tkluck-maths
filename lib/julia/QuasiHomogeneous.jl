module QuasiHomogeneous

monomials_of_grading(total_grading::I, variable_gradings::NTuple{0, I}) where I <: Integer = Channel(ctype=NTuple{0, I}) do ch
end

monomials_of_grading(total_grading::I, variable_gradings::NTuple{1, I}) where I <: Integer = Channel(ctype=NTuple{1, I}) do ch
    if total_grading % variable_gradings[1] == 0
        push!(ch, ( div(total_grading, variable_gradings[1]), ))
    end
end

monomials_of_grading(total_grading::I, variable_gradings::NTuple{2, I}) where I <: Integer = Channel(ctype=NTuple{2, I}) do ch
    for exp1 in 0:div(total_grading, variable_gradings[1])
        remaining1 = total_grading - variable_gradings[1] * exp1
        if remaining1 % variable_gradings[2] == 0
            push!(ch, (exp1, div(remaining1, variable_gradings[2])))
        end
    end
end

monomials_of_grading(total_grading::I, variable_gradings::NTuple{3, I}) where I <: Integer = Channel(ctype=NTuple{3, I}) do ch
    for exp1 in 0:div(total_grading, variable_gradings[1])
        remaining1 = total_grading - variable_gradings[1] * exp1
        for exp2 in 0:div(remaining1, variable_gradings[2])
            remaining2 = remaining1 - variable_gradings[2] * exp2
            if remaining2 % variable_gradings[3] == 0
                push!(ch, (exp1, exp2, div(remaining2, variable_gradings[3])))
            end
        end
    end
end

monomials_of_grading(total_grading::I, variable_gradings::NTuple{4, I}) where I <: Integer = Channel(ctype=NTuple{4, I}) do ch
    for exp1 in 0:div(total_grading, variable_gradings[1])
        remaining1 = total_grading - variable_gradings[1] * exp1
        for exp2 in 0:div(remaining1, variable_gradings[2])
            remaining2 = remaining1 - variable_gradings[2] * exp2
            for exp3 in 0:div(remaining2, variable_gradings[3])
                remaining3 = remaining2 - variable_gradings[3] * exp3
                if remaining3 % variable_gradings[4] == 0
                    push!(ch, (exp1, exp2, exp3, div(remaining3, variable_gradings[4])))
                end
            end
        end
    end
end

monomials_of_grading(total_grading::I, variable_gradings::NTuple{N, I}) where I <: Integer where N = Channel(ctype=NTuple{N, I}) do ch
    for exp in 0:div(total_grading, variable_gradings[1])
        remaining = total_grading - variable_gradings[1] * exp
        for other_gradings in monomials_of_grading(remaining, tuple(variable_gradings[2:end]...))
            push!(ch, tuple(exp, other_gradings...))
        end
    end
end

using PolynomialRings: construct_monomial, formal_coefficients

function generic_quasihomogeneous_polynomial(total_grading::I, variable_gradings::NTuple{N, I}, R, next_coeff::Function) where I <: Integer where N
    monomials = [construct_monomial(R, e) for e in monomials_of_grading(total_grading, variable_gradings)]
    if(length(monomials) == 0)
        return next_coeff() * zero(R) # hack to make sure it has the same type
    else
        return sum(next_coeff()*m for m in monomials)
    end
end

function generic_quasihomogeneous_map(gradings::Array{I}, variable_gradings::NTuple{N,I}, R, next_coeff::Function) where I <: Integer where N

    return [ generic_quasihomogeneous_polynomial(g, variable_gradings, R, next_coeff) for g in gradings ]

end

degrees_of_homogeneous_map(rank::I, highest_free_generator_grading_source::I, highest_free_generator_grading_target::I, total_grading::I) where I <: Integer = Channel(ctype=Array{I,2}) do ch
    for h in monomials_of_grading(highest_free_generator_grading_source, ntuple(i->1, rank-1))
        for v in monomials_of_grading(highest_free_generator_grading_target, ntuple(i->1, rank-1))
            horizontal_gradings = [0; cumsum(collect(h))]
            vertical_gradings   = [0; cumsum(collect(v))]

            push!(ch, I[ total_grading - i + j for i in vertical_gradings, j in horizontal_gradings])
        end
    end
end

degrees_of_matrix_factorizations(rank::I, highest_free_generator_grading_source::I, highest_free_generator_grading_target::I, total_grading::I) where I <: Integer = Channel(ctype=Array{I,2}) do ch
    for t in [div(total_grading,2)] #0:total_grading
        for h in monomials_of_grading(highest_free_generator_grading_source, ntuple(i->1, rank-1))
            for v in monomials_of_grading(highest_free_generator_grading_target, ntuple(i->1, rank-1))
                horizontal_gradings = [0; cumsum(collect(h))]
                vertical_gradings   = [0; cumsum(collect(v))]

                phi = I[ t - i + j for i in vertical_gradings, j in horizontal_gradings]
                psi = I[ total_grading - t - i + j for i in horizontal_gradings, j in vertical_gradings]
                z = fill(-1, size(phi))
                push!(ch, [ z phi; psi z ])
            end
        end
    end
end

generic_matrix_factorizations(rank::I, highest_free_generator_grading_source::I, highest_free_generator_grading_target::I, total_grading::I, variable_gradings::NTuple{N,I}, R, c::Symbol) where I <: Integer where N = Channel() do ch
    for gradings in degrees_of_matrix_factorizations(rank, highest_free_generator_grading_source, highest_free_generator_grading_target, total_grading)
        F = formal_coefficients(R,c)
        next_coeff() = take!(F)
        push!(ch, generic_quasihomogeneous_map(gradings, variable_gradings, R, next_coeff))
    end
end




end

