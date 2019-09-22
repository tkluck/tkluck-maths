using Combinatorics: permutations

using RootSystems

function RootSystems.coefficients_on_simple_roots(Φ::RootSystem, ::Val{:TOP})
    return sum(coefficients_on_simple_roots(Φ, α) for α in positive_roots(Φ))
end

function extend_simple_root_values_to_positive_roots(Φ::RootSystem, simple_root_values)
    @assert length(simple_root_values) == rank(Φ)

    map(enumerate(positive_roots(Φ))) do I
        i, α = I
        sum(c * v for (c, v) in zip(coefficients_on_simple_roots(Φ, α), simple_root_values))
    end
end

@doc """
Try to find an additive assignment of positive values to the positive
roots in `root_system`, which have more than the expected number (2 in
most cases) of ways of assigning these numbers.
"""
function find_indistinguishable_configurations(Φ::RootSystem, searchrange)
    expected_number = dynkin_diagram_automorphism_count(Φ)
    top_coeffs = coefficients_on_simple_roots(Φ, Val(:TOP))

    #for simple_values in ((1, 1, 4, 2, 3),)
    #for simple_values in ((5, 4, 3, 11, 2),)
    #for simple_values in ((7, 10, 9, 3, 2),)
    #for simple_values in ((1, 1, 1, 4, 3, 3),)
    #for simple_values in ((3, 2, 1, 3, 4, 1),)
    while true
        simple_values = ntuple(_ -> rand(searchrange), Val(rank(Φ)))

        all_values = extend_simple_root_values_to_positive_roots(Φ, simple_values)

        max_possible_simple_value = begin
            top_coeffs = coefficients_on_simple_roots(Φ, Val(:TOP))

            sort!(top_coeffs)
            v = sort(all_values)

            min_coeff, coeffs = top_coeffs[1], top_coeffs[end:-1:2]

            min_total_rest = sum(a*b for (a,b) in zip(coeffs, @view v[1:length(coeffs)]))

            (sum(all_values) - min_total_rest)//min_coeff
        end

        subset_that_can_be_assigned_to_simple_root = filter(v -> v <= max_possible_simple_value, all_values)

        possible_assignments = Set()
        sorted_all_values = sort(all_values)
        total_of_values = sum(all_values)

        for other_assignment in permutations(subset_that_can_be_assigned_to_simple_root, rank(Φ))
            sum(i -> other_assignment[i]*top_coeffs[i], eachindex(top_coeffs)) == total_of_values || continue

            other_all_values = extend_simple_root_values_to_positive_roots(Φ, other_assignment)

            if sorted_all_values == sort(other_all_values)
                push!(possible_assignments, other_assignment)
            end
        end

        if length(possible_assignments) > expected_number
            @show possible_assignments
        end
    end
end
