module Polynomials

type Exponent{NumVars}
    e::NTuple{NumVars, Int}
end

function =={NumVars}(a::Exponent{NumVars}, b::Exponent{NumVars})
    return a.e == b.e
end

function +{NumVars}(a::Exponent{NumVars}, b::Exponent{NumVars})
    return Exponent{NumVars}(
        ntuple(Val{NumVars}) do i
            return a.e[i] + b.e[i]
        end
    )
end

function Base.isless{NumVars}(a::Exponent{NumVars}, b::Exponent{NumVars})
    return a.e < b.e
end

immutable _Polynomial{R<:Number, NumVars} <: Number
    coeffs::Vector{ Pair{Exponent{NumVars}, R} }
end

Base.convert{R<:Number, NumVars}(::Type{_Polynomial{R, NumVars}}, c0::R) =
    _Polynomial{R, NumVars}([Pair(Exponent(ntuple(i->0, Val{NumVars})), c0)])

Base.promote_rule{R<:Number,NumVars,S<:Number}(::Type{_Polynomial{R, NumVars}}, ::Type{S}) =
    _Polynomial{promote_type(R,S), NumVars}

function +{R, T, NumVars}(a::_Polynomial{R, NumVars}, b::_Polynomial{T, NumVars})
    S = promote_type(R, T)
    res = Vector{ Pair{Exponent{NumVars}, S} }()

    state_a = start(a.coeffs)
    state_b = start(b.coeffs)
    while !done(a.coeffs, state_a) && !done(b.coeffs, state_b)
        ((exponent_a, coefficient_a), next_state_a) = next(a.coeffs, state_a)
        ((exponent_b, coefficient_b), next_state_b) = next(b.coeffs, state_b)

        if exponent_a < exponent_b
            append!(res, [Pair(exponent_a, coefficient_a)])
            state_a = next_state_a
        elseif exponent_b < exponent_a
            append!(res, [Pair(exponent_b, coefficient_b)])
            state_b = next_state_b
        else
            coeff = coefficient_a + coefficient_b
            if coeff != 0
                append!(res, [Pair(exponent_a, coeff)])
            end
            state_b = next_state_b
            state_a = next_state_a
        end
    end

    append!(res, collect(rest(a.coeffs, state_a)))
    append!(res, collect(rest(b.coeffs, state_b)))

    return _Polynomial{S, NumVars}(res)
end

function  *{R,T,NumVars}(a::_Polynomial{R,NumVars}, b::_Polynomial{T,NumVars})
    S = promote_type(R,T)
    res = Vector{ Pair{Exponent{NumVars}, S} }()

    # this version of Julia doesn't support multiple 'for' in a comprehension
    # so use the 2d syntax and then collect as a workaround. This incurs an
    # extra copy
    summands = collect([
        Pair(exp_a + exp_b, coeff_a * coeff_b)
        for (exp_a, coeff_a) in a.coeffs, (exp_b, coeff_b) in b.coeffs
    ])
    sort!(summands)

    last_exp = Union{}
    for (exponent, coef) in summands
        if exponent == last_exp
            _, cur_coef = res[end]
            res[end] = exponent => cur_coef + coef
        else
            append!(res, [ Pair(exponent, coef) ])
            last_exp = exponent
        end
    end

    return _Polynomial{S, NumVars}(res)
end

#function leading_term{R,NumVars}(p::_Polynomial)
#    return p.coeffs[end]
#end
#
#function groebner_basis{R,NumVars}(polynomials::Vector{_Polynomial{R,NumVars}})
#    transformation = eye(R, length(polynomials))
#
#    exponent, coeff = leading_term(
#
#
#
#
#end

function Base.show{R}(io::IO, p::_Polynomial{R, 1})
    frst = true
    for (e, c) in p.coeffs
        (i,) = e.e
        if !frst
            print(io, " + ")
        else
            frst = false
        end
        print(io, c)
        if i == 1
            print(io, " x")
        elseif i > 1
            print(io, " x^$(i)")
        end
    end
end

x = _Polynomial{Int, 1}([Pair(Exponent{1}((0,)), 0), Pair(Exponent{1}((1,)), 1)])

function random_polynomial()
    res = sum([ rand(1:100) * x^i for i = 1:10 ])
end

function random_matrix()
    return Matrix([ random_polynomial() for i = 1:10, j = 1:10 ])
end

end
