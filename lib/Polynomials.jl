module Polynomials

export groebner_basis, syzygies, minimal_groebner_basis


import Base: +,==,*,//,-,convert,promote_rule,show,cmp,isless,zero

type Exponent{NumVars}
    e::NTuple{NumVars, Int}
end

function =={E <: Exponent}(a::E, b::E)
    return a.e == b.e
end

function +{NumVars}(a::Exponent{NumVars}, b::Exponent{NumVars})
    return Exponent{NumVars}(
        ntuple(Val{NumVars}) do i
            return a.e[i] + b.e[i]
        end
    )
end

function cmp{E <: Exponent}(a::E, b::E)
    # degrevlex
    if(sum(a.e) == sum(b.e))
        return -cmp(a.e, b.e)
    else
        return cmp(sum(a.e), sum(b.e))
    end
end

isless{E <: Exponent}(a::E, b::E) = cmp(a, b)<0

typealias _Monomial{R<:Number, NumVars} Tuple{Exponent{NumVars}, R}

_Monomial{R<:Number,NumVars}(a::Exponent{NumVars}, b::R)= _Monomial((a,b))

coefficient{R,NumVars}(a::_Monomial{R, NumVars})::R = a[2]

function //{R,NumVars,S}(a::_Monomial{R, NumVars}, b::S)
    exponent, coeff = a
    T = promote_type(R,S)
    return _Monomial(exponent, coeff // b)
end

immutable _Polynomial{R<:Number, NumVars} <: Number
    coeffs::Vector{ _Monomial{R, NumVars} }
end

num_variables{R<:Number, NumVars}(::Type{_Polynomial{R,NumVars}}) = NumVars
variables{R<:Number, NumVars}(::Type{_Polynomial{R,NumVars}}) = ntuple(Val{NumVars}) do i
    exponent = ntuple(Val{NumVars}) do j
        i == j ? 1 : 0
    end
    _Polynomial([_Monomial(Exponent(exponent), one(R))])
end

num_variables{R<:Number, NumVars}(p::_Polynomial{R,NumVars}) = num_variables(typeof(p))
variables{R<:Number, NumVars}(p::_Polynomial{R,NumVars}) = variables(typeof(p))

convert{R<:Number, NumVars}(::Type{_Polynomial{R, NumVars}}, c0::R) = (
    c0 != 0
       ? _Polynomial{R, NumVars}([_Monomial(Exponent(ntuple(i->0, Val{NumVars})), c0)])
       : _Polynomial{R, NumVars}([])
)

promote_rule{R<:Number,NumVars,S<:Number}(::Type{_Polynomial{R, NumVars}}, ::Type{S}) =
    _Polynomial{promote_type(R,S), NumVars}

convert{R<:Number, NumVars}(::Type{_Polynomial{R, NumVars}}, c::_Monomial{R,NumVars}) = (
    coefficient(c) != 0
       ?  _Polynomial{R, NumVars}([c])
       : _Polynomial{R, NumVars}([])
)

promote_rule{R<:Number,NumVars}(::Type{_Polynomial{R, NumVars}}, ::Type{_Monomial{R, NumVars}}) =
    _Polynomial{R, NumVars}

+{R,NumVars}(a::_Monomial{R,NumVars},b::_Polynomial{R,NumVars})=+(promote(a,b)...)
*{R,NumVars}(a::_Monomial{R,NumVars},b::_Polynomial{R,NumVars})=*(promote(a,b)...)
-{R,NumVars}(a::_Monomial{R,NumVars},b::_Polynomial{R,NumVars})=-(promote(a,b)...)
+{R,NumVars}(a::_Polynomial{R,NumVars},b::_Monomial{R,NumVars})=+(promote(a,b)...)
*{R,NumVars}(a::_Polynomial{R,NumVars},b::_Monomial{R,NumVars})=*(promote(a,b)...)
-{R,NumVars}(a::_Polynomial{R,NumVars},b::_Monomial{R,NumVars})=-(promote(a,b)...)

iszero{P <: _Polynomial}(p::P)= length(p.coeffs) == 0
iszero{P <: _Polynomial}(p::Vector{P}) = all(iszero(p_i) for p_i in p)
iszero{P <: _Polynomial}(p::Matrix{P}) = all(iszero(p_i) for p_i in p)

function +{R, T, NumVars}(a::_Polynomial{R, NumVars}, b::_Polynomial{T, NumVars})
    S = promote_type(R, T)
    res = Vector{ _Monomial{R, NumVars} }()

    state_a = start(a.coeffs)
    state_b = start(b.coeffs)
    while !done(a.coeffs, state_a) && !done(b.coeffs, state_b)
        ((exponent_a, coefficient_a), next_state_a) = next(a.coeffs, state_a)
        ((exponent_b, coefficient_b), next_state_b) = next(b.coeffs, state_b)

        if exponent_a < exponent_b
            append!(res, [_Monomial(exponent_a, coefficient_a)])
            state_a = next_state_a
        elseif exponent_b < exponent_a
            append!(res, [_Monomial(exponent_b, coefficient_b)])
            state_b = next_state_b
        else
            coeff = coefficient_a + coefficient_b
            if coeff != 0
                append!(res, [_Monomial(exponent_a, coeff)])
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
    res = Vector{ _Monomial{R, NumVars} }()

    # the following seems to be implemented through a very naive version
    # of append! that does a reallocation at every step. So implement
    # it manually below
    #summands = [
    #    _Monomial(exp_a + exp_b, coeff_a * coeff_b)
    #    for (exp_a, coeff_a) in a.coeffs for (exp_b, coeff_b) in b.coeffs
    #]
    summands = Vector{_Monomial{S,NumVars}}(length(a.coeffs) * length(b.coeffs))
    ix = 1
    for (exp_a, coeff_a) in a.coeffs
        for (exp_b, coeff_b) in b.coeffs
            summands[ix] = _Monomial(exp_a + exp_b, coeff_a * coeff_b)
            ix += 1
        end
    end
    assert( ix == length(summands)+1)
    sort!(summands)

    last_exp = Union{}
    for (exponent, coef) in summands
        if exponent == last_exp
            _, cur_coef = res[end]
            res[end] = _Monomial(exponent, cur_coef + coef)
        else
            append!(res, [ _Monomial(exponent, coef) ])
            last_exp = exponent
        end
    end
    res = [m for m in res if coefficient(m) != 0]
    return _Polynomial{S, NumVars}(res)
end

function -{P <: _Polynomial}(f::P)
    return P([_Monomial(exponent, -coeff) for (exponent, coeff) in f.coeffs])
end

function -{P <: _Polynomial}(a::P, b::P)
    return a + -b
end

function =={P <: _Polynomial}(a::P, b::P)
    return a.coeffs == b.coeffs
end

function leading_term{P <: _Polynomial}(p::P)
    if length(p.coeffs) > 0
        return p.coeffs[end]
    else
        throw(ArgumentError("The zero polynomial $( p ) does not have a leading term"))
    end
end

typealias _ModuleElement{P <: _Polynomial} Vector{P}
typealias _HomModuleElement{P <: _Polynomial} Matrix{P}
typealias _AbstractModuleElement{P <: _Polynomial} Union{P, _ModuleElement{P}, _HomModuleElement{P}}
typealias _AbstractModuleElementVector{P <: _Polynomial} Union{AbstractVector{P}, AbstractVector{_ModuleElement{P}}, AbstractVector{_HomModuleElement{P}}}

zero{P <: _Polynomial}(a::AbstractVector{Vector{P}}) = [[0 for _ in a_i] for a_i in a]

type _ModuleMonomial{M <: _Monomial}
    m::M
    pos::Int
end

*{P<:_Polynomial}(a::P, x::_ModuleElement{P})= P[ a*x_i for x_i in x ]
*{R<:Number, NumVars}(x::_ModuleElement{_Polynomial{R, NumVars}}, a::_Monomial{R, NumVars})= convert(_Polynomial{R,NumVars}, a) * x
*{R<:Number, NumVars}(a::_Monomial{R, NumVars}, x::_ModuleElement{_Polynomial{R, NumVars}})= convert(_Polynomial{R,NumVars}, a) * x

*{P<:_Polynomial}(a::P, x::_HomModuleElement{P})= P[ a*x_i for x_i in x ]
*{R<:Number, NumVars}(x::_HomModuleElement{_Polynomial{R, NumVars}}, a::_Monomial{R, NumVars})= convert(_Polynomial{R,NumVars}, a) * x
*{R<:Number, NumVars}(a::_Monomial{R, NumVars}, x::_HomModuleElement{_Polynomial{R, NumVars}})= convert(_Polynomial{R,NumVars}, a) * x

function leading_term{P<:_Polynomial}(a::_ModuleElement{P})
    for (i, f_i) in enumerate(a)
        if !iszero(f_i)
            return _ModuleMonomial(leading_term(f_i), i)
        end
    end
    throw(ArgumentError("The zero element $( a ) does not have a leading term"))
end

function leading_term{P<:_Polynomial}(a::_HomModuleElement{P})
    for (i, f_i) in enumerate(a)
        if !iszero(f_i)
            return _ModuleMonomial(leading_term(f_i), i)
        end
    end
    throw(ArgumentError("The zero element $( a ) does not have a leading term"))
end

function _lcm_multipliers{NumVars}(a::Exponent{NumVars}, b::Exponent{NumVars})
     _lcm = Exponent{NumVars}(
        ntuple(Val{NumVars}) do i
            return max(a.e[i], b.e[i])
        end
    )

    multiplier_a = Exponent{NumVars}(
        ntuple(Val{NumVars}) do i
            return _lcm.e[i] - a.e[i]
        end
    )
    multiplier_b = Exponent{NumVars}(
        ntuple(Val{NumVars}) do i
            return _lcm.e[i] - b.e[i]
        end
    )

    return multiplier_a, multiplier_b
end

function _is_constant{M <: _Monomial}(a::M)
    exponent, coeff = a
    return sum(exponent.e) == 0
end

function _lcm_multipliers{M <: _Monomial}(a::M, b::M)
    exp_a, coeff_a = a
    exp_b, coeff_b = b

    m_exp_a, m_exp_b = _lcm_multipliers(exp_a, exp_b)

    return _Monomial(m_exp_a, coeff_b), _Monomial(m_exp_b, coeff_a)
end

function _monomial_div{M<: _Monomial}(a::M, b::M)::Nullable{M}
    mul_a, mul_b = _lcm_multipliers(a, b)

    if _is_constant(mul_a)
        return mul_b // coefficient(mul_a)
    else
        return nothing
    end
end

function _monomial_div{M<: _Monomial}(a::_ModuleMonomial{M}, b::_ModuleMonomial{M})::Nullable{M}
    if a.pos != b.pos
        return nothing
    else
        return _monomial_div(a.m, b.m)
    end
end

_maybe_lcm_multipliers{M <: _Monomial}(a::M, b::M)::Nullable{Tuple{M,M}} = _lcm_multipliers(a,b)
function _maybe_lcm_multipliers{M <: _Monomial}(a::_ModuleMonomial{M}, b::_ModuleMonomial{M})::Nullable{Tuple{M,M}}
    if a.pos != b.pos
        return nothing
    else
        return _lcm_multipliers(a.m, b.m)
    end
end

function _lead_div_with_remainder{R,NumVars}(f::_Polynomial{R,NumVars}, g::_Polynomial{R,NumVars})::Tuple{Nullable{_Polynomial{R,NumVars}}, _Polynomial{R,NumVars}}
    maybe_factor = _monomial_div(leading_term(f), leading_term(g))

    if isnull(maybe_factor)
        return nothing, f
    else
        factor = get(maybe_factor)
        return factor, f - g*factor
    end
end

_monomials{R, NumVars}(f::_Polynomial{R, NumVars}) = reverse(f.coeffs)

function _monomials{R<:Number, NumVars}(f::_ModuleElement{_Polynomial{R, NumVars}})
    return [
        _ModuleMonomial{_Monomial{R, NumVars}}(m, i)
        for (i, f_i) in enumerate(f)
        for m in _monomials(f_i)
    ]
end

function _monomials{R<:Number, NumVars}(f::_HomModuleElement{_Polynomial{R, NumVars}})
    return [
        _ModuleMonomial{_Monomial{R, NumVars}}(m, i)
        for (i, f_i) in enumerate(f)
        for m in _monomials(f_i)
    ]
end

function _div_with_remainder{P <: _Polynomial}(f::_AbstractModuleElement{P}, g::_AbstractModuleElement{P})::Tuple{Nullable{P}, _AbstractModuleElement{P}}
    if iszero(f)
        return zero(P), f
    else
        for monomial in _monomials(f)
            maybe_factor = _monomial_div(monomial, leading_term(g))
            if !isnull(maybe_factor)
                factor = get(maybe_factor)
                return factor, f - g*factor
            end
        end
        return nothing, f
    end
end

function _reduce{P <: _Polynomial}(f::_AbstractModuleElement{P}, G::_AbstractModuleElementVector{P})
    factors = transpose(zeros(P, length(G)))
    frst = true
    more_loops = false
    f_red = f
    while frst || more_loops
        frst = false
        more_loops = false
        for (i, g) in enumerate(G)
            q, f_red = _div_with_remainder(f_red, g)
            if !isnull(q)
                factors[1, i] += get(q)
                more_loops = true
            end
            if iszero(f_red)
                return f_red, factors
            end
        end
    end

    return f_red, factors
end



function groebner_basis{P <: _Polynomial}(polynomials::_AbstractModuleElementVector{P})

    result = copy(polynomials)
    transformation =[ P[ i==j ? 1 : 0 for i in eachindex(polynomials)] for j in eachindex(polynomials)]


    pairs_to_consider = [
        (i,j) for i in eachindex(result) for j in eachindex(result) if i < j
    ]

    while length(pairs_to_consider) > 0
        (i,j) = pop!(pairs_to_consider)
        a = result[i]
        b = result[j]
        lt_a = leading_term(a)
        lt_b = leading_term(b)

        maybe_multipliers = _maybe_lcm_multipliers(lt_a, lt_b)
        if !isnull(maybe_multipliers)
            m_a, m_b = get(maybe_multipliers)
            S = m_a * a - m_b * b

            # potential speedup: wikipedia says that in all but the 'last steps'
            # (whichever those may be), we can get away with a version of _reduce
            # that only does lead division
            (S_red, factors) = _reduce(S, result)

            factors[1, i] -= m_a
            factors[1, j] += m_b

            if !iszero(S_red)
                new_j = length(result)+1
                append!(pairs_to_consider, [(new_i, new_j) for new_i in eachindex(result)])
                append!(result, [S_red])

                tr = [ sum(-f * transformation[x][y] for (x,f) in enumerate(factors)) for y in eachindex(polynomials) ]
                append!(transformation, [tr])
            end
        end
    end

    flat_tr = [ transformation[x][y] for x=eachindex(result), y=eachindex(polynomials) ]

    return result, flat_tr
end

function minimal_groebner_basis{P <: _Polynomial}(polynomials::_AbstractModuleElementVector{P})

    (basis, transformation) = groebner_basis(polynomials)

    redundant = Set{Int}()
    for i in eachindex(basis)
        #(p_red, factors) = _reduce(basis[i], [b for (j,b) in enumerate(basis) if j!=i && !(j in redundant)])
        if iszero(basis[i])
            push!(redundant, i)
        end
    end

    necessary = [ i for i in eachindex(basis) if !(i in redundant) ]
    return basis[ necessary ], transformation[ necessary, : ]

end

function syzygies{P <: _Polynomial}(polynomials::_AbstractModuleElementVector{P})
    pairs_to_consider = [
        (i,j) for i in eachindex(polynomials) for j in eachindex(polynomials) if i < j
    ]

    result = Vector{_ModuleElement{P}}()

    for (i,j) in pairs_to_consider
        a = polynomials[i]
        b = polynomials[j]
        lt_a = leading_term(a)
        lt_b = leading_term(b)

        maybe_multipliers = _maybe_lcm_multipliers(lt_a, lt_b)
        if !isnull(maybe_multipliers)
            m_a, m_b = get(maybe_multipliers)
            S = m_a * a - m_b * b

            (S_red, syzygy) = _reduce(S, polynomials)
            if !iszero(S_red)
                throw(ArgumentError("syzygies(...) expects a Groebner basis, so S_red = $( S_red ) should be zero"))
            end
            syz_vector= vec(syzygy)
            syz_vector[i] -= m_a
            syz_vector[j] += m_b

            (syz_red, _) = _reduce(syz_vector, result)
            if !iszero(syz_red)
                append!(result, [syz_red])
            end
        end
    end

    (result, _) = minimal_groebner_basis(result)
    flat_result = [ result[x][y] for x=eachindex(result), y=eachindex(polynomials) ]

    return flat_result
end

import Iterators: product

function monomials_not_in_ideal{R <: Number, NumVars}(monomials::Vector{_Polynomial{R, NumVars}})
    vars = variables(_Polynomial{R, NumVars})

    degree_bounds = [0 for _ in vars]
    for monom in monomials
        exponent_tuple = monom.coeffs[1][1].e
        variables_appearing = [(i,e) for (i,e) in enumerate(exponent_tuple) if e != 0]
        if length(variables_appearing) == 1
            (variable_index,variable_exp), = variables_appearing
            if degree_bounds[ variable_index ] == 0
                degree_bounds[ variable_index ] = variable_exp
            elseif variable_exp < degree_bounds[ variable_index ]
                degree_bounds[ variable_index ] = variable_exp
            end
        end
    end
    if any(b == 0 for b in degree_bounds)
        throw(ArgumentError("monomials_not_in_ideal: result is not a finite set"))
    end

    result = eltype(monomials)[]
    for p in product([0:(b-1) for b in degree_bounds]...)
        m = prod(v^p[i] for (i,v) in enumerate(vars))
        if !any(!isnull(_monomial_div(m.coeffs[1], d.coeffs[1])) for d in monomials)
            append!(result, [m])
        end
    end
    return result

end

function show{R}(io::IO, p::_Polynomial{R, 1})
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

function show{R}(io::IO, p::_Polynomial{R, 2})
    frst = true
    if length(p.coeffs) == 0
        show(io, zero(R))
    end
    for (e, c) in p.coeffs
        (i,j) = e.e
        if !frst
            print(io, " + ")
        else
            frst = false
        end
        if c != 1 || i == 0
            print(io, c)
        end
        if i == 1
            print(io, " x")
        elseif i > 1
            print(io, " x^$(i)")
        end
        if j == 1
            print(io, " y")
        elseif j > 1
            print(io, " y^$(j)")
        end
    end
end

function show{R}(io::IO, p::_Polynomial{R, 3})
    frst = true
    if length(p.coeffs) == 0
        show(io, zero(R))
    end
    for (e, c) in p.coeffs
        (i,j,k) = e.e
        if !frst
            print(io, " + ")
        else
            frst = false
        end
        if c != 1 || i == 0
            assert(c != 0)
            print(io, c)
        end
        if i == 1
            print(io, " x")
        elseif i > 1
            print(io, " x^$(i)")
        end
        if j == 1
            print(io, " y")
        elseif j > 1
            print(io, " y^$(j)")
        end
        if k == 1
            print(io, " z")
        elseif k > 1
            print(io, " z^$(k)")
        end
    end
end

x = _Polynomial{Int, 2}([_Monomial(Exponent{2}((1,0)), 1)])
y = _Polynomial{Int, 2}([_Monomial(Exponent{2}((0,1)), 1)])

function random_polynomial()
    res = sum([ rand(0:100) * x^i for i = 0:10 ])
end

function random_matrix()
    return Matrix([ random_polynomial() for i = 1:10, j = 1:10 ])
end

end
