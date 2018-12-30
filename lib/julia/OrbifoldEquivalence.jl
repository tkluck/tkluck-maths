module OrbifoldEquivalence

using Base.Iterators
using SparseArrays: spzeros
using LinearAlgebra: det, diagind, diagm, I

using Primes: factor

using PolynomialRings
using PolynomialRings: allvariablesymbols, basering
using PolynomialRings.NamedPolynomials: unused_variable
using QuasiHomogeneous: monomials_of_grading, find_quasihomogeneous_degrees, quasidegree
using QuasiHomogeneous: gradings, generic_quasihomogeneous_map
using GröbnerSingular: singular_modulo
using MatrixFactorizations: columns

function supertrace(Q::AbstractMatrix)
    n,m = size(Q)
    n == m && n%2 ==0 || throw(ArgumentError("Cannot compute supertrace of $n x $m matrix"))

    firsthalf = diagind(Q)[1:end÷2]
    secondhalf = diagind(Q)[end÷2+1:end]
    return sum(i->Q[i], firsthalf) - sum(i->Q[i], secondhalf)
end

function multivariate_residue(g, f, vars...)
    R = base_extend( eltype(f) )
    G,tr = gröbner_transformation(f)

    # TODO: compute that R/G is finite dimensional; otherwise, this computation
    # does not terminate
    M = spzeros(R, length(vars), length(f))
    for (row,v) in enumerate(vars)
        x = convert(R,v)
        factors, x_red = divrem(x,G)
        while !iszero(x_red)
            x = x^2
            factors, x_red = divrem(x,G)
        end
        M[row,:] = factors*tr
    end

    f_transformed = M * f
    g_transformed = g * det(M)

    term, r = divrem(prod(f_transformed), prod(convert(R,v) for v in vars))
    @assert(iszero(r))

    return coefficient(g_transformed, term, vars...)
end

function quantum_dimension(Q::AbstractMatrix, W, left_vars, right_vars)
    g = supertrace(prod( diff(Q, v) for v in left_vars) * prod( diff(Q, v) for v in right_vars) )
    f = [diff(W,v) for v in left_vars]

    return constant_coefficient(multivariate_residue(g, f, left_vars...), right_vars...)
end


enumerate_admissible_gradings(f::Function, N::Integer, W, vars...) = enumerate_admissible_gradings(f, Val{N}, W, vars...)

function enumerate_admissible_gradings(f::Function, ::Type{Val{N}}, W, vars...) where N
    gr = find_quasihomogeneous_degrees(W, vars...)
    vgr = gradings(gr)
    total_grading = quasidegree(W,gr)
    grading = e -> sum(e.*vgr)

    term_exponents = map(i->i[1], expansion(W, vars...))
    n_terms = length(term_exponents)
    n_vars  = length(vars)

    all_exponents = vcat([collect(t) for t in term_exponents]...)

    admissible_rows = Set{NTuple{N,Int}}()

    for ee in Base.Iterators.product([0:k for k in all_exponents]...)
        divisor_gradings = sort([grading(ee[i:i+n_vars-1]) for i=1:n_vars:length(ee)])

        if N == length(divisor_gradings)
            push!(admissible_rows, ntuple(i->divisor_gradings[i], N))
        else
            throw(AssertionError("Not implemented"))
        end
    end

    for row in admissible_rows
        col = ntuple(i -> total_grading - row[i], N)

        if all( (row .+ (col[k] - col[1])) in admissible_rows for k = 2:N)
            m = [ row[j] + (col[k] - col[1]) for j=1:N, k=1:N ]
            n = [ col[k] + (row[j] - row[1]) for j=1:N, k=1:N ]
            z = fill(-1, size(m))
            res = f( [ z m; n z] )
            if res === :break
                break
            elseif res === :continue
                continue
            else
                throw(RuntimeError("Don't know how to continue after $res"))
            end
        end
    end
end

function equivalence_exists(W, Wvars, V, Vvars, rank)

    R,allvars = polynomial_ring(Wvars..., Vvars...)

    vgr = find_quasihomogeneous_degrees(W - V, Wvars..., Vvars...)

    found = false

    enumerate_admissible_gradings(rank, W - V, Wvars..., Vvars...) do gr
        next_coeff = formal_coefficients(R,:c)
        c1 = next_coeff()
        Q = generic_quasihomogeneous_map(gr, vgr, next_coeff)

        C = flat_coefficients(Q^2 - c1^2*(V-W)*one(Q), Wvars..., Vvars...)

        @assert(!(coefficient((-c1^2).terms[1]) in C))

        qdim1 = quantum_dimension(Q,W,Wvars,Vvars)
        qdim2 = quantum_dimension(Q,V,Vvars,Wvars)

        if iszero(qdim1) || iszero(qdim2)
            @info("Found an admissible grading distribution, but its quantum dimension vanishes identically")
            return :continue
        end

        @info("Found a potentially interesting grading distribution: doing the full computation.")

        # to dense monomials
        converted = to_dense_monomials([qdim1; qdim2; C])
        qdim1 = converted[1]
        qdim2 = converted[2]
        C = converted[3:end]


        CC = groebner_basis(C)

        qdim1_red = rem(qdim1, CC)
        qdim2_red = rem(qdim2, CC)

        if !iszero(qdim1_red) && !iszero(qdim2_red)
            @info("Found one!")
            found = true
            return :break
        else
            @info("Unfortunately, for this grading distribution, no solutions exist with a non-vanising quantum dimension.")
        end
    end
    return found ? true : nothing
end

function is_orbifold_equivalent(W, Wvars, V, Vvars, max_rank=Inf)

    for rank = Base.Iterators.countfrom(2)  # FIXME: generic_quasihomogeneous_map breaks on rank=1
        rank > max_rank && break
        @info("Trying rank=$rank")
        if equivalence_exists(W, Wvars, V, Vvars, rank) == true
            return true
        end
    end
    return nothing
end

function variables_appearing(f)
    vars = allvariablesymbols(typeof(f))
    appears = [false for _ in vars]
    for (p,c) in expansion(f, vars...)
        if !iszero(c)
            appears .|= (!iszero).(p)
        end
    end
    return Symbol[v for v in vars[appears]]
end

function check_orbifold_equivalence_for_some_values(Q, W, V)
    W_vars = variables_appearing(W)
    V_vars = variables_appearing(V)
    lqdim = quantum_dimension(Q, W, W_vars, V_vars)
    rqdim = quantum_dimension(Q, V, V_vars, W_vars)
    I = flat_coefficients(Q^2 - (W-V)*one(Q), W_vars..., V_vars...)
    a0 = unused_variable(I)
    G = gröbner_basis([I; 1 - a0*lqdim*rqdim])
    if rem(1, G) == 0
        return false, lqdim, rqdim, G
    else
        return true, lqdim, rqdim, G
    end
end

function divisors(n)
    fs = factor(n)
    ee = collect(values(fs))
    (prod(b^e for (b,e) in zip(keys(fs), exps)) for exps in product(ntuple(i -> 0:ee[i], length(fs))...))
end

function central_charge_preimage(a, b, k)
    I = promote_type(typeof.((a, b, k))...)
    preimage = Set{Tuple{I, I, I}}()
    n = k == 2 ? 1 : 0
    for l in (0, 1, 2)
        m = l == 2 ? 1 : 0
        RHS = (n - a*b)^2 + (a+b-k)*(m*(a+b-k) - (a*b-n)*l)
        if !iszero(RHS)
            for p::I in divisors(RHS)
                q = RHS ÷ p
                c, r1 = divrem(p + a*b - n, a + b - k)
                d, r2 = divrem(q + a*b - n, a + b - k)
                c, d = (p+a*b-n)/(a+b-k), (q+a*b-n)/(a+b-k)
                if iszero(r1) && iszero(r2) && c > 1 && d > 1
                    push!(preimage, (c,d,l))
                end
            end
        end
    end
    return preimage
end

function central_charge_preimage(a, b, k, x, y)
    P = promote_type(typeof(x), typeof(y))
    preimage = Set{P}()
    for (a′, b′, k′) in central_charge_preimage(a, b, k)
        if k′ == 0
            push!(preimage, x^a′ + y^b′)
        elseif k′ == 1
            push!(preimage, x^a′ + x*y^b′)
        elseif k′ == 2
            push!(preimage, y*x^a′ + x*y^b′)
        else
            throw("Unexpected value: $k′")
        end
    end
    return preimage
end

function central_charge(f::Polynomial, vars::Symbol...)
    gradings = find_quasihomogeneous_degrees(f, vars...)
    d = quasidegree(f, gradings)
    scale = 2//d

    return sum(1 - scale*d for (v, d) in gradings)

end

function matrix_factorization_from_resolution(f, F)
    res = base_extend.(zero(F))
    G, tr = gröbner_transformation(columns(F))
    rhs = f * one(F)
    for i in axes(F,2)
        q, r = divrem(rhs[:,i], G)
        @assert iszero(r)
        res[:,i] = transpose(q * tr)
    end
    return [0I F; res 0I]
end

function herzog_sanders_factorization(f::Polynomial, a::Integer, vars::Symbol...)
    R = typeof(f)
    gr = find_quasihomogeneous_degrees(f, vars...)
    N = quasidegree(f, gr)

    a >= N || throw("herzog_sanders_factorization: need `a` greater than quasidegree of f; quasidegree $f is $N")

    I = [monomials_of_grading(a, gr); f]
    Ω = collect(transpose(I))
    while size(Ω)[1] != size(Ω)[2]
        @show size(Ω)
        Z = diagm(1=>[f for _=1:size(Ω,1)])
        I = singular_modulo(Ω, Z)
        I = xrem.(I, f)
        nonzero_cols = [i for i in axes(I, 2) if !iszero(@view I[:, i])]
        Ω = I[:, nonzero_cols]
    end
    @show size(Ω)
    return matrix_factorization_from_resolution(f, Ω)
end

export supertrace, multivariate_residue, quantum_dimension, equivalence_exists, is_orbifold_equivalent
export check_orbifold_equivalence_for_some_values
export unit_matrix_factorization
export central_charge_preimage, central_charge
export herzog_sanders_factorization


end
