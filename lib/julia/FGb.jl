module FGb

export FGb_with, groebner

using PolynomialRings
using PolynomialRings: construct_monomial
using PolynomialRings.Polynomials: terms
using PolynomialRings.Terms: coefficient, monomial
using PolynomialRings.NamedPolynomials: names

enter_INT() = ccall((:FGb_int_enter_INT, :libfgb), Void, ())
init_urgent(bytes_coeff, bytes_exponent, cargo_cult_drldrl, cargo_cult_100000, cargo_cult_0) = ccall((:FGb_int_init_urgent, :libfgb), Void, (UInt32,UInt32,Cstring,UInt32,UInt32), bytes_coeff, bytes_exponent, cargo_cult_drldrl, cargo_cult_100000, cargo_cult_0)
log_output = Ref{UInt}(0)
init(cargo_cult_1, cargo_cult_1_again, cargo_cult_0, logging_file_handle) = ccall((:FGb_int_init, :libfgb), Void, (UInt32,UInt32,UInt32,Ref{UInt}), cargo_cult_1, cargo_cult_1_again, cargo_cult_0, logging_file_handle)
reset_coeffs(cargo_cult_1) = ccall((:FGb_int_reset_coeffs, :libfgb), Void, (UInt32,), cargo_cult_1)
reset_expos(nvars1, nvars2, varnames) = ccall((:FGb_int_reset_expos, :libfgb), Void, (UInt32,UInt32,Ref{Cstring}), nvars1, nvars2, varnames)
internal_threads(n) = ccall((:FGb_int_internal_threads, :libfgb), Void, (UInt32,), n)
mod_internal_threads(n) = ccall((:FGb_internal_threads, :libfgb), Void, (UInt32,), n)

creat_poly(num_terms) = ccall((:FGb_int_creat_poly, :libfgb), Ptr{Void}, (UInt32,), num_terms)
exp = UInt32[0,1,1,0,0]
set_expos2(poly, term, exp, explen) = ccall((:FGb_int_set_expos2, :libfgb), Void, (Ptr{Void}, UInt32, Ref{UInt32},UInt32), poly, term, exp, explen)
set_coeff_gmp(poly, term, coeff) = ccall((:FGb_int_set_coeff_gmp, :libfgb), Void, (Ptr{Void}, UInt32, Ref{BigInt}), poly,term,coeff)
full_sort_poly2(poly) = ccall((:FGb_int_full_sort_poly2, :libfgb), Void, (Ptr{Void},), poly)

mutable struct SFGB_Comp_Desc
    _compute::UInt32   # in C, an enum
    _nb::Int32
    _force_elim::Int32
    _off::UInt32
    _index::UInt32
    _zone::UInt32
    _memory::UInt32

    _nb2::Int32
    _force_elim2::Int32
    _bk2::UInt32
    _aggressive2::Int32
    _dlim::Int32
    _skip::Int32
    SFGB_Comp_Desc() = begin
        env = new()
        env._compute=1 #FGB_COMPUTE_GBASIS
        env._nb=0
        env._force_elim=0 # if force_elim=1 then return only the result of the elimination
                          # (need to define a monomial ordering DRL(k1,k2) with k2>0 )
        env._off=0        # use to select another subset of prime numbers
        env._index=500000 # This is is the maximal size of the matrices generated by F4
                          # you can increase this value according to your memory
        env._zone=0       # should be 0 
        env._memory=0     # should be 0 
        env
    end
end

groebner(input_basis, input_basis_length, output_basis, _mini, _elim, cputime, _bk0, _step0, _elim0, _env) = ccall((:FGb_int_groebner, :libfgb), UInt32,
       (Ref{Ptr{Void}},UInt32,Ref{Ptr{Void}},UInt32,UInt32,Ref{Float64},UInt32, Int32, UInt32, Ref{SFGB_Comp_Desc}),
       input_basis, input_basis_length, output_basis, _mini, _elim, cputime, _bk0, _step0, _elim0, _env)

nb_terms(p) = ccall((:FGb_int_nb_terms, :libfgb), UInt32, (Ptr{Void},), p)
export_poly_INT_gmp2(n_variables, n_terms, coefficients, exponents, p) = ccall((:FGb_int_export_poly_INT_gmp2, :libfgb), UInt32, (UInt32, UInt32, Ref{Ptr{BigInt}}, Ref{UInt32}, Ptr{Void}), n_variables, n_terms, coefficients, exponents, p)

reset_memory() = ccall((:FGb_int_reset_memory, :libfgb), Void, ())
exit_INT() = ccall((:FGb_int_exit_INT, :libfgb), Void, ())

struct FGbPolynomial{T<:NamedPolynomial}
    p::Ptr{Void}
end

in_FGb = false

function FGb_with(f::Function, ::Type{NP}) where NP<:NamedPolynomial
    global in_FGb
    assert(!in_FGb)

    enter_INT()
    init_urgent(4,2,"DRLDRL",100000,0)
    log_file_handle = Ref{UInt}(0)
    init(1,1,0,log_file_handle)

    reset_coeffs(1)
    T = names(NP)
    varnames = [String(fieldtype(T, i)) for i in 1:nfields(T)]

    reset_expos(length(varnames), 0, varnames)

    internal_threads(1)
    mod_internal_threads(1)

    in_FGb = true
    res = try
        f(p -> convert(FGbPolynomial{NP}, p))
    finally
        in_FGb = false
        reset_memory()
        exit_INT()
    end
    res
end

import Base: convert, show

function show(io::IO, f::FGbPolynomial{T}) where T
    print(io, "FGb(")
    show(io, convert(T,f))
    print(io, ")")
end

function convert(::Type{FGbPolynomial{T}}, f::T) where T<:NamedPolynomial
    p = creat_poly(length(terms(f.p)))
    for (i,t) in enumerate(terms(f.p))
        exp = UInt32[e for e in monomial(t).e]
        set_expos2(p,i-1,exp,length(exp))
        set_coeff_gmp(p,i-1,coefficient(t))
    end
    full_sort_poly2(p)
    FGbPolynomial{T}(p)
end

function convert(::Type{T}, f::FGbPolynomial{T}) where T<:NamedPolynomial
    n_vars = nfields(names(T))
    n_terms = nb_terms(f.p)
    exponents = Vector{UInt32}(n_terms * n_vars)
    coeff_ptrs = Vector{Ptr{BigInt}}(n_terms)

    code = export_poly_INT_gmp2(n_vars, n_terms, coeff_ptrs, exponents, f.p)

    coefficients = map(unsafe_load, coeff_ptrs)

    sum( c * construct_monomial(T, exponents[(1+k*n_vars):((k+1)*n_vars)]) for (k,c) in zip(0:(n_terms-1), coefficients) )
end

function groebner(F::Vector{T}) where T <: FGbPolynomial
    OUTPUTSIZE=100000
    output_basis = Vector{Ptr{Void}}(OUTPUTSIZE)
    env = SFGB_Comp_Desc()
    cputime = Ref{Float64}(0)
    n = groebner(map(f->f.p, F), length(F), output_basis, 1, 0, cputime,0, -1, 0, env)
    [T(output_basis[i]) for i=1:n]
end

end