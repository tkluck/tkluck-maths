using PolynomialRings
using PolynomialRings: terms
using OrbifoldEquivalence
using QuasiHomogeneous

using StatProfilerHTML

const n_variables = 3
const dense_monomials = false
const use_FGb = false          # only with dense monomials!

# within the if block, we define:
# c1 - an unused formal_coefficient(R), to be used for homogenization
# Q - generic matrix with correct gradings
# W, V - the two potentials for which we compute equivalence
# qdim1 , qdim2 - the qdims of Q w.r.t V and W

# -----------------------------
# one variable case: reflexivity
# -----------------------------
if n_variables == 1
    R,(x,y) = polynomial_ring(Rational{BigInt}, :x, :y)

    const d = 6
    W = x^d
    V = y^d

    ch = formal_coefficients(R,:c)
    c1 = take!(ch)

    gr=[-1 1;d-1 -1]
    vgr=(1,1)

    Q = QuasiHomogeneous.generic_quasihomogeneous_map(gr, vgr, R, ch)

    qdim1 = @coefficient quantum_dimension(Q,W,[:x],[:y]) x^0*y^0
    qdim2 = @coefficient quantum_dimension(Q,V,[:y],[:x]) x^0*y^0

# -----------------------------
# two variable case: reflexivity
# -----------------------------
elseif n_variables == 2
    R,(x,y,u,v) = polynomial_ring(Rational{BigInt}, :x, :y, :u, :v)

    const d = 3
    W = x^d - y^d
    V = u^d - v^d

    ch = formal_coefficients(R,:c)
    c1 = take!(ch)

    gr = [-1 -1 d-1 d-1;
          -1 -1 1    1;
          1  d-1 -1 -1;
          1  d-1 -1 -1]

    vgr=(1,1,1,1)

    Q = QuasiHomogeneous.generic_quasihomogeneous_map(gr, vgr, R, ch)

    qdim1 = @coefficient quantum_dimension(Q,W,[:x,:y],[:u,:v]) x^0*y^0*u^0*v^0
    qdim2 = @coefficient quantum_dimension(Q,V,[:u,:v],[:x,:y]) x^0*y^0*u^0*v^0

# -----------------------------
# three variable case: E14 ~ Q10
# -----------------------------
elseif n_variables == 3
    using MatrixFactorizations
    R,(x,y,z,u,v,w) = polynomial_ring(Rational{BigInt}, :x, :y, :z, :u, :v, :w)

    W = x^4 + y^3 + x*z^2
    V = u^4*w + v^3 + w^2

    ch = formal_coefficients(R,:c)
    c1 = take!(ch)

    gr = MatrixFactorizations.E14_Q10_grading()
    vgr= (6,8,9,3,8,12)

    Q = QuasiHomogeneous.generic_quasihomogeneous_map(gr, vgr, R, ch)

    qdim1 = @coefficient quantum_dimension(Q,W,[:x,:y,:z],[:u,:v,:w]) x^0*y^0*z^0*u^0*v^0*w^0
    qdim2 = @coefficient quantum_dimension(Q,V,[:u,:v,:w],[:x,:y,:z]) x^0*y^0*z^0*u^0*v^0*w^0

end

# -----------------------------
# computation
# -----------------------------
C = [coefficient(t) for entry in (Q^2 - c1^2*(V-W)*eye(Int,size(Q)...)) for t in terms(entry.p)]

if dense_monomials
    yada = to_dense_monomials([qdim1; qdim2; C])
    qdim1 = yada[1]
    qdim2 = yada[2]
    C = yada[3:end]
end

if !use_FGb || !dense_monomials
    println(STDERR, "Computing")
    CC = groebner_basis(C, Val{false}, max_degree=max(deg(qdim1),deg(qdim2)))
    println(STDERR, "Done")
else
    using FGb
    CC = FGb_with(eltype(C)) do FGbPolynomial
        println(STDERR, "Converting to FGb representation...")
        G = map(FGbPolynomial, C_dense)
        println(STDERR, "Computing groebner basis...")
        @time gr = groebner(G)
        println(STDERR, "Done")
        map(g -> convert(eltype(C_dense),g), gr)
    end
end

qdim1_red,factors = red(qdim1, CC)
qdim2_red,factors = red(qdim2, CC)
@show qdim1_red
@show qdim2_red
