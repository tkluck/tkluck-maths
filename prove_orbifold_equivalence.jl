using PolynomialRings
using PolynomialRings: terms
using OrbifoldEquivalence
using QuasiHomogeneous

using FGb

using StatProfilerHTML

const n_variables = 1

# within the if block, we define:
# c1 - an unused formal_coefficient(R), to be used for homogenization
# Q - generic matrix with correct gradings
# W, V - the two potentials for which we compute equivalence
# qdim1 , qdim2 - the qdims of Q w.r.t V and W

# -----------------------------
# one variable case: reflexivity
# -----------------------------
if n_variables == 1
    R = @ring! ℚ[x,y]

    const d = 6
    W = x^d
    V = y^d

    @ring! ℤ[c[]]
    c1 = c()

    gr=[-1 1;d-1 -1]
    vgr=[:x=>1,:y=>1]

    Q = QuasiHomogeneous.generic_quasihomogeneous_map(gr, vgr, c)

    qdim1 = quantum_dimension(Q,W,[:x],[:y])
    qdim2 = quantum_dimension(Q,V,[:y],[:x])

# -----------------------------
# two variable case: reflexivity
# -----------------------------
elseif n_variables == 2
    R = @ring! ℚ[x,y,u,v]

    const d = 3
    W = x^d - y^d
    V = u^d - v^d

    @ring! ℤ[c[]]
    c1 = c()

    gr = [-1 -1 d-1 d-1;
          -1 -1 1    1;
          1  d-1 -1 -1;
          1  d-1 -1 -1]

    vgr=[:x=>1,:y=>1,:u=>1,:v=>1]

    Q = QuasiHomogeneous.generic_quasihomogeneous_map(gr, vgr, c)

    qdim1 = quantum_dimension(Q,W,[:x,:y],[:u,:v])
    qdim2 = quantum_dimension(Q,V,[:u,:v],[:x,:y])

# -----------------------------
# three variable case: E14 ~ Q10
# -----------------------------
elseif n_variables == 3
    using MatrixFactorizations
    R = @ring! ℚ[x,y,z,u,v,w]

    W = x^4 + y^3 + x*z^2
    V = u^4*w + v^3 + w^2

    @ring! ℤ[c[]]
    c1 = c()

    gr = MatrixFactorizations.E14_Q10_grading()
    vgr=[:x=>6,:y=>8,:z=>9,:u=>3,:v=>8,:w=>12]

    Q = QuasiHomogeneous.generic_quasihomogeneous_map(gr, vgr, c)

    qdim1 = quantum_dimension(Q,W,[:x,:y,:z],[:u,:v,:w])
    qdim2 = quantum_dimension(Q,V,[:u,:v,:w],[:x,:y,:z])

end

# -----------------------------
# computation
# -----------------------------
C = [coefficient(t) for entry in (Q^2 - c1^2*(V-W)*eye(Int,size(Q)...)) for t in terms(entry)]


CC = groebner_basis(C)

qdim1_red = rem(qdim1, CC)
qdim2_red = rem(qdim2, CC)
@show qdim1_red
@show qdim2_red
