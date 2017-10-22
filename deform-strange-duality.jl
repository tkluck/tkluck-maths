using PolynomialRings; using MFDeformations; using MatrixFactorizations; using OrbifoldEquivalence; using CommutativeAlgebras
using FlintNumbers

@ring ℚ[a,b,c]
Q, f1,f2, g = MatrixFactorizations.E14_Q10_matrix_and_equations()
W = MatrixFactorizations.E14_Q10()
f = f1*f2

F = FlintNumberField(@ring(ℚ[a,b,c])/Ideal(f1,g,c))
QQ = map(@ring(F[x,y,z,u,v,w]), Q)
@assert QQ^2 == W*eye(QQ)

gr = QuasiHomogeneous.find_quasihomogeneous_degrees(W,:x,:y,:z,:u,:v,:w)
T = MFDeformations.graded_implicit_tangent_space(Q->Q^2, QQ, gr)

R = eltype(eltype(T))
QQQ = sum(prod, zip(formal_coefficients(R, :ε), T))
@show @linear_coefficients quantum_dimension(QQQ, W, [:u,:v,:w], [:x,:y,:z]) ε[]
@show @linear_coefficients quantum_dimension(QQQ, W, [:x,:y,:z], [:u,:v,:w]) ε[]
