using PolynomialRings; using MFDeformations; using MatrixFactorizations; using OrbifoldEquivalence; using CommutativeAlgebras
using FlintNumbers

@ring! ℚ[a,b,c]
Q, f1,f2, g = MatrixFactorizations.E14_Q10_matrix_and_equations()
W = MatrixFactorizations.E14_Q10()
f = f1*f2

F = FlintNumberField(@ring(ℚ[a,b,c])/Ideal(f1,g,c))
@ringname F :𝔽
FF = @ring(F[x,y,z,u,v,w])
QQ = map(FF, Q)
@assert QQ^2 == W*eye(QQ)

gr = QuasiHomogeneous.find_quasihomogeneous_degrees(W,:x,:y,:z,:u,:v,:w)
T = MFDeformations.graded_implicit_tangent_space(Q->Q^2, QQ, gr)

R = eltype(eltype(T))
@polyvar ε[]
QQQ = sum(prod, zip(ε[], T))
@show @linear_coefficients quantum_dimension(QQQ, W, [:u,:v,:w], [:x,:y,:z]) ε[]
@show @linear_coefficients quantum_dimension(QQQ, W, [:x,:y,:z], [:u,:v,:w]) ε[]

e = formal_coefficients(eltype(QQ), :e)
g = diagm(e[1:8])
I = eye(QQ)

TG = @linear_coefficients (I+g)*QQ*(I-g) e[]

@ring! ℚ[ε]
Tc, = @linear_coefficients Q(b=b+ε,c=c+2ε) ε

@show TG
@show Tc

for Tg in TG
    @assert iszero( FF.( Q*Tg + Tg*Q ))
end
@assert iszero( FF.( Q*Tc + Tc*Q ))

to_vector, to_polynomial_array = MFDeformations.finite_subspace_conversion([T; TG; [Tc]], :u,:v,:w,:x,:y,:z)

# TG has a linear dependency coming from scalar multiplication ⊆ the maximal torus
@show PolynomialRings.Util.LinAlgUtil.nullspace(hcat(to_vector.(TG)...))

# so TG + Tc has dimension 8. T has dimension 9. As it turns out, T[4] is the
# independent one (the following kernel is zero):
@show PolynomialRings.Util.LinAlgUtil.nullspace(hcat(to_vector.([[T[4]]; TG[1:7]; [Tc]])...))

@show sparse( T[4] )
