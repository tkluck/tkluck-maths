using PolynomialRings: polynomial_ring
using MFDeformations: diff, H1, deformation


using MatrixFactorizations




println("Finding deformations for:")
Q = T2()
display(Q); println()

A, (t,u,v) = polynomial_ring(Int, :t, :u, :v)

Qdef =deformation(Q, :t, :u, :v)
println("Deformation is:")
display(Qdef); println()
println("which squares to:")
display(Qdef^2); println()

