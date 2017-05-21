using PolynomialRings: polynomial_ring
using MFDeformations: diff, H1, deformation


using MatrixFactorizations




println("Finding deformations for:")
Q = T2()
display(Q); println()

symbols = [Symbol("W$i") for i=1:40]
Qdef =deformation(Q, symbols...)
println("Deformation is:")
display(Qdef); println()
println("which squares to:")
display(Qdef^2); println()

