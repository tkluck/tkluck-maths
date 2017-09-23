using PolynomialRings
using MFDeformations: deformation

using MatrixFactorizations




println("Finding deformations for:")
Q = base_extend(T2())
display(Q); println()

symbols = [Symbol("w$i") for i=1:40]
Qdef =deformation(Q, symbols...)
println("Deformation is:")
display(Qdef); println()
println("which squares to:")
display(Qdef^2); println()

