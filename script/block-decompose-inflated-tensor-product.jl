using PolynomialRings
using MatrixFactorizations

@ring! ℚ[x,y,z]

# A = unit_matrix_factorization(x^3, [:x], [:y])
# B = unit_matrix_factorization(y^3, [:y], [:z])
A = [0 x-y; x^2 + x*y + y^2 0]
B = [0 y-z; y^2 + y*z + z^2 0]


# run this line twice! The first time gives a 'world age' related error.
# the second time around, that has been fixed.
Q = ⨶(A,B,y^3,:y)

# ----
# Find a conjugation (sequence of row/col operations) so Q gets
# a block structure.
# ----
TR = Q[1:4,5:8]
BL = Q[5:8,1:4]

step1 = [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1]
step2 = [1 0 0 0; z 1 0 0; -x^2 0 1 0; -x 0 0 1]
step3 = [1 z 1 -x; 0 1 0 0; 0 0 1 0; 0 0 0 1]
step4 = [1 0 0 0; 0 1 0 0; 0 -x + -z 1 0; 0 -1 0 1]
step5 = [1 0 0 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]
step6 = [1 0 0 0; 0 1 0 0; 0 0 1 -z; 0 0 0 1]
step7 = [1 0 0 0; 0 1 0 0; 0 0 1 z; 0 0 0 1]
step8 = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
step9 = [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]

pets1 = [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1]
pets2 = [1 0 0 0; -z 1 0 0; x^2 0 1 0; x 0 0 1]
pets3 = [1 -z -1 x; 0 1 0 0; 0 0 1 0; 0 0 0 1]
pets4 = [1 0 0 0; 0 1 0 0; 0 x + z 1 0; 0 1 0 1]
pets5 = [1 0 0 0; 0 1 0 -1; 0 0 1 0; 0 0 0 1]
pets6 = [1 0 0 0; 0 1 0 0; 0 0 1 z; 0 0 0 1]
pets7 = [1 0 0 0; 0 1 0 0; 0 0 1 -z; 0 0 0 1]
pets8 = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
pets9 = [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]

R = eltype(TR)

TR2 = step8 * step6 * step4 * step2 * step1 * TR * step3 * step5 * pets7 * step9
BL2 = pets9 * step7 * pets5 * pets3 * BL * pets1 * pets2 * pets4 * pets6 * pets8

# Best block structure I could find
QQ = R[zeros(TR) TR2; BL2 zeros(BL)]
