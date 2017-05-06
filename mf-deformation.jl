using Polynomials

x = Polynomials._Polynomial{Int, 3}([(Polynomials.Exponent((1,0,0)),1)])
y = Polynomials._Polynomial{Int, 3}([(Polynomials.Exponent((0,1,0)),1)])
z = Polynomials._Polynomial{Int, 3}([(Polynomials.Exponent((0,0,1)),1)])

println(x*y*z)

d1 = [z   y^2 x^3 0 ; y -x*z 0 x^3 ; x 0 -x*z -y^2 ; 0 x -y z  ]
d0 = [x*z y^2 x^3 0 ; y -z   0 x^3 ; x 0 -z   -y^2 ; 0 x -y x*z]
z = zeros(typeof(x), 4,4)

Q = [ z d1; d0 z ]

display(Q)
println()
display(Q^2)
println()

using Modules

delta(Q, i,j) = begin d=zeros(Q); d[i,j]=one(eltype(Q)); d end
dQ = [ Q*delta(Q,i,j) + delta(Q,i,j)*Q for i in indices(Q,1), j in indices(Q,2) ]

Polynomials.iszero(zero(typeof(x)))
Polynomials.iszero.(dQ)

dQ = HomspaceMorphism([ dQ[i,j][k,l] for i=1:8,j=1:8,k=1:8,l=1:8])

K = kernel(dQ)
display( K )
println()

display( Q*K[1] + K[1]*Q ); println()
display( Q*K[end] + K[end]*Q ); println()


