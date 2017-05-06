using Polynomials
using Polynomials: x, y, _reduce

f = x*y^3 + y
g = x+y
h = y*y

println(f)
println(g)
println(h)
println(_reduce(f, [g]))

println(groebner_basis([f, g]))
println(groebner_basis([f, g, h]))
@time groebner_basis([f, g, h])

(basis, transformation) = groebner_basis([f, g])
println(basis)
println(syzygies(basis))

println( syzygies(basis) * basis )


using Modules

F = ModuleMorphism{typeof(x)}([x 0; 0 y])

println(lift([x, -y], 3*x + x*y))

println(lift(F, [2*x; x*y]))

G = ModuleMorphism{typeof(x)}([x x*y; y y^2])
println(kernel(G))

println(G.m * kernel(G))
