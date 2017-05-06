using Polynomials
using Polynomials: x, y, _reduce

f = x*y^3 + y
g = x+y
h = y*y

println(f)
println(g)
println(h)
println(_reduce(f, [g]))

display(groebner_basis([f, g]))
println()
display(groebner_basis([f, g, h]))
@time groebner_basis([f, g, h])

(basis, transformation) = groebner_basis([f, g])
display(basis)
println()
display(syzygies(basis))
println()
display( syzygies(basis) * basis )
println()


using Modules

F = ModuleMorphism{typeof(x)}([x 0; 0 y])

display(lift([x, -y], 3*x + x*y))
println()

display(lift(F, [2*x; x*y]))
println()

G = ModuleMorphism{typeof(x)}([x x*y; y y^2])
display(kernel(G))
println()

display(G.m * kernel(G))
println()
