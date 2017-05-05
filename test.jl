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


using Modules
using Modules: lift

F = Morphism{typeof(x),2,2}([x 0; 0 y])

println(lift([x, -y], 3*x + x*y))

println(lift(F, [2*x; x*y]))

