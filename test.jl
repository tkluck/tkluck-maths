using Polynomials
using Polynomials: x, y, _reduce

f = x*y^3 + y
g = x+y
h = y*y

println(f)
println(g)
println(h)
println(_reduce(f, [g]))

println(groebner_basis([f, g, h]))
@time groebner_basis([f, g, h])


