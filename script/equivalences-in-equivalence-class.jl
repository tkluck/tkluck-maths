using PolynomialRings
using OrbifoldEquivalence

# alternate groebner basis implementations
# using FGb
#PolynomialRings.Backends.Gröbner.set_default(PolynomialRings.Backends.Gröbner.Buchberger())

@ring! ℤ[x,y,u,v]

potentials = [x^4 + x*y^3, x^3 + x*y^4, x^3*y + x*y^3, x^6 + y^3, x^4 + y^4]

let i = 1 #for i = 1:length(potentials)
    let j = 1 #for j = 1:i
        @show potentials[i], potentials[j]
        @show equivalence_exists(potentials[i], [:x,:y], potentials[j](x=u,y=v), [:u,:v], 2)
    end
end
