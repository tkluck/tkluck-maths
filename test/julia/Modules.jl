using Base.Test
using PolynomialRings: polynomial_ring
using PolynomialRings.Modules: ModuleMorphism, lift, kernel

@testset "Modules" begin

    R,(x,y) = polynomial_ring(Rational{Int}, :x, :y)

    begin
        mapping = [x, -y]
        image = 3*x + x*y
        preimage = lift(mapping, image)
        @test get(preimage) * mapping == [ image ]
    end

    begin
        F = ModuleMorphism{typeof(x)}([x 0; 0 y])
        image = [2*x; x*y]
        F( get(lift(F, image)) ) == image
    end

    begin
        G = ModuleMorphism{typeof(x)}([x x*y; y y^2])
        zz = G.m * kernel(G)
        @test zz == zero(zz)
    end

end
