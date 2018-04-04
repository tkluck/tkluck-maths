using Base.Test
using GröbnerSingular

@testset "Groebner testset" begin
    # will run these tests, but with GröbnerSingular set as the default
    # engine
    include(Pkg.dir("PolynomialRings", "test", "Groebner.jl"))
end
