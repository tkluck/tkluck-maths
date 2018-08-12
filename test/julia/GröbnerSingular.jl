using Test
using PolynomialRings
using GröbnerSingular
GröbnerSingular.enable()

@testset "Groebner testset" begin
    # will run these tests, but with GröbnerSingular set as the default
    # engine
    include(joinpath(dirname(pathof(PolynomialRings)), "../test/Groebner.jl"))
end
