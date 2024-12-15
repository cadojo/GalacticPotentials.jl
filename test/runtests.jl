#
# Unit tests for GalacticPotentials.jl
#

using GalacticPotentials, Test
using ModelingToolkit, Symbolics

using GalacticPotentials

@testset verbose=true "Galactic Potentials" begin
    for name in names(GalacticPotentials)
        !occursin("Potential", "$name") && continue
        occursin("Potentials", "$name") && continue

        @eval field = $name()

        @testset showtiming=true "$name" begin
            N = length(unknowns(field))
            M = length(parameters(field))
            @test field isa ODESystem
            @test calculate_jacobian(complete(field; split = false)) isa AbstractMatrix
        end
    end
end

@testset verbose=true "Milky Way Potentials" begin
    mw = GalacticPotentials.Bovy2014()
    mws = structural_simplify(mw)

    u = randn(6)

    problem = ODEProblem(mws, u, (0.0, 10.0), [])
    du = problem.f(u, problem.p, 0.0)

    @test du isa AbstractVector
    @test du ≉ zeros(length(du))
end