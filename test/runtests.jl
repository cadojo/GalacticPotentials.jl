#
# Unit tests for GalacticPotentials.jl
#

using GalacticPotentials, Test
using ModelingToolkit, Symbolics

using GalacticPotentials: ScalarField

@testset verbose=true "Scalar Fields" begin
    @independent_variables t
    p = @parameters b
    q = @variables x(t) y(t) z(t)

    field = ScalarField(
        (b // 2) * sum(q .^ 2),
        t,
        q,
        p;
        name = :SomeField
    )

    @testset showtiming=true "Constructors" begin
        @test field isa ModelingToolkit.AbstractSystem
    end

    @testset showtiming=true "Calculations" begin
        @test calculate_jacobian(field) isa AbstractMatrix
    end
end

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
    mw = GalacticPotentials.MilkyWay.Bovy2014()
    mws = structural_simplify(mw)

    u = randn(6)

    problem = ODEProblem(mws, u, (0.0, 10.0), [])
    @test problem.f(u, problem.p, 0.0) isa AbstractVector
end