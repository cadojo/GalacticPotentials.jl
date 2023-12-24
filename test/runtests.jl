#
# Unit tests for GalacticPotentials.jl
#

using GalacticPotentials, Test
using ModelingToolkit, Symbolics

using GalacticPotentials: AbstractField, AbstractScalarField, ScalarField

@testset verbose = true "Scalar Fields" begin
    @variables t
    p = @parameters b
    q = @variables x(t) y(t) z(t)

    field = ScalarField(
        (b // 2) * sum(q .^ 2),
        t,
        q,
        p;
        name=:SomeField
    )

    @testset showtiming = true "Constructors" begin
        @test field isa ModelingToolkit.AbstractSystem
    end

    @testset showtiming = true "Calculations" begin
        @test all(calculate_jacobian(field) - calculate_gradient(field) .== 0)
        @test calculate_hessian(field) isa AbstractMatrix
    end

end

@testset verbose = true "Galactic Potentials" begin
    for name in names(GalacticPotentials)
        !occursin("Potential", "$name") && continue
        occursin("Potentials", "$name") && continue

        @eval field = $name()

        @testset showtiming = true "$name" begin
            N = length(states(field))
            M = length(parameters(field))
            @test field isa AbstractField
            @test calculate_gradient(field) isa AbstractVector
            @test ODESystem(field) isa ODESystem
            @test ODEProblem(field, randn(N), (rand(), rand()), randn(M)) isa ODEProblem
        end
    end

end
