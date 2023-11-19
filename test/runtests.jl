#
# Unit tests for GalacticPotentials.jl
#

using GalacticPotentials, Test
using ModelingToolkit, Symbolics

import GalacticPotentials: ScalarSystem

@testset "Scalar Fields" begin
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

    @testset "Constructors" begin
        @test field isa ModelingToolkit.AbstractSystem
    end

    @testset "Calculations" begin
        @test all(calculate_jacobian(field) - calculate_gradient(field) .== 0)
        @test calculate_hessian(field) isa AbstractMatrix
    end


end