#
# Unit tests for GalacticPotentials.jl
#

using GalacticPotentials, Test
using ModelingToolkit, Symbolics

import GalacticPotentials: ScalarSystem

@testset "ScalarSystem Construction" begin
    @variables t
    p = @parameters b
    q = @variables x(t) y(t) z(t)

    @test ScalarSystem(
        (b // 2) * sum(q .^ 2),
        t,
        q,
        p;
        name=:SomeSystem
    ) isa ModelingToolkit.AbstractSystem
end