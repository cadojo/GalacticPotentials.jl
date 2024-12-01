using ModelingToolkit

function ScalarField(value, t, u, p; name, simplify = true, kwargs...)
    @variables Φ(t) [output = true]

    eqs::Vector{Equation} = vcat(
        (Differential(t)^2).(u) .~ -ModelingToolkit.gradient(value, u, simplify = simplify)
    )
    return ODESystem(vcat(Φ ~ value), t; name = name, kwargs...)
end

function PlummerPotential(; name = :PlummerPotential, kwargs...)
    @independent_variables t
    u = @variables x(t) [input = true] y(t) [input = true] z(t) [input = true]
    p = @parameters G m b

    value = -G * m / sqrt(b^2 + x^2 + y^2 + z^2)
    return ScalarField(value, t, u, p; name = name, kwargs...)
end

function MilkyWay(; name = :MilkyWay, alias_equations = false, kwargs...)
    @named disk = PlummerPotential()
    @named bulge = PlummerPotential()
    @named halo = PlummerPotential()

    @independent_variables t
    @variables Φ(t) [output = true]
    D = Differential(t)
    u = @variables x(t) [input = true] y(t) [input = true] z(t) [input = true]
    du = @variables dx(t) [input = true] dy(t) [input = true] dz(t) [input = true]

    aliases = [
        disk.x => x,
        disk.y => y,
        disk.z => z,
        bulge.x => x,
        bulge.y => y,
        bulge.z => z,
        halo.x => x,
        halo.y => y,
        halo.z => z
    ]

    grad(sys) = ModelingToolkit.gradient(sys.Φ, [sys.x, sys.y, sys.z])

    eqs = vcat(
        Φ ~ disk.Φ + bulge.Φ + halo.Φ,
        D.(u) .~ du,
        D.(du) .~ -(grad(disk) .+ grad(bulge) .+ grad(halo))
    )

    if alias_equations
        append!(eqs, [alias.first ~ alias.second for alias in aliases])
    end

    return compose(
        ODESystem(
            eqs, t; name = :MilkyWay, defaults = Dict(aliases)),
        disk, bulge, halo
    )
end

mw = MilkyWay(alias_equations = true)