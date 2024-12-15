#
# Can we use structs to de-duplicate code? This file is experimental!
#

abstract type AbstractPotential end

struct ScalarPotential <: AbstractPotential
    value::Num
    u::Vector{Num}
    du::Union{Nothing, Vector{Num}}
    p::Union{Nothing, Vector{Num}}
    name::Symbol
    ScalarPotential(value, u; du = nothing, p = nothing, name) = new(value, u, du, p, name)
end

function ModelingToolkit.ODESystem(potential::ScalarPotential; simplify = true, kwargs...)
    @variables Φ(t) # t is globally defined as MTK.t_nounits

    ∂V∂u = Symbolics.gradient(potential.value, potential.u, simplify = simplify)

    if isnothing(potential.du)
        u = potential.u
        eqs = [
            D.(D.(potential.u)) .~ -∂V∂u
        ]
    else
        u = vcat(potential.u, potential.du)
        eqs = [
            D.(potential.u) .~ potential.du,
            D.(potential.du) .~ -∂V∂u
        ]
    end

    defaults = (; name = potential.name)
    options = merge(defaults, kwargs)

    if isnothing(potential.p)
        return ODESystem(eqs, t; options...)
    else
        return ODESystem(eqs, t, u, p; options...)
    end
end