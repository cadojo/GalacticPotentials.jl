function ModelingToolkit.ODESystem(field::AbstractScalarField)

    t = ModelingToolkit.get_iv(field)
    u = states(field)
    p = parameters(field)
    Δ = Differential(t)

    name = field.name

    eqs = Δ.(u) .~ calculate_gradient(field)

    return ODESystem(
        eqs, t, u, p; name=name,
    )
end

function SciMLBase.ODEProblem(field::AbstractScalarField, args...; kwargs...)

    return ODEProblem(
        ODESystem(field),
        args...; kwargs...,
    )

end