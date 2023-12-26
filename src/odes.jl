function ModelingToolkit.ODESystem(field::AbstractScalarField)

    t = ModelingToolkit.get_iv(field)
    u = states(field)
    u̇ = let du = map(x -> Symbol(first(split(x, "("))), "Δ" .* string.(u))
        vcat(
            (@variables($(δ)(t)) for δ in du)...
        )
    end

    p = parameters(field)
    Δ = Differential(t)

    name = field.name

    eqs = vcat(
        Δ.(u) .~ u̇,
        Δ.(u̇) .~ -calculate_gradient(field),
    )

    return ODESystem(
        eqs, t, vcat(u, u̇), p; name=name,
    )
end


function SciMLBase.ODEProblem(field::AbstractScalarField, args...; kwargs...)

    return ODEProblem(
        ODESystem(field),
        args...; kwargs...,
    )

end