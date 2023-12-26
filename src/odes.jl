function ModelingToolkit.ODESystem(field::AbstractScalarField; var_map_to_dvs::AbstractDict{Symbol,Symbol}=Dict(:x => :ẋ, :y => :ẏ, :z => :ż))

    t = ModelingToolkit.get_iv(field)
    u = states(field)
    u̇ = let du = map(x -> Symbol(var_map_to_dvs[Symbol(first(split(string(x), "($(Symbolics.value(t)))")))]), u)
        vcat(
            (@variables($(δ)(t)) for δ in du)...
        )
    end

    p = parameters(field)
    Δ = Differential(t)

    name = field.name # TODO: when ModelingToolkit.jl updates, change to get_name

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