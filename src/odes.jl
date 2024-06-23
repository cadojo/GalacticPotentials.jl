"""
Construct an `ODESystem` from an `AbstractScalarField` by taking the gradient of 
the field's value with respect to the state variables. Symbols which describe 
the state derivatives can be provided via the `var_map_to_dvs` keyword argument.
Alternatively, default state derivative symbols are used: :ẋ, :ẏ, and :ż for 
:x, :y, :z states, and :Δ\$(state) otherwise.

## Example

    ODESystem(field; var_map_to_dvs = Dict(:x => :ẋ, :y => :ẏ, :z => :ż))
"""
function ModelingToolkit.ODESystem(field::AbstractScalarField;
        var_map_to_dvs::Union{<:AbstractDict, <:Nothing} = nothing)
    t = ModelingToolkit.get_iv(field)
    u = unknowns(field)

    if isnothing(var_map_to_dvs)
        if string.(collect(u)) == ["x($t)", "y($t)", "z($t)"][CartesianIndices(u)]
            var_map_to_dvs = Dict(:x => :ẋ, :y => :ẏ, :z => :ż)
            du = map(
                x -> Symbol(var_map_to_dvs[Symbol(first(split(
                    string(x), "($(Symbolics.value(t)))")))]),
                u)
        else
            du = map(
                x -> Symbol(:Δ, Symbol(first(split(string(x), "($(Symbolics.value(t)))")))),
                u)
        end
    else
        du = map(
            x -> Symbol(var_map_to_dvs[Symbol(first(split(
                string(x), "($(Symbolics.value(t)))")))]),
            u)
    end

    u̇ = getfield.(
        vcat((@variables($(δ)(t)) for δ in du)...),
        :val
    )

    p = parameters(field)
    Δ = Differential(t)

    name = ModelingToolkit.get_name(field)

    eqs = vcat(
        Δ.(u) .~ u̇,
        Δ.(u̇) .~ -calculate_gradient(field)
    )

    return ODESystem(
        eqs, t, vcat(u, u̇), p; name = name
    )
end

function SciMLBase.ODEProblem(field::AbstractScalarField, args...; kwargs...)
    return ODEProblem(
        complete(ODESystem(field); split = false),
        args...; kwargs...
    )
end