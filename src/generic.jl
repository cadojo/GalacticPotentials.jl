const EMPTY_GRADIENT = Num[]

"""
An abstract supertype for all scalar expressions which behave like systems, such as 
scalar potentials.
"""
abstract type AbstractScalarSystem <: ModelingToolkit.AbstractSystem end

"""
A scalar potential field.
"""
struct ScalarSystem <: AbstractScalarSystem
    expr::Num
    iv::Symbolics.BasicSymbolic{Real}
    states::Vector
    ps::Vector
    jac::Base.RefValue{Any}
    name::Symbol

    function ScalarSystem(expr, iv, states, ps, jac, name; checks::Union{Bool,Int}=true)
        if checks == true || (checks & ModelingToolkit.CheckUnits) > 0
            ModelingToolkit.all_dimensionless([states; ps]) || ModelingToolkit.check_units([expr])
        end
        new(expr, iv, states, ps, jac, name)
    end
end

ModelingToolkit.get_eqs(::ScalarSystem) = Equation[]
ModelingToolkit.get_systems(::ScalarSystem) = ModelingToolkit.AbstractSystem[]

function ScalarSystem(expr, iv, states, ps; name=nothing)

    if isnothing(name)
        throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    end

    jac = Base.RefValue{Any}(EMPTY_GRADIENT)

    return ScalarSystem(expr, Symbolics.value(iv), Symbolics.value.(Symbolics.scalarize(states)), Symbolics.value.(Symbolics.scalarize(ps)), jac, name)

end