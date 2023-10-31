#
# Interfaces for dynamical expressions, like scalar potential fields
#

abstract type AbstractField <: ModelingToolkit.AbstractTimeDependentSystem end
abstract type AbstractScalarField <: AbstractField end

# 
# The implementations below are copied and modified versions of systems found in 
# ModelingToolkit.jl. See the ModelingToolkit LICENSE below for more information.
# 
# 
# The ModelingToolkit.jl package is licensed under the MIT "Expat" License:
#
# > Copyright (c) 2018-22: Yingbo Ma, Christopher Rackauckas, Julia Computing, and
# > contributors
# > 
# > Permission is hereby granted, free of charge, to any person obtaining a copy
# > 
# > of this software and associated documentation files (the "Software"), to deal
# > 
# > in the Software without restriction, including without limitation the rights
# > 
# > to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# > 
# > copies of the Software, and to permit persons to whom the Software is
# > 
# > furnished to do so, subject to the following conditions:
# > 
# > The above copyright notice and this permission notice shall be included in all
# > 
# > copies or substantial portions of the Software.
# > 
# > THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# > 
# > IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# > 
# > FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# > 
# > AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# > 
# > LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# > 
# > OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# > 
# > SOFTWARE.

# The code in `src/structural_transformation/bipartite_tearing/modia_tearing.jl`,
# which is from the [Modia.jl](https://github.com/ModiaSim/Modia.jl) project, is
# licensed as follows:

# MIT License

# Copyright (c) 2017-2018 ModiaSim developers

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

Base.@kwdef struct ScalarField <: AbstractScalarField
    """
    A tag for the system. If two systems have the same tag, then they are
    structurally identical.
    """
    tag::UInt
    """Vector of equations defining the system."""
    value::Num
    """Independent variables."""
    iv::SymbolicUtils.BasicSymbolic{Real}
    """Unknown variables."""
    states::Vector
    """Parameters."""
    ps::Vector
    """Array variables."""
    var_to_name::Any
    """Observed states."""
    observed::Vector{Equation}
    """
    Time gradient. Note: this field will not be defined until
    [`calculate_tgrad`](@ref) is called on the system.
    """
    tgrad::Base.RefValue{Any}
    """
    Jacobian matrix. Note: this field will not be defined until
    [`calculate_jacobian`](@ref) is called on the system.
    """
    jac::Base.RefValue{Any}
    """
    The name of the system.
    """
    name::Symbol
    """
    The internal systems. These are required to have unique names.
    """
    systems::Vector{ScalarField}
    """
    The default values to use when initial conditions and/or
    parameters are not supplied in `ODEProblem`.
    """
    defaults::Dict
    """
    Type of the system.
    """
    connector_type::Any
    """
    Metadata for the system, to be used by downstream packages.
    """
    metadata::Any
    """
    Metadata for MTK GUI.
    """
    gui_metadata::Union{Nothing,ModelingToolkit.GUIMetadata}
    """
    Cache for intermediate tearing state.
    """
    tearing_state::Any
    """
    Substitutions generated by tearing.
    """
    substitutions::Any
    """
    If a model `sys` is complete, then `sys.x` no longer performs namespacing.
    """
    complete::Bool
    """
    The hierarchical parent system before simplification.
    """
    parent::Any

    function ScalarField(tag, value, iv, states, ps, var_to_name, observed, tgrad, jac, name, systems, defaults, connector_type, metadata=nothing, gui_metadata=nothing, tearing_state=nothing, substitutions=nothing, complete=false, parent=nothing; checks::Union{Bool,Int}=true)
        if checks == true || (checks & ModelingToolkit.CheckUnits) > 0
            ModelingToolkit.all_dimensionless([states; ps]) || ModelingToolkit.check_units(value)
        end
        new(tag, value, iv, states, ps, var_to_name, observed, tgrad, jac, name, systems, defaults,
            connector_type, metadata, gui_metadata, tearing_state, substitutions, complete,
            parent)
    end
end

function ScalarField(
    value, iv, states, ps;
    observed=[],
    name=nothing,
    defaults=Dict(),
    systems=ScalarField[],
    connector_type=nothing,
    checks=true,
    metadata=nothing,
    gui_metadata=nothing)

    if isnothing(name)
        throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    end

    # Move things over, but do not touch array expressions
    #
    # # we cannot scalarize in the loop because `eqs` itself might require
    # scalarization
    value = ModelingToolkit.scalarize(value)
    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end

    tgrad = Base.RefValue{Any}(ModelingToolkit.EMPTY_TGRAD)
    jac = Base.RefValue{Any}(ModelingToolkit.EMPTY_TGRAD)
    defaults = ModelingToolkit.todict(defaults)
    defaults = Dict{Any,Any}(value(k) => value(v) for (k, v) in pairs(defaults))

    iv = ModelingToolkit.scalarize(iv)
    states = ModelingToolkit.scalarize(states)
    ps = ModelingToolkit.scalarize(ps)

    iv = ModelingToolkit.value(iv)
    states = ModelingToolkit.value.(states)
    ps = ModelingToolkit.value.(ps)

    var_to_name = Dict()
    ModelingToolkit.process_variables!(var_to_name, defaults, states)
    ModelingToolkit.process_variables!(var_to_name, defaults, ps)
    isempty(observed) || ModelingToolkit.collect_var_to_name!(var_to_name, (eq.lhs for eq in observed))

    ScalarField(Threads.atomic_add!(ModelingToolkit.SYSTEM_COUNT, UInt(1)),
        value, iv, states, ps, var_to_name, observed, tgrad, jac, name, systems, defaults,
        connector_type, metadata, gui_metadata, checks=checks)

end


function ModelingToolkit.calculate_tgrad(sys::AbstractScalarField;
    simplify=false)
    isempty(ModelingToolkit.get_tgrad(sys)[]) || return ModelingToolkit.get_tgrad(sys)[]  # use cached tgrad, if possible

    # We need to remove explicit time dependence on the state because when we
    # have `u(t) * t` we want to have the tgrad to be `u(t)` instead of `u'(t) *
    # t + u(t)`.
    vs = ModelingToolkit.detime_dvs.(sys.value)
    iv = ModelingToolkit.get_iv(sys)
    xs = states(sys)
    rule = Dict(map((x, xt) -> xt => x, ModelingToolkit.detime_dvs.(xs), xs))
    vs = substitute.(vs, Ref(rule))
    tgrad = expand_derivatives.(map(Differential(iv), vs), simplify)
    reverse_rule = Dict(map((x, xt) -> x => xt, ModelingToolkit.detime_dvs.(xs), xs))
    tgrad = Num.(substitute.(tgrad, Ref(reverse_rule)))
    ModelingToolkit.get_tgrad(sys)[] = tgrad
    return tgrad
end

function ModelingToolkit.calculate_jacobian(sys::AbstractScalarField; simplify=false, dvs=states(sys))
    if isequal(dvs, states(sys))
        cache = ModelingToolkit.get_jac(sys)[]
        if cache isa Tuple && cache[2] == simplify
            return cache[1]
        end
    end

    vs = sys.value

    jac = ModelingToolkit.gradient(vs, dvs, simplify=simplify)

    if isequal(dvs, states(sys))
        ModelingToolkit.get_jac(sys)[] = jac, simplify # cache Jacobian
    end

    return jac
end

function ModelingToolkit.generate_jacobian(sys::AbstractScalarField, vs=states(sys), ps=parameters(sys);
    simplify=false, kwargs...)
    jac = ModelingToolkit.calculate_jacobian(sys, simplify=simplify)
    pre = ModelingToolkit.get_preprocess_constants(jac)
    return Symbolics.build_function(jac, vs, ps; postprocess_fbody=pre, kwargs...)
end


function ModelingToolkit.calculate_hessian(sys::AbstractScalarField; sparse=false, simplify=false, dvs=states(sys))
    if isequal(dvs, states(sys))
        cache = ModelingToolkit.get_jac(sys)[]
        if cache isa Tuple && cache[2] == simplify
            return cache[1]
        end
    end

    vs = ModelingToolkit.calculate_jacobian(sys)
    hess = ModelingToolkit.jacobian(vs, dvs, simplify=simplify)

    return hess
end


function ModelingToolkit.generate_hessian(sys::AbstractScalarField, vs=states(sys), ps=parameters(sys);
    sparse=false, simplify=false, kwargs...)
    hess = ModelingToolkit.calculate_hessian(sys, sparse=sparse, simplify=simplify)
    pre = ModelingToolkit.get_preprocess_constants(hess)
    return Symbolics.build_function(hess, vs, ps; postprocess_fbody=pre, kwargs...)
end

function ModelingToolkit.generate_function(sys::AbstractScalarField, dvs=states(sys), ps=parameters(sys);
    kwargs...)
    vs = sys.value
    pre, sol_states = ModelingToolkit.get_substitutions_and_solved_states(sys)

    return Symbolics.build_function([vs], ModelingToolkit.value.(dvs), ModelingToolkit.value.(ps); postprocess_fbody=pre,
        states=sol_states, kwargs...)
end

function ModelingToolkit.jacobian_sparsity(sys::AbstractScalarField)
    ModelingToolkit.jacobian_sparsity([sys.value],
        states(sys))
end

function ModelingToolkit.hessian_sparsity(sys::AbstractScalarField)
    [ModelingToolkit.hessian_sparsity([sys.value],
        states(sys)) for eq in equations(sys)] # TODO: this does not look right
end

function Base.show(io::IO, mime::MIME"text/plain", sys::AbstractField)
    val = sys.value
    vars = states(sys)
    nvars = length(vars)
    if val isa AbstractArray && eltype(eqs) <: Num
        nvs = count(v -> !(v isa Connection), eqs)
        Base.printstyled(io, "Model $(nameof(sys)) with $nvs"; bold=true)
        # nextras = n_extra_equations(sys)
        # if nextras > 0
        #     Base.printstyled(io, "("; bold=true)
        #     Base.printstyled(io, neqs + nextras; bold=true, color=:magenta)
        #     Base.printstyled(io, ") "; bold=true)
        # end
        Base.printstyled(io, "values\n"; bold=true)
    else
        Base.printstyled(io, "Model $(nameof(sys))\n"; bold=true)
    end
    # The reduced equations are usually very long. It's not that useful to print
    # them.
    #Base.print_matrix(io, eqs)
    #println(io)

    rows = first(displaysize(io)) ÷ 5
    limit = get(io, :limit, false)

    Base.printstyled(io, "States ($nvars):"; bold=true)
    nrows = min(nvars, limit ? rows : nvars)
    limited = nrows < length(vars)
    defs = ModelingToolkit.has_defaults(sys) ? ModelingToolkit.defaults(sys) : nothing
    for i in 1:nrows
        s = vars[i]
        print(io, "\n  ", s)

        if defs !== nothing
            val = get(defs, s, nothing)
            if val !== nothing
                print(io, " [defaults to ")
                show(IOContext(io, :compact => true, :limit => true,
                        :displaysize => (1, displaysize(io)[2])), val)
                print(io, "]")
            end
            description = ModelingToolkit.getdescription(s)
            if description !== nothing && description != ""
                print(io, ": ", description)
            end
        end
    end
    limited && print(io, "\n⋮")
    println(io)

    vars = parameters(sys)
    nvars = length(vars)
    Base.printstyled(io, "Parameters ($nvars):"; bold=true)
    nrows = min(nvars, limit ? rows : nvars)
    limited = nrows < length(vars)
    for i in 1:nrows
        s = vars[i]
        print(io, "\n  ", s)

        if defs !== nothing
            val = get(defs, s, nothing)
            if val !== nothing
                print(io, " [defaults to ")
                show(IOContext(io, :compact => true, :limit => true,
                        :displaysize => (1, displaysize(io)[2])), val)
                print(io, "]")
            end
            description = getdescription(s)
            if description !== nothing && description != ""
                print(io, ": ", description)
            end
        end
    end
    limited && print(io, "\n⋮")

    if ModelingToolkit.has_torn_matching(sys) && ModelingToolkit.has_tearing_state(sys)
        # If the system can take a torn matching, then we can initialize a tearing
        # state on it. Do so and get show the structure.
        state = ModelingToolkit.get_tearing_state(sys)
        if state !== nothing
            Base.printstyled(io, "\nIncidence matrix:"; color=:magenta)
            show(io, mime, ModelingToolkit.incidence_matrix(state.structure.graph, Num(Sym{Real}(:×))))
        end
    end
    return nothing
end
