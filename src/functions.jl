#
# Mathematical function registration with ModelingToolkit
#

"""
Computes the lower incomplete gamma function.
"""
function lowergamma(a, x)
    p, q = SpecialFunctions.gamma_inc(a, x)
    return p * gamma(a)
end

"""
Computes the gamma function through `SpecialFunctions.gamma`.

!!! note
    This wrapper is required for Symbolics function registration purposes. 
"""
function gamma(x)
    return SpecialFunctions.gamma(x)
end

@register_symbolic lowergamma(x::Real, y::Real)
@register_symbolic gamma(x::Real)

function Symbolics.derivative(::typeof(gamma), args::NTuple{1, Any}, ::Val{1})
    return ForwardDiff.derivative(gamma, first(args))
end

function Symbolics.derivative(::typeof(lowergamma), args::NTuple{2, Any}, ::Val{1})
    x, y = args
    return ForwardDiff.derivative(x -> lowergamma(x, y), x)
end

function Symbolics.derivative(::typeof(lowergamma), args::NTuple{2, Any}, ::Val{2})
    x, y = args
    return ForwardDiff.derivative(y -> lowergamma(x, y), y)
end