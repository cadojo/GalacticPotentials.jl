#
# Mathematical function registration with ModelingToolkit
#

"""
Computes the lower incomplete gamma function.
"""
function lowergamma(a, x)
    p, q = gamma_inc(a, x)
    return p * gamma(a)
end

@register_symbolic lowergamma(x, y)
@register_symbolic gamma(x)

function ModelingToolkit.derivative(::typeof(gamma), args::NTuple{1, Any}, ::Val{1})
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