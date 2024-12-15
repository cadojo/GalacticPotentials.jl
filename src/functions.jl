#
# Mathematical function registration with ModelingToolkit
#

"""
Computes the lower incomplete gamma function.
"""
function lowergamma(a, x)
    p, q = SpecialFunctions.gamma_inc(a, x)
    return p * SpecialFunctions.gamma(a)
end

@register_symbolic lowergamma(x::AbstractFloat, y::AbstractFloat)
# @register_symbolic SpecialFunctions.gamma(x::AbstractFloat)
@register_symbolic HypergeometricFunctions.M(
    a::AbstractFloat, b::AbstractFloat, x::AbstractFloat)

function Symbolics.derivative(::typeof(lowergamma), args::NTuple{2, Any}, ::Val{1})
    a, x = args
    return (x^a / a) * HypergeometricFunctions.M(a, a + 1, -x) # https://github.com/JuliaMath/HypergeometricFunctions.jl/issues/50#issuecomment-1397363491
end

function Symbolics.derivative(::typeof(lowergamma), args::NTuple{2, Any}, ::Val{2})
    a, x = args
    return -x^(a - 1) * exp(-x) # https://en.wikipedia.org/wiki/Incomplete_gamma_function#Derivatives
end
