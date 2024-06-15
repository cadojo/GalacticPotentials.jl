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