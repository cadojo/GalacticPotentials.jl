"""
All symbolic variables used throughout `GalacticPotentials`.
"""
module Symbols

export t, D, V, x, y, z, ẋ, ẏ, ż, u, u̇

import ModelingToolkit: @variables, t_nounits as t, D_nounits as D

@variables V(t) [output = true]
u = @variables x(t) [input = true] y(t) [input = true] z(t) [input = true]
u̇ = @variables ẋ(t) [input = true] ẏ(t) [input = true] ż(t) [input = true]

end