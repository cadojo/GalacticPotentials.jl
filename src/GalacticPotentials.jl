module GalacticPotentials

using DocStringExtensions
@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(SIGNATURES)

                                         $(DOCSTRING)
                                         """

@template (TYPES, CONSTANTS) = """
                               $(TYPEDEF)

                               $(DOCSTRING)
                               """

using Reexport
@reexport using ModelingToolkit

using Symbolics
using LaTeXStrings
using LinearAlgebra
using ForwardDiff
using FiniteDiff
using SciMLBase
using Memoize

import SpecialFunctions
import HypergeometricFunctions

module Symbols

export t, D, V, x, y, z, ẋ, ẏ, ż, u, u̇

import ModelingToolkit: @variables, t_nounits as t, D_nounits as D

@variables V(t) [output = true]
u = @variables x(t) [input = true] y(t) [input = true] z(t) [input = true]
u̇ = @variables ẋ(t) [input = true] ẏ(t) [input = true] ż(t) [input = true]

end

using .Symbols

export
       ScalarField,
       HarmonicOscillatorPotential,
       HenonHeilesPotential,
       HernquistPotential,
       IsochronePotential,
       JaffePotential,
       KeplerPotential,
       KuzminPotential,
       LogarithmicPotential,
       LongMuraliBarPotential,
       MiyamotoNagaiPotential,
       NFWPotential,
       PlummerPotential,
       PowerLawCutoffPotential,
       SatohPotential,
       StonePotential

include(joinpath(@__DIR__, "gen", "expressions.jl"))
include("functions.jl")
include("potentials.jl")

end # module GalacticPotentials