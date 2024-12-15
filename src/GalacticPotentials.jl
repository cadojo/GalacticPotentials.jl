module GalacticPotentials

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

using DocStringExtensions
@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(SIGNATURES)

                                         $(DOCSTRING)
                                         """

@template (TYPES, CONSTANTS) = """
                               $(TYPEDEF)

                               $(DOCSTRING)
                               """

using ModelingToolkit
using Symbolics
using LaTeXStrings
using LinearAlgebra
using Memoize

import SpecialFunctions
import HypergeometricFunctions

# module files
include("symbols.jl")

# load modules
using .Symbols

# non-module files
include(joinpath(@__DIR__, "gen", "expressions.jl"))
include("functions.jl")
include("primitives.jl")
include("milkyway.jl")

end # module GalacticPotentials