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

using Symbolics, SymbolicUtils
using LaTeXStrings
using LinearAlgebra
using ForwardDiff
using FiniteDiff
using SciMLBase
using Memoize

import SpecialFunctions
import HypergeometricFunctions

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
include("generic.jl")
include("potentials.jl")
include("odes.jl")

end # module GalacticPotentialsj