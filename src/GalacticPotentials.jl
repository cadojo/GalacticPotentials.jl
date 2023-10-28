module GalacticPotentials

export
    GalacticPotential,
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

using Symbolics, SymbolicUtils
using LaTeXStrings
using ModelingToolkit
using LinearAlgebra

include(joinpath(@__DIR__, "gen", "expressions.jl"))
include("generic.jl")

end # module GalacticPotentials