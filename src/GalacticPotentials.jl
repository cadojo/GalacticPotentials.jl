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

using Symbolics, SymbolicUtils
using LaTeXStrings
using LinearAlgebra
using ForwardDiff
using ModelingToolkit
using SciMLBase

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
    # PowerLawCutoffPotential,
    # SatohPotential,
    # StonePotential,
    states,
    parameters,
    ODESystem,
    ODEProblem

include(joinpath(@__DIR__, "gen", "expressions.jl"))
include("generic.jl")
include("potentials.jl")
include("odes.jl")

end # module GalacticPotentials