module GalacticPotentials

export
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

using LaTeXStrings
using ModelingToolkit
using LinearAlgebra

include(joinpath(@__DIR__, "gen", "expressions.jl"))

"""
The potential of a harmonic oscillator.

    $(expressions["HarmonicOscillatorPotential"])

# Extended Help

## References

[1] http://gala.adrian.pw/en/latest/potential/index.html
"""
function HarmonicOscillatorPotential(N::Integer; name=:HarmonicOscillatorPotential)
    @parameters ω[1:N] [description = "The frequencies of the harmonic oscillator"]
    @variables x[1:N] [input = true, description = "The N-dimensional state of the oscillator"]
    @variables V [output = true, description = "The potential of the oscillator."]

    return NonlinearSystem(
        [V ~ sum(1 // 2 * ωᵢ^2 * xᵢ^2 for (ωᵢ, xᵢ) in zip(ω, x))],
        collect(x),
        collect(ω);
        name=name
    )
end

end # module GalacticPotentials
