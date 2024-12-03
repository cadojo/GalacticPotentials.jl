function ScalarField(value, t, u, p; name, gradient = true, simplify = true, kwargs...)
    if gradient
        symbols = [Symbol(:∂Φ∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂Φ∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        eqs = vcat(Φ ~ value, ∂Φ∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    else
        eqs = [Φ ~ value]
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The potential due to a harmonic oscillator.

\$$(LATEX_EXPRESSIONS["HarmonicOscillatorPotential"])\$
"""
function HarmonicOscillatorPotential(
        N::Integer = 1; name = :HarmonicOscillator, kwargs...)
    if N > 1
        @variables (τ(t))[1:N] [input = true]
        @parameters ω[1:N]

        τ = collect(τ)
        ω = collect(ω)

        value = (1 // 2) * ω ⋅ ω * τ ⋅ τ
    elseif N == 1
        @variables τ(t) input=true
        @parameters ω

        value = (1 // 2) * ω^2 * τ^2

        τ = [τ]
        ω = [ω]
    else
        error("Argument `$N` must be greater than zero!")
    end

    return ScalarField(value, t, τ, ω; name = name, kwargs...)
end

"""
The Henon-Heiles potential.

\$$(LATEX_EXPRESSIONS["HenonHeilesPotential"])\$
"""
function HenonHeilesPotential(; name = :HenonHeilesPotential, kwargs...)
    value = x^2 * y + (1 // 2) * x^2 - (1 // 3) * y^3 + (1 // 2)y^2

    return ScalarField(value, t, [x, y], Num[]; name = name, kwargs...)
end

"""
The Hernquist potential.

\$$(LATEX_EXPRESSIONS["HernquistPotential"])\$
"""
function HernquistPotential(; name = :HernquistPotential, kwargs...)
    @parameters G m c

    value = -(G * m) / (c + sqrt(x^2 + y^2 + z^2))
    return ScalarField(value, t, [x, y, z], [G, m, c]; name = name, kwargs...)
end

"""
The Isochrone potential.

\$$(LATEX_EXPRESSIONS["IsochronePotential"])\$
"""
function IsochronePotential(; name = :IsochronePotential, kwargs...)
    @parameters G m b

    value = -(G * m) / (b + sqrt(b^2 + x^2 + y^2 + z^2))
    return ScalarField(value, t, [x, y, z], [G, m, b]; name = name, kwargs...)
end

"""
The Jaffe potential.

\$$(LATEX_EXPRESSIONS["JaffePotential"])\$
"""
function JaffePotential(; name = :JaffePotential, kwargs...)
    @parameters G m c

    value = G * m * log10(
                sqrt(x^2 + y^2 + z^2) / (c + sqrt(x^2 + y^2 + z^2)) / c
            )

    return ScalarField(value, t, [x, y, z], [G, m, c]; name = name, kwargs...)
end

"""
The Kepler potential.

\$$(LATEX_EXPRESSIONS["KeplerPotential"])\$
"""
function KeplerPotential(; name = :KeplerPotential, kwargs...)
    @parameters G m

    value = -G * m / sqrt(x^2 + y^2 + z^2)
    return ScalarField(value, t, [x, y, z], [G, m]; name = name, kwargs...)
end

"""
The Kuzmin potential.

\$$(LATEX_EXPRESSIONS["KuzminPotential"])\$
"""
function KuzminPotential(; name = :KuzminPotential, kwargs...)
    @parameters G m a

    value = -(G * m) / sqrt(x^2 + y^2 + (a + abs(z))^2)
    return ScalarField(value, t, [x, y, z], [G, m, a]; name = name, kwargs...)
end

"""
The logarithmic potential.

\$$(LATEX_EXPRESSIONS["LogarithmicPotential"])\$
"""
function LogarithmicPotential(; name = :LogarithmicPotential, kwargs...)
    @parameters v r q[1:3]

    q = collect(q)

    value = (1 // 2) * v^2 * log10(r^2 + z^2 / q[3]^2 + y^2 / q[2]^2 + x^2 / q[1]^2)

    return ScalarField(value, t, [x, y, z], vcat(v, r, q); name = name, kwargs...)
end

"""
The long Murali-bar potential.

\$$(LATEX_EXPRESSIONS["LongMuraliBarPotential"])\$
"""
function LongMuraliBarPotential(; name = :LongMuraliBarPotential, kwargs...)
    p = @parameters G m a b c α

    value = G * m *
            log10(
                (-a + x * cos(α) + y * sin(α) +
                 sqrt((b + sqrt(c^2 + z^2))^2 + (-x * sin(α) + y * cos(α))^2 +
                      (a - x * cos(α) - y * sin(α))^2)) / (
                a + x * cos(α) + y * sin(α) +
                sqrt((b + sqrt(c^2 + z^2))^2 + (-x * sin(α) + y * cos(α))^2 +
                     (a + x * cos(α) + y * sin(α))^2)
            )
            ) / 2a

    return ScalarField(value, t, [x, y, z], p; name = name, kwargs...)
end

"""
The Miyamoto-Nagai potential.

\$$(LATEX_EXPRESSIONS["MiyamotoNagaiPotential"])\$
"""
function MiyamotoNagaiPotential(; name = :MiyamotoNagaiPotential, kwargs...)
    p = @parameters G m a b

    value = -G * m / sqrt(x^2 + y^2 + (a + sqrt(b^2 + z^2))^2)

    return ScalarField(value, t, [x, y, z], p; name = name, kwargs...)
end

"""
The NFW potential.

\$$(LATEX_EXPRESSIONS["NFWPotential"])\$
"""
function NFWPotential(; name = :NFWPotential, kwargs...)
    p = @parameters G m a b c r

    value = -G * m * log10(1 + sqrt(z^2 / c^2 + y^2 / b^2 + x^2 / a^2) / r) /
            sqrt(z^2 / c^2 + y^2 / b^2 + x^2 / a^2)

    return ScalarField(value, t, [x, y, z], p; name = name, kwargs...)
end

"""
The Plummer potential.

\$$(LATEX_EXPRESSIONS["PlummerPotential"])\$
"""
function PlummerPotential(; name = :PlummerPotential, kwargs...)
    p = @parameters G m b

    value = -G * m / sqrt(b^2 + x^2 + y^2 + z^2)
    return ScalarField(value, t, [x, y, z], p; name = name, kwargs...)
end

"""
The power-law cutoff potential.

!!! warning
    Not yet implemented!

\$$(LATEX_EXPRESSIONS["PowerLawCutoffPotential"])\$
"""
function PowerLawCutoffPotential(; name = :PowerLawCutoffPotential, kwargs...)
    p = @parameters G m a α c

    value = G * α * m * lowergamma(3 // 2 - α // 2, (x^2 + y^2 + z^2) / c^2) /
            (2 * sqrt(x^2 + y^2 + z^2) * SpecialFunctions.gamma(5 // 2 - α // 2)) -
            3 * G * m * lowergamma(3 // 2 - α // 2, (x^2 + y^2 + z^2) / c^2) /
            (2 * sqrt(x^2 + y^2 + z^2) * SpecialFunctions.gamma(5 // 2 - α // 2)) +
            G * m * lowergamma(1 - α // 2, (x^2 + y^2 + z^2) / c^2) /
            (c * SpecialFunctions.gamma(3 // 2 - α // 2))

    return ScalarField(value, t, [x, y, z], p; name = name, kwargs...)
end

"""
The Satoh potential.

!!! warning
    Not yet implemented!

\$$(LATEX_EXPRESSIONS["SatohPotential"])\$
"""
function SatohPotential(; name = :SatohPotential, kwargs...)
    p = @parameters G m a b

    value = -G * m / sqrt(a * (a + 2 * sqrt(b^2 + z^2)) + x^2 + y^2 + z^2)

    return ScalarField(value, t, [x, y, z], p; name = name, kwargs...)
end

"""
The StonePotential potential.

!!! warning
    Not yet implemented!

\$$(LATEX_EXPRESSIONS["StonePotential"])\$
"""
function StonePotential(; name = :StonePotential, kwargs...)
    p = @parameters G m rᵪ rₕ

    value = -2 * G * m *
            (-rᵪ * atan(sqrt(x^2 + y^2 + z^2) / rᵪ) / sqrt(x^2 + y^2 + z^2) +
             rₕ * atan(sqrt(x^2 + y^2 + z^2) / rₕ) / sqrt(x^2 + y^2 + z^2) +
             0.5 * log((rₕ^2 + x^2 + y^2 + z^2) / (rᵪ^2 + x^2 + y^2 + z^2))) /
            (-π * rᵪ + π * rₕ)

    return ScalarField(value, t, [x, y, z], p; name = name, kwargs...)
end

"""
Galactic potentials for our home galaxy: the Milky Way.
"""
module MilkyWay

using GalacticPotentials
using GalacticPotentials.Symbols
using ModelingToolkit
using Memoize

"""
A potential field for the Milky Way galaxy, based off of Dr. Bovy's 2015 paper.
"""
function Bovy2014(; name = :BovyMilkyWayPotential, kwargs...)

    # gala 
    # default_disk = dict(m=68193902782.346756 * u.Msun, a=3.0 * u.kpc, b=280 * u.pc)
    # default_bulge = dict(m=4501365375.06545 * u.Msun, alpha=1.8, r_c=1.9 * u.kpc)
    # default_halo = dict(m=4.3683325e11 * u.Msun, r_s=16 * u.kpc)

    # galpy
    # bp= PowerSphericalPotentialwCutoff(alpha=1.8,rc=1.9/8.,normalize=0.05)
    # mp= MiyamotoNagaiPotential(a=3./8.,b=0.28/8.,normalize=.6)
    # np= NFWPotential(a=16/8.,normalize=.35)
    # MWPotential2014= bp+mp+np

    disk = MiyamotoNagaiPotential(; name = :Disk)
    bulge = PowerLawCutoffPotential(; name = :Bulge)
    halo = NFWPotential(; name = :Halo)

    G = 6.6743e-11

    defaults = [
        disk.a => 3,
        disk.b => 280,
        disk.G => G,
        disk.m => 68193902782.346756,
        bulge.c => 1.9,
        bulge.m => 4501365375.06545,
        bulge.α => 1.8,
        bulge.G => G,
        halo.a => 1.0,
        halo.b => 1.0,
        halo.c => 1.0,
        halo.r => 16,
        halo.m => 4.3683325e11,
        halo.G => G
    ]

    aliases = [
        disk.x => x,
        disk.y => y,
        disk.z => z,
        bulge.x => x,
        bulge.y => y,
        bulge.z => z,
        halo.x => x,
        halo.y => y,
        halo.z => z
    ]

    u = [x, y, z]
    du = [ẋ, ẏ, ż]

    grad(sys) = [sys.∂Φ∂x, sys.∂Φ∂y, sys.∂Φ∂z]

    eqs = vcat(
        Φ ~ disk.Φ + bulge.Φ + halo.Φ,
        D.(u) .~ du,
        D.(du) .~ -(grad(disk) .+ grad(bulge) .+ grad(halo)),
        [alias.first ~ alias.second for alias in aliases]
    )

    return compose(
        ODESystem(
            eqs, t;
            name = name, defaults = Dict(vcat(defaults, aliases))),
        disk, bulge, halo
    )
end

end