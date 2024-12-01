function ScalarField(value, t, u, p; name, simplify = true, kwargs...)
    @variables Φ(t)

    Δ = Differential(t)

    du = map(
        x -> Symbol(:Δ, Symbol(first(split(string(x), "($(Symbolics.value(t)))")))),
        u
    )

    Δu = getfield.(
        vcat((@variables($(δ)(t)) for δ in du)...),
        :val
    )

    eqs = vcat(
        Δ.(u) .~ Δu,
        Δ.(Δu) .~ -ModelingToolkit.gradient(value, u, simplify = simplify)
    )
    return ODESystem(vcat(eqs, Φ ~ value), t, vcat(u, Δu, Φ), p; name = name, kwargs...)
end

"""
The potential due to a harmonic oscillator.

\$$(LATEX_EXPRESSIONS["HarmonicOscillatorPotential"])\$
"""
function HarmonicOscillatorPotential(
        N::Integer = 1; name = :HarmonicOscillator, kwargs...)
    if N > 1
        @independent_variables t
        @variables (τ(t))[1:N]
        @parameters ω[1:N]

        τ = collect(τ)
        ω = collect(ω)

        value = (1 // 2) * ω ⋅ ω * τ ⋅ τ
    elseif N == 1
        @independent_variables t
        @variables τ(t)
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
    @independent_variables t
    @variables x(t) y(t)

    value = x^2 * y + (1 // 2) * x^2 - (1 // 3) * y^3 + (1 // 2)y^2

    return ScalarField(value, t, [x, y], Num[]; name = name, kwargs...)
end

"""
The Hernquist potential.

\$$(LATEX_EXPRESSIONS["HernquistPotential"])\$
"""
function HernquistPotential(; name = :HernquistPotential, kwargs...)
    @independent_variables t
    @variables x(t) y(t) z(t)
    @parameters G m c

    value = -(G * m) / (c + sqrt(x^2 + y^2 + z^2))
    return ScalarField(value, t, [x, y, z], [G, m, c]; name = name, kwargs...)
end

"""
The Isochrone potential.

\$$(LATEX_EXPRESSIONS["IsochronePotential"])\$
"""
function IsochronePotential(; name = :IsochronePotential, kwargs...)
    @independent_variables t
    @variables x(t) y(t) z(t)
    @parameters G m b

    value = -(G * m) / (b + sqrt(b^2 + x^2 + y^2 + z^2))
    return ScalarField(value, t, [x, y, z], [G, m, b]; name = name, kwargs...)
end

"""
The Jaffe potential.

\$$(LATEX_EXPRESSIONS["JaffePotential"])\$
"""
function JaffePotential(; name = :JaffePotential, kwargs...)
    @independent_variables t
    @variables x(t) y(t) z(t)
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
    @independent_variables t
    @variables x(t) y(t) z(t)
    @parameters G m

    value = -G * m / sqrt(x^2 + y^2 + z^2)
    return ScalarField(value, t, [x, y, z], [G, m]; name = name, kwargs...)
end

"""
The Kuzmin potential.

\$$(LATEX_EXPRESSIONS["KuzminPotential"])\$
"""
function KuzminPotential(; name = :KuzminPotential, kwargs...)
    @independent_variables t
    @variables x(t) y(t) z(t)
    @parameters G m a

    value = -(G * m) / sqrt(x^2 + y^2 + (a + abs(z))^2)
    return ScalarField(value, t, [x, y, z], [G, m, a]; name = name, kwargs...)
end

"""
The logarithmic potential.

\$$(LATEX_EXPRESSIONS["LogarithmicPotential"])\$
"""
function LogarithmicPotential(; name = :LogarithmicPotential, kwargs...)
    @independent_variables t
    @variables x(t) y(t) z(t)
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
    @independent_variables t
    u = @variables x(t) y(t) z(t)
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

    return ScalarField(value, t, u, p; name = name, kwargs...)
end

"""
The Miyamoto-Nagai potential.

\$$(LATEX_EXPRESSIONS["MiyamotoNagaiPotential"])\$
"""
function MiyamotoNagaiPotential(; name = :MiyamotoNagaiPotential, kwargs...)
    @independent_variables t
    u = @variables x(t) y(t) z(t)
    p = @parameters G m a b

    value = -G * m / sqrt(x^2 + y^2 + (a + sqrt(b^2 + z^2))^2)

    return ScalarField(value, t, u, p; name = name, kwargs...)
end

"""
The NFW potential.

\$$(LATEX_EXPRESSIONS["NFWPotential"])\$
"""
function NFWPotential(; name = :NFWPotential, kwargs...)
    @independent_variables t
    u = @variables x(t) y(t) z(t)
    p = @parameters G m a b c r

    value = -G * m * log10(1 + sqrt(z^2 / c^2 + y^2 / b^2 + x^2 / a^2) / r) /
            sqrt(z^2 / c^2 + y^2 / b^2 + x^2 / a^2)

    return ScalarField(value, t, u, p; name = name, kwargs...)
end

"""
The Plummer potential.

\$$(LATEX_EXPRESSIONS["PlummerPotential"])\$
"""
function PlummerPotential(; name = :PlummerPotential, kwargs...)
    @independent_variables t
    u = @variables x(t) y(t) z(t)
    p = @parameters G m b

    value = -G * m / sqrt(b^2 + x^2 + y^2 + z^2)
    return ScalarField(value, t, u, p; name = name, kwargs...)
end

"""
The power-law cutoff potential.

!!! warning
    Not yet implemented!

\$$(LATEX_EXPRESSIONS["PowerLawCutoffPotential"])\$
"""
function PowerLawCutoffPotential(; name = :PowerLawCutoffPotential, kwargs...)
    @independent_variables t
    u = @variables x(t) y(t) z(t)
    p = @parameters G m a α c

    value = G * α * m * lowergamma(3 // 2 - α // 2, (x^2 + y^2 + z^2) / c^2) /
            (2 * sqrt(x^2 + y^2 + z^2) * SpecialFunctions.gamma(5 // 2 - α // 2)) -
            3 * G * m * lowergamma(3 // 2 - α // 2, (x^2 + y^2 + z^2) / c^2) /
            (2 * sqrt(x^2 + y^2 + z^2) * SpecialFunctions.gamma(5 // 2 - α // 2)) +
            G * m * lowergamma(1 - α // 2, (x^2 + y^2 + z^2) / c^2) /
            (c * SpecialFunctions.gamma(3 // 2 - α // 2))

    return ScalarField(value, t, u, p; name = name, kwargs...)
end

"""
The Satoh potential.

!!! warning
    Not yet implemented!

\$$(LATEX_EXPRESSIONS["SatohPotential"])\$
"""
function SatohPotential(; name = :SatohPotential, kwargs...)
    @independent_variables t
    u = @variables x(t) y(t) z(t)
    p = @parameters G m a b

    value = -G * m / sqrt(a * (a + 2 * sqrt(b^2 + z^2)) + x^2 + y^2 + z^2)

    return ScalarField(value, t, u, p; name = name, kwargs...)
end

"""
The StonePotential potential.

!!! warning
    Not yet implemented!

\$$(LATEX_EXPRESSIONS["StonePotential"])\$
"""
function StonePotential(; name = :StonePotential, kwargs...)
    @independent_variables t
    u = @variables x(t) y(t) z(t)
    p = @parameters G m rᵪ rₕ

    value = -2 * G * m *
            (-rᵪ * atan(sqrt(x^2 + y^2 + z^2) / rᵪ) / sqrt(x^2 + y^2 + z^2) +
             rₕ * atan(sqrt(x^2 + y^2 + z^2) / rₕ) / sqrt(x^2 + y^2 + z^2) +
             0.5 * log((rₕ^2 + x^2 + y^2 + z^2) / (rᵪ^2 + x^2 + y^2 + z^2))) /
            (-π * rᵪ + π * rₕ)

    return ScalarField(value, t, u, p; name = name, kwargs...)
end

"""
Galactic potentials for our home galaxy: the Milky Way.
"""
module MilkyWay

using GalacticPotentials
using ModelingToolkit
using Memoize

"""
A potential field for the Milky Way galaxy, based off of Dr. Bovy's 2015 paper.
"""
function Bovy2014(; name = :BovyMilkyWayPotential, kwargs...)

    # default_disk = dict(m=68193902782.346756 * u.Msun, a=3.0 * u.kpc, b=280 * u.pc)
    # default_bulge = dict(m=4501365375.06545 * u.Msun, alpha=1.8, r_c=1.9 * u.kpc)
    # default_halo = dict(m=4.3683325e11 * u.Msun, r_s=16 * u.kpc)

    @independent_variables t
    @variables Φ(t) x(t) y(t) z(t) Δx(t) Δy(t) Δz(t)
    disk = MiyamotoNagaiPotential(; name = :Disk)
    bulge = PowerLawCutoffPotential(; name = :Bulge)
    halo = NFWPotential(; name = :Halo)

    Δ = Differential(t)

    aliases = [
        x => disk.x,
        y => disk.y,
        z => disk.z]

    eqs = vcat(
        Φ ~ disk.Φ + bulge.Φ + halo.Φ,
        Δ(Δx) ~ Δ(@namespace(disk.Δx)) + Δ(@namespace(bulge.Δx)) + Δ(@namespace(halo.Δx)),
        Δ(Δy) ~ Δ(@namespace(disk.Δy)) + Δ(@namespace(bulge.Δy)) + Δ(@namespace(halo.Δy)),
        Δ(Δz) ~ Δ(@namespace(disk.Δz)) + Δ(@namespace(bulge.Δz)) + Δ(@namespace(halo.Δz)),
        [alias.first ~ alias.second for alias in aliases]
    )

    return compose(
        ODESystem(eqs, t, [x, y, z, Δx, Δy, Δz, Φ],
            vcat(parameters(disk), parameters(bulge), parameters(halo)); name = :MilkyWay),
        disk, bulge, halo
    )
end

end