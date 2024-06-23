"""
The potential due to a harmonic oscillator.

\$$(LATEX_EXPRESSIONS["HarmonicOscillatorPotential"])\$
"""
@memoize function HarmonicOscillatorPotential(
        N::Integer = 1; name = :HarmonicOscillator, kwargs...)
    if N > 1
        @variables t (x(t))[1:N]
        @parameters ω[1:N]

        x = collect(x)
        ω = collect(ω)

        value = (1 // 2) * ω ⋅ ω * x ⋅ x
    elseif N == 1
        @variables t x(t)
        @parameters ω

        value = (1 // 2) * ω^2 * x^2

        x = [x]
        ω = [ω]
    else
        error("Argument `$N` must be greater than zero!")
    end

    return ScalarField(value, t, x, ω; name = name, kwargs...)
end

"""
The Henon-Heiles potential.

\$$(LATEX_EXPRESSIONS["HenonHeilesPotential"])\$
"""
@memoize function HenonHeilesPotential(; name = :HenonHeilesPotential, kwargs...)
    @variables t x(t) y(t)

    value = x^2 * y + (1 // 2) * x^2 - (1 // 3) * y^3 + (1 // 2)y^2

    return ScalarField(value, t, [x, y], Num[]; name = name, kwargs...)
end

"""
The Hernquist potential.

\$$(LATEX_EXPRESSIONS["HernquistPotential"])\$
"""
@memoize function HernquistPotential(; name = :HernquistPotential, kwargs...)
    @variables t x(t) y(t) z(t)
    @parameters G m c

    value = -(G * m) / (c + sqrt(x^2 + y^2 + z^2))
    return ScalarField(value, t, [x, y, z], [G, m, c]; name = name, kwargs...)
end

"""
The Isochrone potential.

\$$(LATEX_EXPRESSIONS["IsochronePotential"])\$
"""
@memoize function IsochronePotential(; name = :IsochronePotential, kwargs...)
    @variables t x(t) y(t) z(t)
    @parameters G m b

    value = -(G * m) / (b + sqrt(b^2 + x^2 + y^2 + z^2))
    return ScalarField(value, t, [x, y, z], [G, m, b]; name = name, kwargs...)
end

"""
The Jaffe potential.

\$$(LATEX_EXPRESSIONS["JaffePotential"])\$
"""
@memoize function JaffePotential(; name = :JaffePotential, kwargs...)
    @variables t x(t) y(t) z(t)
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
@memoize function KeplerPotential(; name = :KeplerPotential, kwargs...)
    @variables t x(t) y(t) z(t)
    @parameters G m

    value = -G * m / sqrt(x^2 + y^2 + z^2)
    return ScalarField(value, t, [x, y, z], [G, m]; name = name, kwargs...)
end

"""
The Kuzmin potential.

\$$(LATEX_EXPRESSIONS["KuzminPotential"])\$
"""
@memoize function KuzminPotential(; name = :KuzminPotential, kwargs...)
    @variables t x(t) y(t) z(t)
    @parameters G m a

    value = -(G * m) / sqrt(x^2 + y^2 + (a + abs(z))^2)
    return ScalarField(value, t, [x, y, z], [G, m, a]; name = name, kwargs...)
end

"""
The logarithmic potential.

\$$(LATEX_EXPRESSIONS["LogarithmicPotential"])\$
"""
@memoize function LogarithmicPotential(; name = :LogarithmicPotential, kwargs...)
    @variables t x(t) y(t) z(t)
    @parameters v r q[1:3]

    q = collect(q)

    value = (1 // 2) * v^2 * log10(r^2 + z^2 / q[3]^2 + y^2 / q[2]^2 + x^2 / q[1]^2)

    return ScalarField(value, t, [x, y, z], vcat(v, r, q); name = name, kwargs...)
end

"""
The long Murali-bar potential.

\$$(LATEX_EXPRESSIONS["LongMuraliBarPotential"])\$
"""
@memoize function LongMuraliBarPotential(; name = :LongMuraliBarPotential, kwargs...)
    @variables t
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
@memoize function MiyamotoNagaiPotential(; name = :MiyamotoNagaiPotential, kwargs...)
    @variables t
    u = @variables x(t) y(t) z(t)
    p = @parameters G m a b

    value = -G * m / sqrt(x^2 + y^2 + (a + sqrt(b^2 + z^2))^2)

    return ScalarField(value, t, u, p; name = name, kwargs...)
end

"""
The NFW potential.

\$$(LATEX_EXPRESSIONS["NFWPotential"])\$
"""
@memoize function NFWPotential(; name = :NFWPotential, kwargs...)
    @variables t
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
@memoize function PlummerPotential(; name = :PlummerPotential, kwargs...)
    @variables t
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
@memoize function PowerLawCutoffPotential(; name = :PowerLawCutoffPotential, kwargs...)
    @variables t
    u = @variables x(t) y(t) z(t)
    p = @parameters G m a α c

    value = G * α * m * lowergamma(3 // 2 - α // 2, (x^2 + y^2 + z^2) / c^2) /
            (2 * sqrt(x^2 + y^2 + z^2) * gamma(5 // 2 - α // 2)) -
            3 * G * m * lowergamma(3 // 2 - α // 2, (x^2 + y^2 + z^2) / c^2) /
            (2 * sqrt(x^2 + y^2 + z^2) * gamma(5 // 2 - α // 2)) +
            G * m * lowergamma(1 - α // 2, (x^2 + y^2 + z^2) / c^2) /
            (c * gamma(3 // 2 - α // 2))

    return ScalarField(value, t, u, p; name = name, kwargs...)
end

"""
The Satoh potential.

!!! warning
    Not yet implemented!

\$$(LATEX_EXPRESSIONS["SatohPotential"])\$
"""
@memoize function SatohPotential(; name = :SatohPotential, kwargs...)
    @variables t
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
@memoize function StonePotential(; name = :StonePotential, kwargs...)
    @variables t
    u = @variables x(t) y(t) z(t)
    p = @parameters G m rᵪ rₕ

    value = -2 * G * m *
            (-rᵪ * atan(sqrt(x^2 + y^2 + z^2) / rᵪ) / sqrt(x^2 + y^2 + z^2) +
             rₕ * atan(sqrt(x^2 + y^2 + z^2) / rₕ) / sqrt(x^2 + y^2 + z^2) +
             0.5 * log((rₕ^2 + x^2 + y^2 + z^2) / (rᵪ^2 + x^2 + y^2 + z^2))) /
            (-π * rᵪ + π * rₕ)

    return ScalarField(value, t, u, p; name = name, kwargs...)
end