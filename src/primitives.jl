#
# Fundamental potential field equations.
#

"""
The potential due to a harmonic oscillator.

\$$(LATEX_EXPRESSIONS["HarmonicOscillatorPotential"])\$
"""
function HarmonicOscillatorPotential(
        N::Integer = 1; name = :HarmonicOscillator, partials = false,
        gradient = false, simplify = true, kwargs...)
    if N > 1
        @variables (τ(t))[1:N] [input = true]
        @variables (τ̇(t))[1:N] [input = true]

        @parameters ω[1:N]

        u = collect(τ)
        u̇ = collect(τ̇)
        p = collect(ω)

        value = (1 // 2) * ω ⋅ ω * τ ⋅ τ
    elseif N == 1
        @variables τ(t) [input = true]
        @variables τ̇(t) [input = true]

        @parameters ω

        value = (1 // 2) * ω^2 * τ^2

        u = [τ]
        u̇ = [τ̇]
        p = [ω]
    else
        error("Argument `$N` must be greater than zero!")
    end

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The Henon-Heiles potential.

\$$(LATEX_EXPRESSIONS["HenonHeilesPotential"])\$
"""
function HenonHeilesPotential(;
        name = :HenonHeilesPotential, partials = false, gradient = false, simplify = true, kwargs...)
    value = x^2 * y + (1 // 2) * x^2 - (1 // 3) * y^3 + (1 // 2)y^2

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The Hernquist potential.

\$$(LATEX_EXPRESSIONS["HernquistPotential"])\$
"""
function HernquistPotential(;
        name = :HernquistPotential, partials = false, gradient = false, simplify = true, kwargs...)
    @parameters G m c

    value = -(G * m) / (c + sqrt(x^2 + y^2 + z^2))

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The Isochrone potential.

\$$(LATEX_EXPRESSIONS["IsochronePotential"])\$
"""
function IsochronePotential(;
        name = :IsochronePotential, partials = false, gradient = false, simplify = true, kwargs...)
    @parameters G m b

    value = -(G * m) / (b + sqrt(b^2 + x^2 + y^2 + z^2))

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The Jaffe potential.

\$$(LATEX_EXPRESSIONS["JaffePotential"])\$
"""
function JaffePotential(;
        name = :JaffePotential, partials = false, gradient = false, simplify = true, kwargs...)
    @parameters G m c

    value = G * m * log10(
                sqrt(x^2 + y^2 + z^2) / (c + sqrt(x^2 + y^2 + z^2)) / c
            )

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The Kepler potential.

\$$(LATEX_EXPRESSIONS["KeplerPotential"])\$
"""
function KeplerPotential(;
        name = :KeplerPotential, partials = false, gradient = false, simplify = true, kwargs...)
    @parameters G m

    value = -G * m / sqrt(x^2 + y^2 + z^2)

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The Kuzmin potential.

\$$(LATEX_EXPRESSIONS["KuzminPotential"])\$
"""
function KuzminPotential(;
        name = :KuzminPotential, partials = false, gradient = false, simplify = true, kwargs...)
    @parameters G m a

    value = -(G * m) / sqrt(x^2 + y^2 + (a + abs(z))^2)

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The logarithmic potential.

\$$(LATEX_EXPRESSIONS["LogarithmicPotential"])\$
"""
function LogarithmicPotential(;
        name = :LogarithmicPotential, partials = false, gradient = false, simplify = true, kwargs...)
    @parameters v r q[1:3]

    q = collect(q)

    value = (1 // 2) * v^2 * log10(r^2 + z^2 / q[3]^2 + y^2 / q[2]^2 + x^2 / q[1]^2)

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The long Murali-bar potential.

\$$(LATEX_EXPRESSIONS["LongMuraliBarPotential"])\$
"""
function LongMuraliBarPotential(;
        name = :LongMuraliBarPotential, partials = false, gradient = false, simplify = true, kwargs...)
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

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The Miyamoto-Nagai potential.

\$$(LATEX_EXPRESSIONS["MiyamotoNagaiPotential"])\$
"""
function MiyamotoNagaiPotential(;
        name = :MiyamotoNagaiPotential, partials = false, gradient = false, simplify = true, kwargs...)
    p = @parameters G m a b

    value = -G * m / sqrt(x^2 + y^2 + (a + sqrt(b^2 + z^2))^2)

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The NFW potential.

\$$(LATEX_EXPRESSIONS["NFWPotential"])\$
"""
function NFWPotential(; name = :NFWPotential, partials = false,
        gradient = false, simplify = true, kwargs...)
    p = @parameters G m a b c r

    value = -G * m * log10(1 + sqrt(z^2 / c^2 + y^2 / b^2 + x^2 / a^2) / r) /
            sqrt(z^2 / c^2 + y^2 / b^2 + x^2 / a^2)

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The Plummer potential.

\$$(LATEX_EXPRESSIONS["PlummerPotential"])\$
"""
function PlummerPotential(;
        name = :PlummerPotential, partials = false, gradient = false, simplify = true, kwargs...)
    p = @parameters G m b

    value = -G * m / sqrt(b^2 + x^2 + y^2 + z^2)
    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The power-law cutoff potential.

\$$(LATEX_EXPRESSIONS["PowerLawCutoffPotential"])\$
"""
function PowerLawCutoffPotential(;
        name = :PowerLawCutoffPotential, partials = false, gradient = false, simplify = true, kwargs...)
    p = @parameters G m a α c

    value = G * α * m * lowergamma(3 // 2 - α // 2, (x^2 + y^2 + z^2) / c^2) /
            (2 * sqrt(x^2 + y^2 + z^2) * SpecialFunctions.gamma(5 // 2 - α // 2)) -
            3 * G * m * lowergamma(3 // 2 - α // 2, (x^2 + y^2 + z^2) / c^2) /
            (2 * sqrt(x^2 + y^2 + z^2) * SpecialFunctions.gamma(5 // 2 - α // 2)) +
            G * m * lowergamma(1 - α // 2, (x^2 + y^2 + z^2) / c^2) /
            (c * SpecialFunctions.gamma(3 // 2 - α // 2))

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The Satoh potential.

\$$(LATEX_EXPRESSIONS["SatohPotential"])\$
"""
function SatohPotential(;
        name = :SatohPotential, partials = false, gradient = false, simplify = true, kwargs...)
    p = @parameters G m a b

    value = -G * m / sqrt(a * (a + 2 * sqrt(b^2 + z^2)) + x^2 + y^2 + z^2)

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end

"""
The StonePotential potential.

\$$(LATEX_EXPRESSIONS["StonePotential"])\$
"""
function StonePotential(;
        name = :StonePotential, partials = false, gradient = false, simplify = true, kwargs...)
    p = @parameters G m rᵪ rₕ

    value = -2 * G * m *
            (-rᵪ * atan(sqrt(x^2 + y^2 + z^2) / rᵪ) / sqrt(x^2 + y^2 + z^2) +
             rₕ * atan(sqrt(x^2 + y^2 + z^2) / rₕ) / sqrt(x^2 + y^2 + z^2) +
             0.5 * log((rₕ^2 + x^2 + y^2 + z^2) / (rᵪ^2 + x^2 + y^2 + z^2))) /
            (-π * rᵪ + π * rₕ)

    eqs = [V ~ value]

    if partials || gradient
        symbols = [Symbol(:∂V∂, Symbol(first(split(string(x), "($(Symbolics.value(t)))"))))
                   for x in u]

        ∂V∂u = getfield.(
            vcat((@variables($(x)(t)) for x in symbols)...),
            :val
        )

        append!(eqs, ∂V∂u .~ Symbolics.gradient(value, u, simplify = simplify))
    end

    if gradient
        append!(eqs, D.(u) .~ u̇)
        append!(eqs, D.(u̇) .~ -∂V∂u)
    end

    return ODESystem(eqs, t; name = name, kwargs...)
end
