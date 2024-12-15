#
# Milky Way potentials.
#

"""
A potential field for the Milky Way galaxy, based off of Dr. Bovy's 2015 paper.
"""
function Bovy2014(;
        name = :BovyMilkyWayPotential, partials = true, gradient = true, kwargs...)

    # gala 
    # default_disk = dict(m=68193902782.346756 * u.Msun, a=3.0 * u.kpc, b=280 * u.pc)
    # default_bulge = dict(m=4501365375.06545 * u.Msun, alpha=1.8, r_c=1.9 * u.kpc)
    # default_halo = dict(m=4.3683325e11 * u.Msun, r_s=16 * u.kpc)

    # galpy
    # bp= PowerSphericalPotentialwCutoff(alpha=1.8,rc=1.9/8.,normalize=0.05)
    # mp= MiyamotoNagaiPotential(a=3./8.,b=0.28/8.,normalize=.6)
    # np= NFWPotential(a=16/8.,normalize=.35)
    # MWPotential2014= bp+mp+np

    disk = MiyamotoNagaiPotential(; partials = partials, name = :Disk)
    bulge = PowerLawCutoffPotential(; partials = partials, name = :Bulge)
    halo = NFWPotential(; partials = partials, name = :Halo)

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
    u̇ = [ẋ, ẏ, ż]

    value = disk.V + bulge.V + halo.V
    eqs = [V ~ value]

    append!(eqs, [alias.first ~ alias.second for alias in aliases])

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

    return compose(
        ODESystem(
            eqs, t;
            name = name, defaults = Dict(vcat(defaults, aliases))),
        disk, bulge, halo
    )
end
