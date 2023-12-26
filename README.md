[![Tests](https://github.com/cadojo/GalacticPotentials.jl/workflows/Tests/badge.svg)](https://github.com/cadojo/GalacticPotentials.jl/actions?query=workflow%3ATests)
[![Docs](https://github.com/cadojo/GalacticPotentials.jl/workflows/Documentation/badge.svg)](https://cadojo.github.io/GalacticPotentials.jl)

# 🌌 `GalacticPotentials.jl`

_An extension of
[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) which provides
common galactic potentials._

## Installation

Choose one of the two lines below!

```julia
Pkg.add("GalacticPotentials.jl")    # in Julia code
```

```julia
pkg> GalacticPotentials             # in Julia's REPL
```

## Usage

Potentials are defined as subtypes of `ModelingToolkit.jl` models. They can be 
converted to `ODESystem` and `ODEProblem` types to interact with the rest of the
`SciML` ecosystem.

```julia
using Plots
using ModelingToolkit
using DifferentialEquations
using GalacticPotentials

let model = ODESystem(PlummerPotential())

    p = @nonamespace Dict(
        model.G => 6.6743e-20, # field strength (km³ kg⁻¹ s⁻²)
        model.m => 6e31,       # mass (kg)
        model.b => 1e-6        # softening parameter (unitless)
    )

    u0 = @nonamespace Dict(
        model.x => 11e5, 
        model.y => 5e5, 
        model.z => 0,
        model.Δx => 1e3, 
        model.Δy => 1e3, 
        model.Δz => 0
    )

    ts = (0.0, 1e6)

    problem = ODEProblem(model, u0, ts, p)
    solution = solve(problem; reltol=1e-14, abstol=1e-14)

    plot(solution; idxs=(:x,:y), label=:none, dpi = 400, aspect_ratio=:equal)
end
```

![](/img/plummer-orbit.png)

## Credits

This package is [bootstrapped](/gen/gala.jl) off of [`gala`](http://gala.adrian.pw)
and [`galpy`](https://docs.galpy.org), two rich Python packages for galactic dynamics. 
I aim to learn about galactic dynamics by integrating the models within these two popular
Python packages into the [SciML](https://sciml.ai) ecosystem. 

The scalar field symbolic-numerics are copied and modified versions of `AbstractSystem`
interfaces within [`ModelingToolkit.jl`](https://github.com/sciml/ModelingToolkit.jl).
The field implementations in this package are highly unstable; they may change in the 
near future.