[![Tests](https://github.com/cadojo/GalacticPotentials.jl/workflows/UnitTests/badge.svg)](https://github.com/cadojo/GalacticPotentials.jl/actions?query=workflow%3AUnitTests)
[![Docs](https://github.com/cadojo/GalacticPotentials.jl/workflows/Documentation/badge.svg)](https://cadojo.github.io/GalacticPotentials.jl)

# GalacticPotentials.jl

_An extension of
[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) which provides
common galactic potentials._

## Installation

Choose one of the two lines below!

```julia
Pkg.add("https://github.com/cadojo/GalacticPotentials.jl")  # in Julia code
```

```julia
pkg> add https://github.com/cadojo/GalacticPotentials.jl    # in Julia's REPL
```

## Credits

I am [bootstrapping](/gen/gala.jl) this package off of [`gala`](http://gala.adrian.pw)
and [`galpy`](https://docs.galpy.org), two rich Python packages for galactic dynamics. 
I aim to learn about galactic dynamics by integrating the models within these two popular
Python packages into the [SciML](https://sciml.ai) ecosystem.

Thank you the developers of `gala` and `galpy`!