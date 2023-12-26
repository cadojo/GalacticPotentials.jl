# `GalacticPotentials.jl`

_Common models within galactic dynamics!_

!!! warning
    This package is under development, and is pre-alpha.

## Overview

This package extends `ModelingToolkit` to represent common galactic
potentials. All available potentials are shown on the [Reference](reference.md) 
page.

## Installation

Choose one of the two lines below!

```julia
Pkg.add("GalacticPotentials.jl")    # in Julia code
```

```julia
pkg> GalacticPotentials             # in Julia's REPL
```

## Usage

This package is intended to be used alongside `ModelingToolkit.jl` and the rest
of the [`SciML`](https://sciml.ai) ecosystem. Scalar potential fields within 
`gala` and `galpy` --- two popular Python packages for galactic dynamics --- 
were used to bootstrap this package. All available potential fields are shown
on the [Reference](reference.md) page.