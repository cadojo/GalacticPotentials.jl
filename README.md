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

let model = PlummerPotential(gradient=true)

    p = @nonamespace Dict(
        model.G => 6.6743e-20, # field strength (km³ kg⁻¹ s⁻²)
        model.m => 6e31,       # mass (kg)
        model.b => 1e-6        # softening parameter (unitless)
    )

    u0 = @nonamespace Dict(
        model.x => 11e5,
        model.y => 5e5,
        model.z => 0,
        model.ẋ => 1e3,
        model.ẏ => 1e3,
        model.ż => 0
    )

    ts = (0.0, 1e6)

    problem = ODEProblem(model, u0, ts, p)
    solution = solve(problem; reltol=1e-14, abstol=1e-14)

    plot(solution; idxs=(:x,:y), label=:none, dpi = 400, aspect_ratio=:equal)
end
```

![](/docs/src/img/plummer-orbit.png)

## Credits

This package is [bootstrapped](/gen/gala.jl) off of [`gala`](http://gala.adrian.pw)
and [`galpy`](https://docs.galpy.org), two rich Python packages for galactic dynamics.
I aim to learn about galactic dynamics by integrating the models within these two popular
Python packages into the [SciML](https://sciml.ai) ecosystem. The Gala and GalPy
licenses are shown below.

<details>

<summary>Gala License</summary>

```
The MIT License (MIT)

Copyright (c) 2012-2024 Adrian M. Price-Whelan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

</details>

<details>

<summary>GalPy License</summary>

```
Copyright (c) 2010, Jo Bovy
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

</details>