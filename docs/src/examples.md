# Example Usage

First, let's use everyone's favorite toy potential: the Plummer potential.

```julia
julia> using GalacticPotentials

julia> field = PlummerPotential()
Model PlummerPotential
States (3):
  x(t)
  y(t)
  z(t)
Parameters (3):
  G
  m
  b
```

The Plummer potential field equation is shown below.

$$\Phi = - \frac{G m}{\sqrt{b^{2} + x^{2} + y^{2} + z^{2}}}$$

Let's assume some massless particle which exists within this field. How will
the particle move? As [previously](potentials.md) described, we can take the
gradient of the scalar field $\Phi$ with respect to the state variables $x$,
$y$, and $z$ to find the force (per unit mass) applied to the particle at all
points in space.

```julia
julia> system = ODESystem(field)
Model PlummerPotential with 6 equations
States (6):
  x(t)
  y(t)
  z(t)
  ẋ(t)
  ẏ(t)
  ż(t)
Parameters (3):
  G
  m
  b
```

The differential equations which define the `system` variable are shown below.
Note that the gradient of the scalar potential field has been integrated into a
system of first-order differential equations: the expressions for the gradient
are shown in the state equations for $\frac{d \dot{x}}{d t}$,
$\frac{d \dot{y}}{d t}$, and $\frac{d \dot{z}}{d t}$.

$\begin{aligned}
\frac{\mathrm{d} x\left( t \right)}{\mathrm{d}t} =& \textnormal{\.{x}}\left( t \right) \\
\frac{\mathrm{d} y\left( t \right)}{\mathrm{d}t} =& \textnormal{\.{y}}\left( t \right) \\
\frac{\mathrm{d} z\left( t \right)}{\mathrm{d}t} =& \textnormal{\.{z}}\left( t \right) \\
\frac{\mathrm{d} \textnormal{\.{x}}\left( t \right)}{\mathrm{d}t} =& \frac{ - G m x\left( t \right)}{\left( \sqrt{b^{2} + \left( x\left( t \right) \right)^{2} + \left( y\left( t \right) \right)^{2} + \left( z\left( t \right) \right)^{2}} \right)^{3}} \\
\frac{\mathrm{d} \textnormal{\.{y}}\left( t \right)}{\mathrm{d}t} =& \frac{ - G m y\left( t \right)}{\left( \sqrt{b^{2} + \left( x\left( t \right) \right)^{2} + \left( y\left( t \right) \right)^{2} + \left( z\left( t \right) \right)^{2}} \right)^{3}} \\
\frac{\mathrm{d} \textnormal{\.{z}}\left( t \right)}{\mathrm{d}t} =& \frac{ - G m z\left( t \right)}{\left( \sqrt{b^{2} + \left( x\left( t \right) \right)^{2} + \left( y\left( t \right) \right)^{2} + \left( z\left( t \right) \right)^{2}} \right)^{3}}
\end{aligned}$

The `ModelingToolkit.jl` `AbstractSystem` interface methods are defined for
all potential fields within `GalacticPotentials.jl`. Specifically, all fields
in `GalacticPotentials.jl` are subtypes of `AbstractTimeDependentSystem`.
Special subtype and method implementations have been added to
`GalacticPotentials.jl` as needed.

```julia
julia> using ModelingToolkit

julia> G = calculate_gradient(field)
3-element Vector{Num}:
 (-((-G*m) / (sqrt(b^2 + x(t)^2 + y(t)^2 + z(t)^2)^2))*x(t)) / sqrt(b^2 + x(t)^2 + y(t)^2 + z(t)^2)
 (-((-G*m) / (sqrt(b^2 + x(t)^2 + y(t)^2 + z(t)^2)^2))*y(t)) / sqrt(b^2 + x(t)^2 + y(t)^2 + z(t)^2)
 (-((-G*m) / (sqrt(b^2 + x(t)^2 + y(t)^2 + z(t)^2)^2))*z(t)) / sqrt(b^2 + x(t)^2 + y(t)^2 + z(t)^2)

julia> J = calculate_jacobian(system)
6×6 Matrix{Num}:
  … # expression left out of documentation for brevity

julia> f = generate_function(field)
:(function (ˍ₋arg1, ˍ₋arg2)
      begin
          (/)((*)((*)(-1, ˍ₋arg2[1]), ˍ₋arg2[2]), (sqrt)((+)((+)((+)((^)(ˍ₋arg2[3], 2), (^)(ˍ₋arg1[1], 2)), (^)(ˍ₋arg1[2], 2)), (^)(ˍ₋arg1[3], 2))))
      end
  end)

```

Special constructors for `ODESystem` and `ODEProblem` -- two `SciML` types --
are defined for all potential fields within `GalacticPotentials.jl`. The
`ODESystem` constructor was already illustrated above. Let's look at the
`ODEProblem` constructor now.

```julia
julia> problem = let
  p = randn(3)
  u0 = randn(6)
  ts = randn(2)

  ODEProblem(field, u0, ts, p)
end
```

It's generally safer to use _variable maps_ to provide initial conditions for
your `ODEProblem`. Variable maps allow for an arbitrary state vector order; the
`ODEProblem` call above assumes the parameter and state vector orders!

```julia
julia> problem = let model = system

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
end
```

With the `ODEProblem` defined, you can use `OrdinaryDiffEq.jl` or
`DifferentialEquations.jl` to numerically integrate the orbit, and `Plots.jl`
to plot the result! For more information, consult the `SciML`
[documentation](https://docs.sciml.ai), or the `GalacticPotentials.jl`
[Getting Started](index.md) page.
