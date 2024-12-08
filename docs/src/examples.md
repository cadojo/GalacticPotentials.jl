# Example Usage

First, let's use everyone's favorite toy potential: the Plummer potential.

```@example example
using GalacticPotentials

field = PlummerPotential()
```

The Plummer potential field equation is shown below, where $V \triangleq V(u,p,t)$ 
is the scalar potential of the field at any point in space (and time).

$$V = - \frac{G m}{\sqrt{b^{2} + x^{2} + y^{2} + z^{2}}}$$

Let's assume some massless particle which exists within this field. How will
the particle move? As [previously](potentials.md) described, we can take the
gradient of the scalar field $\Phi$ with respect to the state variables $x$,
$y$, and $z$ to find the force (per unit mass) applied to the particle at all
points in space.

```@example example
system = PlummerPotential(gradient=true)
```

Note there are 10 equations! Most of these equations are observables, i.e. 
they can be reconstructed lazily from the solution. We can reduce the system
to 6 unknowns and 6 equations by calling `ModelingToolkit.structural_simplify`.

```@example example
using ModelingToolkit
system = structural_simplify(system)
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
all potential fields within `GalacticPotentials.jl` because they are all 
`ModelingToolkit.ODESystem` instances!

```@example example
J = ModelingToolkit.calculate_jacobian(system)
```

```@example example
f = ModelingToolkit.generate_function(system)
```

As with any `ModelingToolkit.ODESystem`, we can construct a numerical problem
by passing the system to an `ODEProblem` constructor.

```@example example
using OrdinaryDiffEq

problem = let
  p = randn(3)
  u0 = randn(6)
  ts = randn(2)

  ODEProblem(system, u0, ts, p)
end
```

It's generally safer to use _variable maps_ to provide initial conditions for
your `ODEProblem`. Variable maps allow for an arbitrary state vector order; the
`ODEProblem` call above assumes the parameter and state vector orders!

```@example example
problem = let model = system

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

  ODEProblem(model, u0, ts, p)
end
```

With the `ODEProblem` defined, you can use `OrdinaryDiffEq.jl` or
`DifferentialEquations.jl` to numerically integrate the orbit.

```@example example
trajectory = solve(problem; abstol=1e-14, reltol=1e-14)
```

Finally, you can use `Plots.jl` to show the result! For more information, 
consult the `SciML` [documentation](https://docs.sciml.ai), or the 
`GalacticPotentials.jl` [Getting Started](index.md) page.

```@example example
using Plots

figure = plot(
  trajectory, idxs=(:x, :y), 
  label=:none, title="Orbit in the Plummer Potential", 
)
```
