# Gravitational Potentials

Gravitational potential fields provide an alternative (approximate) approach
to large-dimensioned n-body systems. Rather than tracking the orbits of _all_
particles in an n-dimensional system, potential fields allow you to integrate
_one orbit at a time_. The Plummer potential field is a common potential
function used in the field of galactic dynamics; the expression for the
Plummer potential is shown below: $G$ is the gravitational field strength,
$m$ is the mass of the central body, and $b$ is a softening parameter to avoid
infinities when the position of the particle ($x$, $y$, $z$) is near the origin.

$$\Phi = - \frac{G m}{\sqrt{b^{2} + x^{2} + y^{2} + z^{2}}}$$

Given any scalar potential field, the field's gradient provides the strength
and force applied to a body at any position in the field per unit mass. If we
treat the orbiting body as a point mass, then we can use the general ordinary
differential equation below to numerically integrate any orbit along the
potential with state vector $u$, parameter vector $p$, and scalar time $t$.

$$\dot{u} = -\nabla \Phi(u,p,t)$$

Note the generality! This kind of _recipe_ is well suited to tools like
`ModelingToolkit.jl`: you may write your expressions _mathematically_ and
let the `SciML` ecosystem generate fast and non-allocating codes, and
efficiently integrate the dynamics forward (or backward) in time.
`GalacticPotentials.jl` provides mathematical descriptions of common scalar
potential fields used in galactic dynamics, and hooks these descriptions into
`ModelingToolkit.jl` types for ease of use.
