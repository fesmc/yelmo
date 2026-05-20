# Physics

This section describes the continuum equations and numerical methods that Yelmo
uses to solve the coupled momentum balance, mass conservation, and (in time)
thermodynamics of the ice sheet. The goal is to document what the code
actually implements, with the same variable names used in the source so that
the equations here can be traced back to specific routines. The notation
follows [Robinson et al. (2020)](https://doi.org/10.5194/gmd-13-2805-2020)
for mass conservation and
[Robinson, Goldberg, and Lipscomb (2022)](https://doi.org/10.5194/tc-16-689-2022)
for the momentum balance.

## Momentum balance

The ice flow is described as a slow, incompressible, gravity-driven flow of a
power-law (Glen) fluid. The full Stokes problem is never solved. Instead Yelmo
offers three approximations, listed in order of *decreasing* fidelity to
Stokes:

| Approximation | Membrane stresses | Vertical shear | Use case |
|---|---|---|---|
| [DIVA](momentum/diva.md) | yes | yes (closure) | first-class momentum balance, valid from slow interior to streaming flow and shelves |
| [SSA](momentum/ssa.md)   | yes | no (plug flow) | floating shelves, fast-streaming grounded ice |
| [SIA](momentum/sia.md)   | no  | yes (closure) | slow grounded interior, used as a *complement* to DIVA in hybrid mode |

DIVA is the recommended choice in Yelmo and is solved as a single 2D problem
for the depth-averaged horizontal velocity. SSA is recovered from DIVA in the
limit of vanishing vertical shear, and SIA is recovered in the limit of
vanishing membrane stresses and large basal drag. All three formulations share
the same driving stress, the same Glen viscosity, and (where applicable) the
same family of friction laws — see [Parameters](../parameters.md) for the
parameter choices.

The fourth page in this subsection documents how the DIVA / SSA momentum
balance is [discretised and solved](momentum/solvers.md), including the two
solvers currently supported: the **residual** assembler inherited from
Yelmo v1, and the **energy** assembler that minimises a discrete energy
functional.

## Mass conservation

The [continuity equation](mass_conservation/index.md) evolves the ice
thickness $H$ in time given the depth-averaged velocity from the momentum
balance and a set of mass-balance inputs (surface, basal, frontal, calving).
Yelmo offers several discretisations of the resulting hyperbolic transport
problem; the two recommended choices — an explicit donor-cell upwind scheme
and an implicit upwind scheme solved with LIS — are documented in
[Numerical solution](mass_conservation/solvers.md).

## Notation used throughout

- $u, v$: horizontal velocity components $[\mathrm{m\,a^{-1}}]$.
- $\bar u, \bar v$: depth-averaged horizontal velocity (the DIVA/SSA unknowns
  and the advecting field in the continuity equation).
- $u_b, v_b$: basal horizontal velocity.
- $H$: ice thickness $[\mathrm{m}]$; $s$ and $b$: ice surface and basal elevations.
- $z$: vertical Cartesian coordinate, $b \le z \le s$.
- $\rho_i$, $g$: ice density and gravitational acceleration.
- $A(T')$: Glen rate factor, depending on the pressure-corrected temperature $T'$
  and (where relevant) the water content $\omega$. $n = n_\mathrm{glen}$ is the Glen exponent (typically 3).
- $\dot\varepsilon_{ij}$: components of the horizontal strain-rate tensor.
- $\dot\varepsilon_e$: effective strain rate, the second invariant of $\dot\varepsilon_{ij}$.
- $\mu$: effective viscosity; $\bar\mu = \frac{1}{H}\int_b^s \mu \,\mathrm dz$: depth-averaged viscosity.
- $\boldsymbol\tau_d = (\tau_{d,x},\tau_{d,y})$: gravitational driving stress.
- $\boldsymbol\tau_b = (\tau_{b,x},\tau_{b,y})$: basal shear stress.
- $\beta$: basal friction coefficient, defined by $\boldsymbol\tau_b = \beta\,\mathbf u_b$.
- $\dot a$: surface mass balance (`smb`).
- $\dot b_g, \dot b_f$: basal mass balance on grounded and floating ice (`bmb`).
- $\dot c$: calving rate (`cmb`; mass loss, $\dot c \ge 0$).

All horizontal fields live on the Arakawa C-grid (called the **ac-grid** in
Yelmo): scalars at cell centres (`aa`-nodes), $u$ on the right face of each cell
(`acx`-nodes), $v$ on the top face (`acy`-nodes), and cross-derivative
quantities at the cell corners (`ab`-nodes).
