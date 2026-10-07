# Physics

This section describes the continuum equations and numerical methods that Yelmo
uses to solve the coupled momentum balance, mass conservation and
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
offers the following approximations, in order of *decreasing* fidelity to
Stokes:

| Approximation | `ydyn.solver` | Membrane stresses | Vertical shear | Use case |
|---|---|---|---|---|
| [DIVA](momentum/diva.md) | `"diva"` (**default**), `"diva-noslip"` | yes | yes (closure) | valid from slow interior to streaming flow and shelves |
| [SSA](momentum/ssa.md)   | `"ssa"` | yes | no (plug flow) | floating shelves, fast-streaming grounded ice |
| [SIA](momentum/sia.md)   | `"sia"` | no  | yes | slow grounded interior without sliding |
| SIA + SSA | `"hybrid"` | yes | yes (SIA) | depth-averaged velocity = SIA shear + SSA basal velocity |

`"fixed"` keeps the velocity unchanged. DIVA is the recommended choice in Yelmo
and is solved as a single 2D problem for the depth-averaged horizontal
velocity. SSA is recovered from DIVA in the limit of vanishing vertical shear,
and SIA in the limit of vanishing membrane stresses and large basal drag. All
formulations share the same driving stress, the same Glen viscosity and, where
applicable, the same [basal friction](basal-friction.md) laws.

The [Numerical solution](momentum/solvers.md) page documents how the DIVA / SSA
momentum balance is discretised and solved, with two linear-system assemblers:
the **energy** assembler (default, `ydyn.ssa_solver = "energy"`), which
minimises a discrete energy functional, and the **residual** assembler, built
from the strong form of the equations. It also describes the velocity limit.

## Vertical velocity

The [vertical velocity](vertical-velocity.md) $w$ is diagnosed from
incompressibility once $u, v$ are known, anchored at the base by the basal
kinematic boundary condition:

$$
w(z) = \underbrace{\frac{\partial b}{\partial t} + u_b \frac{\partial b}{\partial x} + v_b \frac{\partial b}{\partial y} + \dot b}_{w_b} \; - \int_b^z \left( \frac{\partial u}{\partial x}\bigg|_{z'} + \frac{\partial v}{\partial y}\bigg|_{z'} \right) \mathrm{d}z'
$$

Alongside $w$ (`uz`), Yelmo also forms $w^\star$ (`uz_star`), the
**sigma-relative** advective vertical velocity used by every scalar advection
scheme that takes horizontal gradients at constant $\zeta$ — thermodynamics and
the age/tracer solver. The two are distinct and not interchangeable; see
[Vertical velocity](vertical-velocity.md).

## Mass conservation

The [continuity equation](mass_conservation/index.md) evolves the ice
thickness $H$ in time given the depth-averaged velocity from the momentum
balance and the mass-balance terms (surface, basal, frontal, calving). Yelmo
offers several discretisations of the advection; the two recommended
choices, an implicit upwind scheme solved with LIS (`ytopo.solver = "impl-lis"`,
default) and an explicit donor-cell upwind scheme (`"expl-upwind"`), are
documented in [Numerical solution](mass_conservation/solvers.md).

## Time stepping

The ice thickness is advanced with an adaptive
[predictor–corrector scheme](timestepping.md), whose time step is set by an
estimate of the truncation error and limited by a Courant condition and a cap
on its growth from one step to the next.

## Thermodynamics

The [thermodynamics](thermodynamics.md) solves for the enthalpy of the ice
column (`ytherm.method = "enth"`, default), from which the temperature and the
water content of temperate ice follow. The ice column is coupled to a bedrock
column (in equilibrium with the geothermal heat flux by default), and the basal
boundary condition determines the basal melt or freeze-on. The temperature and
water content feed back into the momentum balance through the Glen rate
factor $A$ and the sub-temperate sliding factor.

## Basal friction

The [basal friction](basal-friction.md) coefficient $\beta$ follows from a
friction law (linear, pseudo-plastic or regularized Coulomb), a bed coefficient
that depends on the till properties and the effective pressure, and a reduction
of sliding below the pressure melting point.

## Calving

By default, the calving front is a [level set](calving.md) that moves with the
ice velocity plus a calving velocity from the von Mises stress law of
Morlighem et al. (2016) (`ycalv.use_lsf = True`,
`calv_flt_method = calv_grnd_method = "vm-m16"`), and front cells are partially
ice covered. Without the level set, calving is a mass-loss rate in the
continuity equation, e.g. from the effective-stress law of Lipscomb et al.
(2019).

## Notation used throughout

- $u, v$: horizontal velocity components $[\mathrm{m\,a^{-1}}]$.
- $\bar u, \bar v$: depth-averaged horizontal velocity (the DIVA/SSA unknowns
  and the advecting field in the continuity equation).
- $u_b, v_b$: basal horizontal velocity.
- $H$: ice thickness $[\mathrm{m}]$; $s$ and $b$: ice surface and basal elevations.
- $z$: vertical Cartesian coordinate, $b \le z \le s$.
- $\rho_i$, $g$: ice density and gravitational acceleration.
- $A(T')$: Glen rate factor, depending on the homologous temperature $T'$
  and optionally the water content $\omega$. $n = n_\mathrm{glen}$ is the Glen exponent (typically 3).
- $\dot\varepsilon_{ij}$: components of the horizontal strain-rate tensor.
- $\dot\varepsilon_e$: effective strain rate, the second invariant of $\dot\varepsilon_{ij}$.
- $\mu$: effective viscosity; $\bar\mu = \frac{1}{H}\int_b^s \mu \,\mathrm dz$: depth-averaged viscosity.
- $\boldsymbol\tau_d = (\tau_{d,x},\tau_{d,y})$: gravitational driving stress.
- $\boldsymbol\tau_b = (\tau_{b,x},\tau_{b,y})$: basal shear stress.
- $\beta$: basal friction coefficient, defined by $\boldsymbol\tau_b = \beta\,\mathbf u_b$.
- $\dot a$: surface mass balance (`smb`).
- $\dot b_g, \dot b_f$: basal mass balance on grounded and floating ice
  (`bmb_grnd`, `bmb_shlf`; combined in `bmb`).
- $\dot f$: frontal mass balance (`fmb`).
- $\dot c$: calving mass balance (`cmb`; negative for mass loss).

All horizontal fields live on the Arakawa C-grid (called the **ac-grid** in
Yelmo): scalars at cell centres (`aa`-nodes), $u$ on the right face of each cell
(`acx`-nodes), $v$ on the top face (`acy`-nodes), and cross-derivative
quantities at the cell corners (`ab`-nodes).
