# SSA — Shallow Shelf Approximation

The SSA is the membrane-stress balance that DIVA reduces to in the limit
of vanishing vertical shear. It treats the ice column as a vertical plug:
the horizontal velocity is independent of $z$, so $\bar{\mathbf u}$
*is* the basal velocity, $\bar{\mathbf u} = \mathbf u_b$. The
implementation is in
[`src/physics/velocity_ssa.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/velocity_ssa.f90).

The SSA is solved by its own Picard loop (`calc_velocity_ssa`) with the same
linear-system assemblers as DIVA (`ydyn.ssa_solver`, see
[Numerical solution](solvers.md)); the difference is in the closure. SSA is
the appropriate model for floating ice shelves (no basal drag, no vertical
shear) and the physical limit in fast-streaming grounded ice. It is used with
`ydyn.solver = "ssa"`, where $\bar{\mathbf u} = \mathbf u_b$, and in the hybrid
mode (`"hybrid"`), where the SSA basal velocity is added to the SIA shear
velocity, $\bar{\mathbf u} = \bar{\mathbf u}^{(\mathrm{SIA})} + \mathbf u_b$.

## Continuum equations

The depth-integrated stress balance has the same shape as the DIVA
equation but with the friction coefficient $\beta$ instead of
$\beta_\mathrm{eff}$:

$$
\frac{\partial}{\partial x}\!\left[\, 2\,\bar\mu\,H\,(2\,\bar\varepsilon_{xx} + \bar\varepsilon_{yy})\,\right]
+ \frac{\partial}{\partial y}\!\left[\, 2\,\bar\mu\,H\,\bar\varepsilon_{xy}\,\right]
- \beta\,\bar u
\;=\; \tau_{d,x},
$$

$$
\frac{\partial}{\partial y}\!\left[\, 2\,\bar\mu\,H\,(2\,\bar\varepsilon_{yy} + \bar\varepsilon_{xx})\,\right]
+ \frac{\partial}{\partial x}\!\left[\, 2\,\bar\mu\,H\,\bar\varepsilon_{xy}\,\right]
- \beta\,\bar v
\;=\; \tau_{d,y},
$$

with depth-averaged viscosity $\bar\mu$, driving stress $\boldsymbol\tau_d$
and friction law $\boldsymbol\tau_b = \beta\,\mathbf u_b$ defined exactly
as in [DIVA](diva.md).

## What is dropped relative to DIVA

The effective strain rate omits the $\frac{\partial u}{\partial z}$ and
$\frac{\partial v}{\partial z}$ terms:

$$
\dot\varepsilon_e^{\,2}
\;=\; \dot\varepsilon_{xx}^{\,2}
   + \dot\varepsilon_{yy}^{\,2}
   + \dot\varepsilon_{xx}\,\dot\varepsilon_{yy}
   + \dot\varepsilon_{xy}^{\,2}
   + \varepsilon_0^{\,2},
$$

so the Glen viscosity $\mu$ depends only on the horizontal strain rates (and,
through $A(T'(z))$, still varies with depth before it is depth-integrated).
There is no F-integral closure: $F_2 \to 0$, so

$$
\beta_\mathrm{eff} \to \beta, \qquad u_b \to \bar u.
$$

The basal stress diagnosed after the solve is simply
$\boldsymbol\tau_b = \beta\,\bar{\mathbf u}$. $\beta$ is the friction
coefficient after the grounding-line scaling, staggering and `beta_min` limit
(see [Basal friction](../basal-friction.md)).

For ice shelves $\beta = 0$ identically and only the membrane terms and
$\boldsymbol\tau_d$ remain — this is the SSA in its purest form.
