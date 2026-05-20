# Mass conservation — continuity equation

Mass conservation evolves the ice thickness $H$ in time given the
depth-averaged horizontal velocity $\bar{\mathbf u} = (\bar u, \bar v)$
from the [momentum balance](../momentum/diva.md) and a set of
**mass-balance inputs** (surface, basal, frontal, calving). The
implementation lives in
[`src/physics/mass_conservation.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/mass_conservation.f90),
with the advection step delegated to
[`src/physics/solver_advection.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/solver_advection.f90)
and the per-timestep driver in
[`src/yelmo_topography.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_topography.f90).

## Continuity equation

Following [Robinson et al. (2020)](https://doi.org/10.5194/gmd-13-2805-2020),
Eq. (2), the ice thickness obeys the vertically integrated
mass-conservation equation

$$
\frac{\partial H}{\partial t}
\;=\; -\,\nabla\!\cdot\!\bigl(H\,\bar{\mathbf u}\bigr)
\;+\; \dot a
\;+\; \dot b_g
\;+\; \dot b_f
\;-\; \dot c,
$$

with

- $H$ the ice thickness;
- $\bar{\mathbf u} = (\bar u, \bar v)$ the depth-averaged horizontal
  velocity;
- $\dot a$ the surface mass balance (`bnd%smb` in the code, units
  $\mathrm{m\,a^{-1}}$ ice equivalent);
- $\dot b_g$ the basal mass balance on grounded ice (`bnd%bmb_grnd`);
- $\dot b_f$ the basal mass balance on floating ice
  (`bnd%bmb_shlf`);
- $\dot c$ the calving rate (positive for mass loss).

All five forcings on the right-hand side are **inputs** to the
mass-conservation step: Yelmo does not generate them — they are
supplied either from a coupled boundary class (`bnd`), from
sub-models (basal hydrology, sub-shelf melt), or from a calving law.
The calving laws themselves are documented separately and are not
covered here.

The continuity equation is written in flux-divergence form so that
mass is conserved by construction up to the boundary fluxes:
the divergence operator $\nabla\!\cdot(H\,\bar{\mathbf u})$ moves
mass from cell to cell without creating or destroying it, and all
source / sink terms appear explicitly on the right-hand side. This
is the property the numerical schemes are required to preserve at
the discrete level — see [Numerical solution](solvers.md).

## Operator-split form used in the code

The Yelmo timestep updates $H$ in an **operator-split** sequence,
applying the advective divergence first and the mass-balance source
terms in turn. Schematically, per (sub-)timestep $\Delta t$:

$$
H^{n+1}
\;=\; H^{n}
\;+\; \Delta t\,\dot H_\mathrm{dyn}
\;+\; \Delta t\,\bigl(\dot a + \dot b_g + \dot b_f - \dot c\bigr),
$$

where $\dot H_\mathrm{dyn} = -\nabla\!\cdot(H\,\bar{\mathbf u})$ is
computed by the advection solver and the mass-balance contributions
are added one by one through the helper
`apply_tendency`, which also enforces non-negativity of $H$ and
clips mass-loss rates to $\dot m \ge -H/\Delta t$. (A few additional
terms — sub-grid discharge `dmb`, relaxation `mb_relax`, and a
boundary residual `mb_resid` — are added in the same loop; they are
diagnostic in most configurations and do not change the continuum
equation above.) The same advection solver is used for tracers such
as the ice area fraction `f_ice` and the level-set field used by
some calving schemes.

The key choice is how the advective tendency $\dot H_\mathrm{dyn}$ is
discretised. Yelmo provides several schemes, selected at runtime by
the namelist parameter `ytopo.solver`; the two recommended choices —
an explicit donor-cell upwind scheme (`expl-upwind`) and an implicit
upwind scheme solved with LIS (`impl-lis`) — are described in
[Numerical solution](solvers.md).
