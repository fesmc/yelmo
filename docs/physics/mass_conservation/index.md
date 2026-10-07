# Mass conservation — continuity equation

Mass conservation evolves the ice thickness $H$ in time given the
depth-averaged horizontal velocity $\bar{\mathbf u} = (\bar u, \bar v)$
from the [momentum balance](../momentum/diva.md) and a set of
**mass-balance terms**. The per-timestep driver is `calc_ytopo_pc` in
[`src/yelmo_topography.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_topography.f90),
the mass-balance helpers are in
[`src/physics/mass_conservation.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/mass_conservation.f90)
and the advection schemes in
[`src/physics/solver_advection.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/solver_advection.f90).

## Continuity equation

Following [Robinson et al. (2020)](https://doi.org/10.5194/gmd-13-2805-2020),
Eq. (2), the ice thickness obeys the vertically integrated
mass-conservation equation

$$
\frac{\partial H}{\partial t}
\;=\; -\,\nabla\!\cdot\!\bigl(H\,\bar{\mathbf u}\bigr)
\;+\; \dot a
\;+\; \dot b
\;+\; \dot f
\;+\; \dot d
\;+\; \dot c,
$$

with

- $\dot a$ the surface mass balance (`smb`, from `bnd%smb`), in
  $\mathrm{m\,a^{-1}}$ ice equivalent;
- $\dot b$ the basal mass balance (`bmb`), combining the grounded basal mass
  balance computed by the [thermodynamics](../thermodynamics.md)
  (`thrm%now%bmb_grnd`) and the sub-shelf melt `bnd%bmb_shlf`;
- $\dot f$ the frontal mass balance at marine ice fronts (`fmb`);
- $\dot d$ the sub-grid discharge (`dmb`);
- $\dot c$ the calving mass balance (`cmb`, negative for mass loss), see
  [Calving](../calving.md).

The continuity equation is written in flux-divergence form, so that the
advection moves mass from cell to cell without creating or destroying it, and
all sources and sinks appear explicitly on the right-hand side.

## Operator-split form used in the code

Each topography step applies the terms one after another, each with
`apply_tendency`, in this order:

1. advection, $-\nabla\cdot(H\bar{\mathbf u})$ (`dHidt_dyn`), with the scheme
   set by `ytopo.solver` (see [Numerical solution](solvers.md));
2. surface mass balance `smb`;
3. basal mass balance `bmb`;
4. frontal mass balance `fmb`;
5. sub-grid discharge `dmb`;
6. calving `cmb`, on the level-set path or the mass-balance path;
7. relaxation towards a reference thickness `mb_relax` (only with
   `ytopo.topo_rel` ≠ 0);
8. boundary and margin corrections `mb_resid`.

The ice fraction `f_ice` is diagnosed from $H$ between the stages
(`update_ice_fraction`); ice thinner than 1 mm (`H_ice_eps`) counts as ice
free. `apply_tendency` adds $\Delta t\,\dot m$ to $H$, sets
negative thickness to zero and resets $\dot m$ to the rate actually applied, so
that the stored rates close the budget. The change made by this limit during
the advection step (mostly ice added where advection would give negative
thickness) is booked as `mb_clip`. The budget of a step is

$$
\frac{\Delta H}{\Delta t} = \dot H_\mathrm{dyn} + \dot m_\mathrm{clip} + \dot m_\mathrm{net} + \dot c,
$$

with `mb_net` = `smb + bmb + fmb + dmb + mb_relax + mb_resid`; the remainder is
written as `mb_err`. The region time series (`yelmo_ts.nc`) include the totals
`mb_relax_tot`, `mb_resid_tot` and `mb_clip_tot`; the global region covers the
whole domain.

The mass-balance rates are prepared by `calc_G_mbal`: melt is limited to the
available ice ($\dot m \ge -H/\Delta t$), ice-free cells cannot melt, and no ice
grows in the open ocean. The surface and basal mass balance are also scaled by
`f_ice` in partially ice-covered cells.

The advective rate is weighted between time levels by the
[predictor–corrector scheme](../timestepping.md), and with
`yelmo.pc_filter_vel = True` (default) $H$ is advected with the mean of the
current and the previous velocity solution. Velocities on faces from a
partially ice-covered cell into an ice-free cell are set to zero
(`set_inactive_margins`), so that partial cells fill before ice spreads
further; with the level set and the subgrid front, faces into cells the front
covers by at least 10 % stay open. This transport velocity
(`calc_transport_velocity`) also sets the Courant limit of the time step (see
[Time stepping](../timestepping.md#transport-velocity-and-courant-limit)).

## Basal mass balance

The grounded basal mass balance (`bmb_grnd`, from the thermodynamics) and the
floating basal mass balance (`bnd%bmb_shlf`) are combined into `bmb` by
`calc_bmb_total` according to `ytopo.bmb_gl_method`, which sets the melt in
partially grounded cells: `"pmp"` (default; partial melt, weighted by the
grounded fraction), `"fmp"` (full melt), `"fcmp"` (flotation criterion),
`"pmpt"` (partial melt over a grounding-zone transition, `gz_Hg0`, `gz_Hg1`)
or `"nmp"` (no melt). With `ytopo.use_bmb = False`, the basal and the frontal
mass balance are not applied.

## Frontal mass balance

The frontal mass balance is selected with `ytopo.fmb_method` (`calc_fmb_total`).
Methods 1–3 act on front cells, i.e. ice cells below flotation height with an
ice-free neighbour:

- 0 (default): prescribed by the boundary field `bnd%fmb_shlf`, wherever it is
  set;
- 1: proportional to the basal melt of the neighbouring ice-free cells
  (Pollard and DeConto, 2012, 2016),
  $$
  \dot f = \dot b_\mathrm{eff}\,\frac{A_f}{A_\mathrm{tot}}\,\theta_f,
  $$
  where $\dot b_\mathrm{eff}$ is the mean `bmb_shlf` of the ice-free
  neighbours, $A_\mathrm{tot} = \Delta x^2$ the cell area, $A_f$ the area of the
  submerged faces (the submerged depth of the cell's ice, times $\Delta x$, for
  each ice-free neighbour), and $\theta_f$ = `ytopo.fmb_scale` (default 1;
  Pollard and DeConto, 2016, suggest 10);
- 2: `bnd%fmb_shlf` scaled by the submerged front area;
- 3: the frontal melt of Rignot et al. (2016) from the ocean thermal forcing
  `bnd%tf_shlf` and the subglacial discharge `bnd%Qd`, scaled by
  `ytopo.fmb_lambda` (ISMIP7 protocol).

The frontal mass balance is lateral: it removes ice at the front without moving
the surface or base of the column.

## Sub-grid discharge

With `ytopo.dmb_method = 1`, grounded ice near the coast loses mass by
sub-grid discharge through outlet glaciers that the grid does not resolve
(`calc_mb_discharge`; parameters `dmb_alpha_max`, `dmb_tau`, `dmb_sigma_ref`,
`dmb_m_d`, `dmb_m_r`). It is off by default.

## Relaxation

With `ytopo.topo_rel` ≠ 0, the thickness is relaxed towards a reference
(`topo_rel_field`: the reference thickness `H_ref`, or the thickness at the
start of the step) with the time scale `topo_rel_tau`, in a subset of cells:
1, floating cells and cells without ice in the reference; 2, as 1 plus the
grounding line; 3, all cells; 4, grounded grounding-zone cells; −1, the
per-cell time scale `bnd%tau_relax`. Relaxation applies only where
`bnd%mask_ice` = `DYNAMIC`. A driver can change `topo_rel` and `topo_rel_tau`
over time with the `&relax` group, with or without the friction optimization
(see [Optimization](../../optimization.md)).

## Boundary and margin corrections

`calc_G_boundaries` books the following corrections as `mb_resid`:

- margin ice with an effective thickness below `ycalv.H_min_flt` (floating,
  default 10 m) or `ycalv.H_min_grnd` (grounded, default 5 m) is removed at the
  rate $H/\tau$, with $\tau$ = `ycalv.H_min_tau` (default 10 yr);
- isolated ice islands are removed, and margin cells thicker than all their
  ice-covered neighbours are reduced to the neighbours' maximum (except subgrid
  front cells);
- the domain boundary conditions and the ice mask `bnd%mask_ice` are applied:
  no ice where `mask_ice` is `NONE`, and the prescribed thickness where it is
  `FIXED`. This is the only place where the mask is imposed: the advection
  moves ice into and out of masked cells, and the change made here is booked,
  so the budget closes at masked cells too.
