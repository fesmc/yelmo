# Numerical solution of the continuity equation

The advective tendency $\dot H_\mathrm{dyn} = -\nabla\!\cdot(H\,\bar{\mathbf u})$
is discretised on the same Arakawa C-grid (the **ac-grid**) used by the
[momentum solvers](../momentum/solvers.md): $H$ lives on aa-nodes (cell
centres) and $\bar{\mathbf u}$ lives on the ac-faces ($\bar u$ on acx,
$\bar v$ on acy). All Yelmo advection schemes evaluate the face fluxes
in **donor-cell upwind** form,

$$
Q^{x}_{i+\tfrac12,j}
\;=\; \bar u_{i+\tfrac12,j}\,H^{\mathrm{up}}_{i+\tfrac12,j},
\qquad
H^{\mathrm{up}}_{i+\tfrac12,j}
\;=\;
\begin{cases}
H_{i,j}   & \text{if } \bar u_{i+\tfrac12,j} > 0,\\
H_{i+1,j} & \text{if } \bar u_{i+\tfrac12,j} < 0,
\end{cases}
$$

with $Q^{y}$ defined analogously on the acy-faces. The corresponding
flux-divergence at cell $(i,j)$ is

$$
\nabla\!\cdot(H\,\bar{\mathbf u})\big|_{i,j}
\;=\; \frac{Q^{x}_{i+\tfrac12,j} - Q^{x}_{i-\tfrac12,j}}{\Delta x}
   + \frac{Q^{y}_{i,j+\tfrac12} - Q^{y}_{i,j-\tfrac12}}{\Delta y}.
$$

The scheme is selected with `ytopo.solver`. The two recommended schemes,
`"impl-lis"` (default) and `"expl-upwind"`, differ only in **how** the upwind
fluxes are evaluated in time. The other options are `"expl"`, `"impl-upwind"`,
`"expl-sico"`, `"impl-sico"`, `"impl-sico-lis"`, `"expl-new"` and `"none"` (no
advection). All schemes are first order in space. The fluxes are
conservative (the flux leaving cell $(i,j)$ across a face is the same flux
entering its neighbour), so the advection does not create or destroy mass,
apart from the flux through the domain edge and, in the explicit schemes, the
rate limiter below. The ice mask `bnd%mask_ice` is not imposed in the
advection: every cell is advected, and `calc_G_boundaries` imposes the mask
afterwards and books the change in `mb_resid` (see
[Mass conservation](index.md)).

The advection solver is called twice per time step, in the predictor and the
corrector of the [time-stepping scheme](../timestepping.md).

## `expl-upwind` — explicit forward-Euler donor-cell upwind

Selected by `ytopo.solver = "expl-upwind"`. The implementation is
`calc_adv2D_expl_upwind` in `solver_advection.f90`; the scheme follows
[Winkelmann et al. (2011)](https://doi.org/10.5194/tc-5-715-2011), Eq. (18).

The flux $Q$ is evaluated using $\bar{\mathbf u}^{\,n}$ and $H^{\,n}$
at the start of the timestep, and the ice thickness is advanced by
explicit forward Euler:

$$
H^{n+1}_{i,j}
\;=\; H^{n}_{i,j}
\;-\; \Delta t\,
\left[
\frac{Q^{x}_{i+\tfrac12,j} - Q^{x}_{i-\tfrac12,j}}{\Delta x}
+ \frac{Q^{y}_{i,j+\tfrac12} - Q^{y}_{i,j-\tfrac12}}{\Delta y}
\right]^{n}
\;+\; \Delta t\,\dot m_{i,j}^{\,n},
$$

where $\dot m$ collects whatever source/sink term is being applied in
the same call (zero in Yelmo, where the mass-balance terms are added
separately). The flux divergence is limited to
$10^{3}\,\mathrm{m\,a^{-1}}$ for safety.

The scheme is first-order accurate in time and space, monotone, and
positivity-preserving for the depth-averaged transport problem when
the CFL condition

$$
\Delta t \;\le\; \min_{i,j}
\frac{1}{|\bar u_{i+\tfrac12,j}|/\Delta x \,+\, |\bar v_{i,j+\tfrac12}|/\Delta y}
$$

is respected. The adaptive time step is limited by this condition with the
Courant number `yelmo.pc_cfl_max` (default 0.5; `yelmo.cfl_max` for
`dt_method = 1`), see [Time stepping](../timestepping.md).

## `impl-lis` — implicit upwind via LIS

Selected by `ytopo.solver = "impl-lis"`. The implementation is
`linear_solver_matrix_advection_csr_2D` in `solver_advection.f90`.

The face fluxes are evaluated with the velocity field $\bar{\mathbf u}^{\,n}$
(frozen at the start of the timestep) but with the unknown thickness
$H^{n+1}$, so the upwind discretisation becomes a sparse linear system
for $H^{n+1}$ instead of an explicit update:

$$
H^{n+1}_{i,j}
\;+\; \frac{\Delta t}{\Delta x\,\Delta y}
\left[
F^{x}_{i+\tfrac12,j}(H^{n+1})
- F^{x}_{i-\tfrac12,j}(H^{n+1})
+ F^{y}_{i,j+\tfrac12}(H^{n+1})
- F^{y}_{i,j-\tfrac12}(H^{n+1})
\right]
\;=\; H^{n}_{i,j}
\;+\; \Delta t\,\dot m_{i,j},
$$

with $F^{x}_{i+\tfrac12,j} = \bar u_{i+\tfrac12,j}\,\Delta y \cdot H^{\mathrm{up}}_{i+\tfrac12,j}$
and analogously for $F^{y}$. Each face contributes a single
off-diagonal coefficient to the row for cell $(i,j)$ — the one
corresponding to its upwind neighbour — and a diagonal contribution
when cell $(i,j)$ itself is the upwind donor. The assembled
operator $A\,H^{n+1} = b$ therefore has at most five non-zeros per
row (centre plus four neighbours), is stored in CSR form, and is solved
each timestep by [LIS](http://www.ssisc.org/lis/) with BiCGSTAB and a
Jacobi preconditioner (tolerance 1e-12, set in `solver_advection.f90`,
not a parameter).

The rows at the domain edge depend on the boundary type of the experiment
(`yelmo.experiment`): zero thickness at the edge (default), zero normal
gradient (`"infinite"`), periodic in both directions, or periodic in one
direction and zero gradient in the other (`"periodic-x"`, `"periodic-y"`).
`"MASK_ICE"` uses zero gradient on all edges. MISMIP3D and TROUGH-F17 have
zero thickness at $x_\mathrm{max}$, zero gradient at the symmetry edge
$x = 0$ and periodic sides in $y$. Border cells with `bnd%mask_ice = NONE` or
`FIXED` do not use the edge condition: they are solved as inner cells with no
flux through the domain edge, so that the flux through their inner faces is
kept and booked in `mb_resid`.

Because the implicit upwind scheme is unconditionally stable for this
linear transport problem, it tolerates larger Courant numbers than
`expl-upwind`. In Yelmo, the time step is set by the predictor–corrector error
estimate and the Courant cap `pc_cfl_max` for both schemes. The price is one
sparse linear solve per call and a slightly more diffusive solution than the
explicit upwind scheme at the same $\Delta t$.

## Choosing between the two

- `impl-lis` (default) is robust where the velocity field has localised
  fast peaks, e.g. at the grounding line and in outlet glaciers.
- `expl-upwind` is fast, simple and fully local, a reasonable choice for
  high-resolution simulations where $\Delta t$ is already small for other
  reasons.

The choice does not change the continuum equation being solved, only the
time integration and the resulting cost / smoothing trade-off.
