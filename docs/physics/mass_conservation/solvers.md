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

The two recommended schemes — selected at runtime via the namelist
parameter `ytopo.solver` — differ only in **how** the upwind fluxes are
evaluated in time. Both are mass-conservative by construction (the flux
leaving cell $(i,j)$ across a face is the same flux entering its
neighbour), so the discrete update of $H$ summed over the whole domain
matches the discrete integral of the source/sink terms exactly, up to
boundary outflows.

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
the same call (zero in the standard dynamic-only call from the
`mass_conservation.f90` driver, since the mass-balance terms are added
separately). A hard rate limiter $|\dot H| \le 10^{3}\,\mathrm{m\,a^{-1}}$
is applied for safety.

The scheme is first-order accurate in time and space, monotone, and
positivity-preserving for the depth-averaged transport problem when
the CFL condition

$$
\Delta t \;\le\; \min_{i,j}
\frac{1}{|\bar u_{i+\tfrac12,j}|/\Delta x \,+\, |\bar v_{i,j+\tfrac12}|/\Delta y}
$$

is respected. Yelmo enforces this via an adaptive timestepper at the
top level (see [Parameters](../../parameters.md)), so the user does not
need to tune the advection step manually — but the scheme is sensitive
to under-resolved velocity peaks at the grounding line.

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
each timestep by [LIS](http://www.ssisc.org/lis/). The default solver
configuration is BiCG with an ILU preconditioner, configurable at
runtime via `ytopo.adv_lis_opt`.

Boundary rows are set per face by the `boundaries` flag: `"zero"`
imposes $H = 0$ at the edge (Dirichlet), `"infinite"` zeroes the
normal derivative (Neumann outflow), and `"periodic"` wraps the
stencil to the opposite edge. The mask values `mask = 0` and
`mask = -1` further allow per-cell pinning of $H$ to zero or to its
previous value.

Because the implicit upwind scheme is unconditionally stable for this
linear transport problem, it tolerates larger timesteps than
`expl-upwind` — typically the limiting factor becomes the *physical*
adaptive timestep used by Yelmo (driven by the momentum balance and
the mass-balance terms) rather than a CFL constraint on the
advection. The price is one sparse linear solve per timestep and a
slightly more diffusive solution than the explicit upwind scheme at
the same $\Delta t$.

## Choosing between the two

- `expl-upwind` is fast, simple, fully local, and the natural default
  for high-resolution simulations where $\Delta t$ is already small for
  other reasons (e.g. fast streaming flow).
- `impl-lis` is the right choice when the velocity field has localised
  fast peaks that would force `expl-upwind` into very small timesteps —
  it lets the global adaptive $\Delta t$ be set by physics elsewhere in
  the model.

Both schemes are mass-conservative to round-off; the choice does not
change the continuum equation being solved, only the time-stepping and
the resulting cost / smoothing trade-off.
