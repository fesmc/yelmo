# Time stepping

A call `yelmo_update(dom,time)` advances the model from its current time to
`time` in one or more internal time steps
([`src/yelmo_ice.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_ice.f90),
[`src/yelmo_timesteps.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_timesteps.f90)).
The internal time step is chosen with `yelmo.dt_method`:

| `dt_method` | Description |
|---|---|
| 0 | No internal time step: one step over the whole interval |
| 1 | Adaptive, limited by the Courant number `cfl_max` (default 0.1) of the [transport velocity](#transport-velocity-and-courant-limit), and reduced by a factor 0.05 when a checkerboard pattern appears in $\partial H/\partial t$ |
| 2 | Adaptive, set by the predictor–corrector error estimate (**default**) |

The first step after initialisation is always `dt_min` long. All parameters on
this page are in the `&yelmo` group.

## Predictor–corrector step

The ice thickness is always advanced with a predictor–corrector scheme, and the
difference between predictor and corrector gives an estimate of the local
truncation error; `dt_method` only sets how the time step is chosen. The step
is (`yelmo_update`):

1. **Predictor**: topography step from $H^n$ to $H^\ast$
   (`calc_ytopo_pc(..., "predictor")`).
2. **Dynamics**: velocity solution on the predicted geometry $H^\ast$
   (`calc_ydyn`).
3. **Corrector**: topography step from $H^n$ to $H^{n+1}$ with the new velocity
   (`calc_ytopo_pc(..., "corrector")`).
4. **Error estimate**: truncation error $\tau$ [m a$^{-1}$] from $H^{n+1}$ and
   $H^\ast$, and its norm $\eta$ (`pc_eta`, below).
5. **Redo**: if $\eta >$ `pc_tol` (default 1 a$^{-1}$) and $\Delta t >$ `dt_min`,
   the step is rejected and repeated with
   $\Delta t \cdot 0.7/(1 + (\eta - \mathrm{pc\_tol})/10)$ (at least `dt_min`).
   At most `pc_n_redo` (default 5) attempts are made; the last one is accepted.
6. **Other components**: material (`calc_ymat`), tracers (`calc_ytrc`),
   thermodynamics (`calc_ytherm`) and hydrology (`calc_yhyd`) are updated.
7. **Advance**: the topography is set to the new thickness. With
   `pc_use_H_pred = True` (default) this is the predictor $H^\ast$, which is
   consistent with the velocity solution; the corrector is then used only for
   the error estimate.

The schemes (`pc_method`) are:

- `"AB-SAM"` (default): Adams–Bashforth predictor and semi-implicit
  Adams–Moulton corrector (Cheng et al., 2017). The predictor advective rate is
  $(1 + \zeta/2)\,f(H^n) - (\zeta/2)\,f(H^{n-1})$, with
  $\zeta = \Delta t_n/\Delta t_{n-1}$, and the corrector rate is the mean of
  $f(H^\ast)$ and $f(H^n)$. The truncation error is
  $\tau = \zeta\,(H^{n+1} - H^\ast)/((3\zeta + 3)\,\Delta t)$.
- `"FE-SBE"`: forward Euler predictor and semi-backward Euler corrector (first
  order). It is also used for the first step.
- `"HEUN"`: Heun's method.

With `pc_filter_vel = True` (default), the ice thickness is advected with the
mean of the current and the previous velocity solution. The velocity fields
themselves are not filtered.

### Error norm

The error norm is the root mean square of the relative truncation error,

$$
\eta = \mathrm{RMS}\left(\frac{|\tau|}{1\,\mathrm m + 0.01\,H}\right),
$$

over the cells of the predictor–corrector mask (`set_pc_mask`). The mask
excludes ice thinner than `pc_eta_H_min` (default 10 m) or slower than
`pc_eta_u_min` (default 0), cells with a cell of `f_ice` < 1 in their 3×3
neighbourhood (including ice-free cells), floating and grounding-line cells, and isolated points with large
$\tau$. The fraction `pc_eta_trim` (default 0) of cells with the largest errors
is left out.

### Next time step

The next time step follows from a PI controller (`pc_controller`, default
`"PI42"`; also `"H312b"`, `"H312PID"`, `"H321PID"`, `"PID1"`). For PI42,

$$
\Delta t_{n+1} = \Delta t_n
\left(\frac{\varepsilon}{\eta_n}\right)^{k_I + k_P}
\left(\frac{\varepsilon}{\eta_{n-1}}\right)^{-k_P},
$$

with $k_I = 2/(5k)$, $k_P = 1/(5k)$, $k$ the order of the scheme, and the target
error $\varepsilon$ = `pc_eps` (default 0.02 a$^{-1}$, which must not exceed
`pc_tol`). The step is then limited:

- to at most `pc_rho_max` (default 2) times the previous step. The controller
  itself has no upper bound: after a step in which the ice hardly moves,
  $\eta_n$ is tiny and PI42 can ask for a step 10$^3$ times longer. The step
  ratio also sets the weights of the AB-SAM predictor, which then extrapolates
  far beyond the last two rates;
- by the Courant number $C$ = `pc_cfl_max` (default 0.5) of the
  [transport velocity](#transport-velocity-and-courant-limit);
- to the interval [`dt_min`, remaining time], with `dt_min` = 0.1 a by default.
  A step longer than half of the remaining time, but shorter than it, is set to
  half of the remaining time, to avoid a very short last step.

The call steps until `time` is reached, however many steps that takes.
`par_load` stops the model unless `pc_rho_max` > 1, `pc_cfl_max` and `cfl_max`
are in (0, 1], and `pc_eta_trim` is in [0, 0.5).

### Transport velocity and Courant limit

The Courant limit is computed from the velocity that advects $H$ in the
predictor, the transport velocity (`calc_transport_velocity`,
[`src/yelmo_topography.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_topography.f90)).
The predictor, the corrector and the Courant limit all call this routine, so
the limit is set by exactly the faces that move ice. On the cell faces
(ac-nodes) it is

1. the depth-averaged velocity $\bar{u}$, or with `pc_filter_vel = True`
   (default) the mean of the current and the previous velocity solution;
2. set to zero on faces that carry no ice (`set_inactive_margins`): faces
   between a partially ice-covered cell ($f_\mathrm{ice} < 1$) and an ice-free
   cell that may not fill. With the level set and a subgrid front
   (`ycalv.use_lsf = True`, `ytopo.front_subgrid` ≠ `"none"`), an ice-free cell
   may fill if the front covers at least `A_FRONT_MIN` = 10 % of it (level-set
   area fraction, `calc_lsf_area_fraction`); otherwise no ice-free cell may
   fill. A partial cell thus fills before ice flows beyond it.

The time step is then limited to

$$
\Delta t \le \min_{i,j} \frac{C}{u_c/\Delta x + v_c/\Delta y + 0.1/\Delta x},
$$

with $u_c$ and $v_c$ the largest speeds of the transport velocity on the faces
of each cell (`calc_adv2D_timestep1`), $C$ = `pc_cfl_max` (`dt_method` = 2) or
`cfl_max` (`dt_method` = 1). The limit is evaluated once per step, from
$f_\mathrm{ice}$, the level set and the velocity at the start of the step, i.e.
with the faces of the predictor. A face that opens during the step (the cell
beyond starts to fill) is not in the limit; the error estimate $\eta$ checks
the step instead.

Closed faces are left out because they can be fast: at a partial front cell,
the speed on the face towards the ocean is the extrapolated front velocity
(e.g. 5 km a$^{-1}$ at Rink Isbræ in GRL-8KM), which limited the step to
0.7 a although no ice crosses the face. With the transport velocity the
Courant limit is a backstop and the error estimate sets the step in the
initMIP runs (1 kyr, mean $\Delta t$: GRL-8KM 0.73 → 1.15 a, GRL-16KM
1.37 → 3.2 a, ANT-16KM 1.2 → 2.4 a). With `front_subgrid = "none"`,
$f_\mathrm{ice}$ is 0 or 1 and no face is closed.

There is no time-step limit from the 3D advection. The horizontal advection in
the thermodynamics is sub-cycled internally instead (see
[Thermodynamics](thermodynamics.md)).

The time steps, `pc_eta` and solver diagnostics are written to `timesteps.nc`
with `log_timestep = True`, and model speed and numerical metrics to
`yelmo_metrics.nc` with `write_metrics = True`.

## Instability checks

At the end of each step, `yelmo_check_kill` stops the model if

- $H$ reaches $10^4$ m, or the depth-averaged speed reaches
  2 `ydyn.ssa_vel_max`;
- $H$, $\bar u$ or $T$ are not finite;
- the mean of the last three values of `pc_eta` exceeds 10 `pc_tol`;
- the last 50 steps of the current `yelmo_update` call (or all of them, if
  fewer) were at `dt_min` (only for `dt_min` ≤ 0.01 a).

Non-finite boundary fields also stop the model. Before stopping, the state is
written to `yelmo_killed.nc`. `disable_kill = True` switches these checks off.

## Reference

Cheng, G., Lötstedt, P., and von Sydow, L. (2017). Accurate and stable time
stepping in ice sheet modeling. J. Comput. Phys., 329, 29–47.
