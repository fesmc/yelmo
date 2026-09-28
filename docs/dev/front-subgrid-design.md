# Subgrid ice-front treatment: design proposal

Status: design agreed 2026-09-28 (decisions in section 8). Nothing implemented yet.

## 1. Problem

Ice that is advected into an ice-free cell at a front is spread over the
whole cell. With the default `ytopo.margin_flt_subgrid = False`, `f_ice` is
binary, so the front is a row of thin cells next to thick ones. This shows up
in three ways in the 2026-09 tests (ANT-16KM, GRL-8KM initmip):

- **Calving does not set the front.** With `vm-l19`, calving demand is ~7400
  km3/yr but only ~700 km3/yr can be applied. About 400 of ~900 floating front
  cells are emptied every step, mostly by the thin-ice and tongue rules
  (`apply_calving_rate_thin`, `calc_calving_rate_tongues`), and the excess is
  lost at the clip in `apply_tendency`. The front sits where thin cells are
  deleted, which is why the front position still depends on dt (~1% floating
  area between Courant 0.5 and 0.8).
- **The SSA converges poorly.** GRL-8KM with LSF and `H_min_tau = 10` yr needs
  11.7 Picard iterations per step (p90 = 20) instead of ~5. The largest
  residuals sit in thin floating front cells (H = 30–40 m next to 0–5 m, speeds
  of 0.7–2.4 km/yr).
- **Time-step dependent removal.** Removing thin cells every step made the
  removal per year scale with the number of steps. This is fixed with
  `H_min_tau`, but the thin cells it acts on are still unphysical.

## 2. What Yelmo does now

Two regimes, both unsatisfactory:

- **`margin_flt_subgrid = False` (all current configs).** `f_ice` is 0 or 1.
  Thin front cells are full cells in the momentum balance, with their true
  (thin) thickness.
- **`margin_flt_subgrid = True`.** `f_ice = H/H_eff` at floating margins, with
  `H_eff` the mean thickness of all ice-covered neighbours
  (`calc_ice_fraction`, topography.f90:372). Partial cells are then:
  - blocked from exporting into ice-free cells (`set_inactive_margins`,
    velocity_general.f90:1635), and fill until f = 1;
  - **invisible to the momentum balance**: `set_ssa_masks` only solves faces
    with an f = 1 neighbour, `calc_ice_front` codes partial cells as ice-free,
    the lateral BC is applied at the full/partial face with full-cell values,
    viscosity and beta are f = 1 only.

  This is essentially the Lipscomb et al. (2019) "inactive calving-front cell"
  scheme. A cell switches from inactive to active when it becomes full, which
  changes the local force balance abruptly. This is the likely reason turning
  it on made the model less stable.

Other relevant facts:

- `tpo%now%H_eff` exists (yelmo_defs.f90:302) but is never computed.
  `calc_H_eff` rebuilds H/f on the fly; with `set_frac_zero` it treats partial
  cells as ice-free (`z_srf`, `H_grnd`).
- SMB/BMB are not scaled by `f_ice` anywhere (`calc_G_mbal`).
- There is no path that returns excess calving or overshooting ice upstream.
- `calc_ice_fraction_lsf` returns a binary mask unless subgrid is on, so
  `f_ice_method = "lsf"` alone changes nothing (confirmed: GRL-8KM identical).

## 3. What CISM does now

From `~/models/CISM` (HEAD 71291a46), `which_ho_calving_front = 1|2`. CISM
removed the inactive-cell scheme in 2022 (c3633cde: "flickering and
oscillations when a cell changed from inactive to active").

- Front cells are marine cells (option 1: floating; option 2: also
  marine-grounded) with an ocean edge neighbour.
- `H_eff` = max over interior edge neighbours of min(H, H_flotation), minus
  `dthck_dx_cf`·distance; diagonal neighbours if no edge neighbour. Floor
  `thck_effective_min` (50 m), cap at flotation. Option 2 limits the effective
  surface slope and height.
- Partial if (H_nbr − H)/dist > `dthck_dx_cf`. `a_eff = min(H/H_eff, 1)`.
- **Partial cells are always active in the velocity solve**, with `H_eff`, the
  surface from `H_eff`, and the ocean-pressure BC on the partial cell's ocean
  face. Because of the floor, the active set never flickers. `f_ground` uses
  the true H.
- Transport is not masked. Afterwards, ice in cells beyond the front is
  returned to the upstream cells it came from (split by inflow); leftover
  calves.
- SMB/BMB in partial cells are scaled by `a_eff`.
- Calving is a lateral volume flux: dH = Cr·dt·H_eff·L_cf/(dx·dy), with the
  front length from the number of ocean faces. Demand that exceeds the column
  is carried to upstream calvable neighbours. The front advances only when
  H > H_eff (excess moved to downstream ocean neighbours).
- Land-terminating margins: no H_eff; cells below `thklim` inactive, surface
  gradient ramped over `thck_gradient_ramp` to avoid on/off jumps.

## 4. Proposal for Yelmo

The core idea is CISM's: a partial front cell is always part of the momentum
balance with an effective thickness, and only its area fraction changes. The
design below maps that onto Yelmo's C-grid.

### 4.1 Effective thickness and area fraction

- New option `ytopo.front_subgrid = "none" | "floating" | "marine"`.
  `margin_flt_subgrid` is retired (reading it is an error). "none" keeps binary
  `f_ice` (current behaviour); the target setting is `"marine"`.
- Front-eligible cells: floating (`"floating"`) or floating + grounded below
  sea level (`"marine"`), with H > 0 and at least one ocean face (ice-free,
  `z_bed < z_sl`).
- `H_eff` from upstream **full** edge neighbours: max of min(H, H_flot) minus
  `ytopo.front_dHdx`·dx; diagonal neighbours (distance √2·dx) if there is no
  full edge neighbour. Floor `ytopo.front_H_eff_min` (default 50 m), cap at
  flotation for floating cells. For marine-grounded cells, limit the effective
  surface as CISM option 2: `z_srf_eff` ≤ `z_srf` of the neighbour + 0.001·dx
  and ≤ `z_srf` + 25 m (if a limit applies, recompute `H_eff` and treat the
  cell as full).
- `a_eff = min(H/H_eff, 1)`; stored as `tpo%now%H_eff` (finally filled) and
  used as `f_ice` at fronts. Everywhere else `H_eff = H`, `f_ice = 1` (or 0).
- `f_ice_method` stays as the way `a_eff` is obtained: `"upstream"` (above) or
  `"lsf"` (geometric fraction from the level set, with `H_eff` still from the
  neighbours so that H/`a_eff` is bounded).
- `calc_H_eff` reads the stored field. The `set_frac_zero` call sites
  (`calc_z_srf_max`, `calc_H_grnd`, `scale_beta_gl_zstar`) are reviewed one by
  one: surface and base use `H_eff` in partial cells; `f_grnd` from the true H
  (as CISM).

### 4.2 Momentum balance

- `H_ice_dyn = H_eff`, `f_ice_dyn = 1` in partial front cells (in
  `calc_ytopo_diagnostic`, as the new default branch; `slab`/`slab-ext` are
  unchanged).
- Everything in the dynamics uses the `_dyn` fields: `set_ssa_masks` (now uses
  `f_ice`), `calc_lateral_bc_stress_2D` (now uses `H_ice`), driving stress,
  viscosity, beta, N (`hydro_calc_N` already takes `f_ice_dyn`).
- `calc_ice_front` is computed from `f_ice_dyn`, so the partial cell is the
  front cell and the lateral BC sits on its ocean face with `H_eff`.
- `fill_partial_ice_cells` / `fill_strain_2D_partial` become unnecessary for
  front cells (they are solved), but stay for other partial cells.

### 4.3 Transport and front advance

Yelmo already has what CISM had to emulate: `set_inactive_margins` zeroes the
ac-velocity on partial→ice-free faces before advection, so a partial cell
cannot export into the ocean and no ice is ever beyond the front. With
upwind/implicit advection a zero face velocity gives a zero flux, so the
negative-thickness problem CISM hit with edge masks does not arise here.

- Keep `set_inactive_margins`, with "partial" defined by `a_eff < 1`.
- Front advance, as CISM `advance_calving_front`: after calving, in a front
  cell with H > `H_eff` the excess plus a small amount (0.1 m) is moved to its
  ocean edge neighbours, split by outward flux. The new cell is a thin partial
  cell. The front therefore advances within the step, without a lag that
  depends on dt.
- The dynamics velocity on a closed face is kept for the momentum balance and
  diagnostics; only the advective flux is zero.

### 4.4 Mass balance

- `calc_G_mbal` takes `a_eff` and scales SMB/BMB in partial cells (the
  per-area rate applies to the covered area only).
- `calc_fmb_total` (methods 1–3) uses the stored `H_eff` for submerged depth
  and the same front-cell definition.

### 4.5 Calving (mass-balance path)

- Rate as lateral volume flux: dH = Cr·dt·`H_eff`·n_ocean_faces·dx/(dx·dy)
  (CISM's `cf_length` for 1/2/3 faces as an option) in front cells. The laws
  (`vm-l19`, `eigen`, `threshold`) supply only Cr.
- If dH > H, the column is removed and the remainder is applied to the
  upstream calvable neighbours (split by their inflow), within the same step
  (CISM `apply_calving_dthck`). This replaces the lost excess at the
  `apply_tendency` clip.
- `apply_calving_rate_thin` and `calc_calving_rate_tongues` are retired for
  this scheme (they exist to clean up thin cells that no longer arise).

### 4.6 Removal rules

- `calc_G_boundaries`: `H_min_flt`/`H_min_grnd` compare `H_eff` (not H/f);
  the "H_eff above neighbour maximum" cap is replaced by the `H_eff` cap in
  4.1.
- `calc_G_remove_fractional_ice`: a partial cell with no full edge or diagonal
  neighbour is an iceberg and is removed (at rate H/`H_min_tau`).

### 4.7 LSF path

- The level set sets the front position; it moves by (u + c)·dt with its own
  substepping, so a retreat of more than one cell per step needs no carry-over.
- `a_lsf` = area fraction of the cell behind the level-set front (marching
  squares, `calc_ice_fraction_lsf`), `H_eff` as in 4.1.
- Thickness follows the level set: in a front cell, ice above `a_lsf`·`H_eff`
  is calved (CISM subgrid calving mask, H/H_eff = 1 − mask; a cell with
  `a_lsf` below a small threshold is emptied). Cells beyond the front
  (`a_lsf` = 0) are emptied as now. If `a_lsf`·`H_eff` > H, the cell keeps its
  ice and fills by advection.
- `a_eff = H/H_eff` after this step, so `f_ice`, the momentum balance and the
  mass balance are the same as in the mass-balance path.

### 4.8 Not in scope here

- Thermodynamics in partial cells: keep the neighbour fill
  (yelmo_thermodynamics.f90:470).
- Land-terminating margins (Greenland): a surface-gradient ramp near a
  minimum thickness, like CISM's `thck_gradient_ramp`. Separate change.

## 5. New and changed parameters

| Parameter | Group | Default | Replaces |
|---|---|---|---|
| `front_subgrid` | ytopo | `"none"` (target `"marine"`) | `margin_flt_subgrid` (retired) |
| `front_H_eff_min` | ytopo | 50 m | – |
| `front_dHdx` | ytopo | 0.0 | – |
| `f_ice_method` | ytopo | `"upstream"` | (kept) |

`calv_thin`, `Hc_ref_thin` become unused with `front_subgrid /= "none"`.

## 6. Implementation plan (commits)

1. `H_eff`/`a_eff` computation and storage (`front_subgrid`, new
   `calc_ice_fraction` variant, `calc_H_eff` reads the field; remove the dead
   `H_lim` branch). No behaviour change with `"none"`.
2. Dynamics geometry: `H_ice_dyn`/`f_ice_dyn` for partial cells, masks and
   lateral BC from the `_dyn` fields, `calc_ice_front` from `f_ice_dyn`.
3. Transport: `set_inactive_margins` on `a_eff`, front advance with excess
   redistribution.
4. Mass balance scaling by `a_eff`; `calc_fmb_total` front definition.
5. Calving as volume flux with upstream carry-over; retire thin/tongue rules
   under the new scheme.
6. Removal rules (`calc_G_boundaries`, fractional/iceberg removal).
7. LSF path: thickness trimmed to the level-set area fraction.

Each commit keeps `front_subgrid = "none"` bit-identical to the current dev.

## 7. Tests

- Benchmarks with `"none"`: bit-identical (EISMINT-moving, TROUGH, CalvingMIP,
  MISMIP3D, symmetry check).
- CalvingMIP exp1/exp2 and MISMIP3D/TROUGH with `"marine"`: symmetry, front
  position vs protocol.
- ANT-16KM and GRL-8KM initmip, 1 ka, Courant 0.5 vs 0.8: Picard iterations,
  calving demand vs applied (the 2026-09 diagnostic), floating-area
  sensitivity, cost.

## 8. Decisions (2026-09-28)

1. All marine fronts: `"marine"` (floating + marine-grounded), with CISM's
   surface limits. Land-terminating margins stay a separate change (4.8).
2. `front_H_eff_min` = 50 m.
3. Excess calving is carried upstream within the step. With LSF it is
   implicit in the level-set motion.
4. Front advance as CISM: excess above `H_eff` moved to the ocean neighbours
   after calving.
5. `margin_flt_subgrid` is retired.
6. The LSF path is part of this effort.
