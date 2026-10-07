# Courant cap on transporting faces: design proposal

Status: proposal (2026-10-07), not implemented, design not yet confirmed.
Branch `cfl-transport` (off dev b7416276).

## 1. Problem

The Courant cap of the adaptive timestep (`yelmo.pc_cfl_max`, `set_adaptive_timestep_pc`
via `calc_adv2D_timestep1`; also `set_adaptive_timestep`) uses `dyn%now%ux_bar`/`uy_bar`
on every face, i.e. the max face speed around each cell. The thickness advection does not
use all of these faces: `calc_G_advec_simple` first calls `set_inactive_margins`, which
zeroes the velocity on faces between a partial cell (`f_ice` < 1) and an ice-free cell
that may not fill (`a_front` < `A_FRONT_MIN`). No mass crosses those faces, but their
speed still limits dt.

Case (GRL-8KM initmip from present day, 1 kyr, dev after merge 25250a38):

- Rink Isbræ, front cell (62,184) (71.7 N, 51.7 W): partial (`f_ice` ≈ 0.8, H ≈ 760 m,
  bed −847 m). Its ocean face to (61,184) moves at ~5 km/yr (cell centre ~2.8 km/yr,
  inflow face ~0.3 km/yr); (61,184) cannot fill (`a_lsf` = 0), so the face is closed for
  transport.
- That face sets the Courant cap: dt 0.73 yr on average vs 0.97 yr in the run before
  25250a38, where the front had advanced into the land-walled bend and the glacier
  stagnated (~250 m/yr). The cap binds in all steps (dt_now = dt_pi, which includes the cap).
- Earlier baseline (dev 39af9511): the binding cells were also partial front cells next
  to ice-free ocean (e.g. (131,110), f 0.45, face 3.2 km/yr).

## 2. Proposal

- New routine in `yelmo_topography`, e.g.

      subroutine calc_transport_velocity(ux_t,uy_t,tpo,dyn,bnd)

  copies `dyn%now%ux_bar`/`uy_bar` and applies
  `set_inactive_margins(ux_t,uy_t,tpo%now%f_ice,tpo%par%boundaries,a_front)`, with
  `a_front` from `calc_lsf_area_fraction(a_front,tpo%now%lsf,tpo%now%H_ice,bnd%z_bed,bnd%z_sl,...)`
  when `use_lsf` and `front_subgrid /= "none"` (absent otherwise). This is the same face
  rule as `calc_ytopo_pc` uses for the advection (a_front logic there: lines with
  `if (allocated(a_front)) call calc_lsf_area_fraction(...)`).
- `yelmo_ice` (adaptive timestep block, ~line 170): pass `ux_t`/`uy_t` to
  `set_adaptive_timestep` and `set_adaptive_timestep_pc` instead of `dyn%now%ux_bar`/`uy_bar`.
- Only faces that the advection closes are excluded. No thickness threshold (an open face
  with a thin upwind cell still counts).

Expected effects:

- `front_subgrid = "none"`: `f_ice` is binary, `set_inactive_margins` closes nothing
  (faces from full cells stay open), so benchmarks without the subgrid front should be
  bit-identical.
- Land margins are unaffected (full cells).
- Partial front cells: the ocean face no longer caps dt unless the cell beyond can fill.

## 3. Open points to confirm

1. Exclude closed faces only (proposal) vs also faces with a negligible upwind thickness.
2. Use the current velocity (as now) or the advected one (`pc_filter_vel`: mean of current
   and previous solution) for the cap. Proposal: current, as now.
3. The face rule depends on `f_ice`/`a_front` at the start of the step; the predictor can
   open a face (cell beyond starts to fill). The PC error (`pc_eta`) should catch that;
   check redo counts.
4. One routine for both uses (advection + cap) would avoid duplicating the `a_front`
   logic; `calc_G_advec_simple` currently applies `set_inactive_margins` internally.

## 4. Tests

Levante baseline clone with dev code: `/work/ba1442/robinson/models/yelmo-fl/merged`
(outputs `output/fl/<case>`: rst, grl16, trough8, eis). Scripts in
`/work/ba1442/robinson/models/yelmo-fl/`: `fl_build.sh <clone>`, `fl_submit.sh <clone> <case>`
(cases rst, grl8, ant16, grl16, eis, halfar, trough8, maskice, mismip3d, calv1, homc, homf),
`fl_summary.py <clone> <cases>` (vs `base`; edit the base path to `merged`),
`dtcent.py`, `dtlim.py` (Courant-min cells from the restart). Delete the clones when done.

- GRL-8 1 kyr from PD: dt should rise towards ~1 yr; check the Rink front and volume.
- GRL-8 restart case (`rst`), GRL-16, ANT-16 1 kyr: dt, `ssa_lim_n`, `iter_redo`, V/A.
- Benchmarks: expected bit-identical where `front_subgrid = "none"` or no partial cells;
  EISMINT symmetry.
