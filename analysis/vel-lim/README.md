# Velocity limit tests (TROUGH-F17)

Scripts used to test `ydyn.ssa_vel_lim_method = "drag"` against the old
`"clip"` in TROUGH-F17 (2026-10-02). The method and the results are
described in `docs/physics/momentum/solvers.md` (section "Velocity limit").

The runs were made on Levante in a throwaway clone
(`/work/ba1442/robinson/models/yelmo-vellim`, branch `vel-lim-drag`), with
4 OpenMP threads on the shared queue. The fetch scripts copy the output
from there into `./data` (not tracked). The plots go to `./plots` (not
tracked).

## Runs

Common: `runme -s -r -q shared --omp 4 -w 12:00:00 -e trough -n par/yelmo_TROUGH-F17.nml -p yelmo.log_timestep=True ...`
(at the time `par/yelmo_TROUGH-F17.nml` had `"clip"` and `ssa_vel_max = 5000`).

0–8 kyr, 2D output every 200 yr (`vl_fetch.sh`, `vl_analyse.jl`):

| run | settings |
|---|---|
| `dev_clip` | dev dfb9fdff executable (`-e libyelmo/bin/yelmo_trough.x.dev`) |
| `br_clip` | branch, `ssa_vel_lim_method = "clip"` |
| `br_drag` | branch, `"drag"`, `ssa_vel_max = 5000` |
| `br_drag6k`, `br_drag8k`, `br_drag10k` | branch, `"drag"`, `ssa_vel_max = 6000, 8000, 10000` |

0–6 kyr, 2D output every 10 yr (`ctrl.dt2D_out=10`; `vl_ts_fetch.sh`, `vl_centreline.jl`):

| run | settings |
|---|---|
| `ts_clip` | `"clip"`, `ssa_vel_max = 5000` |
| `ts_drag5k` ... `ts_drag50k` | `"drag"`, `ssa_vel_max = 5000, 6000, 8000, 10000, 15000, 20000, 30000, 50000` |

`vl_ts_fetch.sh` extracts the centreline (`yc = 0`) with `ncks` on
Levante before copying.

## Scripts

- `vl_analyse.jl`: solver statistics per activation window (steps, steps
  at `dt_min`, Picard iterations, solves at `ssa_iter_max`, linear-solver
  failures, pc_eta), bit-identity of `br_clip` and `dev_clip`, maximum
  speeds and surge timing and volume.
- `vl_centreline.jl [x_km]`: time series of `uxy_bar` and `H_ice` on the
  centreline at `x_km` (default 300). Writes
  `plots/vl_centreline_x<x_km>.png` (copied to
  `docs/img/vel-lim-trough-x300.png`).

Usage (with an environment that has NCDatasets and CairoMakie, e.g.
`--project=..` after `Pkg.instantiate()`):

```bash
./vl_fetch.sh && julia --project=.. vl_analyse.jl
./vl_ts_fetch.sh && julia --project=.. vl_centreline.jl 300
```
