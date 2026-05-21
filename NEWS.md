# Yelmo release/tag notes

## v2.0 (2026-05-21)

First release under the new `fesmc/yelmo` repository home. v2.0 collects the substantial development that has happened since the v1.0 release described in the 2020 GMD model description paper.

### Dynamics & velocity solvers

- **DIVA solver added and now default.** Other solvers available: `"sia"`, `"ssa"`, `"hybrid"`. An `"l1l2"` solver also exists but is broken, not recommended, and will not be developed further.
- **Two SSA discretizations**: `ssa_solver = "residual"` (Yelmo's traditional FD formulation) and `ssa_solver = "energy"` (new energy-functional FE form, as in Yelmo.jl). Files `solver_ssa_ac.f90` and `solver_ssa_ac_energy.f90` replace the old `solver_ssa_sico5.F90`.
- New `velocity_general.f90` / `velocity_ssa.f90` / `velocity_diva.f90` consolidate shared velocity-solver machinery versus solver-specific routines, and a generic `solver_linear.F90` wraps Lis.
- New `&ydyn` knobs: `uz_method`, `visc_method`/`visc_const`, `eps_0` (strain-rate regularization), `scale_T`/`T_frz` (cold-ice friction scaling), `ssa_lat_bc`, plus separate Lis option strings per SSA formulation.

### Topography, calving & geometry

- **Level-set calving front** (`src/physics/calving/lsf_module.f90`) with two modes (`snap` and `redist` Sussman/Osher), controlled by a new `&ycalv` group split out from `&ytopo`.
- Calving rewritten into a dedicated `src/physics/calving/` package (`calving_aa.f90`, `calving_ac.f90`) with new laws: `vm-l19`, `vm-m16` (lsf-only), `stress-b12`, separately controlled for floating vs grounded fronts.
- **Sub-grid discharge mass balance** (`src/physics/discharge.f90`, new `dmb_*` parameters in `&ytopo`).
- **Frontal mass melt** (`fmb_method`, `fmb_scale`) for shelf-front melting separate from basal melt.
- Grounding-zone basal melt: new `bmb_gl_method` (`fcmp`/`fmp`/`nmp`/`pmp`/`pmpt`) with grounding-zone penetration controls (`gz_Hg0`, `gz_Hg1`) for the latter.
- New geometry controls: `dHdt_dyn_lim`, `grad_lim_zb`, `margin_flt_subgrid`, `f_ice_method` (upstream vs LSF area fraction), `topo_rel_field`.

### Thermodynamics

- **Active bedrock layer**: new `&ytherm` rock block (`rock_method`, `nzr_aa`, `H_rock`, `cp_rock`, `kt_rock`, `zeta_scale_rock`) for lithospheric heat conduction.
- **Quadrature-based basal heating**: `qb_method = 2` (quadrature) added alongside the simple staggered estimate.

### Time stepping

- **Improved PC controller**: new options `pc_use_H_pred`, `pc_filter_vel`, `pc_corr_vel`, `pc_n_redo` (limits repeated retries on a single step) and `disable_kill` (let the model keep going through instabilities for diagnostics). Default flow now uses predicted thickness (`pc_use_H_pred = True`) + velocity filtering (`pc_filter_vel = True`). Error metric `pc_eta` now scaled by both an absolute and relative tolerance parameter (currently hard coded).

### Boundary conditions, masking, coupling

- New `bnd%mask_ice` integer mask framework, replacing the older `ice_allowed`/`tau_relax==0` patchwork. Facilitates regional modeling, which is now fully supported.
- **C API for external coupling**: `src/yelmo_c_api.f90` exposes Yelmo to non-Fortran callers.
- Restart machinery: per-field opt-in (`restart_z_bed`, `restart_H_ice`, `restart_relax`) and explicit relax-from-restart-to-input.

### Parameters & configuration

- `&ycalv`, `&ytill`, `&yneff` groups split out of `&ydyn`/`&ytopo` for cleaner separation.
- Top-level `&yelmo` block now supports overriding which namelist groups feed each subcomponent (`nml_ytopo`, `nml_ycalv`, etc.) — lets one config file drive multiple Yelmo instances with different physics.
- `phys_const = "Earth"` makes the constants set explicit/swappable.
- Default `&ymat` enhancement factors raised from `enh_shear = enh_stream = 2.0` to `3.0`; `de_max` raised from `0.5` to `2.0`.

### Testing & tooling

- Many new unit tests in development: `test_levelset.f90`, `test_roots.f90`, `test_ssa_energy*.f90`, `test_variables.f90`, plus driver-level test cases for `calving`, `ismiphom`, `mask_ice`, `slab` - not all working. Symcheck (`symcheck.jl`, `README_symcheck.md`) for symmetry-based regression checking.
- `&opt` namelist group (basal-friction optimization driver) integrated.
- New `runme`/`.runme/` workflow tooling, `scripts/` directory, batch helpers (`run_calvingmip.sh`, `run_unit_tests.sh`).
- Host configuration now driven by `config.py`.
- External shared utilities (`ncio`, `nml`, `interp1D`, `gaussian_filter`) moved out to the `fesm-utils` package and consumed as a dependency.

### Documentation

- New Quarto-based documentation site under `docs/`, published at https://fesmc.github.io/yelmo/. Replaces the previous https://palma-ice.github.io/yelmo-docs site.
- Per-component variable references (`yelmo-variables-{ybound,ydata,ydyn,ymat,ytherm,ytopo}.md`) plus dedicated pages for physics, I/O, optimization, remapping, benchmarks, and HPC notes.

### Repository

- Repo migrated from `palma-ice/yelmo` to `fesmc/yelmo`; v2.0 marks the first release under the new home.

## v1.15 (2026-01-19)

- Use of Gaussian Quadrature module from fesm-utils for calculating the Jacobian of velocity (strain-rate tensor), vertical velocity, dynamic viscosity (DIVA, SSA), basal friction, and other quantities.
- Added `uz_lim` to vertical velocity (improves stability for edge cases).
- Added new parameter `ytopo.dHdt_dyn_lim` to be able to limit rate of change of ice thickness due to dynamics can be (c
an help with stability).
- Implementation of switches to test different staggering methods (simple staggering versus Gaussian Quadrature, etc.): `ytherm.qb_method` and `ydyn.uz_method`.
- Implementation of LSF for calving, including CalvMIP test cases.
- Separation of calving parameters from topo parameters in namelist groups.
- Converted all further instances of get_neighbor_indices to get_neighbor_indices_bc_codes. Overall change led to spee
dup of 10% on a 16km Greeland run.
- OpenMP improvements means significant speedups are now possible for high-resolution runs.
- calc_bmb_total: bug fix; removed all traces of grounded_melt parameter, which was no longer used, and also removed optional argument mask_pd.
- Introduced `ydyn.scale_T` and `ydyn.T_frz` to control a linear reduction in friction until cf_ref in the case that ice is frozen at the base. This should make basal velocities more consistent with expectations, even when background friction is artificially low.
