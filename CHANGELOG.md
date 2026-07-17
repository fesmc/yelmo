# Yelmo changelog

## v2.3.1 (2026-07-17)

- **Fix ISMIPHOM aborting at startup.** `par/yelmo_ISMIPHOM.nml` still set
  `ydyn.solver = "l1l2"` after the L1L2 solver was deleted in v2.3, and
  `yelmo_check_enum` stops hard on an unknown value. Switched to `diva`; note this
  changes ISMIP-HOM benchmark results, as DIVA is a different approximation.
  Stale `l1l2` doc comments removed across `par/` and `input/yelmo_defaults.nml`
  (no parameter values changed).
- **Refresh the bundled `yelmo-config` defaults/enums snapshot**, which predated the
  v2.3 `&ytrc` refactor — installs without a checkout served pre-refactor parameters.
- **Fix `yelmo-config snapshot` writing to the installed package** instead of the
  checkout, which is how the snapshot went stale unnoticed. Added a drift check
  pinning the committed snapshot to `input/yelmo_defaults.nml`.

## v2.3 (2026-07-15)

> **Note:** ISMIPHOM aborts at startup in this release (`ydyn.solver = "l1l2"`
> outlived the L1L2 solver). Fixed in v2.3.1 — use that instead. Other
> configurations are unaffected.

Enthalpy is now the default thermodynamics solver, a passive-tracer subsystem is
added, and the horizontal thermal advection is rewritten. This release also folds
in the fesm-utils/coords compatibility shift and a large dynamics/topography
bug-fix audit.

### Thermodynamics: enthalpy solver by default

- **`ytherm.method = "enth"` is now the default**, replacing the temperature solver.
  The enthalpy formulation is validated against the Kleiner et al. (2015) Experiment A
  and B polythermal benchmarks (vendored data + standalone column driver) and against
  EISMINT-2 A/F, where enth now matches the temp solver.
- Numerous enthalpy fixes en route to a robust 2D solver: basal mass balance in
  `calc_enth_column`, `Q_ice_b` sign, CR-dependence of basal melt at a melting base,
  cross-CTS conduction with the cold-side conductivity (fixes 2D margin NaNs),
  `Q_b`/`Q_lith` unit conversion to `J m-2 a-1` (fixes a cold-base bias in 2D), and a
  Péclet-hybrid upwinding scheme for vertical enthalpy advection.
- Dead poly-enthalpy solver and unused thermodynamics helpers removed.

### Horizontal thermal advection

- **Conservative flux-form van Leer / MUSCL horizontal advection** with a flux-limited
  2nd-order scheme (`advecxy_order`, default 2) and adaptive sub-cycling for stability.
- **Floating-base thermal coupling to the ocean** ("conducting freezing base") — fixes
  shelf-front speckle and N–S asymmetry in the trough. A reflection-axis selector was
  added to the offline `test_symmetry` gate.
- New `H_ice_thin` thin-ice threshold exposed as a parameter; `uz_star` used
  consistently for age/enhancement vertical advection in sigma coordinates.

### Passive-tracer subsystem (`%trc`)

- New `ytrc_class` passive-tracer framework with three parallel backends
  (`euler` / `tracer` / `elsa`), including restart support and staggered-velocity
  API. Elsa layer stacks and Lagrangian particle clouds are harmonized onto a
  gridded deposition-time field; isochrones are keyed by deposition time (`time_iso`),
  and transient gridded tracer stats (`trc_count`, `trc_depth_iso`) are emitted.
- `elsa`/`tracer` registered as nested yelmo dependencies in `configme`, driven via
  nested defaults files.

### Dynamics, topography & mass conservation bug-fix audit (2026-07)

A systematic audit fixed ~25 verified bugs across the solvers, including:

- **Velocity / basal drag**: Picard L2 residual (`norm_method=1`) missing `sqrt` and
  per-direction masks; per-node `cbn` in power-plastic friction; grounded-end test in
  GL flux weighting; last-column coverage in `scale_beta_gl_fraction`; zero-gradient
  x-border beta BC; `taud_gl_method=-1` now a hard error stop.
- **SIA**: down-neighbour stress terms use `k-1`; OMP `private` clause corrected.
- **Deformation / rate factor**: `T_prime` clamped to the homologous cap `T0`; periodic
  BC writes the last column; `dyy` margin stencil j-index fix.
- **Time stepping**: RK23 advances with the 2nd-order solution; RK4 truncation error
  made cumulative; `pc_k=1` for FE-SBE; adaptive `dt` no longer rounds to zero / int
  overflow.
- **Topography / calving / mass conservation**: `pi->phi` typo in the fraction-above-zero
  check; inner-column copy for the infinite BC; error-stops on uninitialized `calv_now`
  (Eigen calving) and unimplemented acy-node GL flux; BC initialized before neighbour-index
  lookups in the LSF calving routines.
- **Staggering / smoothing**: `stagger_aa_ab_ice` averages over iced corners only;
  `smooth_gauss_2D` halo/corner fixes; unified periodic wrap offset in `set_boundaries_2D_aa`.

### fesm-utils / coords compatibility shift

**Yelmo now requires `fesm-utils` dev at `3f415cc` (2026-06-26) or later.** Building
against an older fesm-utils will fail to compile. fesm-utils folded its standalone
`coordinates` library into `utils/src/coords/` and moved several symbols:

- `mv` / `TOL` / `TOL_UNDERFLOW` moved from module `precision` to a new `constants` module.
- `nc_read_interp` moved from `mapping_scrip` to `ncio_interp`, and its mapping argument
  changed from `mps=(map_scrip_class)` to `map=(map_class)`.

`yelmo_io` was migrated to the unified `map_class` / `map_read` path (restart interpolation).
A new `restart_interp_gen` switch selects how the conservative restart map is built:
`"cdo"` (load a pre-generated SCRIP map; default, unchanged behaviour) or `"coords"`
(generate the weights in-package via the coords library, no cdo dependency and no map file).
The domain grid now uses the fesm-utils/coords `grid_class` directly.
`FastHydrology` is now built as its own dependency (its compile rule moved into yelmo);
the bundled `FastHydrology` and `FastIsostasy` libraries need no source changes but must be
rebuilt against the same fesm-utils.

## v2.2 (2026-06-18)

Parameter-load validation, a `yelmo-config` CLI, and the `&yhyd` basal-hydrology
namelist migration.

- **Parameter validation at load time**: enum / range / ordering checks across every
  `*_par_load`, file-existence checks, and a loud failure when `domain`/`grid_name` are
  still `"None"`. Includes an enum-check audit that fixed `bmb_gl_method`, `ssa_lat_bc`,
  the `ytopo` solver, `calv_flt_method`/`cf_min`, and therm/material/calv_grnd dispatch
  mismatches; dead `calvingmip` and paleo-shear params dropped.
- **`yelmo-config` tool**: a parameter-management CLI with a bundled defaults+enums
  snapshot (works without a checkout), a self-update command, and curated constraints.
- **Basal hydrology → `&yhyd`**: namelists migrated and canonicalized; FastHydrology
  `mdot` fed in SI m/s; FastHydrology compile rule moved into yelmo as its own dependency.
- Config reorganization: all parameter files consolidated into `par/` (from `par-gmd/`),
  obsolete configs moved to a `legacy/` subfolder.

## v2.1.3 (2026-06-15)

- **ISMIP7 retreat law** for marine-terminating glaciers, implemented as a calving law
  (merge of `calving-blasco-dev`, PR #3). Working Greenland frontal-retreat configuration;
  additional calving diagnostic output fields for 2D analysis.

## v2.1 (2026-06-11)

FastHydrology coupling and the integer `mask_ice` convention.

- **FastHydrology basal hydrology integrated**: build wiring, run each step into
  `dom%hyd`, and I/O + restart surfacing. New `&fhyd` namelist block. The old `&yneff`
  group and the `thrm%H_w` bucket are retired; `neff_method = 6` plumbs `hyd%N` into
  `dyn%N_eff`. Synced with the FastHydrology `W_til`/`W_til_max`/`overflow` rename.
- **`mask_ice` convention** switched to `{0 = NONE, 1 = FIXED, 2 = DYNAMIC}` via named
  constants (`MASK_ICE_NONE`/`FIXED`/`DYNAMIC`).
- **Build system**: `config/common.mk` splits out dependency wiring for `configme`
  (`FASTHYDROROOT`, FFTW); `runme` is now a pip-installed package rather than a bundled
  script.
- Docs restructured to lead with a quick start; `initmip` default `time_end` lowered from
  20 kyr to 1 kyr.

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
