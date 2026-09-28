# Yelmo changelog

## Unreleased (dev)

Mostly fixes from a whole-source audit (2026-09). Most default runs change a
little. MISMIP3D and DIVA runs change more.

### Changes that affect existing par files

- **`&yhyd marine_rho_sw` and `k24_latent_heat_water` are no longer read under Yelmo.**
  FastHydrology takes `rho_sw` and `L_ice` from the physical-constants record Yelmo
  hands it, so these two keys are inert in a coupled run (they still work for
  FastHydrology's standalone drivers). The keys stay in `yelmo_defaults.nml`, so no
  par file needs editing; a value set there is simply not the source of truth.

- **`ydyn.pc_corr_vel` removed** (it was unused). `nml_validate` stops on unknown
  parameters, so delete it from external par files.
- **`ydyn.ssa_lat_bc = "floating"` now means floating fronts only.** The producer
  and consumer of the ice-front mask disagreed on its codes, so `"floating"` acted as
  `"marine"` (floating and grounded marine fronts). The codes are now shared
  (`MASK_FRNT_*` in `yelmo_defs`). The default and all shipped par files use
  `"marine"`, which keeps previous results. External par files that set `"floating"`
  change behaviour at grounded marine fronts: switch to `"marine"` to keep it.
- **Benchmark physical constants follow the published protocols**
  (`input/yelmo_phys_const.nml`). EISMINT uses ρ_ice = 910. MISMIP3D uses g = 9.8,
  ρ_sw = 1000 and a 365-day year. There are new `MISMIPplus`, `ISMIPHOM` and
  `CALVINGMIP` groups. MISMIP3D changes a lot: Stnd at 16.1 ka has x_gl 400→500 km
  and volume +33%. EISMINT H changes by +0.27% and ISMIP-HOM by about −1%.
- **CalvingMIP uses constant N_eff = 1 Pa again** (`bkt_N_closure = 0`,
  `const_N = 1`). Since the `neff_method` retirement it had been mapped to the
  overburden closure with `cf_ref = 1e4`, which froze the bed.
- **`k24_eta_w` in the par files is now SI** (5.70e-11 Pa s). The old value was a
  stale per-year one, about 3e7 times too large.
- **Removed parameters:** `yelmo.cfl_diff_max` (it was unused, and so was the
  placeholder `dt_diff` output) and `ydyn.cb_sia` (it was read, but its code block
  was empty). Delete them from external par files.
- **`ytopo.margin_flt_subgrid` is replaced by `ytopo.front_subgrid`** ("none",
  "floating" or "marine"; default "none" = the previous `False`), with
  `front_H_eff_min` (50 m) and `front_dHdx` (0). Front cells get an effective
  thickness `H_eff` from their thickest interior neighbour, following the CISM
  subgrid calving front, and `f_ice = H/H_eff`; `tpo%now%H_eff` is now filled.
  Replace `margin_flt_subgrid` in external par files.
- **Partial front cells take part in the dynamics** (step 2 of
  docs/dev/front-subgrid-design.md). The geometry for the surface, gradients,
  front mask and velocity solver is the active ice column: every ice-covered
  cell, with `H_ice_dyn = H_eff` and `f_ice_dyn = 1`. The front boundary
  condition is on the partial cell's ocean face. `H_grnd` uses the actual
  thickness. With `front_subgrid = "none"` results are unchanged.
- **The ice front advances as in CISM** (step 3): after calving, a front cell
  with more ice than `H_eff` passes the excess (plus 0.1 m) to its ice-free
  ocean neighbours, split by the outward velocity. Partial front cells still
  export nothing into the ocean (`set_inactive_margins`). Only with
  `front_subgrid /= "none"` and mass-balance calving.
- **Mass balance and calving at subgrid fronts** (steps 4–5, with
  `front_subgrid /= "none"`). SMB and BMB in partial cells act on the covered
  area only (scaled by `f_ice`). Floating calving demand (`vm-l19`, `eigen`,
  `threshold`) is scaled by the front length (1, √2, 2 for 1, 2, ≥3 ocean faces),
  and demand beyond a front cell's ice is taken from its upstream neighbours,
  split by inflow (CISM `apply_calving_dthck`); before, it was lost at the
  clip. The thin-ice and tongue calving rules are not used with the subgrid
  front.
- **Removal rules at subgrid fronts** (step 6). `H_min_flt`/`H_min_grnd` compare
  the stored `H_eff`. The cap of margin cells at their thickest neighbour is
  not applied to subgrid front cells (the front advance handles them; land
  margins keep it). A partial cell is removed as an iceberg when it has no
  full edge or diagonal neighbour (before: edge only).
- **Subgrid fronts with the level set** (step 7, `use_lsf` with
  `front_subgrid /= "none"`). Front-cell thickness follows the level set: cells
  with less than 10% of their area behind the front are emptied, and front cells
  (also those touching the ocean at a corner) hold at most `a_lsf`·`H_eff`
  (CISM subgrid calving mask), repeated up to 3 times. `f_ice = H/H_eff` as in
  the mass-balance path. `ytopo.f_ice_method` is removed (it had no effect with
  `"none"`); `calc_ice_fraction_lsf` becomes `calc_lsf_area_fraction`.
- **`ydyn.ssa_lat_bc = "slab"` and `"slab-ext"` are removed** (and
  `extend_floating_slab`). They were only set in the ISMIP-HOM and SLAB-S06
  pars, whose periodic domains are fully ice-covered, so they had no effect;
  those pars now use "marine".
- **New parameters:**
  - `yelmo.log_mb_check` (default false) prints a global mass-budget check every
    step. It replaces the hard-coded `check_mb`.
  - `ytopo.slope_bg_x` and `ytopo.slope_bg_y` (default 0) add a uniform background
    slope to the surface and bed gradients. They are for periodic domains whose
    geometry is tilted; the tilt itself is not in `z_srf`/`z_bed`.
  - `yelmo.pc_cfl_max` (default 0.5) is the Courant-number cap on the
    predictor-corrector timestep. It was hard-coded to 0.5, but the half-step
    rule in `limit_adaptive_timestep` also acted on this cap and halved it, so the
    effective value was 0.25. The rule now only acts on the time left in the
    call. The effective cap is therefore now 0.5 (was 0.25), which changes
    results. At 0.8 or 1.0, ANT-16/ANT-8 gain ice systematically (ANT-16 at 2 ka
    with 0.8: +40e3 km3, in fast grounded ice and at grounding lines) and GRL-8
    fits observations worse with 1.0.
  - `yelmo.pc_eta_H_min` (default 10 m, was hard-coded), `yelmo.pc_eta_u_min`
    (default 0) and `yelmo.pc_eta_trim` (default 0) define which points enter the
    predictor-corrector error norm: thinner or slower ice is left out, and
    `pc_eta_trim` drops that fraction of points with the largest errors, so a few
    flickering cells cannot set the timestep alone. Defaults keep previous results.
  - `ycalv.H_min_tau` (default 10 yr): margin ice thinner than `H_min_flt` /
    `H_min_grnd` and isolated partial cells are removed at the rate H/`H_min_tau`
    (all of it when dt ≥ `H_min_tau`). They were removed completely every step,
    so the removal per year grew with the number of steps: in ANT-16KM initmip
    (`H_min_flt = 75` m) it was the largest sink at the ice front, and Courant
    0.8 instead of 0.5 gave +40e3 km3 after 2 ka (now +3e3 km3 after 1 ka).
    `H_min_tau = 0` gives the previous behaviour; benchmarks are unchanged.
  - `ytopo.dHdt_dyn_lim` is removed, with the tendency limit in
    `apply_tendency`. It clipped the dynamic thickness change at ±100 m/yr
    cell by cell, which does not conserve mass. In ANT-16KM/GRL-8KM initmip it
    acted every step for 1 kyr, at fast marine fronts and grounding lines, not
    only at the start; without it the runs are as stable, with the same steps
    and SSA iterations (ANT-16 200 yr: +11e3 km3 ice, more calving).
  - `ycalv.tau_ice` is split into `tau_ice_flt` and `tau_ice_grnd` (both
    250 kPa): the `vm-m16` ice strength for floating and marine-grounded fronts
    (Morlighem et al., 2016 use separate values).
  - `ycalv.H_min_flt` default and initmip value 10 m (was 75 m). With
    `H_min_tau = 10` yr, ANT-16KM initmip at 1 ka: +23e3 km3 ice, +1.2% floating
    area, and ~1.6× faster.
- **TROUGH-F17 and MISMIP3D use `pc_eps = 1e-2`** (was 1.0). With 1.0 the
  controller let dt reach 5 yr during fast flank sliding, where a lateral mode
  grew about 1e5-fold from single-precision round-off: TROUGH (8 km) was up to
  12 m mirror-asymmetric at t = 200 yr and 6 m at 3.2 ka, now ≤ 2e-3 m, at the
  same cost. TROUGH volume at 5 ka +1.1%; MISMIP3D x_gl 520.17 → 520.08 km.
- **`yelmo.pc_filter_vel` (default true) changed meaning.** The velocity
  solution is no longer filtered. The thickness update is advected with the mean
  of the current and previous solutions, a true two-step mean (before, it was a
  running average over all past steps). All velocity outputs are now mutually
  consistent. TROUGH mean H changes by +1.4%, mostly from a shift in surge timing.
- **Periodic benchmark grids:**
  - TROUGH-F17, MISMIP+ and MISMIP3D use a y-grid centred on y = 0 with
    `ny = ly/dx` rows, so the period is exactly ly. Before, both walls were rows,
    which made the channel one cell too wide. The run stops if ly/dx is not an
    integer.
  - SLAB-S06 and RAYMOND are now periodic in both directions for every component,
    and carry their bed tilt as `slope_bg_x`.
  - TROUGH (8 km, 5 ka) volume per unit width changes by −3.7%; MISMIP+ by about −0.5%.
- **Benchmark drivers:**
  - MISMIP+ runs again.
  - SLAB-S06 uses constant N_eff = 1 Pa and no thermal scaling of c_bed. Before,
    the flow was about 1e-6 m/yr.
  - The MISMIP driver honours `ctrl.time_end > 0` (≤ 0 runs to the protocol end)
    and reports `x_gl` on the centreline and a new `x_gl_edge` at the channel wall.
  - The TROUGH and MISMIP drivers initialise the level set from the ice thickness.
    Before, all marine ice calved at t = 0 when `use_lsf` was on.
  - The TROUGH driver writes a restart.
  - The `yelmo_slab.x` program and the `make slab` target are retired. SLAB-S06
    runs through `-e trough`.
- **Output file names follow the yelmox/CLIMBER-X convention:** `yelmo2D.nc` →
  `yelmo.nc`, `yelmo1D.nc` → `yelmo_ts.nc`, regional `yelmo1D_<name>.nc` →
  `yelmo_ts_<name>.nc`, and initmip `yelmo2Dsm.nc` → `yelmo_sm.nc`. This applies
  to all test drivers, the default region file names, and the analysis scripts.
- **Requires FastHydrology dev ≥ `905a81d`**. It adds `hydro_calc_N`,
  `hydro_init_state` taking `H_ice`, and optional `periodic_x`/`periodic_y` in
  `hydro_init`.

### Answer-changing fixes

- **The hydrology converts its per-year inputs with the domain's own year.**
  `yelmo_hydrology` divided `bmb_w`, `uxy_b` and `A_glen_b` by FastHydrology's
  module-level `SEC_PER_YEAR` (3.1556926e7), and FastHydrology used the same
  constant for `dt_sec` and for the `bkt_till_rate` m/a → m/s conversion, so
  `&Earth sec_year` never reached the hydrology at all. All four now use
  `bnd%c%sec_year`, which `yhyd_par_load` also passes to `hydro_init`. Bitwise
  identical for every group whose `sec_year` is the CF/UDUNITS year (EISMINT,
  500 yr: all 246 restart and 111 output variables unchanged). **MISMIP3D**
  declares 3.1536e7 and runs the bucket (`method_til = 1`), so its SI till rate
  and `dt_sec` move by 0.066%; over 500 yr of RF that reaches exactly one field,
  the diagnostic `hyd_dW_til_dt` (6.6e-4 relative), with `W_til`, `N`, `p_w` and
  every ice field unchanged -- the bucket has no source there, since
  `bmb_grnd = 0`.
- **Basal water and freshwater flux:** the ρ_ice/ρ_w ratio was inverted in the basal
  water predictor and the freshwater flux.
- **N_eff is evaluated on the current dynamics geometry.** Since `neff_method` was
  retired, N_eff came from hydrology computed on geometry 1–2 steps old. Newly iced
  or newly grounded cells therefore had N = 0 (free sliding). `calc_ydyn_neff` now
  calls FastHydrology's state-free `hydro_calc_N` on `H_ice_dyn`. K24 and external N
  still use `hyd%now%N`. FastHydrology now starts N at overburden, so K24 no longer
  starts from N = 0. TROUGH mean H changes by −1.4%, and the "beta appears to be
  zero" warnings in CalvingMIP and EISMINT with SSA/DIVA are gone.
- **DIVA:** β_eff is now computed on ac-nodes from staggered β and F2. It was
  staggered from aa-nodes, which could reverse u_b, and no-slip gave 1/F2 = Inf.
  This affects every DIVA run: TROUGH uxy +1.4%.
- **Vertical grid:** the quad3D strain-rate `dzx`/`dzy` used the wrong vertical
  faces. The enthalpy conductivity and internal-melt layer thickness used stale
  `nz_ac = nz_aa-1` indexing. The sign of the depth term in the shelf-base freezing
  point was wrong. Enthalpy `Q_ice_b` is now output in mW m-2.
- **Stress and strain:** the 2D stress used the previous viscosity. Strain rates
  are reset at partial and ice-free cells.
- **Ice fronts facing land:** a grounded-below-sea-level ("marine") or floating
  front cell got the ocean-front stress condition on all its ice-free faces, even
  where the ice-free neighbour is bedrock above sea level (nunataks, fjord walls,
  holes). With `ssa_lat_bc = "marine"` this pushed thick trough ice into such
  holes and made them flicker between empty and refilled. The ice-free side of a
  front is now marked ocean (`MASK_FRNT_ICE_FREE`, −1) or land
  (`MASK_FRNT_ICE_FREE_LAND`, −2), and a face to land is treated as a front
  grounded above sea level, with no water back-pressure. GRL-8KM (200 yr,
  `pc_eps = 0.01`): the flickering cells stay ice-free and the run takes 1226
  instead of 1433 steps; the remaining error is dominated by one fast marine
  outlet front. Benchmarks are unchanged.
- **Partially ice-covered cells** get the 2D strain rates and `visc_bar` of their
  fully ice-covered neighbours (they were zero). Before, `calc_eps_eff` and
  `calc_tau_eff` patched this separately, and the `vm-m16` calving law saw zero
  stress, so partial front cells never calved. This changes von Mises and eigen
  calving runs (Antarctica 32 km, 1 ka, `vm-l19`: shelf volume −0.6%).
- **Hydrology coupling:** K24 now receives `uxy_b` and `A_glen` in SI units.
  FastHydrology's ρ_ice, ρ_w, ρ_sw and g now come from Yelmo's domain constants.
- **Eulerian tracer advection** read values already updated earlier in the same
  sweep. EISMINT age asymmetry is now 2.4 yr (was 78 yr).
- **Periodic boundaries:** there is now one true-wrap convention everywhere (period
  n, no halo; `periodic-x` wraps in x and is infinite in y). Beta staggering, the
  halo routines and the impl-lis advection builder follow it. The thermodynamics
  now solves the periodic edge columns instead of copying them. The SSA masks, the
  predictor-corrector mask, the CFL check and the f_grnd_acx/acy edges wrap too.
  FastHydrology is told which directions are periodic, so it no longer overwrites
  their rims. Non-periodic runs are unchanged. The ISMIP-HOM shift error drops from
  0.22 to 1e-6, and ISMIP-HOM uses `slope_bg_x` instead of a tilted geometry that
  jumped at the wrap. MISMIP+ is now mirror-symmetric to 1e-3 m.
- **Mass conservation:** the 10% `mb_resid` overshoot is applied only where ice is
  removed. The margin thickness limit can no longer add ice at fractional cells.
- **Linear solver:** an all-zero RHS returns x = 0. Before, Lis returned NaN, which
  was clamped to +u_max. `limit_vel` now lets NaN through, so it is caught by
  `yelmo_check_kill`.

### Non-default options

- **Calving:**
  - The ISMIP7 retreat now acts along the front normal (−∇lsf), so a stagnant
    front also retreats.
  - It uses thermal forcing relative to the local seawater freezing point (new
    `calc_T_freeze_sw`) and the true water depth.
  - Negative subglacial discharge is clipped to 0 instead of giving NaN.
  - Eigencalving is zero unless both eigenvalues are positive.
- **LSF:**
  - The level set is advected by its own advective-form upwind solver
    (∂φ/∂t + w·∇φ = 0, subcycled for CFL), independent of `ytopo.solver`. Before,
    it reused the flux-form thickness solver, which added a −φ∇·w term.
  - The ice velocity is extended into the ocean from ice faces, taking the nearest
    source on either side, so the result no longer depends on sweep order.
    CalvingMIP stays mirror-symmetric (exp1 V +0.13%).
  - The level set has its own all-dynamic mask and is no longer held at 0 where ice
    is not allowed. Ocean cells are pinned to lsf = +1.
  - The `LSFsnap` time counter no longer overflows.
- **SSA:** fixes for `visc_method = 0` and `2`, and for the energy solver with the
  "mask" boundary condition.
- **Grounding line:**
  - Fixes for `taud_gl_method = 2, 3` and `beta_gl_stag = 4` remove a N–S
    asymmetry in TROUGH.
  - `beta_gl_stag = 4` now integrates velocity rather than flux (Gladstone et al.
    2010). The MISMIP3D RF hysteresis gap goes from 464 to 365 km.
- **K24 hydrology:** the latent heat is taken from Yelmo's `L_ice`.
- **Discharge:** `dmb_method = 1` no longer scales `dist_grline` by dx a second time.
- **`ytopo.margin2nd`:** the one-sided margin gradient was twice too large and
  failed the EISMINT symmetry check. It now passes (Linf/Hmax 2e-6).
- **OpenMP:** fixed races on `cb_ref_now`, `is_margin` and `bmb_int`.

### Diagnostics

- `qq_gl_acx/acy` now hold the ice flux across the grounding line [m3/a]; they
  were allocated but never set. `qq_acx/acy` use the upwind thickness, as the
  advection does (was the mean of the two cells).
- Regional budgets integrate fluxes over the region, weight by the true projected
  cell area and `f_ice`, and now close: TROUGH `cmb` was always 0. Regional values
  on polar stereographic grids change by 3–5%.
- `mb_err` is now a true mass-balance residual (it was `dHidt_dyn`), and
  `tot_dHidt_dyn` is added to the budget.
- NH-* predefined grids: y0 = −5400 km (was +5400 km).
- Domains with no regions file default to `regions = 0` (Greenland keeps 1.3).
- New `write_metrics` option writing `yelmo_metrics.nc`.
- The restart read no longer writes a debug `yelmo_restart_init.nc`.
- The SSA "beta appears to be zero" warning only counts rows that are actually
  solved. MISMIP3D printed it 176 times from boundary rows.
- The mass-budget check (`log_mb_check`) reports the absolute residual and the
  residual relative to the gross throughput. The old percentage error blew up at
  equilibrium.
- Domains without a dedicated case mark borders as fixed only in non-periodic
  directions. Eurasia stops with an error if the regions file has no Eurasia codes.

### Other

- C API: read-only getters for `dta%pd%uxy_s` and `dta%pd%H_grnd`, and setters for
  `hyd_N`/`hyd_W_til`. `yhyd.is_external` lets a host model own N_eff.
- `hyd%now%q` is written to restarts. MISMIP3D is handled in
  `ybound_define_mask_ice`.
- Removed dead code:
  - `calc_adv3D_timestep*`, `dt_adv3D`, `index_north`/`south` and
    `calc_diff2D_timestep`.
  - Unused staggering, boundary, regularisation and extrapolation helpers.
  - The unused SIA basal-velocity routines, `ydyn_set_borders`,
    `update_ssa_mask_convergence` and `grounding_line_flux.f90`.
- Test drivers: output timing uses 64-bit integers.
- New docs page on numerical precision (why the symmetry check needs double
  precision for DIVA, and a plan for double-precision internals).

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
