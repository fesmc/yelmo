# Yelmo changelog

## Unreleased (dev)

Mostly fixes from a whole-source audit (2026-09). Most default runs change a
little. MISMIP3D and DIVA runs change more.

### Changes that affect existing par files

- **Optimization methods as strings.** `opt.opt_cf` and `opt.opt_tf` are no
  longer logical: `opt_cf = "none" | "L21"` (`optimize_cb_ref`), `opt_tf =
  "none" | "L21" | "L21-points"` (`optimize_tf_corr_basin`, one correction per
  basin, or `optimize_tf_corr`, point by point). Other values stop the model in
  `optimize_par_load`. initmip: `True` -> `"L21"`, `False` -> `"none"`; `opt_tf`
  must stay `"none"` (no ocean model).

- **Relaxation separate from the optimization.** The topography relaxation of a
  spin-up moves out of `&opt` into its own group `&relax` (`relax_params`,
  `relax_par_load`, `relax_update` in `libs/ice_optimization.f90`): `opt.rel_tau1`,
  `rel_tau2`, `rel_time1`, `rel_time2`, `rel_m` -> `relax.tau1`, `tau2`, `time1`,
  `time2`, `m`, plus `relax.topo_rel` (the `ytopo.topo_rel` mode while active,
  was 4 fixed). After `time2` the par-file `ytopo.topo_rel` and `topo_rel_tau`
  apply again (was `topo_rel = 0`). initmip: `equil_method = "relax"` uses
  `&relax` too (was `topo_rel = 2`, 50 yr until `ctrl.time_equil`);
  `ctrl.time_equil` is removed.

- **Mass budget at `mask_ice` cells.** The implicit advection (`impl-lis`) imposed
  H = 0 in `MASK_ICE_NONE` cells and H = H in `MASK_ICE_FIXED` cells inside the solve,
  so ice that flowed into masked cells disappeared, and ice leaving fixed cells appeared,
  without any budget term (ISLAND4: ~20 % of the initial volume over 5 kyr). Masked cells
  are now advected like other cells (on a non-periodic domain border with no flux through
  the domain edge), and `calc_G_boundaries` alone enforces the mask and books the change
  in `mb_resid`. The global region (`yelmo_ts.nc`) covers the whole domain, and the time
  series has `mb_resid_tot`, `mb_clip_tot` and `mb_relax_tot`, so the budget closes from
  the time series. Results change only where ice reaches masked cells (MASK_ICE: < 0.03 m;
  ISLAND4: 4e-7 in volume); the other benchmarks are bit-identical.

- **New `yelmo.pc_rho_max`** (default 2): the pc adaptive timestep grows by at most this
  factor per step (`set_adaptive_timestep_pc`). The controllers had no upper bound: after
  a nearly static step, `pc_eta` reached its 1e-8 floor and PI42 asked for ~4000 times the
  last step. In TROUGH-F17 (4 km, `dt_min` = 1e-3) dt jumped from 1e-3 to 2.5 a, the
  AB-SAM predictor extrapolated with weights ~±1250 and the run blew up. Only steps where
  the controller asks for more than `pc_rho_max` change (TROUGH-F17 8 km: 5 steps in the
  first 8 a; GRL-16: 9 steps in 1 kyr, up to 2.3x).
- **`yelmo.pc_tol` 5 -> 1** in the defaults and in the par files that had 5. A step is
  redone if `pc_eta` > `pc_tol`; with `pc_eps` = 0.01-0.02, a bad step with `pc_eta`
  ~0.4 was accepted before. The instability kill (mean `pc_eta` > 10 `pc_tol`) now acts at 10.
- **`ydyn.slide_T`, `gamma_T`, `lambda_min` replaced by `ydyn.frz_scale`, `frz_efold`,
  `frz_min`** (`calc_f_slide`, `calc_beta`). β was divided by `f_slide`, which scales
  the sliding speed by `f_slide**(1/q)`: with q = 1/3 (TROUGH) `gamma_T` = 1 K was a
  0.33 K e-fold and `lambda_min` = 1e-6 a speed floor of 1e-18. β is now multiplied by
  `f_slide**(-q)` (q = `beta_q`, 1 for `beta_method` 0, 1), with
  `f_slide = frz_min + (1-frz_min)*exp(T_prime_b/frz_efold)`, so `frz_efold` is the
  e-folding temperature of the sliding speed for any friction law. `f_slide = 1` also at
  partially floating cells (`f_grnd < 1`),
  which inherited the sub-shelf base temperature (T'_b ~ -1.9 K) and froze, and where the
  bed is wet (`hyd_W > 0` or `hyd_W_til > 0`). Defaults and initmip: `frz_scale = True`,
  `frz_efold = 3`, `frz_min = 1e-3` (was effectively 1 K and 1e-6 for the linear law); a
  1 K speed e-fold was too sharp at 4 km in TROUGH-F17 (purges at the velocity limit,
  irregular cycles, 1.4-1.8x cost). Off in the benchmarks. Replace the old keys in
  external par files (`nml_validate` stops).
- **New `ytherm.gl_temperate`** (default True; False in COLUMN-SLAB and FRONT-SLAB). Holds the base of fully grounded
  cells next to floating ice or open ocean at the pressure melting point (ocean-wetted
  bed), in `calc_enth_column` and `calc_temp_column`; freeze-on there is not limited by
  the capacity rule. Otherwise newly grounded cells keep the sub-shelf base temperature
  (T'_b ~ -1.9 K) and count as frozen.
- **Mirror-symmetric one-sided strain rates at ice fronts** (`calc_jacobian_vel_3D_uxyterms`,
  `jvel%dxx`/`dyy`). At a front with ice on the low-index side, the second-order
  one-sided stencil on the faces i, i-1, i-2 tested `f_ice` of cell i-2 instead of
  cell i-1 (the cell between faces i-2 and i-1), unlike the mirror case and
  `calc_strain_rate_horizontal_2D`. The strain rates, viscosity and principal stresses
  near fronts were not mirror symmetric; in the ISLAND4 benchmark this switched the
  calving rate at single front cells (in double precision, symmetry error 3e-6 → <1e-13
  over 100 yr).

- **Ice thinner than 1 mm is ice free for the subgrid front scheme** (`H_ice_eps` in
  `calc_ice_fraction` and `calc_front_cells`, was `H_ice > 0`). A cell holding a round-off
  amount of ice became a front cell with the full reference thickness (`H_eff >=
  front_H_eff_min`), so the front force jumped by one cell on round-off and broke the
  symmetry of the ISLAND4 benchmark (CISM uses `thck > eps11` for the same masks, in double
  precision). CalvingMIP Exp1, MISMIP+ Ice0, MISMIP3D Stnd and A2: unchanged (volume
  differences <= 2e-5).

- **No shear stress at ice-margin corners in the SSA solvers** (`stagger_visc_aa_ab`,
  both assemblers). The corner viscosity was the mean over the ice-covered cells around
  the corner, so corners on a calving front coupled the front faces to the u = 0 faces of
  the ice-free cells: a drag on the velocity along the front. It is now the mean over the
  four cells when all are fully ice covered, and zero otherwise (traction-free margin).
  Found with benchmark A2 (radial floating shelf, docs/dev/benchmark-protocol): rms error
  44 % (energy) / 15 % (residual) → 0.07 %. CalvingMIP Exp1 (25 km, 10 kyr): grounding-line
  radius 573 → 539 km, axis-to-diagonal spread 39 → 28 km (no orientation trend), volume
  −15 %. TROUGH-F17: volume −1.5 %, max speed 971 → 781 m/yr. MISMIP+ Ice0 and MISMIP3D
  Stnd unchanged (straight fronts).

- **`ytherm.use_strain_sia` replaced by `ytherm.strain_heating = "full" | "sia" | "none"`**
  (default `"full"`, same as `use_strain_sia = False`). `"none"` switches strain heating
  off, which the analytic thermodynamics benchmarks need. Par files using
  `use_strain_sia` must be updated (all files in `par/` are).
- **`yelmo.pc_eps` default 1.0 → 0.02** (input/yelmo_defaults.nml, par/yelmo_initmip.nml).
  With the RMS pc norm, pc_eta stays at 1e-3 - 5e-2 in GRL/ANT runs, so pc_eps >= 0.2 never
  limited dt; 0.02 removes the 8-km outlet checkerboard (GRL-8: 0 persistent cells) and lets
  GRL-8/GRL-4 run where pc_eps 1 was killed. Benchmark par files (<= 1e-2) are unchanged.
- **`ymat.de_max` default 2 → 100 a⁻¹, removed from the par files except the trough ones**
  (`input/yelmo_defaults.nml`). The cap on the effective strain rate dates from a less stable
  version. It only enters strain heating and the material viscosity and stresses (`mat%now`),
  not the DIVA/SSA viscosity, but it can reduce strain heating in fast-stream shear margins
  (1–2 a⁻¹). TROUGH-F17 at 4 km with `de_max` = 100 vs 0.5: surge peaks 11.6–22 vs 11–26 km/yr.

- **`yhyd.bkt_floating_mode` default 1 → 0** (input/yelmo_defaults.nml and all par files).
  MARGIN_FILL (1) saturated W_til on grounded cells next to floating ice, so newly grounded
  ice near the grounding line kept N ≈ 0.03–0.5 P0 for centuries. In TROUGH-F17 this
  ungrounded the trough flanks near the front and moved the grounding line ~80 km upstream
  of PISM. ZERO (0) only zeroes W_til on floating cells; the till of newly grounded ice then
  refreezes or drains, as in PISM.
- **TROUGH-F17: `ssa_vel_max = 5e4` m/yr** (par/yelmo_TROUGH-F17.nml, was 1e4), so that
  the surge peak (about 24 000 m/yr at 4 km) is not set by the limit.
- **Smooth velocity limit is the default** (`ydyn.ssa_vel_lim_method = "drag"`,
  `ssa_vel_max = 1e4` m/yr in the defaults and all par files, was a per-component
  clip at 5000 m/yr). A drag τ_c·x², x = (s − 0.8·u_max)/(0.2·u_max), acts on all free
  faces above 0.8·u_max (see below); it is Newton-linearised in the matrix and does not
  enter τ_b or the frictional heating. The speed settles near 0.8–0.85·u_max. New parameter `ssa_vel_lim_tau` (1e5 Pa, drag at u_max).
  `"clip"` is kept. TROUGH-F17 activations (clip vs drag at 5000 m/yr, 0–8 kyr):
  9957 steps / 1161 at dt_min → 4802 / 5, Picard at `ssa_iter_max` in about 50 % of
  activation solves → 0 %. Runs whose speed stays below 4000 m/yr (grounded) and
  5000 m/yr (floating) are unchanged; faster runs change. Details:
  `docs/physics/momentum/solvers.md`, scripts in `analysis/vel-lim/`.
- **The velocity-limit drag acts on all free faces** (grounded, floating and ice-front
  faces, `ssa_mask` 1–4, no f_grnd weight), with half weight at front faces as for
  the friction. It is passed to both assemblers as a separate linearised term (k, r)
  instead of through copies of β and τ_d, which the energy assembler ignored at front
  faces. With `ssa_solver="residual"`, lateral-bc front faces are clipped at
  `ssa_vel_max` instead. Grounded-only drag let front faces run away: ANT-32 killed
  at t = 0.1 yr, GRL-8 (Helheim cliff) at 89 yr, GRL-4 (Jakobshavn) at 41 yr.
- **pc error norm switch** (`pc_norm_L8` in yelmo_timesteps.f90, hard-coded `.FALSE.`):
  the L8 norm of the scaled pc error is kept next to the RMS (default, unchanged results).
  L8 removes the outlet 2Δx checkerboard at pc_eps ~0.03 but costs more steps; see the
  comment there for the 2026-10-02 tests.
- **Velocity-limit diagnostic**: `ssa_lim_n` in timesteps.nc (faces where the limit acts
  after the last Picard iteration: drag onset 0.8·u_max, or clipped), and a log line per
  `yelmo_update` call with the number of affected steps and the maximum face count.
- **`yelmo_check_kill` velocity limit is 2·`ssa_vel_max`** (was a fixed 1e4 m/yr);
  `ssa_vel_max` must be > 0.

- **Capacity basal boundary condition is the default** (`ytherm.basal_bc_method =
  "capacity"`); the till-water rule `"wtil"` still works but is deprecated. A grounded
  base is held at the pressure melting point if the freezing it needs is at most the
  freeze-on capacity C of the bed water; otherwise all of C is frozen and the base
  cools under the flux condition. `ytherm.cap_source = "auto"` (default) takes C from
  the hydrology model with water transport on (`"hyd"`, e.g. K24's routed inflow) and
  from the bucket's `W_til` otherwise (`"till"`). Methods other than `"enth"` keep
  `"wtil"`. GRL-16KM initmip, 1 kyr: about 3% of grounded cells change from held at
  T_pmp to the flux condition with K24; within 0.2% without transport.
- **Hydrology heat in the basal balance under either rule**: Q_diss + Q_sens from the
  hydrology model enter `bmb_grnd` and the frozen-bed flux condition also with
  `"wtil"`. Zero without a model that computes them.
- **K24 source from terms**: `hydro_update` gives K24 the geothermal heat, the heat
  conducted into the ice (`Q_ice_b`), the drained englacial water (new field
  `thrm%now%melt_int`, as `i_eb`) and the sliding fields; `bmb_grnd` drives the bucket
  only. `yhyd.k24_long_coupling_water` is replaced by `k24_coupling_length_kamb86`
  (twice the old value; external par files need the rename, `nml_validate` stops),
  and K24 gains the FastHydrology.jl routing, fill, friction and sliding-law options.
- New output: `bmb_grnd_star`, `bc_b`, `bmb_clamp`, `melt_int` (ytherm), `hyd_C_frz`,
  `hyd_Q_diss`; C-API getters `thrm_bmb_grnd_star`, `thrm_bc_b`, `thrm_bmb_clamp`,
  `thrm_melt_int` and setters `hyd_C_frz`, `hyd_Q_diss`, `hyd_Q_sens` for an external
  hydrology model.
- **Capacity rule: a cold base must be warmed before it is held at T_pmp** (#12,
  `calc_enth_column`). `bmb_star` did not check that the base was at T_pmp, so a cold
  base with a little water jumped to T_pmp in one step for free (the base node has no
  volume). `bmb_star` now includes the heat to warm the base to T_pmp within the step,
  `rho_ice*max(enth_pmp - enth_b, 0)*dz/dt`, over the same `dz` as the conductive flux.
  Zero for a base at T_pmp. A 0.01 K cold-base threshold was tried first; T_pmp drifts
  past it with thinning ice and it broke D4 symmetry on ISLAND4.
- **New `yhyd.k24_ub_hook`** (default True; #13): with a hydrology whose N depends on
  u_b (`hydro_N_responds_to_ub`, K24 transport), DIVA recomputes N
  (`hydro_N_from_ub`) and c_bed after each iteration, and converges only once c_bed
  has settled. N computed once per step from the previous u_b alternated between two
  states (GRL-16 with K24: 3215 cells; 95 with the hook and `k24_sliding_law = 4`).
  Other hydrologies are unchanged. Requires FastHydrology dev (dd21a65).
- **C API `yelmo_set_neff_callback(cb, tag, alias)`** (#14): a host (e.g. Julia K24 via
  YelmoMirror) supplies N from |u_b| inside the DIVA iteration, through the same hook;
  a null pointer unregisters. Rebuild `libyelmo_c_api.so`.

- **`ytherm.qb_method` renumbered**: 1 = faces, 2 = faces to quadrature nodes
  (default), 3 = simple stagger (was 1), 4 = quadrature (was 2). A par file with
  `qb_method = 1` or `2` now selects an energy-consistent method; set 3 or 4 to keep
  the former one (see Answer-changing fixes).
- **Level-set calving is the default** (`ycalv.use_lsf = True`, `calv_flt_method` and
  `calv_grnd_method = "vm-m16"`) in the defaults and in `par/yelmo_initmip.nml`.
  The mass-balance path (`vm-l19` for floating ice, no grounded calving) left thick
  grounded marine cliffs whose front faces ran at `ssa_vel_max` and set the
  timestep. 1 kyr on 16 threads: ANT-8KM 27.5 → 16.7 min, GRL-16KM 1.1 → 0.4 min.
  Benchmarks with a protocol calving law (`kill-pos` in MISMIP3D, MISMIP+ and
  TROUGH) and MASK_ICE keep the mass-balance path; all other shipped par files
  already set `use_lsf = True`.
- **`par/yelmo_Antarctica.nml` removed.** No in-repo driver could run it (its
  `&ctrl` is the yelmox layout). Antarctica runs use `par/yelmo_initmip.nml` with
  `ctrl.set_nm = "set_ant_pd"` or `"set_ant_lgm"`, `yelmo.domain = "Antarctica"`
  and `yelmo.grid_name = "ANT-32KM"` (or 16/8 km).
- **`ycalv.tau_ice_grnd = 1 MPa`** (was 250 kPa, equal to `tau_ice_flt`): the
  ice strength of `vm-m16` at grounded marine fronts, in line with the ~1 MPa yield
  strength of grounded cliffs (Bassis and Walker, 2012) and common `vm-m16` practice.
  `tau_ice_flt` stays 250 kPa. Not calibrated.
- **`ytherm.cp_rock` replaced by `ytherm.rhoc_rock`** (2.0e6 J m-3 K-1), the
  volumetric heat capacity of the bedrock, equal to the former ρ_rock·cp_rock;
  `rho_rock` left `input/yelmo_phys_const.nml`. The `enth_rock` output and restart
  field are removed. Results are bit-identical. Replace `cp_rock` in external par
  files (`nml_validate` stops).
- **`&Earth sec_year` = 31556926 s** (365.2422 d, the CF/UDUNITS year; was
  31536000 s) in `input/yelmo_phys_const.nml`, as in the other groups except
  MISMIP3D. Changes every real-domain run (per-year/per-second conversions, e.g.
  `kt` in J a-1 m-1 K-1, the hydrology). yelmox and climber-x-input track this file.
- **Physical constants come from fesm-utils `phys_constants`.** `yelmo_init` takes an
  optional constants record from a coupled host; standalone runs read
  `input/yelmo_phys_const.nml` as before. The file's groups changed: `rho_a` is now
  `rho_asth`, and `cp_ice`, `cp_w`, `cp_ocn` and `area_seasurf` are added. Copies of
  this file (yelmox, climber-x-input) must follow. The prognostic path is
  bit-identical.
- **`opt.use_yelmo_cf_min` and `opt.opt_cf_min` removed** (`libs/ice_optimization`).
  They were read but never used: the `cf_ref` floor of the optimisation is
  `ytill.cf_min`. Remove them from `&opt` in par files.
- **`ytopo.front_subgrid = "marine"` by default and in all par files** (was
  `"none"` except Antarctica and initmip). Runs without marine ice (EISMINT,
  HALFAR, ISMIP-HOM, SLAB-S06, MASK_ICE) are unchanged; TROUGH, MISMIP3D, MISMIP+
  and CalvingMIP now use subgrid fronts. Set `"none"` for the binary front.
- **`ydyn.scale_T`/`T_frz` replaced by `ydyn.slide_T`, `gamma_T`, `lambda_min`.**
  `scale_T` raised `c_bed` of a frozen bed only up to `cf_ref*N_eff`, so frozen
  beds still slid (0.2–1 m/a in ANT-32km initmip, 10–100 m/a with low `cf_ref`).
  Now β is divided by `f_slide = max(lambda_min, exp(T_prime_b/gamma_T))` after
  the friction law, for any `beta_method` (new output `f_slide`). On in the
  defaults and initmip (`gamma_T=1`, `lambda_min=1e-6`), off in the
  benchmarks. Replace the old keys in external par files (`nml_validate` stops).
  Benchmarks equal the old `scale_T=0`; where `scale_T=1` stiffened a cold bed
  they change (MISMIP3D Stnd: H up to 7.7 m; TROUGH: < 1 cm).
- **`&yhyd marine_rho_sw` and `k24_latent_heat_water` are no longer read under Yelmo.**
  FastHydrology takes `rho_sw` and `L_ice` from the physical-constants record Yelmo
  hands it, so these two keys are inert in a coupled run (they still work for
  FastHydrology's standalone drivers). The keys stay in `yelmo_defaults.nml`, so no
  par file needs editing; a value set there is simply not the source of truth.

- **`yelmo.pc_corr_vel` removed** (it was unused). `nml_validate` stops on unknown
  parameters, so delete it from external par files.
- **`ydyn.ssa_lat_bc = "floating"` now means floating fronts only.** The producer
  and consumer of the ice-front mask disagreed on its codes, so `"floating"` acted as
  `"marine"` (floating and grounded marine fronts). The codes are now shared
  (`MASK_FRNT_*` in `yelmo_defs`). The default and all shipped par files except
  MASK_ICE (`"all"`) use `"marine"`, which keeps previous results. External par
  files that set `"floating"` change behaviour at grounded marine fronts: switch to
  `"marine"` to keep it.
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
  "floating" or "marine"; "none" = the previous `False`; default now "marine",
  see above), with `front_H_eff_min` (50 m) and `front_dHdx` (0). Front cells get an effective
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
  (also those touching the ocean at a corner) hold at most `a_lsf`·`H_ref`
  (CISM subgrid calving mask; `H_ref` from the interior neighbours, one pass, see
  "Level-set subgrid fronts" below). `f_ice = H/H_eff` as in
  the mass-balance path. `ytopo.f_ice_method` is removed (it had no effect with
  `"none"`); `calc_ice_fraction_lsf` becomes `calc_lsf_area_fraction`.
- **`ydyn.ssa_lat_bc = "slab"` and `"slab-ext"` are removed** (and
  `extend_floating_slab`). They were only set in the ISMIP-HOM and SLAB-S06
  pars, whose periodic domains are fully ice-covered, so they had no effect;
  those pars now use "marine".
- **New, removed and changed parameters:**
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
  - Performance (ANT-8KM, 16 threads, 100 yr: 31.5 -> 21.9 min): the Picard
    viscosity relaxation runs as an OpenMP loop (was serial), and the thickness
    advection uses the Jacobi preconditioner instead of ILU, whose row-by-row
    teardown in LIS was serial. Results are identical in single precision.
    Further serial work moved into OpenMP loops: advection solver BiCGSTAB
    (BiCG's transposed product was serial), border types in the SSA assembly
    evaluated once (not per cell), power-law basal friction, thermodynamic
    horizontal advection, `calc_ice_fraction`/`calc_front_cells`, the DIVA
    viscosity copy and the LIS array copies.
  - **SSA energy solver fixed and made the default** (`ydyn.ssa_solver =
    "energy"`, CG tolerance 1e-2 in all pars and defaults). Its calving-front
    rows applied the front force twice: the front-face driving stress
    (`taud*dx*dy`, taken across the front) equals the boundary work
    `taul_int*dy`, and both were on the RHS. Shelves spread too fast
    (MISMIP3D Stnd grounding line 440 km instead of 520 km, TROUGH -26% ice).
    Now only the boundary work is applied; MISMIP3D (520.07 vs 520.08 km) and
    TROUGH (+0.3% ice) match the residual solver. CalvingMIP differs only in
    the ragged front ring (weak vs strong front condition). Tolerance 1e-2 and
    1e-4 give the same results.
  - **SSA energy assembly rewritten element by element** from the energy
    density (4x4 local Hessians of cell and corner terms, mapped to free,
    Dirichlet, tied or ghost unknowns). Free-slip edges that ice reaches are
    folded into their inner unknown (T^T K T) instead of a two-entry
    constraint row, so K is symmetric for every boundary type (before:
    non-symmetric with `infinite`, `mask`, `periodic-x/-y` and the MISMIP3D/
    TROUGH right edge once ice reached it). Rows are assembled in parallel.
    Benchmarks, ANT-16/ANT-8 and GRL-8 are unchanged to round-off (none has
    ice at a free-slip edge). The LIS copy passes only the stored entries.
  - `ytopo.dHdt_dyn_lim` is removed, with the tendency limit in
    `apply_tendency`. It clipped the dynamic thickness change at ±100 m/yr
    cell by cell, which does not conserve mass. In ANT-16KM/GRL-8KM initmip it
    acted every step for 1 kyr, at fast marine fronts and grounding lines, not
    only at the start; without it the runs are as stable, with the same steps
    and SSA iterations (ANT-16 200 yr: +11e3 km3 ice, more calving).
  - `ycalv.tau_ice` is split into `tau_ice_flt` (250 kPa) and `tau_ice_grnd`
    (now 1 MPa, see above): the `vm-m16` ice strength for floating and marine-grounded fronts
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
- **Requires fesm-utils dev ≥ `842b22a`**: `phys_constants` (db17be3),
  `tstep_update` advancing on every call (78a9f91), ncio a7df7c5, and the `act`
  argument of the gaussian-quadrature node routines (`gq*_to_nodes_acx/acy` and
  `gq2D/gq3D_to_nodes_aa`, 842b22a).
- **`ytopo.margin2nd` removed.** Taken at the face, the one-sided quadratic reduces
  to the plain difference, which is what ran (the option was off in all par files).
  Delete the key from external par files (`nml_validate` stops).
- **Requires FastHydrology dev ≥ `905a81d`**. It adds `hydro_calc_N`,
  `hydro_init_state` taking `H_ice`, and optional `periodic_x`/`periodic_y` in
  `hydro_init`.
- **`ytopo.surf_gl_method` and `ydyn.ssa_beta_max` removed.** Both were read but
  never used (the surface is always `calc_z_srf_max`). Remove them from external
  par files (`nml_validate` stops).
- **`yelmo.experiment` and `ytherm.enth_cp_method` are validated.** An unknown
  value now stops the run. Before, a typo silently gave `"zeros"` boundaries
  (`experiment`) or the constant heat capacity (`enth_cp_method`). Valid
  experiments: `None`, `EISMINT`, `MISMIP3D`, `MISMIP+`, `TROUGH-F17`, `SLAB`,
  `ISMIPHOM`, `slab`, `periodic`, `periodic-xy`, `periodic-x`, `infinite`,
  `MASK_ICE`.

### Answer-changing fixes

- **`yelmo_update` steps until `time` is reached** (`yelmo_ice.f90`). The loop ran at most
  `nstep = ceiling((time-time_now)/dt_min)` steps, so steps shorter than `dt_min` (the last
  steps to reach `time`, `dt_min` rounded down in `limit_adaptive_timestep`) could make the
  call return before `time`. The loop now exits when `time` is reached; `dt_save` is replaced
  by running counters, and the `dt_min` kill check is unchanged. EISMINT, TROUGH-F17:
  unchanged.
- **Courant limit of the time step on the transport velocity** (`calc_transport_velocity`,
  `yelmo_update`). The limit used `ux_bar`/`uy_bar` on all faces, including faces that
  the thickness advection closes (partial cell next to an ice-free cell that may not
  fill) and without the `pc_filter_vel` mean. At partial front cells the closed ocean face
  carries the extrapolated front speed (Rink Isbræ, GRL-8KM: 5 km/yr) and set dt although
  no ice crosses it. The predictor, corrector and Courant limit now use the same transport
  velocity. InitMIP 1 kyr, mean dt: GRL-8KM 0.73 -> 1.15 a, GRL-16KM 1.37 -> 3.2 a,
  ANT-16KM 1.2 -> 2.4 a, no redos, volume within 2e-4. EISMINT, HALFAR, MASK_ICE,
  MISMIP3D, ISMIP-HOM, TROUGH-F17 and CalvingMIP Exp1: unchanged (only the `dt_adv`
  diagnostic in the last two). See docs/physics/timestepping.md.
- **`calv_flt_method = "kill"` / `"kill-pos"` act after the front advance**
  (`calc_ytopo_calving`). The kill ran before `calc_G_front_advance`
  (`front_subgrid`), so ice the advance pushed into the kill region survived the
  step and was removed in the next one: the first cell beyond the front filled
  and emptied on alternate steps, and the SSA front moved by one cell each step
  (TROUGH-F17: SSA-active points 12805/12837, first Picard change 1400-2600 m/yr).
  Now no floating ice remains in the kill region after each step.
- **Ice-free land does not count as ice in the level-set area fraction**
  (`calc_lsf_area_fraction`). The level set is held at −1 on land, and these values
  entered `a_lsf` of the neighbouring cells as ice. An ocean cell beyond the front
  (lsf > 0) with land on three sides had `a_lsf` ≈ 0.4, so it stayed open for
  filling and the front cell next to it spread ice into it. A few millimetres
  there made a 50-m (`front_H_eff_min`) column in the SSA that closed the front
  cell's ocean face: GRL-8KM (108,55) moved sideways at ~8 km/yr (velocity-limit
  drag active in 86 % of steps, dt 0.43 yr).
- **Front cells with `f_ice` < `A_FRONT_MIN` (0.1) are ice-free in the momentum
  balance** (`calc_ytopo_diagnostic`: `H_ice_dyn` = `f_ice_dyn` = 0). They keep
  their ice and fill by transport. Every cell with ice was active, and partial cells
  take part with `H_eff` ≥ `front_H_eff_min`, so a few millimetres of ice (e.g. land
  ice spilling into the ocean, slivers next to a front cell) were 50-m columns that
  closed the ocean face of the neighbouring front cell. With `front_subgrid = "none"`
  only cells with H ≤ 1 mm change (`f_ice` = 0 there). Both calving paths. With both
  changes above (1 kyr): GRL-16 velocity-limit drag active in 2 % of steps (was 12 %),
  GRL-8 restart case dt 0.43 → 0.76 yr; GRL-8 from PD dt 0.97 → 0.73 yr, since Rink
  Isbræ no longer advances into a land-walled fjord bend and flows at ~2.8 km/yr.
  Benchmarks: EISMINT/HALFAR round-off (≤ 1 mm cells), MISMIP3D transient ≤ 3 m
  (final GL unchanged). TROUGH-F17: without the kill-after-advance change above,
  one inactive cell beyond the front gained snow (0 → 68 m in 2 kyr); with it, the
  cell stays empty.
- **Robin temperature profile uses `const_kt` and `const_cp`** (`define_temp_robin_3D`,
  methods `"robin"` and `"robin-cold"`). It used `kt` and `cp` from the current `T_ice`,
  which in `yelmo_init_state` is still 0 K: k = 9.83 W m-1 K-1 and c = 146 J kg-1 K-1.
  The basal gradient was then G/9.83 instead of G/k, and the profile far too
  diffusive (TROUGH-F17, H = 500 m, G = 70 mW m-2, a = 0.3 m/yr: T_b = -16.5 °C,
  should be about -10 °C).
- **Gaps in the topography files are filled** (`yelmo_init_topo`, `ydata_load`):
  where a dataset has missing values (e.g. outside its coverage, as in the ISMIP7
  obs files), there is no ice, the bed comes from the nearest valid cell (fesm-utils
  `fill_nearest`), the surface from the bed and the ice thickness (sea level 0), and
  `z_bed_sd` is 0. Before, the bed was -9999 there, and `grad_lim_zb` pulled the
  valid bed down from those cliffs up to 85 cells into the domain (ISMIP7 GRL-8KM:
  21,000 cells, by up to 7 km); the initial `z_srf` and `z_bed_sd` read the raw fill
  values. New `ydata_fill_topo_gaps`. Files without gaps are unchanged.
- **Basal frictional heating: `ytherm.qb_method = 2` by default.** New options 1
  ("faces", PR #11) and 2 ("faces to quadrature nodes") form the friction work
  `|taub_acx*ux_b|` and `|taub_acy*uy_b|` on the C-grid faces, where both factors live.
  1 averages the two faces of each cell to the aa-node; 2 brings both terms to the
  quadrature points and sums them there (a 1-2-1 smoothing of 1 across each face). The
  heat of a face next to a cell that is not fully ice covered (`Q_b = 0`) goes to the
  fully covered cells, so the domain total is the work done by basal friction in the
  momentum balance (GRL-16: 1.0000 with 2, 0.995 with 1). The former options 1 (simple
  stagger) and 2 (quadrature) are now 3 and 4; they multiply the magnitudes of
  separately interpolated vectors and do not (4: −10.5 % GRL-16; ISMIP7 spin-ups
  −1.8 % GrIS 8 km, +6.8 % AIS 16 km). GRL-16, 1 kyr from a spun-up state, friction
  fixed: grounded basal melt +12 % with 2 (+10 % with 1) vs 4; temperate area and
  volume change at the noise level. `ytherm_par_load` stops on a `qb_method` other
  than 1-4 (before, an unknown value left `Q_b` unchanged).
- **Thermodynamics on the dynamic column** (review MAT-3). The enthalpy solve,
  `T_pmp`, `T_shlf`, the SIA strain heating, the Robin/linear profiles and the
  bedrock shelf temperature use `H_ice_dyn` (the column of `uz_star`) instead of
  `H_ice/f_ice` or `H_ice`. `f_ice == 1` still decides where the column is solved.
  Changes front and `H_eff`-floor cells (1511 in ANT-16, 353 in GRL-16): basal T′ at
  ocean fronts +0.01 K (ANT) and +0.05 K (GRL); TROUGH max|ΔH| 0.1 m; CalvingMIP,
  MISMIP3D and EISMINT keep their ice thickness. ANT/GRL-16 300 yr changes are at
  the noise level of a 1e-4 perturbation of `enh_shear`.
- **One front classification for the level-set trim and `f_ice`** (review TPO-4).
  `calc_front_cells` with the level-set area fraction makes eligible cells cut by
  the front (`a_lsf` < 1) that touch the ocean only at a corner front cells, in both
  `calc_G_lsf_front` and `calc_ice_fraction`. These cells are now partial (`f_ice` ≈
  `a_lsf`, `H_ice_dyn = H_ref`) instead of thin full cells. 187 such cells in ANT-16,
  37 in GRL-16, 60 in CalvingMIP exp1 (10 ka). ANT-16 300 yr: floating area −0.3 %,
  calving +6.5 %; GRL-16: floating area +17 %, calving +21 %; CalvingMIP exp1:
  area −0.17 %, equivalent front radius 753.5 → 752.9 km. Level-set runs only.
- **Exact `f_ice == 1` test for full ice cover** (review OMP-6) in the DIVA
  viscosity, `calc_F_integral`, `calc_visc_eff_int`, the SSA masks and assemblers,
  `calc_ice_front` and `gen_mask_bed` (was `is_equal`, tolerance 1e-5). No cell
  with 1−1e-5 < `f_ice` < 1 occurred; all validation runs are bit-identical.
- **Half drag at mask-4 front faces** (review DYN-7). Front faces treated as inner
  SSA (grounded fronts, faces to ice-free land) get ½β like mask-3 faces, in the
  residual and the energy assembler. GRL-16 300 yr (~650 such faces): volume
  −4.5e-4, floating area +5.6 %, grounding-line flux −4.7 %; ANT-16: floating area
  +0.14 %, grounding-line flux +2.2 %; CalvingMIP exp1 at 1 ka (advancing grounded
  margin): volume −0.2 %; no change at 10 ka, EISMINT (SIA) unchanged.

- **Sub-temperate sliding and basal drag at grounding lines and margins**
  (review 2026-10-01). β is divided by `f_slide` on aa-nodes in `calc_beta`,
  before staggering, like any other spatial variation of friction
  (`scale_beta_slide` is removed). Grounding-line faces of frozen grounded cells
  get the grounded β/`f_slide` (before, `f_slide` was averaged with the floating
  side's 1, so at most twice the drag), and interior frozen/temperate faces take the
  mean of β/`f_slide` instead of β over the mean of `f_slide`. The friction laws
  with quadrature (`beta_method` 1, 2, 3) map `c_bed` to the quadrature points with
  corner means over grounded ice only; before, the zero `c_bed` of floating and
  ice-free neighbours reduced β by 12–44 % at grounding-line and margin cells.
  MISMIP3D with `beta_method = 2` now has its grounding line at 520.2 km, as with
  `beta_method = 4` (was one cell short, 499.7 km). ANT-16 200 yr: grounding-line
  flux -13 %, floating area -1.6 %, calving -4.5 %; GRL-16 floating area -4.2 %;
  CalvingMIP exp1 volume +3.6 %. Calibrated `cb_ref` fields should be
  re-optimised.
- **Front geometry and calving stress at subgrid fronts** (review 2026-10-01).
  Full front cells holding more ice than their reference `H_eff` keep their own
  column in the dynamics (`H_ice_dyn = max(H_eff, H_ice)`; 15 such cells in GRL-16,
  43 in ANT-16, up to 550 m of excess). The material viscosity, and so the calving
  stress, is masked with the geometry of the last velocity solution (new restart
  fields `H_ice_solv`, `f_ice_solv`; older restarts fall back to `H_ice_dyn`,
  `f_ice_dyn`) instead of the start-of-step geometry, and cells that received ice
  since that solution take their neighbours' principal stresses before `vm-m16`
  and `vm-l19`. Before, newly filled front cells had zero stress, so no stress
  calving, for one to two steps. ANT-16 200 yr: calving +7%, floating area -0.3%;
  GRL-16: floating area -1.9%; volume change <= 2e-5. Benchmarks without stress
  calving (TROUGH, CalvingMIP exp1, EISMINT) keep their ice thickness.
- **HALFAR uses the adaptive predictor-corrector timestep** (`dt_method = 2` in
  `par/yelmo_HALFAR.nml`). The fixed 1-yr explicit step exceeded the SIA stability
  limit at dx = 2 km (0.71 yr at t = 0): the run was killed at ~15 yr ("velocity too
  fast") once `dHdt_dyn_lim` was removed, and before that it finished with ~100 m rms
  error and +30% volume. Now 12 m rms after 200 yr (dx = 2 km), converging with dx.
- **Vertical velocity at ice fronts** (review 2026-10-01). The sigma-coordinate
  transform (`uz`, `uz_star`, Jacobian corrections) uses the gradients of the ice
  column: at faces between ice and ice-free ocean (bed below sea level) the jump
  to sea level is replaced by the gradient of the adjacent interior face; margins
  to ice-free land keep their slope. Before, half of the front cliff entered the
  basal `uz` and `uz_star` of full front cells: surface `uz_star + smb` was +6.95
  m/yr at TROUGH-8 fronts and +9.3 m/yr at MISMIP3D fronts (now ~0), and front
  columns were too warm (TROUGH -10.6 -> -11.8 °C). The uz quadrature, the cross
  sigma terms of the Jacobian and the velocity-dependent beta (`beta_method` 2, 3)
  only use velocity faces with a solution. Ice thickness: TROUGH <= 0.5 m,
  CalvingMIP 3 m at 10 ka, MISMIP3D and EISMINT unchanged.
- **Calving-front stress and strain at subgrid fronts** (review 2026-09-29).
  With `front_subgrid /= "none"`, `calc_ymat` built the viscosity and the calving
  stress on `f_ice` instead of the dynamics' `f_ice_dyn`/`H_ice_dyn`: partial front
  cells took a neighbour's viscosity, or none (zero stress, no calving: 22 such
  cells in GRL-16, 60 in GRL-8). Now on the active geometry; the partial-cell fills
  are removed. Calving +10%, volume change <= 0.02% (200 yr).
- **Front strain rates no longer use the zero velocity of ocean faces.** Cross
  derivatives (`dxy`/`dyx`) and the quadrature corner means of the strain-rate
  tensor and the DIVA/SSA viscosity only use velocity faces next to ice. A translating
  slab now has zero strain at straight and 45-degree fronts (was up to 3V/(16dx),
  `make front_strain`). Benchmarks change by <= 0.6%.
- **Vertical velocity from the kinematic rates of the ice column; no uz clamps.**
  The basal rate was `dzsdt - dHidt`, which at subgrid fronts mixed the `H_eff`
  surface with the true thickness and counted calving, front advance and removals
  as vertical motion: basal `uz` of 3.5-5.5 km/yr in the first step of GRL-16/ANT-16
  and 0.3-1.6 km/yr spikes later, hidden by the ±10 m/yr clamps on `uz`/`uz_star`.
  Those clamps also cut the physical w in fast outlets (8% of the ice flux in
  GRL-16), which entered the enthalpy advection as a spurious sigma velocity.
  Now `dzbdt_kin`/`dzsdt_kin` come from the vertical column change (advection, smb,
  bmb, relaxation) and the bedrock/sea-level rates between `yelmo_update` calls
  (floating columns float), zero in partial front cells; the clamps are removed.
- **Energy SSA: half the basal drag at calving-front faces**, the ice half of the
  face's control area. Grounded marine fronts are 2.2-2.5x faster than before in
  GRL-16/ANT-16/GRL-8; floating fronts are unchanged.
- **Level-set subgrid fronts** (`use_lsf` with `front_subgrid /= "none"`).
  The trim of front cells to `a_lsf*H_eff` used an `H_eff` that depends on the
  cell's own thickness (surface limit, no-neighbour fallback), so each step lowered
  the next target: a 1000 m marine cliff at `a_lsf=0.7` ended at 247 m and floating,
  1-2-cell-wide tongues thinned away under a stationary front. The trim now uses a
  reference from the interior neighbours (`calc_front_H_ref`), and cells without one
  are not trimmed. `a_lsf` includes the cell centre (continuous for 1-cell-wide
  tongues, which were emptied at once). Eligible cells are emptied by area
  (`a_lsf < 0.1`) instead of by the centre `lsf > 0`, ice flows from partial cells
  into ice-free cells the level set covers, and only ice-free cells without an
  ice-covered edge neighbour are reset to ocean in the level set (`make lsf_front`).
  TROUGH (4 km, zero calving): front speed 0.61-0.73 u -> 0.88-0.94 u, ice removed
  35-38% -> 7-8% of the inflow. ANT-16 (200 yr): trim-caused grounded-to-floating
  switches 328 -> 17, volume +0.05%. Runs without the level set are unchanged.
- **Level set crosses the coast.** Grounded faces with a mean bed above mean sea
  level had `cr = -u` (no level-set motion), which also closed many faces between a
  land and a marine cell, so the land value could not reach the first marine cell.
  With a subgrid front that cell filled with ice while `lsf > 0`, and the value
  spread along grid-aligned flow lines (CalvingMIP: near-zero `lsf` under the ice
  on the axes, a one-cell ice-free channel on each axis in exp2). Now only faces
  between two land cells are closed.
  CalvingMIP (marine): no `lsf > 0` under interior ice (was 84 cells), exp2 axis
  front 742 km at 1 kyr (was 601 km), exp1 volume +0.5%, ice crosses the coast on
  the axis at 1400 yr (was 4200). ANT-16 (200 yr): volume +0.001%, calving -5%;
  re-advance after a retreat to the coast regains 23% more area. Runs without the
  level set are unchanged.
- **Mass balance on the ice fraction after transport.** smb, bmb, fmb and dmb
  used `f_ice` from before the advection step, so a cell that had just received ice
  (f_ice = 0) took the full per-area melt: with the level set, ~100 ice-free cells
  behind the ANT-16 front melted all inflow every step (bmb_shlf = -2 m/yr on the
  whole cell, ~23 km3/yr at t=200) and never filled.

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
- **Partially ice-covered cells** carry their own 2D strain rates and `visc_bar`
  (they were zero), as part of the active ice geometry (see "Calving-front stress
  and strain at subgrid fronts" above; a first version copied them from fully
  ice-covered neighbours). Before, `calc_eps_eff` and `calc_tau_eff` patched this
  separately, and the `vm-m16` calving law saw zero stress, so partial front cells
  never calved.
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
- **SSA masks: floating ice ends at a wall across ice-free land.** A face between
  floating ice and ice-free land was solved as an inner SSA face, driven by the slope
  from the bare land surface with no friction, and ran at `ssa_vel_max` (the special
  case meant to set it to zero was never reached). These faces set the Courant
  timestep. They now have zero velocity. ANT-16KM / GRL-16KM, 200 yr: faces at the
  velocity limit 7 → 0 / 49 → 0, timesteps 198 → 52 / 261 → 73, volume -0.02% /
  +0.04%, calving -9% / -11%. TROUGH, CalvingMIP and EISMINT are bit-identical.
- **An imposed beta is used as given.** With `beta_method = -1` (beta on aa-nodes) or
  `beta_gl_stag = -1` (beta on ac-nodes) the friction coefficient is no longer
  modified: no grounding-line scaling, no zero under floating ice, no `beta_min`
  limit and no `f_slide` scaling. Before, an imposed beta was divided by `f_slide`
  again in every Picard iteration, so it kept growing over a frozen bed. With
  `beta_method = -1` the face values are still staggered from the imposed beta as
  set by `beta_gl_stag`. ISMIP-HOM C and F (`f_slide = 1`) are bit-identical.
- **Level set (`lsf_method = "snap"`): two free cells on each side of the front.**
  The snap left one free cell, so the cell ahead of the front was held at +1 until
  the front cell changed sign, and the front moved at only 0.87 of the advection
  speed for small Courant numbers (0.92 at C = 0.13), depending on the timestep. With
  two free cells it moves at 0.96-0.98. Every level-set run with a moving front
  changes (fronts up to 12% faster, in advance and retreat). TROUGH with zero
  calving: ice removed by the front trim 7.3% → 3.8% of the inflow. CalvingMIP exp2:
  maximum retreat 0.92 → 0.98 of the prescribed one; exp1: steady front 5 km further
  out, volume +1.3%. ANT-16KM, 200 yr: floating area +1%, volume within 0.02%. Runs
  without the level set are bit-identical.
- **`yelmo_set_time` also sets the tracer and hydrology clocks and `bnd%time_n`.**
  After `yelmo_update_equil` the first step skipped tracer advection and hydrology
  (dt = 0) and had zero bedrock and sea-level rates.
- **`calc_strain_rate_horizontal_2D` (DIVA/SSA viscosity):** the one-sided `dvdy`
  had no first-order fallback on the last row (it kept the centred difference across
  the ice-free face), and the ±2 neighbours were not BC-aware, so periodic seams used
  first order. Changes results only at periodic seams (TROUGH, MISMIP3D) and on the
  last rows.

- **Level set in periodic domains:** the advection and redistance of the level set
  wrap in periodic directions, as the snap and the area fraction do; "infinite"
  replaces only zero borders. Changes only runs with ice at a periodic seam (none of
  the shipped benchmarks).
- **Iceberg rule after the front advance** (mass-balance calving, subgrid fronts): a
  front cell left by the advance at `H_eff` − 0.1 m counts as full, so the cell it
  fed at a convex corner is no longer removed as an iceberg. TROUGH, MISMIP3D and
  MASK_ICE are bit-identical.

### Non-default options

- New `ytrc.elsa_restart` (default `True`): `False` starts elsa's layers fresh also
  on a restart, e.g. for a transient run on another clock than the spin-up it
  restarts from. The other backends restart as before.
- **Imposed `beta_acx/acy` (`beta_gl_stag = -1`)** is no longer raised to `beta_min`
  or overwritten at the domain borders.
- **New `yelmo.experiment = "periodic-y"`**: periodic in y (true wrap, period ny) and
  infinite in x, the transpose of `periodic-x` in all solvers. A periodic-y strip
  reproduces the transposed periodic-x strip to round-off (B2 flowline, 32 km, 500 yr).

- **Frontal melt `ytopo.fmb_method = 3`** (Rignot et al., 2016, ISMIP7 protocol)
  with the new boundary field `bnd%tf_shlf` (thermal forcing) and the subglacial
  discharge `bnd%Qd`, scaled by the new `ytopo.fmb_lambda` (default 1). The rate
  (m/yr, negative discharge clipped to 0) comes from `calc_melt_rate_rignot16`,
  shared with the ISMIP7 retreat (which changes only at round-off: its constants
  were single precision).
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
- **OpenMP:** fixed races on `cb_ref_now`, `is_margin` and `bmb_int`.

### Diagnostics

- Level-set calving diagnostics (`calv_rate_flt/grnd`, `cmb_flt_x/y`, `cmb_grnd_x/y`,
  `cr_acx/acy`) are saved with the predictor/corrector fields, so the output shows the
  same stage as `lsf` and `cmb` (they mixed predictor and corrector). `calv_rate_*`
  are built from `cr_acx/cr_acy`, the face rates the level set uses (law chosen by
  `f_grnd_acx/acy`). Model state is unchanged.
- New `mb_clip` (ytopo): the clip of negative thickness after transport. It was
  booked in `dHidt_dyn`, which is now pure transport; `mb_err` and the
  `log_mb_check` line include it. ANT-16/GRL-16: 0.16/0.13 km3/yr, 0.2%/2% of calving.
- New `uz_srf_err` (ydyn): `uz_star` at the surface + smb on fully ice-covered cells
  (0 if the kinematic rates and uz agree). ANT-16: median |.| 0.03 (grounded) and
  0.08 m/yr (floating), 95th percentile 0.7 and 1.1 m/yr; GRL-16: 0.08 and 0.9,
  1.6 and 4.6 m/yr.
- `taub` uses the friction of the SSA matrix: `beta_min` at grounded faces with zero
  β (`beta_eff` for DIVA) is set once before the solve, not inside the assemblers,
  where `taub` stayed 0. Results are bit-identical.
- With `ytill.method = -1` (external `cb_ref`), `calc_ydyn` stops if `cb_ref < 0`
  under grounded ice, or if it is 0 everywhere under grounded ice (not set by the
  driver). Local zeros are allowed (SLAB-S06 free-slip centreline). `yelmo_init_state`
  no longer sets an all-zero `cb_ref` to 1; its high-beta fallback applies only to
  `beta_method = -1`.

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

- `yelmo-config` checks the remaining numeric par_load checks: `yelmo.dt_min`,
  `pc_rho_max`, `ydyn.ssa_vel_max`, `neff_nxi`, `ytherm.qb_method`, `advecxy_order`,
  `advecxy_cfl`, `advecxy_nmax`, `cap_W_floor`, `cap_eps`.
- `ytherm.basal_bc_method` and `cap_source` are checked with `yelmo_check_enum` (so
  `yelmo-config` knows them); `yelmo.dt_min` must be > 0 (was: not 0).
- Docs: C API reference (docs/c-api.md); pages synced with dev.
- New module `yelmo_remapping` (`src/yelmo_remapping.f90`, re-exported by `use yelmo`):
  `yelmo_remap` (2D: horizontal with a coords map; 3D: vertical interpolation onto a
  Yelmo axis, then horizontal), `yelmo_load_map` (map from a file's axes; replaces
  `yelmo_restart_load_map` in the restart reader) and `yelmo_read_remap` (read and map a
  2D or 3D field from NetCDF). See docs/remapping.md.
- `yelmo_init_state(..., thrm_method="prescribed", T_ice=T_ice)`: initial ice temperature
  from an external field on the Yelmo grid, capped at the pressure melting point, with
  the consistent enthalpy.
- `yelmo_opt.x` (`tests/yelmo_opt.f90`, `make opt`, runme alias `opt`) removed:
  basal friction optimization is the initmip spin-up option
  `ctrl.equil_method = "opt"`.
- `restart_interpolated` compares the restart grid spacing in m (restart `xc` is in
  km); before, every interpolated restart counted as coarser than the model grid.
  With `yelmo.restart_z_bed = True`, a restart from a finer grid now uses the
  interpolated restart bedrock.
- `make clean` also cleans elsa and tracer (like FastHydrology), so switching
  between `openmp=0` and `openmp=1` no longer links stale sub-library objects.
- New public `yelmo_restart_init(dom, filename, time)`: the restart branch of
  `yelmo_init_state` (restart fields, active predictor-corrector, passive-tracer
  backends from their sidecars). Coupled drivers that restart from their own
  bundles call it, so elsa and tracer are initialized after a restart too.
- All fatal error paths in `src/` and the `tests/` drivers end with `error stop 1`
  (was `stop`, exit status 0), so a failed run is reported as FAILED by SLURM,
  including the restart writers when an io table lists an unknown variable.
- The SSA/DIVA Picard loop stops (`error stop 1`) when the velocity on a solved
  face or the residual is not finite, reporting the iteration and the first bad
  face and writing `ssa_check.nc`. It replaces the unreachable "strange case" test
  (change > 9999 m/yr).
- The `RALSTON` predictor-corrector branches are removed (`yelmo.pc_method`
  already accepted only `FE-SBE`, `AB-SAM` and `HEUN`).
- ISMIP-HOM Experiment F (`yelmo_ismiphom`, `ctrl.experiment = "EXPF1"` no slip,
  `"EXPF2"` slip ratio 1; Pattyn et al., 2008): a 1000 m slab on a 3° slope over a
  Gaussian bed bump relaxes to steady state with zero SMB, n = 1 and
  A = 2.140373e-7 Pa⁻¹ a⁻¹. The driver sets the 100 km domain, n, A, β and lifts
  `ydyn.taud_lim` (the driving stress is 4.7e5 Pa). Steady state after ~1000 yr
  (`time_end = 2000`); along the central flowline Yelmo (DIVA or hybrid) is within
  0.4-0.6 m (surface) and 0.2-0.5 m/a (speed) of the full-Stokes results. Before,
  `"EXPF"` was a placeholder that ran without ice.
- runme OpenMP jobs run one thread per physical core (`--hint=nomultithread`,
  `OMP_PLACES=cores`) with `KMP_BLOCKTIME=0`. On Levante `shared`, 16 threads were
  placed on 8 cores with 2 SMT threads each; ANT-8KM 200 yr: 20.5 -> 15.8 min,
  results identical. Jobs are charged for the full cores.
- C API: read-only getters for `dta%pd%uxy_s` and `dta%pd%H_grnd`, and setters for
  `hyd_N`/`hyd_W_til`. `yhyd.bkt_N_closure = -1` lets a host model own N_eff
  (set with `yelmo_set_var2D("hyd_N")`).
- C API: getters for `tpo_dzsdt_kin`, `tpo_dzbdt_kin`, `tpo_dHidt_vert`,
  `tpo_calv_rate_flt` and `tpo_calv_rate_grnd`. The `dyn_f_slide` setter is removed
  (`f_slide` is recomputed in every `calc_ydyn`).
- The restart carries the previous-call bedrock and sea level of
  `ybound_update_rates` (`z_bed_n`, `z_sl_n`, `bnd_time_n`, `bnd_rates_init`), so a
  continued run has the straight run's `dz_bed_dt`/`dz_sl_dt` on its first step (they
  were zero). It is restored only when the start time equals the restart time; a
  restart used as a state at another time, and old restarts, give zero rates as before.
- The restart read of the 3D `enh_bnd` uses 3D start/count.
- `hyd%now%q` is written to restarts. MISMIP3D is handled in
  `ybound_define_mask_ice`.
- Removed dead code:
  - `calc_adv3D_timestep*`, `dt_adv3D`, `index_north`/`south` and
    `calc_diff2D_timestep`.
  - Unused staggering, boundary, regularisation and extrapolation helpers.
  - The unused SIA basal-velocity routines, `ydyn_set_borders`,
    `update_ssa_mask_convergence` and `grounding_line_flux.f90`.
- Test drivers: output timing uses 64-bit integers.
- Kill check: NaN and Inf are found with a bit test on the exponent (`is_finite`
  in `yelmo_tools`), which is kept under `-Ofast` and catches single cells; the
  `maxval(abs(x-x))` form missed isolated NaNs and never caught Inf. Non-finite
  forcing (`bnd` fields) stops the run at the start of `yelmo_update`. A kill ends
  with `error stop` (exit status 1, SLURM reports FAILED).
- The 12 heaviest ice-only OpenMP loops (enthalpy columns, DIVA/SSA viscosity, uz,
  strain rates, beta) use `schedule(dynamic,64)`: the static split left some
  threads with up to twice the mean number of ice cells. Time in `yelmo_update`
  -8.6% / -6.5% (ANT-8KM, 16 / 32 threads) and -8.4% / -4.9% (GRL-8KM); results
  bit-identical.
- A redone timestep restores only topography and dynamics (`tpo`, `dyn`), the
  components the predictor-corrector modifies, instead of a copy of the whole model
  every step; no copy with `pc_n_redo = 1`. `yelmo_update_equil` saves only the
  parameters. The `update_others_pc` option (compile-time, off) is removed.
- Linear solver status is no longer ignored: LIS returned success at its iteration
  limit or on breakdown. The `ssa:` log line now ends with the linear iterations and
  status (`| lin  187 MAXITER`); a `yelmo_update` call that had such solves (SSA or
  impl-lis thickness advection) prints one summary line; the timestep log
  (`log_timestep`) has `ssa_lin_iter`, `ssa_lin_fail`, `adv_lin_iter`, `adv_lin_fail`.
- Picard loop of the SSA/DIVA solver: the effective viscosity uses 2D strain arrays
  (unused 3D quadrature branch removed), and the serial passes (beta staggering,
  `beta_eff`, basal velocity, relaxation, convergence norms, LIS vector transfer) run
  in parallel. Main loop -9% (ANT-8KM) and -18% (GRL-8KM) at 16 threads; results
  bit-identical.
  The convergence norm is summed per row and the rows serially, so the Picard
  residual no longer depends on the number of threads (it did by 1e-15..1e-13
  relative, without changing a convergence decision in the tested runs).
- elsa coupling for elsa v3.0.0 (fesmc/yelmo#10): Yelmo passes its native arrays to
  `elsa_init`/`elsa_update` (no double-precision copies of the 3D velocities every
  step), and stops at start-up if `use_elsa = True` and `ytrc.time_end` is not later
  than the start time. `ytrc.time_end` defaults to 0.0 (present day). elsa v3.0.0
  applies the time-mean forcing over each coupling period; its restart sidecars
  (`*_elsa.nc`) written by earlier versions no longer load.
- **Calving diagnostics `cmb_flt` and `cmb_grnd`** are now the applied calving (all
  removal at the front) in floating (`f_grnd = 0`) and grounded cells, on both calving
  paths, so `cmb_flt + cmb_grnd = cmb` and their regional sums in `yelmo_ts.nc` are
  mass fluxes. With the level set they held the magnitude of the front calving speed
  at every cell, and the regional "mass balance" was that speed times the cell area.
  The front speeds are new fields `calv_rate_flt` and `calv_rate_grnd` [m/yr], set at
  front cells (ice with an ice-free ocean edge neighbour) on the level-set path.
- **Driver time loops** use fesm-utils `tstep_update`, which now advances on every
  call. The first loop pass was a zero-length
  step at `time_init`: it re-ran the model update at dt = 0, wrote `time_init` twice
  in some outputs, and in `yelmo_initmip` applied one extra `optimize_cb_ref` update
  at every start and restart. `yelmo_initmip` now writes output and restarts at the
  top of the loop (`time_init` on the first pass, `time_end` on the last), reads the
  timeline with `tstep_init(ts, path_par, "ctrl", dtt)`, writes restarts on the new
  `trst` output schedule (always at `time_end`), and has `ctrl.restart_mode`:
  `"state"` (start at `time_init` from a restart state, as before) or `"continue"`
  (stop unless `time_init` equals the restart file's time). Cold starts change
  slightly; a continuation now reproduces the straight run. New
  `yelmo_regions_write_init` creates the regional output files without writing a
  record, so `yelmo_initmip` no longer writes its `time_init` regional record twice.
- `speed_tpo` in the timestep log (`log_timestep`) is now the speed of the topography
  step (predictor + corrector + advance); it was always zero.
- `yelmo_trough` and `yelmo_mismip` keep the restart ice thickness
  (`restart_H_ice=True`); the analytic initial thickness overwrote it.
- **Restart interpolation without cdo:** `restart_interp_gen` (compile-time, `yelmo_io.f90`)
  is now `"coords"`: a restart on another grid builds the conservative weights in-package
  (cached in `maps/`), so no cdo and no pre-generated SCRIP map are needed.
- `make regridding` writes `yelmo_test_regridding.x` (it overwrote `yelmo_ismiphom.x`);
  `make usage` lists all targets; runme alias `mask_ice`.
- The docs variable tables are copied from `input/yelmo-variables-*.md` when the site is
  rendered (Quarto pre-render), so they cannot drift from the code; new ytrc page.
- Output metadata: `pc_eta`/`eta_avg` are the pc error norm in 1/yr (was labelled m/yr,
  "maximum"); units of `enth` (J kg-1), `Q_strn`, `kt`, `Q_rock` (mW m-2) and `advecxy`
  corrected in the ytherm table.
- New docs page on numerical precision (why the symmetry check needs double
  precision for DIVA, and a plan for double-precision internals).
- New `yelmo.mask_border` sets the ice mask on the domain border: `"auto"` (default,
  by domain as before), `"none"` (no ice), `"fixed"` (thickness prescribed) or
  `"dynamic"` (left as the domain mask defines it). `"none"` and `"fixed"` skip
  periodic directions. `ybound_define_mask_ice` builds the mask in two parts (where
  ice is allowed in the domain, then the border) and takes `mask_border` as a new
  argument.
- `yelmo_init_grid(grd, grid0)` builds the Yelmo grid from a coords `grid_class`
  defined elsewhere (e.g. a coupler's grid description), with the axes converted to
  meters. Use it before `yelmo_init(..., grid_def="none")`.
- `yelmo_init` takes boundary fields from a coupled driver that owns the domain
  definition: optional `regions`, `basins`, `mask_ice` (where ice is allowed;
  `mask_border` is applied on top) and `topo_pd`, `topo_init` (new
  `ytopo_input_class`: `H_ice`, `z_bed`, optional `z_bed_sd`, `z_srf`). Each one
  replaces the matching file read; the processing that follows (`z_bed_f_sd`,
  englacial lakes, smoothing, `grad_lim_zb`, `init_topo_state`, references) is the
  same. Without them `yelmo_init` is unchanged. A file path (e.g.
  `yelmo_data.pd_topo_path`) is only checked when its file is read.

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
A new compile-time `restart_interp_gen` switch (module parameter in `yelmo_io.f90`) selects
how the conservative restart map is built: `"cdo"` (load a pre-generated SCRIP map; default
in this release, unchanged behaviour) or `"coords"`
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
