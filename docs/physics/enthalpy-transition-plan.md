# Enthalpy transition plan (ytherm)

Design document for switching Yelmo's production thermodynamics from the
temperature-based solver (`method="temp"`) to the enthalpy-based solver
(`method="enth"`) for v2.0.

Status: **default flipped to `enth`** on branch `therm-dev` (`input/yelmo_defaults.nml`).
Phases 0–2 complete (all 1D benchmarks pass); Phase 3 (2D/3D): the 2D margin-column
robustness crash is fixed (CTS diffusivity) and a geothermal-flux units bug that
made the enthalpy base spuriously cold is fixed; enth now matches temp on
EISMINT-2 A and F, and horizontal advection was upgraded to a flux-limited
2nd-order scheme (`advecxy_order`, default 2) that reduces the EISMINT spoke
asymmetry ~23% while staying TVD/margin-safe (see log). **Post-flip validation
still pending: a Greenland (initmip, 16 km) and an Antarctica spin-up**, enth vs
temp, on real geometry.
Author: thermodynamics review, 2026-07.

## Progress log

- **P4 (early cleanup, done):** removed dead `ice_enthalpy_poly.f90` + test, and
  unused helpers `calc_Q_bedrock`, `calc_hires_cell`, `interp_bilin_pt`.
- **P0 (done):** standalone driver `tests/test_enthalpy.f90` (`make enthalpy`)
  with NetCDF output; Kleiner (2015) data vendored to `tests/data/Kleiner2015/`.
- **P1 (done):** finished `calc_enth_column` — removed the upstream `stop`, added
  flux-based `calc_bmb_grounded_enth`, made `Q_ice_b` swappable, and **fixed a
  sign bug**: the enth path computed `Q_ice_b = +kt·(T₂−T₁)/dz`, the opposite of
  `calc_temp_column`, so basal melt was wrong. Cold-limit equivalence (T2) now
  passes (max|ΔT| = 9.6×10⁻⁴ K) and the Kleiner-A transient shows a physical
  melt/refreeze cycle. Magnitude vs the Kleiner reference is Phase 2.
- **P2 Exp A (done):** validated against Kleiner (2015) Experiment A. Corrected
  the protocol (warm phase −10 °C, β = 7.9×10⁻⁸, free basal-water accumulation)
  and **fixed a real CR-dependence bug**: the CTS diffusivity treatment choked
  the basal interface with the temperate K₀ whenever the basal node was
  temperate, making the diagnosed melt collapse at the nominal cr = 0.1. The K₀
  choke now applies only to a genuine temperate *layer* (k_cts ≥ 2); a melting
  base (k_cts = 1) conducts to the Dirichlet-pmp base with cold-ice κ. Basal melt
  is now **cr-independent** (matching the three reference models) and converges
  to the analytic solution (Eq A14): nz = 401 → warm 2.32 vs 2.33, cold −2.01 vs
  −2.03 mm/a; the phase-III transient tracks Eq A15 to ~1%. Driver test **T3**
  added; re-vendored the correct −10 °C analytic reference (the previous file was
  a −5 °C variant).
- **P2 Exp B (done):** validated against Kleiner (2015) Experiment B (steady
  polythermal slab): 200 m, 4° incline, no geothermal flux, strain heating
  `Ψ = 2A(ρg sinγ)⁴(H−z)⁴`, constant downward advection `vz = −0.2 m/a`, surface
  −3 °C, β = 0, L = 3.35×10⁵. The cold region matches the analytic `T(z)`
  (Appendix A2) to RMS 0.06 °C; the temperate layer converges to the analytic as
  the conductivity ratio `cr → 0` (the K₀→0 limit, matching the paper's CR sweep
  10⁻¹→10⁻⁵). At cr = 10⁻⁴, nz = 201: **CTS = 19.5 m** (analytic 19.0),
  **base ω = 0.0205** (analytic 0.0207). Driver test **T4** added. Follow-up: an
  intermediate cr ≈ 10⁻² shows a coarse-grid instability (stable at nz = 401);
  the accurate/stable regime is cr ≤ 10⁻³.

All three standalone benchmarks (**T2 cold-limit, T3 Exp A, T4 Exp B**) now pass.

- **P3 (in progress, blocked):** the first-ever 2D enthalpy run (EISMINT-2 EXPA,
  `method="enth"`) crashes at ~13 ka with a NaN at the base of a thin margin
  column (temp runs the full 100 ka). Also fixed a blocker along the way: the
  shipped benchmark par files still listed `till_rate`/`H_w_max` under `&ytherm`
  (moved to `&fhyd`), which `nml_validate` rejected — removed from 5 configs.
  Two enthalpy failure modes identified in thin, strain-heated, polythermal
  margin columns:
  1. **Advection-dominated** oscillation — **fixed** by Péclet-hybrid upwinding
     of the vertical advection (Spalding 1972, as the paper prescribes).
  2. **Near-singular temperate block** — with `enth_cr = 10⁻³` the temperate-layer
     nodes (κ = cr·κ_cold) blow up to ~−10⁷ J/kg while the cold layer stays fine.
     The temperate block becomes an almost-isolated Neumann problem (near-zero
     flux at the base via `enth(1)=enth(2)` and at the CTS top via the K₀ choke)
     with an internal strain-heating source, so it grows unboundedly. Exp B
     survives because its strong downward advection couples the block; the thin,
     weakly-advected EISMINT margin column (uz ≈ 0.003) does not.

  **Three candidate fixes tried, all still crash at ~13 ka** (each kept the 1D
  benchmarks passing):
  - **Diffusivity floor / raise `enth_cr`** (the approved Option A): cr = 10⁻²,
    5×10⁻², 10⁻¹ all still crash (cr = 10⁻¹ crashes *earlier*, at 9 ka — the
    failure is non-monotonic in cr, confirming it is conditioning, not magnitude).
    Only cr = 1.0 (no temperate reduction at all) is stable.
  - **Double-precision Thomas solve** (single precision was a suspect, given
    enth ~ 5×10⁵ with O(10) differences): still crashes → the discretized system
    is genuinely near-singular, not just rounding.
  - **Pin the temperate base at the pmp** (Dirichlet instead of `enth(1)=enth(2)`):
    still crashes → the singularity is in the block interior / CTS-top coupling,
    not only the base.

  **Assessment:** 2D robustness of the enthalpy solver at thin polythermal
  margins is a genuine, multi-faceted problem (the Kleiner benchmarks are all
  1-D columns; robust 3-D enthalpy is years of engineering in models like PISM).
  Likely also involves the very large horizontal enthalpy advection at margins
  (advecxy ~ −30, i.e. cp× the temperature-solver value). This is the gating item
  before the enth default flip and needs a dedicated effort — candidates: adopt
  Aschwanden et al. (2012)'s exact CTS/temperate discretization (PISM's, proven in
  3-D); limit/clip the enthalpy solution or the horizontal enth advection at
  margins; or special robust handling for thin margin columns.

- **P3 (2D crash FIXED):** the near-singular temperate block was caused by the
  **CTS diffusivity choke**, not the base or advection. The cold↔temperate
  interface was forced to ~2·K0 (harmonic mean of κ_temp, κ_cold) plus an
  explicit `kappa_a = kappa(k-1)` override at `k_cts+1`; with the zero-flux
  temperate base this isolated the temperate sub-block into a near-singular
  Neumann island with an internal strain source. **Fix:** conduct across the CTS
  with the **cold-side (upper-node) κ** (Blatter & Greve 2015, Eq. 25; the
  icetemp reference, which explicitly rejects the harmonic mean at the CTS and
  had this half-written in a commented-out `kappa_b = kappa(k+1)` line). The
  block now drains its strain heat upward by cold-ice conduction toward the CTS;
  both CTS rows use the same cold κ, conserving the interface flux. Result: **2D
  EISMINT-A `method="enth"` runs the full 100 ka** (was NaN at 13.24 ka), dome
  H0 = 4008 m (analytic ~3990). 1D benchmarks: T2 PASS, T3 PASS, T4 PASS at
  nz=201/401, cr=1e-4 (CTS 20.4 m vs analytic 19.0, within the 2 m gate — ~1.4 m
  less accurate than the choke, the accuracy/robustness trade being worth it; the
  advection-dominated Exp B is nearly insensitive to the CTS conductive flux).
  Commit `276215c5`.

- **P3 (EISMINT-2 A+F comparison + geothermal-flux bug FIXED):** the
  enth-vs-temp comparison (`analysis/compare_enth_temp.jl`, EISMINT-2 A and F,
  100 ka) first exposed a large enth cold-base bias: summit basal homologous
  temperature enth −32.6 °C vs temp −16.0 (EISMINT-A), a near-isothermal cold
  divide column, with a thicker dome and strongly amplified grid-aligned basal-
  temperature spokes (T′_b azimuthal-asymmetry RMS enth 6.6 vs temp 1.4 in F).
  **Root cause:** `calc_enth_column` used the basal/bedrock heat fluxes
  (`Q_b`, `Q_lith`) raw in [mW m⁻²], while `calc_temp_column` converts them to
  [J m⁻² a⁻¹] (`×1e-3·sec_year`). Both 2D-driver callers pass raw [mW m⁻²], so
  the enthalpy base got ~3.15×10⁴× too little geothermal heat. (The Kleiner 1D
  tests hid it: the standalone driver pre-converted `Q_lith` and used `Q_b=0`.)
  The Péclet upwinding was ruled out first (centered advection gave the same
  −32.6). **Fix:** convert `Q_b`/`Q_lith` inside `calc_enth_column` like the
  temp path, and pass raw `Q_rock` from the standalone driver — unifying the
  units contract ([mW m⁻²] in). Commit `e7cde320`. **Result:** enth now matches
  temp on both experiments — summit T′_b −16.0 vs −16.0 (A), −29.7 vs −29.7 (F);
  domes, volumes, melt fractions agree; the spoke asymmetry collapsed to temp's
  level (A 0.28 vs 0.30; F 1.46 vs 1.36). enth and temp now differ from the pure-
  SIA no-slip EISMINT reference identically (that offset is the Yelmo hybrid+
  sliding config, not the solver). All three 1D benchmarks (T2/T3/T4) unchanged.
  **Remaining before the default flip:** an Antarctica spin-up, enth vs temp
  (Phase 3 proper).

---

## 1. Motivation

Yelmo carries a full enthalpy formulation alongside the temperature solver, but
production has always used `temp` because `enth` "seemed less trustworthy". This
document records what the enthalpy path actually needs, why it has never run, and
a phased plan (with standalone validation) to make it the default.

## 2. Current state of `ytherm`

### 2.1 Dispatch

`calc_ytherm` ([`src/yelmo_thermodynamics.f90`](../../src/yelmo_thermodynamics.f90))
selects on `thrm%par%method ∈ {enth, temp, robin, robin-cold, linear, fixed}`.

- `enth` and `temp` share one driver, `calc_ytherm_enthalpy_3D`, which loops
  columns and calls `calc_enth_column` or `calc_temp_column`.
- `temp` is `enth` with `enth_cr=1.0` and `omega_max=0.0` forced at par-load
  (so temperate ice has cold-ice conductivity and zero water content).
- Default is `method="temp"` (`input/yelmo_defaults.nml`); **every** shipped par
  file uses `temp` or `fixed`. No shipped configuration uses `enth`.

The enthalpy *scaffolding* (horizontal advection of `enth`, `convert_to/from_enthalpy`,
CTS index, omega diagnosis, per-column enth output) is present and structurally
exercised. Only the per-column solve is incomplete.

### 2.2 The blocking defect

`calc_enth_column` ([`src/physics/ice_enthalpy.f90`](../../src/physics/ice_enthalpy.f90))
terminates the run:

```fortran
! Calculate basal mass balance
write(*,*) "calc_enth_column:: routine needs to be updated with Q_rock etc."
stop
```

`git blame` places this `stop` in the **original v1.15 import** — it is unfinished
upstream code, not a regression. So `method="enth"` has never been runnable in this
repository. The intended helper `calc_bmb_grounded_enth` does not exist in
`thermodynamics.f90` (only referenced in a comment and by the dead poly module).

### 2.3 Reference implementation exists

The testbed [`github.com/alex-robinson/icetemp`](https://github.com/alex-robinson/icetemp)
(the ancestor of this code) contains a **complete, working** `calc_enth_column` and
`calc_bmb_grounded_enth`, and ships **Kleiner et al. (2015)** benchmark data. This is
our reference for finishing the Yelmo routine and for validation (see §5).

## 3. Gap list (temp vs enth)

| # | Item | `temp` | `enth` | Action |
|---|------|--------|--------|--------|
| 1 | Basal mass balance | `calc_bmb_grounded` (works) | **`stop`; no bmb** | Port `calc_bmb_grounded_enth` |
| 2 | Basal heat flux `Q_ice_b` | temperature gradient | temp gradient (enth gradient commented "Problematic") | Decide + validate discretization |
| 3 | CTS treatment | none | post-solve `enth(1)=enth(2)` fixup; zero-flux BC; one-sided diffusivity jump at `k_cts+1` | Justify or replace vs benchmark |
| 4 | Basal BC units/sign | `-(Q_b+Q_rock)/kt` | `(Q_b+Q_lith)/kt*cp`, pre-converted `Q_lith` | Reconcile so cold-limit ≡ temp |
| 5 | Surface BC | raw `min(T_srf,T0)` | `min(T_srf,T0)*cp(nz)` | Reconcile |

**Correctness anchor:** in the cold limit (`omega=0`, `enth_cr=1`) the enthalpy
solver must reduce *exactly* to the temperature solver. This equivalence is the
cheapest, strongest regression test and is currently untested.

### 3.1 Reference physics for item 1 (flux-based bmb)

From the icetemp testbed (`thermodynamics.f90`):

```fortran
net_enth = enth_b - enth_pmp_b            ! basal enthalpy above pmp
Q_net    = Q_b + Q_ice_b + Q_geo_now      ! net basal energy flux [J a-1 m-2]
bmb_grnd = -Q_net / (rho_ice*L_ice - net_enth)
```

This is the **flux-based** form (chosen for consistency with `calc_bmb_grounded`),
with an enthalpy correction `- net_enth` in the denominator accounting for latent
heat already stored as basal water. The Yelmo port must adapt the interface to the
current boundary set (`Q_rock` bedrock flux, `W_til` predictor, `f_grnd` weighting).
The design should keep the bmb closure **swappable** (flux-based default, leaving
room to experiment with an enthalpy-state form).

## 4. Standalone test harness

A single-column driver already exists — `tests/test_icetemp.f90`, target
`make icetemp` — that builds an EISMINT column, computes the Robin analytical
solution, time-steps a solver, and writes NetCDF. **It is bit-rotted**: it calls
`calc_enth_column`/`calc_temp_column` with old, shorter signatures and will not
compile against the current routines. `tests/test_icetemp_poly.f90` targets the
dead poly solver.

Plan for the harness (Phase 0):

1. Update calls to current signatures; drop the poly test.
2. Add automated pass/fail assertions:
   - **T1** temp-solver vs Robin steady state (existing comparison, made quantitative).
   - **T2** enth ≡ temp in the cold limit (bit-for-bit or within tol).
   - **T3** Kleiner (2015) **Experiment A** — transient basal melt rate vs
     `data/Kleiner2015/Kleiner2015_EXPA_Fig2-IIIa-melt.txt`.
   - **T4** Kleiner (2015) **Experiment B** — steady polythermal column: CTS
     position, `T(z)`, `omega(z)`, basal enthalpy vs
     `Kleiner2015_EXPB_analytic_nz401_z.dat` and `..._enth_b.txt`.
3. Reference data vendored under `tests/data/Kleiner2015/`.

## 5. Kleiner et al. (2015) benchmark

Two experiments with published/analytic references (vendored from icetemp):

- **Experiment A** — parabolic surface-temperature forcing over a cold slab;
  tests the cold↔temperate base transition and basal melt onset. Reference:
  basal melt-rate time series.
- **Experiment B** — steady polythermal column with a temperate basal layer;
  tests CTS position and water content. Reference: analytic profile at `nz=401`
  and basal enthalpy series (`T_ref = 173.15 K`).

These target gap items 2–3 directly (CTS discretization, enthalpy-gradient flux).

## 6. Phased plan

**Phase 0 — Revive the standalone harness.** Update `test_icetemp`, vendor Kleiner
data, add tests T1–T4. De-risks everything downstream with a fast build/run loop.

**Phase 1 — Finish the enthalpy column solver.** Port `calc_bmb_grounded_enth`
(flux-based, swappable), remove the `stop`, resolve the `Q_ice_b` discretization,
reconcile basal/surface BC units and sign so cold-limit ≡ temp (test T2).

**Phase 2 — Validate CTS / polythermal.** Justify or replace the CTS fixups against
Kleiner B (T4); confirm energy conservation.

**Phase 3 — 2D/3D validation.** EISMINT (`make benchmarks`) and an Antarctica
spin-up, `enth` vs `temp`: basal temperate area, bmb, thickness. Only then consider
changing the default.

**Phase 4 — Cleanup (lands early, incrementally).**
- Delete dead `src/physics/ice_enthalpy_poly.f90` and `tests/test_icetemp_poly.f90`
  (never compiled; would not compile — uses nonexistent `calc_bmb_grounded_enth`).
- Remove dead helpers `calc_Q_bedrock`, `calc_hires_cell` (zero call sites).
- Hoist magic numbers (`enth_ref=273.15*2009.0`, `H_ice_thin=10.0`) to constants.
- Remove stale comments/signature cruft (`W_til` "kept for signature stability" in
  `define_temp_robin`; "remerge into icetemp" module note).
- Track the OpenMP-in-bedrock "leads to NaNs" note and the gated `ajr symtest` block.

## 7. Open technical questions

- **Q_ice_b discretization** (gap 2): testbed uses enthalpy-gradient
  `kappa*rho_ice*(enth(2)-enth(1))/dz`; Yelmo reverted to temperature gradient and
  flagged the enthalpy form "Problematic". Kleiner B should decide which is correct
  and stable.
- **CTS boundary condition** (gap 3): are the zero-flux BC and one-sided diffusivity
  jump physical, or stability band-aids? Validate against Kleiner B.
- **bmb closure form**: flux-based confirmed as default; keep an interface seam for
  an alternative enthalpy-state closure.
