# Yelmo benchmark protocol (draft)

*Status: draft for iteration, 2026-10-05. Nothing here is implemented yet unless marked "exists".*

## Purpose

This protocol defines a ladder of benchmark tests that verifies Yelmo in the configuration it is used in for production runs. The EISMINT experiments are based on the shallow-ice approximation and use artificial forcing, while MISMIP3D and MISMIP+ depend strongly on the treatment of lateral walls. Neither is a good test of the full momentum solver with land, marine and floating margins together. Here we define tests that are quick to run, detect symmetry and conservation errors directly, and end with the full production physics on a closed domain.

## Design principles

1. **No lateral boundaries.** Every domain is either a closed island surrounded by open ocean, with a calving mask that keeps ice away from the domain border, or periodic in y. The result therefore does not depend on wall conditions.
2. **Symmetry built in.** Geometry and forcing are invariant under the D4 group of the square grid (rotations by 90° and reflections about the axes and diagonals) wherever possible. Any asymmetry in the solution is then a model error. A transposed grid (x↔y) must reproduce the transposed solution to round-off, which detects i/j indexing errors.
3. **Production defaults as the base.** All tests use the default parameter file `input/yelmo_defaults.nml`, which holds the production physics (momentum solver, enthalpy, friction law, grounding-line treatment, front subgrid scheme, calving). The parameter file of each test lists only the settings it changes: geometry, forcing and the one ingredient it isolates. The ladder thus ends at the production setup and not at a test-only configuration.
4. **Short by design.** Tests with a known steady state start from it and measure drift. Symmetry errors grow exponentially, so they appear within about 1 kyr. Diagnostic tests need a single velocity solve.
5. **One ingredient per step.** Each tier adds one component, so that a failure points to its cause.
6. **Isothermal and coupled pairs.** Every prognostic test runs isothermally (constant rate factor) and with coupled enthalpy. A failure only in the coupled run points to the thermodynamics or the thermomechanical coupling, as found for EISMINT EXPF with the enthalpy solver.
7. **Inputs are generated externally.** Geometry, forcing, initial states and analytic reference fields are generated in Julia before the run and read as input. The Julia code follows the conventions of IceSheetBenchmarks.jl, where the tests will move once they are established. The protocol thus does not require the analytic solutions (e.g., IceColumnSolutions.jl) to be implemented in Yelmo.

## Overview

| Tier | Test | Isolates | Reference | Cost |
|---|---|---|---|---|
| 0 | Column tests (exists) | Vertical thermodynamics | IceColumnSolutions.jl, Robin (1955) | s |
| A | A1 Slab (exists) | Grounded momentum balance, friction | Analytic (Schoof, 2006) | s |
| A | A2 Radial floating shelf | Floating momentum balance, calving-front stress condition, staircase front | Analytic | s |
| A | A3 ISLAND4 diagnostic | Full solver on land, marine and floating ice | D4 symmetry to round-off | s |
| A | A4 3D thermodynamics strip | 3D enthalpy solver, vertical velocity, transients | IceColumnSolutions.jl | s–min |
| B | B1 ISLAND4-L land dome | Land margin, sliding, thermomechanical coupling | Symmetry, mass conservation | min |
| B | B2 Marine flowline | Grounding-line dynamics | Schoof (2007), Tsai et al. (2015) | min |
| B | B3 CalvingMIP circular Exp1 (exists) | Floating front, front advection | Front fixed at r = 750 km, symmetry | min |
| C | C0 ISLAND4 spin-up | Full production setup | Run once per version | h |
| C | C1 Control from C0 | Drift, symmetry, conservation | Zero drift | min |
| C | C2 Perturbations from C0 | Response of all margin types | Resolution convergence, regression | min |
| C | C3 Grid tests | Grid orientation and indexing | Transposed and rotated runs | min |
| D | GRL, ANT at 32 km | Realistic geometry | Regression, budget closure, timing | min |

## Common setup

### Grids

The grid must be symmetric about the origin, with nx = ny and x_i = −x_{nx+1−i}. The origin can lie on a cell center (nx odd) or on a cell corner (nx even). Both are used in C3. Standard resolutions for the ISLAND4 domain (x, y ∈ [−800, 800] km) are 32, 16 and 8 km.

### Metrics

Every test writes the same diagnostics, which the scorecard script reads.

- **Symmetry error.** For each element g of D4 and each field f, e_g(f) = max|f − g·f| / max|f|. Scalar fields are mapped directly. Velocity components change sign and swap between the x- and y-faces of the C-grid, as described in [`tests/README_symcheck.md`](../../../tests/README_symcheck.md). Fields checked: H_ice, ux, uy, enthalpy (all layers), f_grnd. For the ISLAND4 domain we report the maximum over all eight group elements and its time evolution.
- **Mass budget residual.** r_M = [ΔV − ∫(SMB + BMB − calving) dt] / V over each output interval, from the Yelmo budget terms.
- **Integrated state.** Ice volume V, grounded area A_g, floating area A_f, and the grounding-line and calving-front radii along the axes and diagonals.
- **Solver statistics.** Picard iterations per solve, mean and minimum time step, wall time.

### Scorecard

`scripts/scorecard.jl` (to be written) reads each run, computes the metrics and prints one pass/fail line per test. The symmetry operators generalize `tests/symcheck.jl` from mirror reflections to the full D4 group.

### Pass criteria

Thresholds are set empirically after the first round of runs and recorded here. The provisional targets are e_g < 1e-10 for diagnostic tests, e_g growing by less than a factor of 10 per kyr in prognostic tests, and |r_M| < 1e-6 per output interval.

## Tier 0: column tests (exists)

`tests/test_icetemp.f90` (experiment `ics`) and `tests/test_enthalpy.f90` (`robin-column`) referee the column solvers against IceColumnSolutions.jl (Moreno-Parada et al., 2024) and the Robin (1955) solution. The driver `analysis/compare_column_analytic.jl` measures the order of accuracy over an nz sweep. These tests cover cold columns only. Temperate and polythermal columns are checked by `robin-column` and indirectly by the coupled runs in tiers B and C.

## Tier A: diagnostic tests

### A1 Slab (exists)

Tilted slab on a periodic domain (Schoof, 2006), as documented in `docs/benchmarks.md`. It tests the grounded momentum balance and the friction law against an analytic solution.

### A2 Radial floating shelf

A2 tests the floating momentum balance and the calving-front stress condition against an exact solution. The setup is a circular shelf of constant thickness H and radius R_s on a Cartesian grid. The front is therefore a staircase, which is the configuration found in production runs.

For constant H and constant rate factor A, the exact solution is isotropic spreading, u = ε̇ x and v = ε̇ y. The membrane stresses are uniform, so the interior momentum balance holds trivially, and the front condition sets ε̇. With ε̇_xx = ε̇_yy = ε̇, the effective strain rate is ε̇_e = √3 ε̇ and the front condition can be written as

$$6\,\eta\,\dot\varepsilon = S, \qquad S = \tfrac{1}{2}\,\rho_i g H \left(1 - \frac{\rho_i}{\rho_w}\right), \qquad \eta = \tfrac{1}{2} A^{-1/n} \dot\varepsilon_e^{(1-n)/n},$$

which gives

$$\dot\varepsilon = A\, S^n\, 3^{-(n+1)/2},$$

where ρ_i and ρ_w are the densities of ice and seawater and n is the flow-law exponent. For n = 3, ε̇ = A S³/9. For comparison, the one-dimensional unconfined shelf gives A S³/8 (Weertman, 1957). *This derivation needs to be checked against the Yelmo definitions of η and ε̇_e before use.*

A shelf that is floating everywhere has a nullspace (rigid translation and rotation). We remove it with a small uniform basal friction β_reg. This modifies the solution by a relative amount of order β_reg R_s² / (η H), so β_reg is chosen to make this term smaller than about 1e-4 while keeping the system well conditioned.

Proposed values: H = 400 m, R_s = 300 km, domain ±400 km, dx = 20, 10 and 5 km, n = 3, A constant. Metrics: relative error of u and v against the analytic solution in the interior and at the front, and D4 symmetry error.

### A3 ISLAND4 diagnostic

A3 runs a single velocity solve with the full production solver on the ISLAND4 geometry (see below), with prescribed thickness and enthalpy. The thickness is a Vialov profile, H(r) = H_0 [1 − (r/R_i)^{(n+1)/n}]^{n/(2n+2)}, with H_0 = 3500 m and R_i = 650 km. Ice floats wherever the flotation criterion is met. Prognostic runs (B1, C0) start from the same profile without the floating part (`init = :vialov_grounded`), so that the shelves form during the run. This produces grounded ice on land and below sea level, a grounding line and small shelves in the embayments. The enthalpy field is the stationary column solution of IceColumnSolutions.jl in each column, computed from the local H, SMB, T_srf and Q_geo. It is written to an input file.

The test passes if all D4 symmetry errors are at round-off level. The same solve on ISLAND4-R (rotated by 45°) gives the dependence on grid orientation, measured by the integrated speed and the grounding-line flux.

### A4 3D thermodynamics strip

A4 tests the 3D enthalpy solver against exact transient and stationary solutions. The column tests in Tier 0 do not reach the 3D code path, the staggering of the vertical velocity or the coupling to the bedrock layer. A4 is designed so that every column satisfies the assumptions of the analytic solution exactly: linear vertical velocity, no horizontal advection of enthalpy, cold base.

The domain is a strip with x ∈ [−L_x, L_x] and y ∈ [0, L_y]. The parameters H, SMB, T_srf and Q_geo vary along y only. The velocity is prescribed (dynamics fixed, thickness fixed) as

$$u = c(y)\,x, \qquad v = 0, \qquad c(y) = \frac{\mathrm{SMB}(y)}{H(y)},$$

with zero basal mass balance. The horizontal divergence equals SMB/H, so continuity gives the linear vertical velocity u_z = −SMB ζ. The thickness is steady, since ∇·(uH) = SMB. Enthalpy varies along y only and v = 0, so horizontal advection is exactly zero, also in the discrete equations. Each row of the strip is thus an independent column case, and a single run covers a range of Péclet numbers and basal heat fluxes. Strain heating is switched off (`ytherm.strain_heating = "none"`). The parameters are chosen so that every base stays below the pressure-melting point.

- **A4a Stationary.** Initialize with the stationary analytic solution and run 5 kyr. The drift must be small. A second run starts from a linear profile and is compared with the analytic solution once it is stationary.
- **A4b Transient.** Start from the stationary state of parameter set P1 and switch to P2 (T_srf − 10 K, or SMB × 2). Compare with the transient solution (Moreno-Parada et al., 2024, Appendix A) at t = 1, 5, 10 and 50 kyr.
- **A4c Strain heating (optional).** With a constant rate factor, plug flow gives depth-uniform strain heating, which the stationary solution includes through the Brinkman number. This tests the strain-heating term.

Each parameter set occupies three identical rows, and only the middle row is compared, in the columns at least two cells from the x-borders. The default vertical velocity (`ydyn.uz_method = 3`) averages the divergence over neighbouring cells, with weights 1/4, 1/2, 1/4 across rows, so a single row per set would mix the parameter sets of adjacent rows (see Status). The strip is symmetric under x → −x, and a transposed run tests the y-direction.

## ISLAND4 domain

ISLAND4 is a closed island with land-terminating, marine-terminating and floating margins, and it is invariant under the D4 group of the grid. It consists of a radial base profile with four V-shaped troughs cut along the diagonals. The bed elevation is defined as

$$z_b(x,y) = B_c - (B_c - B_l)\frac{r^2}{R_0^2} - \sum_{k=1}^{4} D(\xi_k)\, T(\xi_k,\eta_k),$$

$$T(\xi,\eta) = \tfrac{1}{2}\left[1 + \tanh\!\left(\frac{(\xi - r_h)\tan(\alpha/2) - \sqrt{\eta^2 + \ell^2}}{\ell}\right)\right], \qquad D(\xi) = D_0 + B_{od} \exp\!\left(-\frac{(\xi - r_{od})^2}{w_{od}^2}\right),$$

where (ξ_k, η_k) are the coordinates along and across the k-th diagonal, and T is a smoothed indicator of the trough. Its walls are straight lines from the apex at ξ = r_h with opening angle α, smoothed over the width ℓ. The parameters are given in the table below.

| Parameter | Value | Meaning |
|---|---|---|
| B_c | 900 m | Bed elevation at the center |
| B_l | −2000 m | Bed elevation of the base profile at R_0 |
| R_0 | 1000 km | Radius scale of the base profile |
| r_h | 250 km | Radius of the trough head |
| α | 120° | Opening angle of the troughs |
| ℓ | 50 km | Width of the trough walls and head |
| D_0 | 1500 m | Trough depth below the base profile |
| B_od | 0 or 700 m | Depth of the overdeepening (0: off) |
| r_od | 350 km | Radius of the deepest point of the overdeepening |
| w_od | 60 km | Half-width of the overdeepening |

The island forms a four-pointed star (Fig. 1a). The coast reaches r ≈ 556 km along the axes and r ≈ 278 km along the diagonals. Ice in the troughs grounds below sea level and feeds the shelves. With the control forcing, the ice also grounds beyond the coast of the peninsulas, so a continuous shelf surrounds the island out to r_lim. Land-terminating margins are therefore tested in ISLAND4-L (B1). The wide opening angle gives smooth, concave coastlines, so that the geometry does not favor the onset of instabilities. With B_od = 0, the trough beds deepen monotonically seaward, so the steady state is unique, which a regression test needs. Setting B_od > 0 adds an overdeepening with a retrograde bed between r ≈ 360 and 448 km (Fig. 1b, d), which tests marine ice-sheet instability. A single parameter thus switches the instability test on and off.

Rotating the troughs by 45° (ISLAND4-R) places them along the grid axes (Fig. 1c). The base profile and the forcing are radial, so ISLAND4-R has the same solution as ISLAND4 in the continuum, rotated by 45°. Differences between the two runs measure the dependence of the discrete solution on grid orientation.

![ISLAND4 geometry. (a) Bed elevation with the coastline (black), the equilibrium line r_ela (red dashed) and the calving mask r_lim (black dotted). (b) The overdeepened variant (B_od = 700 m). (c) ISLAND4-R, rotated by 45°. (d) Bed profiles along the axis and the diagonal (lines in a, b) and the SMB profile (red, right axis). (e) Bed across the trough at r = 450 km, without (solid) and with (dashed) overdeepening.](figures/island4.png)

### Forcing

The forcing is radial or depends on the local state only, so it preserves the D4 symmetry. Fixed fields are generated in Julia with the geometry. Forcing that depends on the evolving ice geometry is computed during the run by a Fortran function in the driver. Each online function has a Julia counterpart, so that the forcing written by the model can be checked against the specification.

- **Surface mass balance (fixed).** SMB(r) = SMB_0 (1 − r/r_ela), with SMB_0 = 0.5 m a⁻¹ and r_ela = 650 km (450 km for ISLAND4-L). For a linear profile, the net SMB of a disc of radius R vanishes at R = 1.5 r_ela. With r_ela = 650 km, this radius lies beyond r_lim, so the marine margins are set by melt and calving, and the SMB on the shelves stays between about −0.1 and +0.1 m a⁻¹. ISLAND4-L keeps r_ela = 450 km, since its land margin needs ablation inside r_lim. We chose an SMB that does not depend on surface elevation, which avoids the elevation feedback and multiple steady states.
- **Surface temperature (online).** T_srf = T_sl − Γ z_s, with T_sl = −10 °C and Γ = 8 K km⁻¹. The temperature follows the evolving surface, which keeps it consistent with the ice geometry. Its feedback on the dynamics acts through the thermodynamics only and is weak.
- **Geothermal heat flux (fixed).** Q_geo = 50 mW m⁻², uniform.
- **Sub-shelf melt (online).** The MISMIP+ Ice1 parameterization (Asay-Davis et al., 2016), m = Ω tanh(H_c/H_c0) max(z_0 − z_d, 0), with Ω = 0.01 a⁻¹, H_c0 = 75 m and z_0 = −200 m, where z_d is the depth of the ice base and H_c the water-column thickness. We chose it because it is a community standard, and the tanh term reduces the melt smoothly toward the grounding line. The MISMIP+ values (Ω = 0.2 a⁻¹, z_0 = −100 m) are sized for the ice flux of MISMIP+, and on ISLAND4 they remove about 30 times more ice than the island supplies (about 3e10 m³ a⁻¹ per trough with r_ela = 450 km). The values used here give no melt below shelves thinner than about 225 m, about 1.5 m a⁻¹ below 400 m of ice and about 7 m a⁻¹ at a grounding line 1 km deep.
- **Calving.** The calving law of the default parameter file, together with a fixed calving mask that removes all ice at r ≥ r_lim = 750 km. The mask keeps ice away from the domain border.
- **Friction.** The friction law of the default parameter file with spatially uniform parameters.

The values of r_ela, Ω and z_0 come from a sweep of 5-kyr runs at 32 km (2026-10-05). With the original values (r_ela = 450 km and the MISMIP+ melt), the shelves were about 50–65 m thick on average and ended inside r_lim. With the values above, the shelves reach r_lim, with a mean thickness of about 150 m and about 90 m at the front. Weaker melt or more accumulation increases this only slightly. The shelf thickness is limited by the geometry: beyond the coast the shelves are unconfined and spread under their own weight (strain rate ∝ H^n), so they thin to a few hundred meters within one or two cells of the grounding line. Raising B_l to −840 m extends the peninsulas to about 660 km, but the troughs then form 90° sectors between narrow spurs and give little lateral support (mean shelf thickness about 165 m, about 110 m at the front). We therefore kept B_l = −2000 m. Thin unconfined shelves are the physical solution of this geometry, and they still test the floating momentum balance, calving at r_lim and sub-shelf melt.

## Tier B: one margin type at a time

### B1 ISLAND4-L land dome

ISLAND4-L is ISLAND4 with the bed raised by 2500 m, so that all ice inside r_lim rests on land. The bed keeps its D4 structure, so the flow is not radial and the transposed and rotated tests remain meaningful. The ice starts from the Vialov profile of A3 with the column enthalpy from IceColumnSolutions.jl, and the run lasts 5 kyr with the production solver and sliding.

The run is not in thermal equilibrium after 5 kyr. However, symmetry and conservation do not require equilibrium, and the transient margin adjustment exercises the land margin. Metrics: growth of the symmetry error over time (isothermal and coupled), mass budget residual, and comparison with ISLAND4-L-R.

### B2 Marine flowline

B2 tests grounding-line dynamics in a strip that is periodic in y, with a small number of cells across the flow (about 3–5). It follows the MISMIP experiment 1 setup (linear prograde bed) for several values of A. The periodic strip is a true one-dimensional problem, so there are no lateral walls.

The reference is the boundary-layer flux of Schoof (2007) for a power-law friction law, and of Tsai et al. (2015) for a Coulomb friction law. The run starts from the semi-analytic steady profile, computed in Julia by integrating the inland profile upstream from the grounding line given by the flux formula. This shortens the run to the adjustment time, about 2–5 kyr, instead of a full spin-up from zero. The test reports the final grounding-line position against the reference, and dx_g/dt over the run, at dx = 1, 2, 4 and 8 km. The boundary-layer formulas are asymptotic and agree with converged SSA solutions to within a few percent, so the tolerance must allow for this.

### B3 CalvingMIP circular Exp1 (exists)

The driver `tests/yelmo_calving.f90` runs CalvingMIP Experiment 1 on the circular domain, where the calving rate equals the front velocity so that the front stays at r = 750 km. Metrics: front radius along the axes and diagonals, and the symmetry error from `tests/symcheck.jl`.

## Tier C: integrated ISLAND4

Tier C runs the production setup on ISLAND4. It is the main target of the protocol.

- **C0 Spin-up.** The run starts from the Vialov profile and column enthalpy of A3 and continues until the drift criterion is met (|dV/dt|/V below a threshold per kyr, to be set). The thermal adjustment takes about 10⁴–10⁵ yr. C0 runs once per model version and resolution (32, 16 and 8 km) on Levante. The restarts are stored and used by C1–C3.
- **C1 Control.** 1 kyr from the C0 restart without changes. Metrics: drift in V, A_g and A_f, symmetry error and mass budget residual.
- **C2 Perturbations.** 500 yr from the C0 restart.
  - C2a Ocean: Ω × 2.
  - C2b Shelf removal: all floating ice removed at each time step (as in ABUMIP).
  - C2c Surface: SMB − 0.1 m a⁻¹ everywhere (r_ela moves to 520 km).
  - C2d Reversibility (optional): C2a followed by 500 yr of control forcing.

  The references are resolution convergence (32, 16, 8 km) and a regression envelope against the previous model version. The run length of 500 yr will be shortened if the signals become clear earlier.
- **C3 Grid tests.**
  - C3a Transposed grid: must reproduce the transposed C1 solution to round-off.
  - C3b ISLAND4-R: C0 and C1 on the rotated geometry. The differences in V, A_g, A_f and the grounding-line radii (axes and diagonals exchanged) measure the dependence on grid orientation.
  - C3c Half-cell shift: origin on a cell corner instead of a cell center.

Isothermal versions of C1 and C2 are not part of the standard suite. They are run to diagnose a failure of the coupled runs.

## Tier D: realistic domains

Greenland and Antarctica at 32 km, 100 yr from a production restart (as in initMIP). Metrics: mass budget residual, time-step statistics, wall time and regression of V, A_g and A_f against the previous version. This tier checks that the production configuration still runs. It does not verify the physics.

## Thermodynamics in each tier

| Tier | Isothermal | Coupled |
|---|---|---|
| 0 | – | Column tests |
| A | A1, A2 | A3 (prescribed enthalpy), A4 |
| B | B1, B2, B3 against their references | B1, B2, B3 for symmetry, conservation and drift |
| C | C1, C2 for diagnosis only | C0–C3 |
| D | – | All |

The analytic references of A1, A2 and B2 assume a constant rate factor, so they apply to the isothermal runs only. The coupled runs are evaluated against symmetry, conservation and drift.

## Implementation

### Driver

All tests in tiers A–C run with one Fortran driver, `tests/yelmo_bench.f90`. The driver contains no geometry. It reads a fixture file generated in Julia, sets the boundary fields and the initial state from it, and runs Yelmo with the time stepping of the parameter file. Forcing that depends on the evolving state (the surface temperature and sub-shelf melt of ISLAND4) is computed by functions in a small Fortran module, selected in the parameter file. A test is therefore defined by its fixture and its parameter file only. This avoids code duplication between tests and gives a single executable for CI testing. The column tests (Tier 0) and the realistic domains (Tier D, `yelmo_initmip`) keep their own drivers.

The fixture holds the grid (xc, yc), the boundary fields (z_bed, z_sl, smb_ref, T_srf, Q_geo, bmb_shlf, T_shlf, H_sed, calv_mask) and the initial ice thickness, with the field names of IceSheetBenchmarks.jl. Optional fields are read when present: a 3D enthalpy or temperature field (A3, A4, B1, C0), a prescribed 3D velocity field for `ydyn.solver = "fixed"` (A4) and a basal friction field (A2). Perturbation experiments (C2) restart from C0 with a different fixture.

The existing slab (A1) and CalvingMIP Exp1 (B3) tests move into the same driver once it works.

### Fixtures

Each test is a Julia type that subtypes `AbstractBenchmark` from IceSheetBenchmarks.jl and implements `state(b, t)` and `write_fixture!(b, path)`, e.g., `Island4Benchmark(:ctrl; dx_km = 16, B_od = 0, rot = 0)`. Helper functions for the geometry (`island4_bed`) and the forcing are exported, as for the existing benchmarks. These include the Julia counterparts of the online Fortran forcing functions, which the scorecard applies to the model state to check the forcing written by the model. The code lives in `tests/bench/` as a Julia project that depends on IceSheetBenchmarks.jl. The types can thus move to IceSheetBenchmarks.jl without changes.

### Running

```bash
make bench
julia --project=tests/bench -e 'import Pkg; Pkg.instantiate()'
julia --project=tests/bench tests/bench/make_fixture.jl input/bench/island4-16km.nc island4 dx_km=16
runme -r -e bench -n par/yelmo_bench_ISLAND4.nml -o output/bench/a3-16km
julia --project=tests/bench tests/bench/check_symmetry.jl output/bench/a3-16km
```

Fixtures are written to `input/bench/` (not tracked), which runme links into every run directory. Keyword arguments of the benchmark constructor are passed as `key=value` (e.g., `B_od=700`, `rot=45`, `land=true`, `exp=smb`). With `time_end = time_init`, the run is diagnostic (A3). At the start of a fresh run, the driver prints the difference between its online forcing and the fixture values computed by the Julia counterparts. Partially ice-covered front cells are excluded from this check, since Yelmo defines their surface elevation with the effective front thickness.

### Status

The ISLAND4 fixture holds the initial ice temperature (`T_ice` on 51 uniform ζ levels) from the stationary column solution of IceColumnSolutions.jl, with w0 = max(SMB, 0) for grounded columns and a linear profile from T_shlf for floating columns. The driver maps it onto the Yelmo levels with `yelmo_read_remap` and initializes with `thrm_method = "prescribed"`.

First A3 run (2026-10-05, ISLAND4, 16 km, Vialov thickness, Robin temperature):

- H_ice and z_srf are exactly D4-symmetric, and the velocity components are symmetric to ~1e-7, which is the precision of the single-precision output.
- The online T_srf and bmb_shlf agree with the Julia counterparts to ~3e-5 K and ~2e-4 m a⁻¹ (100 front cells excluded).
- The initial SSA solve needs 21 Picard iterations, one more than the default `ydyn.ssa_iter_max = 20`.

First A2 run (2026-10-05, 10 km, H = 400 m, R_s = 300 km, A = 1e-18, β_reg = 1e-3, Picard converged to 1e-5):

- The velocity is D4-symmetric to output precision, but does not reproduce the analytical solution. With the default energy assembler, the interior speed is ~13–33 % of the exact value (rms error 44 %); with `ssa_solver = "residual"` the rms error is ~15 %. DIVA and SSA give the same result.
- The viscosity agrees with Glen's law for the simulated strain rates, so A and n are applied correctly. The membrane stress N_rr diagnosed from the output is ~1.2 H S near the front and ~0.5 H S at the centre, which implies a distributed resistance of ~190 Pa inside the shelf (β_reg u is ~0.1 Pa).
- The result does not change with `front_subgrid`, with a linear-solver tolerance of 1e-10, or with `ssa_lat_bc` ("all", "floating", "marine", "none") for the energy assembler.
- Cause: the corner (shear) viscosity on the ice margin coupled the front faces to the u = 0 faces of the ice-free cells, a drag on the velocity along the front. With zero corner viscosity on the margin (dev eb2879de), both assemblers match the analytical solution to ~0.07 % (rms).

Prognostic symmetry (2026-10-05, ISLAND4, 32 km, 20–100 yr). With the default solver tolerances and the Vialov initial state including floating ice, the symmetry error reached 2.6 % after 20 yr with a driver step `dtt = 10` (1e-5 with `dtt = 1`). The diagnosis found three contributions:

- The default tolerances (Picard `ssa_iter_conv = 1e-2`, linear solver `-tol 1.0e-2`) determine the velocity only to ~1 %, so a single solve can turn a round-off asymmetry of 1e-6 into 1e-3, with symmetric masks and inputs.
- The Vialov initial state puts up to ~3 km of floating ice in the troughs, with speeds at the velocity limit and fast-moving fronts, which amplify round-off.
- A cell holding a round-off amount of ice became a front cell with the full reference thickness, so the front force jumped by one cell (fixed with `H_ice_eps = 1 mm`, as CISM's `eps11`).

The ISLAND4 parameter file therefore uses tight solver tolerances (Picard 1e-5, linear 1e-8), and prognostic runs start from `init = :vialov_grounded` (Vialov profile where grounded, no initial floating ice; A3 keeps `init = vialov`). With these, the symmetry error stays at ~1e-5 (H) and ~7e-5 (velocity) over 100 yr with `dtt = 10`, which is round-off growth in single precision.

A double-precision build (docs/numerics-precision.md) separates round-off from real asymmetry: any error above about 1e-12 points to the code. It showed a deterministic asymmetry of about 1e-5 after 100 yr, which was traced through the level-set calving and the principal stress to the 3D strain rates. At a front with ice on the low-index side, the one-sided stencil in `calc_jacobian_vel_3D_uxyterms` tested `f_ice` of cell i−2 instead of i−1, while the mirror case tests i+1 (fixed in dev 987f3a8c). The strain rates near fronts were thus not mirror-symmetric, which switched the calving rate of single front cells.

Symmetry after the fix (2026-10-05, ISLAND4, 32 km, 1 kyr, `dtt = 10`, albedo):

| Run | 100 yr | 300 yr | 500 yr | 1000 yr |
|---|---|---|---|---|
| Double, tight tolerances | 3e-14 / 1e-12 | 1e-13 / 6e-13 | 2e-12 / 4e-12 | 6e-13 / 1e-11 |
| Double, default tolerances | 7e-15 / 7e-14 | 3e-10 / 6e-10 | 4e-11 / 6e-11 | 1e-8 / 7e-8 |
| Single, tight tolerances | 4e-6 / 4e-5 | 9e-7 / 1e-5 | 6e-3 / 3e-3 | 3e-3 / 0.9 |
| Single, default tolerances | 4e-7 / 6e-6 | 2e-5 / 1e-4 | 1e-3 / 7e-3 | 3e-2 / 1.0 |

The values are the maximum D4 errors of H and velocity. In double precision the solution stays symmetric to round-off over 1 kyr. In single precision the error stays at round-off level for about 300 yr and then jumps at isolated cells, where round-off differences switch a threshold (front cells, calving) in the thin outer shelf. These runs used the original forcing. With the tuned forcing (thicker shelves ending at r_lim), the single-precision error stays below about 2e-6 (H) and 4e-6 (velocity) over 1 kyr with tight tolerances, which is round-off level. With default tolerances it stays bounded at about 1e-3, the precision to which these tolerances determine the velocity. The single-precision breakdown was thus caused by the thin fronts of the original forcing. The H_ice_eps floor and the Jacobian fix change the final states of CalvingMIP Exp1, MISMIP+, MISMIP3D, TROUGH-F17 and A2 by less than 3e-5 in volume, with identical grounded and floating areas.

A4 (2026-10-05, 10 km, 7 × 38 cells, 12 parameter sets with Pe ≈ 0.5–33, T_srf = −40 °C, Q_geo = 50 mW m⁻²):

- With a single row per parameter set, the default vertical velocity (`uz_method = 3`) gave a surface uz of up to 5 times −SMB next to a jump in SMB/H between rows, and temperature errors of ~7 K after 5 kyr. Methods 2 and 3 average the divergence over the neighbouring cells; method 1 uses the face differences of the thickness equation and is exact here. In realistic flow the averaging is an advantage: method 3 gives a 1.5–4 times smaller surface mismatch `uz_srf_err` and a 2–5 times smoother uz than method 1 in Greenland (16 km, 1 kyr), ISLAND4 (16 km, 1 kyr) and TROUGH-F17 (4 km, 5 kyr). Both methods leave a mismatch of several m a⁻¹ at grounding lines and outlet margins, since neither is the discrete flux divergence of the thickness equation. The method also changes the ISLAND4 state after 1 kyr considerably (temperate basal fraction 0.17 with method 1 and 0.08 with method 3, volume difference 9 %), which shows the thermomechanical sensitivity of the early spin-up. We kept method 3 and use three rows per parameter set.
- Stationary test (A4a), 50 kyr from the analytic profile: the maximum error is 0.45, 0.14 and 0.03 K for nz = 10, 20 and 40, i.e., second-order convergence to the analytic solution, and the columns of a row agree to output precision (~1e-4 K).
- Transient tests (A4b), 50 kyr from the stationary P1 profile: IceColumnSolutions.jl used opposite sign conventions for Pe in its stationary and transient solutions (inherited from Moreno-Parada et al., 2024), and the wrong sign of the heat source for Pe ≠ 0; both are fixed in fesmc/IceColumnSolutions.jl#7, with finite-difference tests. With the corrected reference, SMB × 2 gives a maximum error of 0.44 K at nz = 10 at all times. The step T_srf − 10 K gives 0.3–2.3 K after 1 kyr and less than 0.9 K from 10 kyr on, dominated by the surface boundary layer, which the vertical grid resolves poorly at early times; the errors converge with nz (row maxima at 5 kyr of 1.4, 0.65 and 0.25 K for nz = 10, 20 and 40). The step itself cannot be compared at t = 0, and the 10-mode series is not converged at 1 kyr for the highest Péclet number (Pe = 33).

### Model changes

- `ytherm.strain_heating = "full" | "sia" | "none"` replaces `ytherm.use_strain_sia` (done). A4 needs strain heating switched off.
- `ydyn.solver = "fixed"` keeps the velocity prescribed by the driver, and the vertical velocity is still computed from continuity (exists).

## Open items

1. **Ensembles with overrides-only parameter files.** `runme -p` can only change parameters that appear in the parameter file, so a parameter taken from the defaults must be listed explicitly before it can be varied (a runme extension is being explored).
2. **Output precision.** The single-precision output limits the symmetry check to ~1e-7. A round-off check needs double-precision output of the checked fields.
3. **Flux-consistent vertical velocity.** A vertical velocity computed from the discrete layer fluxes of the thickness equation would be exact in A4 and remove the mismatch at grounding lines and margins (design: [`docs/dev/uz-flux-consistent.md`](../uz-flux-consistent.md)).
4. **ISLAND4 spin-up.** The 5-kyr tuning runs are not in equilibrium (volume still falls by about 3 % per kyr from the Vialov start), and the bed is temperate only in the troughs. The forcing must be confirmed by the first C0 run.
5. **Pass thresholds.** To be set after the first round of runs. Single precision with tight tolerances stays at about 1e-6 over 1 kyr on ISLAND4 at 32 km.

## Scripts and figures

| File | Purpose |
|---|---|
| `scripts/island4_map.jl` | Figure 1: ISLAND4 geometry and forcing profiles |
| `scripts/scorecard.jl` | Metrics and pass/fail summary (to be written) |

The scripts use the Julia environment in `scripts/`:

```bash
julia --project=docs/dev/benchmark-protocol/scripts -e 'import Pkg; Pkg.instantiate()'
julia --project=docs/dev/benchmark-protocol/scripts docs/dev/benchmark-protocol/scripts/island4_map.jl
```

## References

- Asay-Davis, X. S., et al. (2016). Experimental design for three interrelated marine ice sheet and ocean model intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1). *Geosci. Model Dev.*, 9, 2471–2497.
- Moreno-Parada, D., Robinson, A., Montoya, M., and Alvarez-Solas, J. (2024). Analytical solutions for the advective–diffusive ice column in the presence of strain heating. *The Cryosphere*, 18, 4215–4232.
- Robin, G. de Q. (1955). Ice movement and temperature distribution in glaciers and ice sheets. *J. Glaciol.*, 2, 523–532.
- Schoof, C. (2006). A variational approach to ice stream flow. *J. Fluid Mech.*, 556, 227–251.
- Schoof, C. (2007). Ice sheet grounding line dynamics: Steady states, stability, and hysteresis. *J. Geophys. Res.*, 112, F03S28.
- Tsai, V. C., Stewart, A. L., and Thompson, A. F. (2015). Marine ice-sheet profiles and stability under Coulomb basal conditions. *J. Glaciol.*, 61, 205–215.
- Weertman, J. (1957). Deformation of floating ice shelves. *J. Glaciol.*, 3, 38–42.
