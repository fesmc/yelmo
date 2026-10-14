# Enthalpy vs temperature solver: residual discretization difference

Follow-up investigation, deferred to a dedicated session. Companion to
[enthalpy-transition-plan.md](enthalpy-transition-plan.md).

## Status

After the 2026-07 enth fixes (cold-floor clamp, `cp_ref`/A1, integral/A2,
`dzsdt`→`uz`), the enthalpy solver reproduces the temperature solver on
initmip-grl-16km to **~0.46 K RMS** basal temperature. The remaining
difference is a genuine discretization inconsistency, not a bug that a small
patch removes. This note records what is established so a future session does
not repeat the dead ends.

## What is established

1. **It is a real inconsistency, not truncation.** Refining the vertical grid
   (`nz` 51 → 201 → 501) leaves the temp↔enth basal-temperature gap essentially
   fixed (~0.69 K at `Q_geo = 70 mW m⁻²` in the `robin-column` test). The two
   solvers **converge to different limits**.

2. **It appears only with variable material properties.** With *constant*
   `cp` and `kt`, temp ≡ enth ≡ Robin analytic to 0.008 K. The divergence
   switches on only when `cp(T) = 146.3 + 7.253 T` and `kt(T)` vary through the
   column. It grows with the basal temperature / vertical gradient (i.e. with
   the amount of property variation).

3. **It is not the harmonic-mean interface indexing.** `calc_temp_column_internal`
   evaluates the interface conductivity at `zeta_ac(k)`/`zeta_ac(k+1)` (consistent
   with `dzeta_a(k)`/`dzeta_b(k)` from `calc_dzeta_terms`), while
   `calc_enth_column_internal` uses `zeta_ac(k-1)`/`zeta_ac(k)`. Aligning enth to
   temp's indexing was **tested and rejected**: it did not close the gap and it
   worsened the Kleiner-B CTS (20.8 → 23.2 m vs analytic 19.0 m). The original
   enth indexing is the empirically-validated one; reverted, not committed.

4. **Root cause is the T-scheme vs E-scheme difference.** temp discretizes the
   temperature equation (diffuse `T` with `kt`); enth discretizes the enthalpy
   equation (diffuse `E = ∫cp dT` with `κ = kt/(ρ cp)`). These are discretely
   equivalent **only** when `cp` is constant. With variable `cp`, the nonlinear
   `E(T)` combined with the harmonic mean of the compound `kt/cp` coefficient
   makes the two discrete operators genuinely different, so they reach different
   steady states. A2 (integral enthalpy) removes the constant-`cp` *definition*
   error but not this operator-level difference.

5. **"Which is more correct" is currently undetermined.** The Robin analytic and
   the Kleiner benchmarks all assume constant `cp`, so none of the existing
   references can referee the variable-property case where temp and enth disagree.

## Can IceColumnSolutions.jl referee it?

<https://github.com/fesmc/IceColumnSolutions.jl> — analytic solutions for the 1D
advection–diffusion column (Moreno-Parada et al., 2024, Appendix A):
`solve_stationary` (steady θ(ζ)) and `solve(par, ts)` (transient θ(ξ,τ) via an
eigenfunction/eigenvalue expansion), parameterized by dimensionless `Pe`
(vertical advection), `Br` (strain heating), `Λ` (horizontal advection), `γ`/`β'`
(geothermal / surface insulation).

- **Helps a lot** for the regimes where temp and enth *should* agree. It is far
  broader than the hand-rolled steady Robin used so far: transient, strong
  advection (`Pe`), strain heating (`Br`), horizontal advection (`Λ`). Use it to
  (a) verify both Yelmo solvers reproduce the exact solution in the
  constant-property limit across this parameter space, and (b) measure each
  scheme's order of accuracy / dispersion. A scheme that is demonstrably more
  accurate there is the more trustworthy one.
- **Not sufficient on its own.** The eigenfunction method needs a *linear,
  constant-coefficient* equation, so it assumes constant `cp`/`kt`. It therefore
  cannot directly referee the variable-`cp(T)`/`kt(T)` divergence — which is
  exactly the disagreement. For that, a separate **high-resolution numerical
  reference** (an independent fine solver, or MMS with a manufactured
  variable-property solution) is needed; IceColumnSolutions.jl can validate that
  reference in the constant-property limit.

## Proposed plan for the new session

1. Wire IceColumnSolutions.jl into the standalone harness (or a Julia driver) and
   validate `calc_temp_column` and `calc_enth_column` against `solve`/`solve_stationary`
   across (Pe, Br, Λ, transient) at constant properties — confirm both are
   correct and quantify their convergence order.
2. Build a variable-property referee: either MMS (manufactured solution with
   `cp(T)`, `kt(T)` — add the analytic source term) or a heavily-refined
   independent numerical solution. Determine which of temp / enth converges to it.
3. If one scheme is demonstrably wrong for variable properties, fix that one; if
   both are "reasonable" discretizations of different-but-valid forms, decide
   whether to unify them onto a single shared discrete operator (so temp and enth
   agree by construction) — noting the cost of revalidating the production temp
   solver against EISMINT.
4. Re-run initmip-grl temp-vs-enth to confirm the gap closes.

## Landing state (2026-07-08, branch `therm-dev`)

Committed and stable: cold-floor clamp, A1 (`cp_ref`), A2 (`enth_cp_method`,
default `integral`), `dzsdt`→`uz` fix, `analysis/compare_grl.jl`. enth ≈ temp to
0.46 K RMS on initmip-grl. This residual is the subject of the investigation
above and does **not** block using the enth solver.

## Session update (2026-07-08): IceColumnSolutions.jl bench (plan step 1)

Step 1 of the plan is implemented and the constant-property regime is validated.

- **Harness.** The standalone column driver `tests/test_icetemp.f90` was
  resurrected against the current API (constants moved into `ybound_const_class`;
  new `calc_temp_column`/`calc_enth_column`/`calc_T_pmp`/bedrock signatures;
  dropped the removed `calc_temp_robin_column` and active-bedrock
  `calc_temp_column_bedrock`). A new `"ics"` experiment was added: constant
  cp/kt, linear vertical velocity, no strain heating, cold column, geothermal
  flux prescribed at the base, uniform ζ grid. It is driven via CLI args
  (`nz cr solver experiment [smb] [Qgeo]`).
- **Referee.** `analysis/compare_column_analytic.jl` orchestrates the harness and
  compares the steady T(ζ) to `IceColumnSolutions.solve_stationary`, over an nz
  refinement sweep and a (Pe, γ) grid. Env in `analysis/Project.toml`
  (devs IceColumnSolutions).
- **Result.** Both solvers reproduce the exact advection–diffusion solution to
  **0.003–0.06 K RMS** across Pe ∈ [0, 22], γ from Q_geo ∈ {15,30,45} mW m⁻².
  **temp ≈ enth everywhere** (differences ~1–2 %), confirming they agree at
  constant properties against a far broader exact solution than the hand-rolled
  Robin. In the pure-diffusion limit (Pe=0) enth is ~2× less accurate than temp
  (e.g. 6.0e-3 vs 2.8e-3 K at Q_geo=30) — a faint but consistent signal.
- **Convergence order could NOT be measured.** Yelmo is built single precision
  (`wp = sp`); the error floors at the sp round-off level (~0.05–0.1 K), so RMS is
  flat vs nz (order ≈ 0.0–0.3) across the whole sweep. Measuring the true order
  needs a double-precision toolchain.
- **Why dp is not a quick flip.** Setting `wp = dp` cascades into the *external*
  `fesm-utils` library (shared by all models): its public routines are locked to
  `wp = sp`, so a dp Yelmo hits interface type mismatches. `fesm-utils/src/subgrid.f90`
  was made precision-generic (generic interface → `_dp` worker + `_sp` wrapper;
  committed/pushed on `fesm-utils:dev`), but `gaussian_quadrature`, `derivatives`,
  `distances`, … need the same treatment (dozens of routines). The cheaper route
  for a future dp bench is to build a separate dp `fesm-utils` variant and link
  Yelmo's dp build against it, rather than genericizing every routine.
- **IceColumnSolutions gotcha (filed).** Its `w0` sign convention is inverted vs
  its docstring: accumulation is `Pe > 0` (real-`erf` branch); the documented
  "negative = downward" lands on an ill-conditioned imaginary-`erf` path and
  gives unphysical results. The referee therefore maps `w0 = +smb`. See
  fesmc/IceColumnSolutions.jl#6 (proposes adopting negative = downward, matching
  Yelmo, and fixing the sign so downward is the well-conditioned branch).

**Still open:** step 2 (variable-property MMS / fine-solver referee) — the actual
`cp(T)` divergence — is untouched and remains the crux. Step 1's analytic bench
can only validate that reference in the constant-property limit.
