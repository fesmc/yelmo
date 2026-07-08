# Yelmo dev-branch bug-fix plan (2026-07)

Derived from the 4-agent audit of dyn / tpo / mat / timestepping. Thermodynamics
excluded (handled on `therm-dev` / enthalpy transition).

## Decisions taken

- **dzsdt (audit #1): IN SCOPE — important.** `dzsdt` is computed at
  `yelmo_topography.f90:409` from a `z_srf` not yet updated for the new `H_ice`
  (refreshed only later in `calc_ytopo_diagnostic`, L418), so it is ~0 where
  `calc_ydyn`→`calc_uz_3D*` consume it, corrupting vertical velocity every step.
  Fix: refresh `z_srf` from current `H_ice` before computing rates, in the
  predictor and advance paths. **Coordinate with therm-dev** (vertical-velocity /
  advection rework) to avoid a double-patch.
- **L1L2 solver (audit #5): DELETE ENTIRELY** — obsolete, broken, superseded by
  DIVA. Removes the membrane-stress bug and the L1L2 periodic-BC bug outright.
- **mask_adv → retire, use `bnd%mask_ice` directly.** `mask_adv` is hardwired to
  `1` everywhere (`yelmo_ice.f90:915`) and carries no information. The advection
  solver's internal convention `{-1=fixed, 0=none, 1=dynamic}` is migrated to the
  `MASK_ICE_*` convention `{NONE=0, FIXED=1, DYNAMIC=2}` (yelmo_defs.f90:50-52),
  `bnd%mask_ice` is threaded through `calc_G_advec_simple`/rk/lsf, and the
  `mask_adv` field is removed (type, alloc/dealloc, io output, c_api). This
  subsumes the `mask==-1` `b_value=H` fix. **Must be atomic** (see WS-GLUE).
- **T_prime rate-factor clamp (audit #6): OPEN — see §Open q.1.** Recommend ceiling
  `T0` (or remove); currently `min(T_prime,T_pmp)`. Not scheduled until decided.
- All other verified findings: proceed.

## Parallelization strategy

Each workstream owns a **disjoint set of files** → parallel agents in separate
worktrees, no merge conflicts. Hot "glue" files + the mask-migration + PC-history
+ dzsdt are consolidated into a **single Wave-2 integration workstream** (they are
interdependent and share files), run after the Wave-1 leaf workstreams merge.

Worktree discipline (CLAUDE.md): `git worktree add .claude/worktrees/<ws> dev`,
build-verify there, commit, merge to `dev` sequentially. No merges while batch
experiments run on `dev`.

### Wave 1 — leaf files, fully parallel (no shared files)

| WS | Files owned | Fixes |
|----|-------------|-------|
| **WS-SIA** | `physics/velocity_sia.f90` | #3 remove `k` from OMP `private` (L178); #4 `tau_*z_n_dn` first index → `k-1` (L187,191) |
| **WS-BDRAG** | `physics/basal_dragging.f90` | #7 power-plastic uses `cbn` not `c_bed(i,j)` (L986); #12 `infinite/mask` BC → zero-gradient in x (L507-508); scale_beta_gl_fraction x-loop `1,nx-1`→`1,nx` (L1147); gl_flux_weight grounded-end test (L1699-1718) |
| **WS-DEFORM** | `physics/deformation.f90` | #10 dyy margin check `f_ice(i,jm2)` + add missing `else` (L756,764-772); #11 `calc_visc_int` periodic writes col `nx` not `nx-1` (L364-372); strain-reset for ice-free (low) |
| **WS-TOOLS** | `yelmo_tools.f90` | #9 `stagger_aa_ab_ice` invert to iced-corner avg (L428) **+ re-verify `stagger_ab_aa_ice` L520 self-composition still correct after the flip**; #15 periodic `set_boundaries_2D_aa` single convention (L1198-1208); #16 `smooth_gauss_2D` reverse right/top halo + fill corners (L1670-1679) |
| **WS-VGEN** | `physics/velocity_general.f90` | #8 Picard `norm_method=1`: add `sqrt`, accumulate x/y under own masks (L1977-2002); `taud_gl_method=-1` no-op (decide: honor `beta_gl_stag` or `error stop`); remove `if(j==6) write` debug (L1287,1351) |
| **WS-MASSCALV** | `physics/mass_conservation.f90`, `physics/calving/calving_ac.f90`, `physics/grounding_line_flux.f90`, `physics/topography.f90` | #13 `infinite` border `H(nx,:)=H(nx-1,:)` (mass_cons L712); #14 uninitialized `BC` in 3 live LSF routines + 2 unused (calving_ac); #17 `pi`→`phi` typo (topography L2430); `qq_gl_acy` empty loop (implement mirror or `error stop`); `calv_now` uninit stub (delete or port); mb-limiter diagnostic (mass_cons L168) |
| **WS-TIME** | `yelmo_timesteps.f90`, `physics/runge_kutta.f90` | floor-rounding can zero dt + int overflow (L794-799); `pc_k=1` for FE-SBE + docstring (L432); rk23 advance with 2nd-order soln + labels (L106-125); rk4 tau `d2` cumulative + drop nudge (L279-300) |

*Note:* `physics/solver_advection.f90` is intentionally **not** in Wave 1 — its BC
fixes (#18 periodic-x, #19 expl) are bundled into WS-GLUE with the mask migration
so the file has a single owner and the convention change lands atomically.

### Wave 2 — integration (single agent, hot/glue files) — after Wave 1 merges

| WS | Files owned | Fixes |
|----|-------------|-------|
| **WS-GLUE** | `yelmo_dynamics.f90`, `yelmo_ice.f90`, `yelmo_topography.f90`, `yelmo_material.f90`, `yelmo_defs.f90`, `physics/solver_advection.f90`, `yelmo_io.f90`, `yelmo_c_api.f90`, `config/Makefile_yelmo.mk`, **delete** `physics/velocity_l1l2.f90` | See below |

WS-GLUE tasks (one coherent change set, keep dispatch + timestepping consistent):

1. **L1L2 deletion**: rm `physics/velocity_l1l2.f90`; rm `use`/dispatch (L195-198)/
   `calc_ydyn_l1l2` (L562-648)/option-string (L1001) in `yelmo_dynamics.f90`; rm 3
   lines in `config/Makefile_yelmo.mk` (120,156,250).
2. **mask migration (atomic)**: in `solver_advection.f90`, replace the
   `{-1,0,1}` mask logic with `{MASK_ICE_NONE=0 → zero row, MASK_ICE_FIXED=1 →
   b_value=H row, MASK_ICE_DYNAMIC=2 → inner discretization}` across **all** solver
   variants; thread `bnd%mask_ice` through `calc_G_advec_simple`/`rk*_2D_step`/
   `LSFupdate` in place of `tpo%now%mask_adv`; delete the `mask_adv` field from
   the ytopo type + alloc/dealloc/init, `yelmo_io` output, `yelmo_c_api`. Also fold
   in #18 (impl-lis `periodic-x` matrix BC) and #19 (`expl` periodic
   mass-conservation: pass BC into the ab-node stagger). Subsumes the `mask==-1`
   `b_value` fix.
3. **dzsdt (#1)**: in `yelmo_topography.f90`, refresh `z_srf` from current `H_ice`
   before computing `dzsdt` (L409) in the predictor + advance paths.
4. **PC history (#2)**: store raw predictor `dHidt_now` as a new field; corrector
   uses `β₃·f*_{n+1} + β₄·f_n`; shift raw `f_n→f_{n-1}` at advance
   (`yelmo_topography.f90` L103,123,150,340-403).
5. **Bootstrap-tau consistency**: branch the tau call on `pc_active` like the beta
   call (`yelmo_ice.f90` L244-333).
6. **Deallocs / init**: `pc%dmb` (`yelmo_topography.f90` L1791); `strn2D%f_shear`
   init + dealloc (`yelmo_material.f90` L401,445,504).
7. **mask_grz / dist_margin (#Tier2)**: compute distances (km) before zone
   classification; compute `dist_margin` (`yelmo_topography.f90` L1061).

## Merge order & gating

1. Wave-1 branches build-verify independently, merge to `dev` in any order
   (disjoint files ⇒ no conflicts). Build after each merge.
2. WS-GLUE branches off `dev` after all Wave-1 merges, rebuilds, merges last.
3. Between merges: no running experiments on `dev`.

## Verification

- **Per fix**: clean build (`make yelmo-static` or host target).
- **Answer-preserving** (dead code, deallocs, diagnostics, guards, debug removal):
  expect bit-identical results on a reference run.
- **Answer-changing** (Picard `sqrt`, SIA `k-1`, boundary fixes, dzsdt, PC history
  #2, mask migration): each needs a before/after benchmark. **dzsdt** and **PC
  history #2** are the highest-risk items — each its own commit + EISMINT/MISMIP
  diff review.
- **Transpose-symmetry regression test** (proposed): run e.g. EISMINT-moving on a
  transposed grid, assert `H(x,y)==H(y,x)ᵀ` to tolerance. Would have caught #10,
  #11, #12, #13 automatically and protects the v2.0 refactors.

## Open questions

1. **T_prime ceiling** (`deformation.f90:415`). `T′ = T_ice − T_pmp + T0` is the
   Greve-Blatter homologous temperature in **Kelvin**; the `+T0` is the homologous
   shift, not a °C→K conversion. The physical cap `T_ice ≤ T_pmp` maps to
   `T′ ≤ T0` (not `T′ ≤ T_pmp`), so `min(T_prime,T_pmp)` cuts ATT by ~1.6× for
   near-temperate ice. Recommend ceiling `T0` (or remove; the EISMINT variant has
   no clamp). **Decision needed before scheduling.**
2. **dzsdt vs therm-dev**: confirm this fix lands here and is coordinated with the
   vertical-velocity work on therm-dev (no double-patch).
3. **mask_adv retirement**: OK to drop the `mask_adv` output variable + c_api
   entry (restart/output configs referencing it will lose it)?
