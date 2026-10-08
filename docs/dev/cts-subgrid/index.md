# Sub-grid CTS in the enthalpy column solver

Status: branch `cts-subgrid-30c974` (from dev e5ef9fc0), 2026-10-08. Column tests done;
3D checks (EISMINT EXPF symmetry, initmip GRL 16 km 1 kyr) pending on albedo. Not merged.

Commits:

1. 094a8c35 `therm: enth vertical diffusion with split sensible/latent face flux`
2. 62567bdc `therm: H_cts from temperate-side water-content extrapolation`
3. 62a46423 `tests: kleiner-b checks a steady CTS, tighter tolerance`

## 1. Problem

`test_enthalpy.x kleiner-b enth 201 1e-4` (Kleiner et al., 2015, Exp B) failed on albedo
with a final CTS height of 21.78 m (analytic 19.0 m, tolerance 2 m), and passed locally with 20.80 m.
The CTS was not steady. It jumped between three nodes every time step (period 3 steps =
6 yr at dt = 2 yr; Fig. 1a). The 500-yr output aliased this into an apparent slow
cycle, so the final snapshot (and PASS/FAIL) depended on the phase. The mean CTS was
about 1.8 m too high at nz = 201. The minmod vertical advection (2961c361) only changed
the phase.

## 2. Cause

The face diffusivity was set per node (cold `Kc` or temperate `K0 = cr*Kc`) and applied to
the full enthalpy difference. At the CTS face the cold diffusivity was used
(`kappa_b = kappa(k_cts+1)`). The enthalpy jump across that face includes the water content
`omega*L` of the temperate node, so cold-ice conduction carried latent heat into the cold
ice. In one step the top temperate node froze, then the next one down did the same, and
then the cold ice warmed back up. The cold side ended up flat at about -21 J/kg below E_pmp instead of
reaching E_pmp with zero gradient, as it does at a melting CTS (Fig. 1b).

## 3. Scheme

The flux on each face is split into a sensible and a latent part:

    F = Kc_f * d(E_s)/dz + K0_f * d(E_l)/dz,   E_s = min(E, E_pmp),   E_l = max(E - E_pmp, 0)

with `Kc_f` the harmonic mean of `kt/(rho*cp)` (whatever the phase) and `K0_f = cr*Kc_f`.
Cold ice conducts its temperature, temperate ice diffuses its water content, and across
the CTS cold ice conducts to the pressure melting point of the temperate node. The flux
is continuous in E, so the CTS settles between nodes.

- Linear system: each node's phase is taken from the start-of-step enthalpy (Dirichlet
  base and surface: from their imposed values). A node then contributes `c*x + d` to a
  face flux (`calc_split_coeffs`; `x = E - enth_ref`): cold `c = Kc`, `d = 0`;
  temperate `c = K0`, `d = (Kc - K0)*(E_pmp - enth_ref)`. Still one tridiagonal solve per step.
- Picard iteration on the phase was tested and not adopted. One solve already gives a
  consistent phase in 99.8 % of the steps, and at nz = 801 / cr = 1e-5 the iteration cycled
  without changing the result.
- Removed: `calc_enth_diffusivity`, the CTS-face override, the melting-base override
  (`k_cts == 1`; it was also applied to a cold base, because `get_cts_index` returns 1
  there), the floating-base `kappa_a = kappa(1)` override and `kappa_aa(1)` cold override.
  The split covers all of these: a melting base (Dirichlet at E_pmp) gets `Kc*(E_2 - E_pmp)`;
  a floating base at T_shlf below temperate ice gets `Kc*(E_pmp,2 - E_1)` plus the small
  latent part, without the latent leak of the old override.
- A cold column gives the same results as before, bit for bit (the products are ordered
  as before: `wp` is single precision).

`H_cts` (output only) used linear interpolation of `E - E_pmp` between the top temperate
node and the cold node above. Because the cold side approaches E_pmp with zero gradient,
this puts the CTS up to one layer too high. It now uses the lower of that and the
extrapolation of the temperate water content to zero (nodes `k_cts-1`, `k_cts`, not the
base node). A cold base gives 0. Previously it could give a spurious height there.

## 4. Results

![Kleiner Exp B, old vs new: (a) CTS every step, (b) profile near the CTS, (c) CTS vs nz, (d) basal water content vs nz](figures/cts_kleiner_b.png)

Exp B, cr = 1e-4, dt = 2 yr, 50 kyr; H_cts over the last 1 kyr, every step:

| nz | old H_cts mean [min, max] | new H_cts | new, linear interp. | old base omega | new base omega |
|---|---|---|---|---|---|
| 51  | 25.28 [23.66, 27.49] | 20.79 | 23.61 | 0.0180 | 0.0181 |
| 101 | 22.86 [19.99, 26.20] | 19.42 | 20.00 | 0.0194 | 0.0194 |
| 201 | 20.87 [19.93, 22.28] | 19.38 | 19.99 | 0.0200 | 0.0200 |
| 401 | 20.04 [18.99, 21.24] | 19.18 | 19.50 | 0.0204 | 0.0204 |
| 801 | 19.37 [18.00, 21.01] | 19.02 | 19.24 | 0.0206 | 0.0205 |

Analytic: CTS 19.0 m, base omega 0.0207. The new CTS is constant over the last 1 kyr at
all nz and for cr = 1e-3, 1e-4, 1e-5. Between cr = 1e-4 and 1e-5 the CTS changes by
< 0.01 m; cr = 1e-3 puts it up to 0.6 m higher (nz = 101). Basal water content is unchanged.

Other column tests (`test_enthalpy.x`, default nz = 51 unless noted):

| Test | Before | After |
|---|---|---|
| cold-limit, max abs(T_temp - T_enth) | 1.22e-4 K | 1.22e-4 K (identical) |
| kleiner-a, warm / cold a_b [mm/a] | 2.292 / -1.984 | 2.292 / -1.984 |
| kleiner-a-cap | 2.292 / -1.984 | 2.292 / -1.984 |
| thin-margin | all ok | all ok |
| robin-column, enth-A1 - Robin, Q_geo = 50 | -0.484 K | -0.519 K |
| robin-column, enth-A2 - temp, Q_geo = 50 | 0.417 K | 0.357 K |

The Robin changes come from the lowest face: the old code used `kappa(2)` there for
any base (the melting-base override), whereas the new code uses the harmonic mean of
nodes 1 and 2, like every other face.

The kleiner-b check is now: CTS within 0.5 m, base omega within 5 %, and an H_cts range
below 5 cm over the last 1 kyr (sampled every step). It passes for nz >= 201. At nz = 101 base omega
is 6 % low.

## 5. Expected 3D impact (to check)

Changes only where a CTS or a temperate node next to cold ice exists:

- No latent heat is conducted into cold ice at the CTS, so expect a lower CTS and less warming of
  the cold ice just above it. More water stays in the temperate layer and drains to the
  bed through `omega_max` (`melt_int`).
- Floating bases below temperate ice refreeze only through the sensible flux, so shelf-base
  refreezing of temperate ice slows down.
- Cold columns: identical.

Checks (albedo): EISMINT EXPF symmetry (20-kyr mean gate) and initmip GRL 16 km 1 kyr,
compared with dev.

## Reproduce

The figure data came from a copy of `tests/test_enthalpy.f90` that writes every step after
49 ka (`kleiner-b enth <nz> 1e-4`), built against dev e5ef9fc0 (old) and this branch (new),
output saved as `<dir>/{old,new}/kb_nz<nz>.nc`:

```bash
julia --project=docs/dev/cts-subgrid/scripts docs/dev/cts-subgrid/scripts/cts_plots.jl <dir> docs/dev/cts-subgrid/figures/cts_kleiner_b.png
```
