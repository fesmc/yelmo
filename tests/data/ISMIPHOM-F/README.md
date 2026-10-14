# ISMIP-HOM Experiment F reference results

Steady-state results of the models that took part in ISMIP-HOM Experiment F,
from the supplement (`tc-2007-0019-sp2.zip`, `ismip_all/`) of

Pattyn, F., Perichon, L., Aschwanden, A., et al. (2008), *Benchmark experiments
for higher-order and full-Stokes ice sheet models (ISMIP-HOM)*, The Cryosphere,
2, 95–108, https://doi.org/10.5194/tc-2-95-2008 (CC BY 3.0).

## Protocol (Sect. 3.6, Table 2)

Slab of mean thickness H0 = 1000 m on a mean slope of 3°, over a Gaussian bed
bump of amplitude 100 m and width σ = 10 km, periodic domain of L = 100 km,
zero surface mass balance, run to steady state. Linear rheology (n = 1) with
A = 2.140373e-7 Pa⁻¹ a⁻¹, so the unperturbed surface speed is 100 m/a.
F1: no slip (slip ratio c = 0); F2: c = 1, β² = 1/(c A H0).

## Files

`<model>f000.txt` (F1) and `<model>f001.txt` (F2) for 7 models: cma1 and oga1
(full Stokes), cma2, fpa1, fsa1, mbr1, mtk1. Columns:

| 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|
| x | y | z_s perturbation (m) | v_x (m/a) | v_y (m/a) | v_z (m/a) |

x and y are in km (= x/H0), relative to the bump centre, in [-50, 50], except
cma1 and cma2, which use normalised coordinates in [0, 1] (bump at 0.5). Grids
differ between models; mtk1 coordinates carry ~1e-3 km jitter and fpa1 has no
row at y = 0 (nearest rows y = ±1.28 km).

Used by `tests/ismiphom_f.jl` (comparison figure in `docs/benchmarks.md`).
