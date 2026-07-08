# Kleiner et al. (2015) enthalpy benchmark data

Reference data for validating Yelmo's enthalpy thermodynamics solver, vendored
from [github.com/alex-robinson/icetemp](https://github.com/alex-robinson/icetemp)
(`data/Kleiner2015/`).

Reference: Kleiner, T., Rückamp, M., Bondzio, J. H., and Humbert, A. (2015),
*Enthalpy benchmark experiments for numerical ice sheet models*,
The Cryosphere, 9, 217–228, https://doi.org/10.5194/tc-9-217-2015.

## Experiment A protocol (Table A1)

Parallel-sided slab, H = 1000 m, no flow, no strain heating. Geothermal flux
q_geo = 0.042 W m⁻². Surface-temperature phases: I (0–100 ka) −30 °C, II
(100–150 ka) **−10 °C**, III (150–300 ka) −30 °C. Basal water accumulates freely
(no drainage, no cap). Constants: ρ_i = 910, c_i = 2009, k_i = 2.1, L = 3.34×10⁵,
β = 7.9×10⁻⁸ K Pa⁻¹.

## Files

- `Kleiner2015_EXPA_analytic_ab_Tw-10.txt` — Experiment A, **analytic** basal
  melt rate a_b(t) for the published −10 °C protocol, from the paper's Eq. A13–A15
  (Fourier series, n = 25). `a_b > 0` = melt, `< 0` = freeze-on. Steady values:
  warm = 2.334 mm a⁻¹, cold = −2.027 mm a⁻¹. Generated from the paper's own
  analytical solution — this is the authoritative reference.

  Note: the previously-vendored `Kleiner2015_EXPA_Fig2-IIIa-melt.txt` (~3.10 mm a⁻¹
  warm plateau) corresponded to a **−5 °C** warm phase, not the published Exp A,
  and was removed to avoid confusion.
- `Kleiner2015_EXPB_analytic_nz401_z.dat` — Experiment B, analytic steady-state
  column at nz=401: columns `z enth T_ice omega`. Steady polythermal column;
  tests CTS position and water content.
- `Kleiner2015_EXPB_extra-series_enth_b.txt` — Experiment B, basal enthalpy
  `enth` vs `time` (header notes `T_ref = 173.15 K`).

Used by the standalone column test driver `tests/test_enthalpy.f90` (tests T3/T4;
see `docs/physics/enthalpy-transition-plan.md`).
