# Kleiner et al. (2015) enthalpy benchmark data

Reference data for validating Yelmo's enthalpy thermodynamics solver, vendored
from [github.com/alex-robinson/icetemp](https://github.com/alex-robinson/icetemp)
(`data/Kleiner2015/`).

Reference: Kleiner, T., Rückamp, M., Bondzio, J. H., and Humbert, A. (2015),
*Enthalpy benchmark experiments for numerical ice sheet models*,
The Cryosphere, 9, 217–228, https://doi.org/10.5194/tc-9-217-2015.

## Files

- `Kleiner2015_EXPA_Fig2-IIIa-melt.txt` — Experiment A, basal melt rate `ab`
  vs `time` (Fig. 2, III a). Transient response of a cold slab to a parabolic
  surface-temperature forcing; tests the cold↔temperate basal transition.
- `Kleiner2015_EXPB_analytic_nz401_z.dat` — Experiment B, analytic steady-state
  column at nz=401: columns `z enth T_ice omega`. Steady polythermal column;
  tests CTS position and water content.
- `Kleiner2015_EXPB_extra-series_enth_b.txt` — Experiment B, basal enthalpy
  `enth` vs `time` (header notes `T_ref = 173.15 K`).

Used by the standalone column test driver `tests/test_icetemp.f90` (tests T3/T4;
see `docs/physics/enthalpy-transition-plan.md`).
