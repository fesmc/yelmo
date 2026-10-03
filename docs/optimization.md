# Basal friction optimization

Yelmo can adjust the basal friction coefficient `cb_ref` during a simulation so
that the simulated ice thickness approaches an observed thickness. The method
follows Lipscomb et al. (2021) and is implemented in `optimize_cb_ref`
([`libs/ice_optimization.f90`](https://github.com/fesmc/yelmo/blob/main/libs/ice_optimization.f90)).
It is used by two programs:

- `yelmo_initmip.x` with `ctrl.equil_method = "opt"` and the `&opt` group of
  `par/yelmo_initmip.nml`. This is the standard route; see the initmip
  [benchmarks](benchmarks.md).
- `yelmo_opt.x` (`tests/yelmo_opt.f90`, `make opt`, runme alias `opt`), a
  stand-alone optimization program with `ctrl.opt_method = "L21"` and the
  `&opt_L21` group. No parameter file for it is included in `par/`.

## Update of `cb_ref`

At every outer time step `dt`, `cb_ref` is updated in each cell where the
simulated speed is non-zero and the observed ice is grounded:

$$
\frac{\mathrm d c}{\mathrm d t} = -\frac{c}{H_0}\left(\frac{H - H_\mathrm{obs}}{\tau_c} + 2\,\frac{\partial H}{\partial t} + \frac{f_\mathrm{tgt}}{\tau_c}\ln\frac{c}{c_\mathrm{tgt}}\right),
$$

with $c$ = `cb_ref`. The thickness error and the thickness change are taken
from the upstream neighbours, weighted by the flow direction. The second term
damps oscillations. The third term, with $f_\mathrm{tgt} = 0.05$ `H0`, pulls
`cb_ref` towards the target field `cb_tgt` (initmip only; zero for
`H0 = 0`). Where there is no observed velocity, the upstream thickness error is
the mean of the two upstream neighbours. The parameters are:

| Parameter | Meaning |
|---|---|
| `tau_c` | [yr] Time scale of the adjustment (default 500 yr) |
| `H0` | [m] Thickness scale of the adjustment. With `H0 = 0` (initmip default), $H_0 = \max(H_\mathrm{obs}, 50\ \mathrm m)$ in each cell. Negative values should not be used |
| `fill_method` | Value in cells where the observed ice is floating or absent: `"cf_min"` or `"target"` (`cb_tgt`). Grounded cells with zero simulated speed get `cf_min` |
| `cf_min`, `cf_max` | Limits of `cb_ref`. In initmip these are `ytill.cf_min` and `ytill.cf_ref` |

`sigma_err` and `sigma_vel` are read, but the smoothing of the thickness error
they control is switched off in `optimize_cb_ref`.

## Relaxation of the ice shelves

The optimization works best when the floating ice and the grounding zone are
first held close to the observed thickness and then gradually released. The
relaxation time scale `tau` (`tpo%par%topo_rel_tau`) is `rel_tau1` until
`rel_time1`, increases to `rel_tau2` at `rel_time2` as

$$
\tau = \tau_1 + (\tau_2 - \tau_1)\left(\frac{t - t_1}{t_2 - t_1}\right)^m,
$$

with $m$ = `rel_m`, and relaxation is switched off after `rel_time2`. initmip
relaxes the grounding zone (`ytopo.topo_rel = 4`), with times counted from
`ctrl.time_init`. `yelmo_opt.x` relaxes the floating ice and the grounding line
(`topo_rel = 2`), with `rel_time1` and `rel_time2` given as model times.

## initmip

In `&opt` (times are counted from `ctrl.time_init`, also after a restart):

- `opt_cf`: switch for the friction optimization, active between
  `cf_time_init` and `cf_time_end` (default 0–15 kyr);
- `cf_init`: initial `cb_ref` everywhere without a restart (must be positive);
- `tau_c`, `H0`, `fill_method`, `rel_tau1`, `rel_tau2`, `rel_time1`,
  `rel_time2`, `rel_m` as above.

The optimization of the ocean thermal forcing (`opt_tf`) is not implemented
and stops the model.

## yelmo_opt.x

`yelmo_opt.x` reads `&ctrl` (`opt_method`, `cb_ref_init_method`, `sigma_err`,
`sigma_vel`, `cf_min`, `cf_max`, `bmb_shlf_const`, `dT_ann`, `z_sl`) and
`&opt_L21` (`time_init`, `time_end`, `rel_tau1`, `rel_tau2`, `rel_time1`,
`rel_time2`, `tau_c`, `H0`), in addition to the Yelmo groups. It requires
`ytill.method = -1`, so that `cb_ref` is not recomputed by Yelmo. The outer time
step is 5 yr and `fill_method = "cf_min"`. The forcing is the present-day
surface temperature (plus `dT_ann`) and surface mass balance, with constant
`bmb_shlf`, sea level zero and $Q_\mathrm{geo}$ = 50 mW m$^{-2}$.

`cb_ref` is initialised according to `cb_ref_init_method`: `"guess"` (from the
driving stress and the observed velocity, `guess_cb_ref`), `"restart"` (from the
restart file) or `"none"` (0.2 everywhere). Without a restart, the program
first runs the model for 20 years with DIVA to smooth the input topography,
which then becomes the target thickness of the optimization, and then for
20 kyr with fixed topography to equilibrate the thermodynamics. This state is
written to `yelmo_restart_init.nc`. Setting `yelmo.restart` to this file skips
the spin-up in later runs; the target thickness is then the observed one. The optimization then runs from `time_init` to
`time_end`, and the final state is written to `yelmo_restart.nc`.

```bash
make opt
runme -r -e opt -o output/opt-test -n <par file>
```

## Reference

Lipscomb, W. H., Leguy, G. R., Jourdain, N. C., Asay-Davis, X., Seroussi, H.,
and Nowicki, S. (2021). ISMIP6-based projections of ocean-forced Antarctic ice
sheet evolution using the Community Ice Sheet Model. The Cryosphere, 15,
633–661. <https://doi.org/10.5194/tc-15-633-2021>
