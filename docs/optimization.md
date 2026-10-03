# Basal friction optimization

Yelmo can adjust the basal friction coefficient `cb_ref` during a simulation so
that the simulated ice thickness approaches an observed thickness. The method
follows Lipscomb et al. (2021) and is implemented in `optimize_cb_ref`
([`libs/ice_optimization.f90`](https://github.com/fesmc/yelmo/blob/main/libs/ice_optimization.f90)).
It is a spin-up option of `yelmo_initmip.x`, with `ctrl.equil_method = "opt"`
and the `&opt` group of `par/yelmo_initmip.nml` (see the initmip
[benchmarks](benchmarks.md)). initmip sets `ytill.method = -1`, so that
`cb_ref` is not recomputed by Yelmo.

## Update of `cb_ref`

At every outer time step `dt`, `cb_ref` is updated in each cell where the
simulated speed is non-zero and the observed ice is grounded:

$$
\frac{\mathrm d c}{\mathrm d t} = -\frac{c}{H_0}\left(\frac{H - H_\mathrm{obs}}{\tau_c} + 2\,\frac{\partial H}{\partial t} + \frac{f_\mathrm{tgt}}{\tau_c}\ln\frac{c}{c_\mathrm{tgt}}\right),
$$

with $c$ = `cb_ref`. The thickness error and the thickness change are taken
from the upstream neighbours, weighted by the flow direction. The second term
damps oscillations. The third term, with $f_\mathrm{tgt} = 0.05$ `H0`, pulls
`cb_ref` towards the target field `cb_tgt` (zero for
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
relaxes the grounding zone (`ytopo.topo_rel = 4`).

## Parameters

In `&opt` (times are counted from `ctrl.time_init`, also after a restart):

- `opt_cf`: switch for the friction optimization, active between
  `cf_time_init` and `cf_time_end` (default 0–15 kyr);
- `cf_init`: initial `cb_ref` everywhere without a restart (must be positive);
- `tau_c`, `H0`, `fill_method`, `rel_tau1`, `rel_tau2`, `rel_time1`,
  `rel_time2`, `rel_m` as above.

The optimization of the ocean thermal forcing (`opt_tf`) is not implemented
and stops the model.

## Reference

Lipscomb, W. H., Leguy, G. R., Jourdain, N. C., Asay-Davis, X., Seroussi, H.,
and Nowicki, S. (2021). ISMIP6-based projections of ocean-forced Antarctic ice
sheet evolution using the Community Ice Sheet Model. The Cryosphere, 15,
633–661. <https://doi.org/10.5194/tc-15-633-2021>
