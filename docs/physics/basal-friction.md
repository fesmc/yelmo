# Basal friction

The basal shear stress of grounded ice is written as

$$
\boldsymbol\tau_b = \beta\,\mathbf u_b,
$$

with the friction coefficient $\beta$ [Pa a m$^{-1}$] (`beta`) computed by
`calc_beta` in
[`src/physics/basal_dragging.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/basal_dragging.f90).
$\beta$ enters the momentum balance directly (SSA) or through the effective
drag $\beta_\mathrm{eff}$ (DIVA, see [DIVA](momentum/diva.md)). For floating
ice $\beta = 0$.

## Friction laws

$\beta$ is the product of the bed coefficient $c_b$ [Pa] (`c_bed`) and a
function of the basal speed. The law is selected with `ydyn.beta_method`:

| `beta_method` | Law | $\beta$ |
|---|---|---|
| -1 | Set externally | — |
| 0 | Constant | `beta_const` |
| 1 | Linear (**default**) | $c_b/u_0$ |
| 2 | Pseudo-plastic power law (Bueler and van Pelt, 2015) | $c_b\,(\lvert u_b\rvert/u_0)^q / \lvert u_b\rvert$ |
| 3 | Regularized Coulomb (Joughin et al., 2019) | $c_b\,\left(\lvert u_b\rvert/(\lvert u_b\rvert + u_0)\right)^q / \lvert u_b\rvert$ |
| 4, 5 | As 2 and 3, evaluated at the cell centre | |

with $q$ = `beta_q` (default 1) and $u_0$ = `beta_u0` (default 100 m a$^{-1}$).
The basal speed is regularized as $\sqrt{\lvert u_b\rvert^2 + u_\mathrm{min}^2}$, with
$u_\mathrm{min} = 10^{-3}$ m a$^{-1}$. For
methods 1–3, $\beta$ is averaged over Gaussian quadrature points of the cell, with
$c_b$ and $u_b$ interpolated to the points. With $q = 1$, method 2 is the linear
law. For $q \to 0$, methods 2 and 3 tend to a plastic bed with
$\lvert\tau_b\rvert \to c_b$; for method 3, $\lvert\tau_b\rvert \le c_b$ for
any speed.

After the friction law, $\beta$ is modified in this order:

1. divided by the sub-temperate sliding factor $f_\mathrm{slide}$ (below);
2. scaled near the grounding line (`beta_gl_scale`): 0, multiplied by
   `beta_gl_f` at the grounding line; 1, reduced linearly to zero as the
   thickness above flotation decreases from `H_grnd_lim` to 0; 2, scaled by
   the thickness above flotation (Zstar, Gladstone et al., 2017); 3, multiplied
   by the grounded fraction;
3. set to zero for floating ice and limited to at least `beta_min` (default
   100 Pa a m$^{-1}$) for grounded ice.

Steps 1–3 are skipped for `beta_method = -1`. $\beta$ is then staggered to the
velocity faces (`beta_gl_stag`): 0, mean of the two cells; 1 (default), the
upstream (grounded) value at grounding-line faces; 2, the downstream value;
3, a mix of the grounded and floating values weighted by the square of the
grounded fraction of the face; 4, weighted by the grounded
share of the velocity across the face (Gladstone et al., 2010).

### Frozen-bed sliding

With `ydyn.frz_scale = True` (default), sliding is reduced where the bed is
below the pressure melting point. After the friction law, $\beta$ is multiplied
by $f_\mathrm{slide}^{-q}$, with

$$
f_\mathrm{slide} = f_\mathrm{min} + (1 - f_\mathrm{min})\exp(T'_b/\gamma_T),
$$

where $T'_b \le 0$ is the basal homologous temperature, $\gamma_T$ is
`ydyn.frz_efold` and $f_\mathrm{min}$ is `ydyn.frz_min` (e.g., Fowler, 1986;
Hindmarsh and Le Meur, 2001). $q$ is `beta_q`, or 1 for `beta_method` 0 and 1.
At a given basal stress this scales the sliding speed by $f_\mathrm{slide}$ for
the linear and power-plastic laws, so $\gamma_T$ is the e-folding temperature of
the sliding speed for any friction law. For the regularized Coulomb law this
holds for $u_b \ll u_0$; for $u_b \gg u_0$ the yield stress becomes
$c_b f_\mathrm{slide}^{-q}$. $f_\mathrm{slide} = 1$ where the base is temperate,
not grounded, in contact with the ocean (partially floating, or next to floating
ice) or wet (`hyd_W` > 0 or `hyd_W_til` > 0). It is output as `f_slide`.

## Bed coefficient

The bed coefficient is

$$
c_b = c_{b,\mathrm{ref}}\,N_\mathrm{eff}
\qquad\text{or}\qquad
c_b = \tan\!\left(\frac{\pi}{180}\,c_{b,\mathrm{ref}}\right)N_\mathrm{eff},
$$

with the effective pressure $N_\mathrm{eff}$ (`N_eff`) and the till property
$c_{b,\mathrm{ref}}$ (`cb_ref`). With `ytill.is_angle = True`, `cb_ref` is a till
friction angle in degrees (Bueler and van Pelt, 2015), otherwise a
dimensionless factor. `cb_ref` represents the properties of the bed (roughness,
till) and is independent of $N_\mathrm{eff}$.

With `ytill.method = 1` (default), `cb_ref` is computed from the bedrock
elevation relative to present-day sea level (e.g., Winkelmann et al., 2011),
selected with `ytill.scale_zb`:

- 0: `cb_ref = cf_ref` everywhere;
- 1 (default): `cb_ref = max(cf_min, cf_ref·λ)`, with
  $\lambda = (z_b - z_0)/(z_1 - z_0)$ limited to $[0, 1]$;
- 2: `cb_ref = max(cf_min, cf_ref·λ)`, with
  $\lambda = \min(1, \exp((z_b - z_1)/(z_1 - z_0)))$, which is $e^{-1}$ at
  $z_b = z_0$.

Here $z_0$, $z_1$ are `ytill.z0`, `ytill.z1` (defaults −300 m, 200 m) and
`cf_ref`, `cf_min` default to 0.8 and 0.1. With `n_sd > 1` (default 10), the
scaling is averaged over Gaussian-weighted samples of $z_b \pm \sigma_{z_b}$
(`z_bed_sd`), to account for sub-grid bed roughness.

The sediment thickness can reduce `cb_ref` further (`ytill.scale_sed`), with a
linear ramp $\lambda_s$ from 0 at `sed_min` to 1 at `sed_max`:

- 0 (default): no sediment scaling;
- 1: `cb_ref = min(cb_ref, cf_min·λ_s + cf_ref·(1 − λ_s))`;
- 2: `cb_ref` multiplied by $1 - (1 - f_\mathrm{sed})\lambda_s$, with
  $f_\mathrm{sed}$ = `f_sed`, and limited to at least `cf_min`;
- 3: as 2, without the lower limit.

With `ytill.method = -1`, `cb_ref` is set externally, e.g. by the
[basal friction optimization](../optimization.md).

## Effective pressure

$N_\mathrm{eff}$ comes from the basal hydrology (FastHydrology, group `&yhyd`).
With the till-water bucket and no water transport (`yhyd.method_transport = 0`),
the closure is chosen with `yhyd.bkt_N_closure`:

| `bkt_N_closure` | $N_\mathrm{eff}$ |
|---|---|
| -1 | Set externally (host model, via `yelmo_set_var2D("hyd_N")`) |
| 0 | Constant, `yhyd.const_N` (grounded ice) |
| 1 | Overburden pressure, $\rho_i g H$ |
| 2 | Marine closure (Leguy et al., 2014), exponent `yhyd.marine_p` |
| 3 | Till closure (Bueler and van Pelt, 2015), parameters `yhyd.till_*` (**default**) |

With water transport (`method_transport = 1`), $N_\mathrm{eff}$ is computed by
the transport model. With `ydyn.neff_nxi` > 0, the cell value of
$N_\mathrm{eff}$ is the average of the hydrology's $N$ interpolated to sub-grid
points (1: the four Gauss points, > 1: an `neff_nxi` × `neff_nxi` grid; default 0,
no interpolation).

## References

- Bueler, E. and van Pelt, W. (2015). Mass-conserving subglacial hydrology in
  the Parallel Ice Sheet Model version 0.6. Geosci. Model Dev., 8, 1613–1635.
- Gladstone, R. M., Payne, A. J., and Cornford, S. L. (2010). Parameterising the
  grounding line in flow-line ice sheet models. The Cryosphere, 4, 605–619.
- Gladstone, R. M., et al. (2017). Marine ice sheet model performance depends on
  basal sliding physics and sub-shelf melting. The Cryosphere, 11, 319–329.
- Joughin, I., Smith, B. E., and Schoof, C. G. (2019). Regularized Coulomb
  friction laws for ice sheet sliding. Geophys. Res. Lett., 46, 4764–4771.
- Leguy, G. R., Asay-Davis, X. S., and Lipscomb, W. H. (2014). Parameterization of
  basal friction near grounding lines in a one-dimensional ice sheet model. The
  Cryosphere, 8, 1239–1259.
- Winkelmann, R., et al. (2011). The Potsdam Parallel Ice Sheet Model (PISM-PIK),
  Part 1: Model description. The Cryosphere, 5, 715–726.
