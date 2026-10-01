
# Parameters

Here important parameter choices pertinent to running
**Yelmo** will be documented. Each section will
outline a specific parameter or set of related parameters.
The author of each section and the date last updated
will apear in the heading, to maintain traceability
in the documentation (since code usually changes over time).

**This is a work in progress!**

::: {.callout-tip}
To browse, build, compare and validate parameter files from the command line, see
[`yelmo-config`](yelmo-config.md).
:::

## Basal friction

Yelmo includes the representation of several friction laws that all take the form:

$$
\tau_b = -\beta u_b
$$

where $\beta$ is composed of a coefficient $c_b$ and potentially another contribution that depends on $u_b$ too:

$$
\beta = c_b f(u_b)
$$

In Yelmo, the field $c_b$ is defined by the variable `c_bed` and has units of Pa. The term $f(u_b)$ is not output in the model, but it contributes with units of yr m$^{-1}$, so $\beta$ finally has units of Pa yr m$^{-1}$. When multiplied with $u_b$, we arrive at $\tau_b$ with units of Pa.

Yelmo calculates $c_b$ (`c_bed`) internally as either:

$$
c_b = c_{\rm b,ref} * N_{\rm eff}
$$

or

$$
c_b = {\rm tan}(c_{\rm b,ref}) * N_{\rm eff}
$$

This is controlled by the user option `ytill.is_angle`. If `ytill.is_angle=True`, then $c_{\rm b,ref}$ (variable `cb_ref` in the code) is considered as an angle and the latter formulation above is used, following e.g., Bueler and van Pelt (2015). If `ytill.is_angle=False`, then `cb_ref` is used as a scalar field directly. In both cases, this field represents the till or basal properties (roughness, etc.) that are rather independent from how the effective pressure $N_{\rm eff}$ (variable `N_eff`) may be defined.

With the variables formulated as above, it is possible to consider `cb_ref` as a tunable field that can be adjusted to improve model performance on a given domain. This can be achieved, for example, by performing optimization via the `ice_optimization` module, which adjusts `cb_ref` as a function of the mismatch of the simulated ice thickness with a target field. Also, `cb_ref` can either be optimized as a scalar field itself, or as an angle that is input to ${\rm tan}(c_{\rm b,ref})$ above.

Another possibility is to tune `cb_ref` as a function of other model or boundary variables. The most common approach is to tune it as as function of the bedrock elevation relative to present-day sea level (e.g., Winkelmann et al., 2011). In Yelmo, this is controlled by the parameter choices in the `ytill` section, and in particular the parameter `ytill.scale_zb` (0: none, 1: linear, 2: exponential). When `ytill.scale_zb=0`, no scaling function is applied and then `cb_ref=ytill.cf_ref` everywhere. When `ytill.scale_zb=1`, a linear scaling is applied so that `cb_ref` goes from `ytill.cf_min` to `ytill.cf_ref` for bedrock elevations between `ytill.z0` and `ytill.z1` (saturating otherwise). Finally, if `ytill.scale_zb=2`, an exponential decay function is applied, such that `cb_ref=ytill.cf_ref` for `z_bed >= ytill.z1`, and decays following a curve that reaches ~30% of its value at `z_bed=ytill.z0`. Finally, all values are limited to a minimum value of `ytill.cf_min`.

### Sub-temperate sliding

The friction laws above cap the basal stress at a value set by $c_b$, so a
frozen bed still slides wherever $\tau_d$ is large compared to $\beta$. With
`ydyn.slide_T=True`, sliding is reduced below the pressure melting point by
dividing $\beta$ by a sliding factor

$$
f_{\rm slide} = {\rm max}\left(\lambda_{\rm min}, \exp(T'_b/\gamma_T)\right)
$$

where $T'_b \le 0$ is the basal homologous temperature, $\gamma_T$ is
`ydyn.gamma_T` [K] and $\lambda_{\rm min}$ is `ydyn.lambda_min` (e.g., Fowler, 1986;
Hindmarsh and Le Meur, 2001). $f_{\rm slide}=1$ where the base is temperate
or not grounded. $\beta$ is divided by $f_{\rm slide}$ (output as `f_slide`) at
the cell centres, after the friction law and before the grounding-line scaling
and the staggering to the velocity nodes, for any `beta_method` except an
imposed $\beta$ (`beta_method=-1`). The velocity nodes across a frozen/temperate
transition get the mean of $\beta/f_{\rm slide}$. For DIVA, $\beta \to \infty$ tends to the
no-slip limit $\beta_{\rm eff}=1/F_2$, so small $\lambda_{\rm min}$ is safe.

## Effective pressure

Effective pressure (`N_eff`, $N_{\rm eff}$) in Yelmo is currently only used in the basal friction formulation as shown above. It provides a mechanism to alter the basal friction as a function of the state of the ice sheet, which is separate from $c_{\rm b,ref}$ (`cb_ref`), which represents the properties of the bed beneath the ice sheet. `N_eff` comes from the basal hydrology (FastHydrology, group `&yhyd`). With the till bucket (`yhyd.method_til=1`, no transport), the closure is chosen with `yhyd.bkt_N_closure`:

```bash
yhyd.bkt_N_closure = [-1,0,1,2,3]
-1: Set N externally (host model, via yelmo_set_var2D("hyd_N")), not modified internally.
 0: Impose a constant value, N_eff = yhyd.const_N
 1: Impose the overburden pressure, N_eff = rho_ice*g*H_ice
 2: Marine closure following Leguy et al. (2014), exponent yhyd.marine_p
 3: Till pressure following Bueler and van Pelt (2015), parameters yhyd.till_*
```

The former `&yneff` group (`yneff.method`) was retired in v2.1.
