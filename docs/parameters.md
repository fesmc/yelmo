
# Parameters

Here important parameter choices pertinent to running
**Yelmo** will be documented. Each section will
outline a specific parameter or set of related parameters.
The author of each section and the date last updated
will apear in the heading, to maintain traceability
in the documentation (since code usually changes over time).

**This is a work in progress!**

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

Another possibility is to tune `cb_ref` as a function of other model or boundary variables. The most common approach is to tune it as as function of the bedrock elevation relative to present-day sea level (e.g., Winkelmann et al., 2011). In Yelmo, this is controlled by the parameter choices in the `ytill` section, and in particular the parameter `ytill.scale=['none','lin','exp']`. When `ytill.scale='none'`, no scaling function is applied and then `cb_ref=ytill.cf_ref` everywhere. When `ytill.scale='lin'`, a linear scaling is applied so that `cb_ref` goes from `ytill.cf_min` to `ytill.cb_ref` for bedrock elevations between `ytill.z0` and `ytill.z1` (saturating otherwise). Finally, if `ytill.scale='exp'`, an exponential decay function is applied, such that `cb_ref=ytill.cf_ref` for `z_bed >= ytill.z1`, and decays following a curve that reaches ~30% of its value at `z_bed=ytill.z0`. Finally, all values are limited to a minimum value of `ytill.cf_min`.

## Effective pressure

Effective pressure (`N_eff`, $N_{\rm eff}$) in Yelmo is currently only used in the basal friction formulation as shown above. It provides a mechanism to alter the basal friction as a function of the state of the ice sheet, which is separate from $c_{\rm b,ref}$ (`cb_ref`), which represents the properties of the bed beneath the ice sheet. The calculation of `N_eff` can be done with several methods:

```bash
yneff.method = [-1,0,1,2,3]
-1: Set N_eff external to Yelmo, do not modify it internally.
 0: Impose a constant value, N_eff = yneff.const
 1: Impose the overburden pressure, N_eff = rho_ice*g*H_ice
 2: Calculate N_eff following the Leguy formulation
 3: Calculate N_eff as till pressure following Bueler and van Pelt (2015). 
 4: Calculate N_eff as a 'two-valued' function scaled by f_pmp using yneff.delta.
 ```
