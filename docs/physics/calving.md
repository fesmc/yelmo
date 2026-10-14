# Calving

Yelmo has two calving paths. By default (`ycalv.use_lsf = True`), the calving
front is a level set that moves with the ice velocity plus a calving velocity.
With `use_lsf = False`, calving is a mass-balance rate applied to front cells.
The allowed calving laws depend on the path. Frontal melt is a separate term
(see [Mass conservation](mass_conservation/index.md#frontal-mass-balance)).

## Default: level-set front with von Mises calving

The defaults (`input/yelmo_defaults.nml`) use a level-set front
(`ycalv.use_lsf = True`), the von Mises law `vm-m16` for floating and
marine-grounded fronts (`calv_flt_method = calv_grnd_method = "vm-m16"`) and
the subgrid front (`ytopo.front_subgrid = "marine"`). This path is computed in
`calc_ytopo_calving_lsf` (`src/yelmo_topography.f90`).

### Calving rate

`vm-m16` follows Morlighem et al. (2016, Eq. 4). On each velocity face the
calving rate is

$$
c = u \, \max\left(0, \frac{\tau_1}{\tau_{\rm ice}}\right)
$$

directed against the ice velocity $u$, where $\tau_1$ is the first principal
value of the depth-averaged deviatoric stress (mean of the two cells, or the
ice cell's value at a face between ice and an ice-free cell) and
$\tau_{\rm ice}$ is the ice strength: `ycalv.tau_ice_flt` (default 250 kPa) on
floating faces (`f_grnd_ac = 0`) and `ycalv.tau_ice_grnd` (default 1 MPa) on
the other marine faces. The floating law is used on faces with a grounded
fraction of zero, the grounded law on the other faces. On faces between two
land cells (bed at or above sea level), $c = -u$, so the level set does not
move there. The front is stationary where
$\tau_1 = \tau_{\rm ice}$.

The stresses come from the last velocity solution. A cell that received ice
since then takes the mean of its edge neighbours that were part of that
solution (`fill_stress_new_ice`).

Other laws on this path:

- `"zero"` / `"none"`: no calving;
- `"equil"`: the front is held in place ($c = -u$);
- `"threshold"`: the front retreats where the ice is thinner than
  `Hc_ref_flt` (floating) or `Hc_ref_grnd` (grounded);
- `"exp1"`–`"exp5"` (floating): the CalvingMIP experiments;
- `"ismip7"` (grounded): frontal retreat from the frontal melt of Rignot et
  al. (2016), with the subglacial discharge `bnd%Qd` and the thermal forcing
  `bnd%T_shlf` minus the seawater freezing point at the water depth (Jenkins,
  1991; ISMIP7 protocol). It does not use `bnd%tf_shlf`.

### Level-set front

The level set $\varphi$ is negative in the ice domain and positive in the
ocean; its zero contour is the calving front (`src/physics/calving/lsf_module.f90`).
It is advected with the front velocity $w = u + c$ (with $c$ against $u$),

$$
\frac{\partial \varphi}{\partial t} + w \cdot \nabla \varphi = 0,
$$

with first-order upwinding on the faces and explicit sub-steps (`LSFupdate`).
$w$ is extended into the ocean from the faces next to ice, and $\varphi$ is
kept in $[-1, 1]$. Cells with the bed at or above sea level are set to
$\varphi = -1$. With `ycalv.lsf_method = "snap"` (default), cells more than
two cells from the front are reset to $\pm 1$ (with `dt_lsf` > 0, the level set
is also reset to $\pm 1$ from the ice mask every `dt_lsf` years); `"redist"`
uses Sussman/Osher redistancing instead (`lsf_redist_n_iter` iterations).
Marine cells where `bnd%mask_ice` is `NONE` are set to $\varphi = 1$.

Ice in cells with $\varphi > 0$ is removed in one step, except in cells that
can be partial front cells (below). Ice-free marine cells behind the front
without an ice-covered edge neighbour are returned to the ocean ($\varphi = 1$),
except with `calv_flt_method = "equil"`.

### Subgrid front

With `ytopo.front_subgrid = "marine"`, ice cells with the bed below sea level
(`"floating"`: floating cells only) can be partially ice covered. A front cell
(such a cell with an ice-free ocean edge neighbour) has an effective thickness
`H_eff`, taken from its thickest interior edge neighbour (or diagonal
neighbour, if there is none) minus `ytopo.front_dHdx`·distance (at least
`ytopo.front_H_eff_min`, with a limit on the rise of the effective surface),
and the ice fraction `f_ice = H_ice/H_eff`. Front cells without an interior
neighbour keep `H_eff = H_ice`, and so do front cells entirely behind the
level-set front (`a_lsf` = 1). Ice thinner than 1 mm (`H_ice_eps`) is ice free
(`f_ice` = 0) and is never a front cell.

With the level set, the thickness of these cells follows the front
(`calc_G_lsf_front`): `a_lsf`, the area fraction of the cell behind the zero
contour, is computed from the cell's centre, edge and corner values of
$\varphi$. Ice-free land neighbours, where $\varphi$ is held at −1, take the
cell's own value, so they do not count as ice. Cells with `a_lsf` < 0.1 are emptied. With the level set, eligible
cells cut by the front (`a_lsf` < 1) that touch the ocean only at a corner are
front cells too; the trim and `f_ice` use this one front/interior
classification. Front cells with `a_lsf` < 1 are trimmed to `a_lsf`·`H_ref`,
with `H_ref` the reference thickness from the interior neighbours, so `f_ice`
is approximately `a_lsf`. In the momentum balance, partial front cells use
`H_eff` as their thickness and the front boundary condition is applied on
their ocean faces. Cells with `f_ice` < `A_FRONT_MIN` = 0.1 keep their ice and
fill by transport, but are ice free in the momentum balance until they reach
0.1 (`H_ice_dyn` = 0); otherwise a film of a few millimetres would be a full
`H_eff` column and turn the ocean face of the neighbouring front cell into an
interior face.

## Mass-balance calving path

With `use_lsf = False`, the calving laws give a calving mass balance `cmb` in
front cells:

- floating ice (`calv_flt_method`):
  - `"zero"` / `"none"`: no calving;
  - `"threshold"`: ice thinner than a critical thickness calves with the time
    scale `calv_tau`. The critical thickness changes from `Hc_ref_flt` to
    `Hc_deep` as the bed deepens from `zb_deep_0` to `zb_deep_1` (bed smoothed
    with `zb_sigma`);
  - `"vm-l19"`: von Mises effective stress (Lipscomb et al., 2019), below;
  - `"eigen"`: eigen calving (Levermann et al., 2012), scaling factor `k2`;
  - `"kill"`: all floating ice is removed;
  - `"kill-pos"`: floating ice is removed where `bnd%calv_mask` is set;
  - both kill methods act at the end of the calving step, after the front
    advance, so no floating ice remains in the kill region after each step;
- grounded ice (`calv_grnd_method`): `"zero"` / `"none"` or `"stress-b12"`
  (Bassis and Walker, 2012). Any grounded law also adds calving of grounded
  ice where the sub-grid bed roughness `z_bed_sd` is large, rising from zero at
  `sd_min` to `calv_grnd_max` at `sd_max`.

With the subgrid front (`ytopo.front_subgrid` ≠ `"none"`), the calving demand
is applied along the front and the front can advance into neighbouring cells.
Without it, thin ice is calved at the rate `calv_thin` below `Hc_ref_thin`
(`vm-l19`, `eigen`), and thin floating tongues are removed.

### Lipscomb et al. (2019)

$$
c = k_\tau\, \tau_{\rm ec}
$$

where $k_\tau$ [m yr$^{-1}$ Pa$^{-1}$] is an empirical constant (`kt_ref`,
default 0.0025, changing to `kt_deep` as the bed deepens from `zb_deep_0` to
`zb_deep_1`) and $\tau_{\rm ec}$ [Pa] is the effective calving stress,

$$
\tau_{\rm ec}^2 = \max(\tau_1,0)^2 + \omega_2 \max(\tau_2,0)^2.
$$

$\tau_1$ and $\tau_2$ are the eigenvalues of the depth-averaged horizontal
deviatoric stress tensor and $\omega_2$ is an empirical weighting constant
(`ycalv.w2`, default 0; Lipscomb et al., 2019, use 25). The calving rate is
applied to floating cells with an ice-free ocean neighbour, as the mass
balance $-\min(H\,c/\Delta x, 2000\ \mathrm{m\,a^{-1}})$.

The stress tensor is $\tau_{ij} = 2\bar\mu\,\dot\varepsilon_{ij}$
(`calc_stress_tensor_2D`), and its eigenvalues are the roots of

$$
\lambda^2 - (\tau_{xx} + \tau_{yy})\,\lambda + \tau_{xx}\tau_{yy} - \tau_{xy}^2 = 0
$$

(`calc_2D_eigen_values`).

## References

- Bassis, J. N. and Walker, C. C. (2012). Upper and lower limits on the
  stability of calving glaciers from the yield strength envelope of ice.
  Proc. R. Soc. A, 468, 913–931.
- Jenkins, A. (1991). A one-dimensional model of ice shelf–ocean interaction.
  J. Geophys. Res., 96(C11), 20671–20677.
- Levermann, A., et al. (2012). Kinematic first-order calving law implies
  potential for abrupt ice-shelf retreat. The Cryosphere, 6, 273–286.
- Lipscomb, W. H., et al. (2019). Description and evaluation of the Community
  Ice Sheet Model (CISM) v2.1. Geosci. Model Dev., 12, 387–424.
- Morlighem, M., et al. (2016). Modeling of Store Gletscher's calving dynamics,
  West Greenland, in response to ocean thermal forcing. Geophys. Res. Lett.,
  43, 2659–2666.
- Rignot, E., et al. (2016). Modeling of ocean-induced ice melt rates of five
  West Greenland glaciers over the past two decades. Geophys. Res. Lett., 43,
  6374–6382.
