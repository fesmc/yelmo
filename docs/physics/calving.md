# Calving schemes

Here is a summary of calving schemes.

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
the other marine faces. Faces between two land cells (bed at or above sea level)
get $c = u$, so the level set does not move there. The front is stationary where
$\tau_1 = \tau_{\rm ice}$.

The stresses come from the last velocity solution. A cell that received ice
since then takes the mean of its edge neighbours that were part of that
solution (`fill_stress_new_ice`).

Other laws on this path: `"zero"`, `"equil"` (front held in place),
`"threshold"`, the CalvingMIP laws `"exp1"`–`"exp5"` (floating) and
`"ismip7"` (grounded).

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
two cells from the front are reset to $\pm 1$; `"redist"` uses Sussman/Osher
redistancing instead.

Ice in cells with $\varphi > 0$ is removed in one step, except in cells that
can be partial front cells (below). Ice-free marine cells behind the front
without an ice-covered edge neighbour are returned to the ocean ($\varphi = 1$).

### Subgrid front

With `ytopo.front_subgrid = "marine"`, ice cells with the bed below sea level
(`"floating"`: floating cells only) can be partially ice covered. A front cell
(such a cell with an ice-free ocean edge neighbour) has an effective thickness
`H_eff`, taken from its thickest interior neighbour minus
`ytopo.front_dHdx`·distance (at least `ytopo.front_H_eff_min`, with a limit on
the rise of the effective surface), and the ice fraction
`f_ice = H_ice/H_eff`. Front cells without an interior neighbour keep
`H_eff = H_ice`.

With the level set, the thickness of these cells follows the front
(`calc_G_lsf_front`): `a_lsf`, the area fraction of the cell behind the zero
contour, is computed from the cell's centre, edge and corner values of
$\varphi$. Cells with `a_lsf` < 0.1 are emptied, and front cells (also those
touching the ocean at a corner) are trimmed to `a_lsf`·`H_ref`, with `H_ref`
the reference thickness from the interior neighbours. `f_ice` is then
approximately `a_lsf`. In the momentum balance, partial front cells use
`H_eff` as their thickness and the front boundary condition is applied on
their ocean faces.

## Lipscomb et al. (2019)

Used on the mass-balance calving path (`use_lsf = False`,
`calv_flt_method = "vm-l19"`).

$$
c = k_\tau \tau_{\rm ec}
$$

where $k_\tau$ (m yr$^{-1}$ Pa$^{-1}$) is an empirical constant and $\tau_{\rm ec}$ (Pa) is the effective calving stress, which is defined by:

$$
\tau_{\rm ec}^2 = \max(\tau_1,0)^2 + \omega_2 \max(\tau_2,0)^2
$$

$\tau_1$ and $\tau_2$ are the eigenvalues of the 2D horizontal deviatoric stress tensor and $\omega_2$ is an empirical weighting constant (`ycalv.w2`).

The eigenvalues $\tau_1$ and $\tau_2$ are calculated from the depth-averaged (2D) stress tensor $\tau_{\rm ij}$ as follows. Given the stress tensor components $\tau_{\rm xx}$, $\tau_{\rm yy}$ and $\tau_{\rm xy}$, we can solve for the real roots $\lambda$ of the tensor from the quadratic equation:

$$
a \lambda^2 + b \lambda + c = 0
$$

where

$$
a = 1.0 \\
b = -(\tau_{\rm xx} + \tau_{\rm yy}) \\
c = \tau_{\rm xx}*\tau_{\rm yy} - \tau_{\rm xy}^2
$$

glissade_velo_higher.F90:

```fortran
tau_xz(k,i,j) = tau_xz(k,i,j) + efvs_qp * du_dz            ! 2 * efvs * eps_xz
tau_yz(k,i,j) = tau_yz(k,i,j) + efvs_qp * dv_dz            ! 2 * efvs * eps_yz
tau_xx(k,i,j) = tau_xx(k,i,j) + 2.d0 * efvs_qp * du_dx     ! 2 * efvs * eps_xx
tau_yy(k,i,j) = tau_yy(k,i,j) + 2.d0 * efvs_qp * dv_dy     ! 2 * efvs * eps_yy
tau_xy(k,i,j) = tau_xy(k,i,j) + efvs_qp * (dv_dx + du_dy)  ! 2 * efvs * eps_xy
```
