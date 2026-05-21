# Calving schemes

Here is a summary of calving schemes.

## Lipscomb et al. (2019)

$$
c = k_\tau \tau_{\rm ec}
$$

where $k_\tau$ (m yr$^{-1}$ Pa$^{-1}$) is an empirical constant and $\tau_{\rm ec}$ (Pa) is the effective calving stress, which is defined by:

$$
\tau_{\rm ec}^2 = \max(\tau_1,0)^2 + \omega_2 \max(\tau_2,0)^2
$$

$\tau_1$ and $\tau_2$ are the eigenvalues of the 2D horizontal deviatoric stress tensor and $\omega_2$ is an empirical weighting constant. For partially ice-covered grid cells (with $f_{\rm ice} < 1$), these stresses are taken from the upstream neighbor.

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
