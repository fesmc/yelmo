# DIVA — Depth-Integrated Viscosity Approximation

DIVA is the first-class momentum balance solved in Yelmo. It is a 2D
(depth-integrated) problem for the horizontal velocity, augmented with a
1D vertical closure that recovers the full 3D velocity field. The formulation
follows [Goldberg (2011)](https://doi.org/10.1017/S0022143000004585),
[Arthern et al. (2015)](https://doi.org/10.1002/2014JF003239),
[Lipscomb et al. (2019)](https://doi.org/10.5194/gmd-12-3199-2019),
and [Robinson et al. (2022)](https://doi.org/10.5194/tc-16-689-2022), whose
notation we use here.

DIVA is valid uniformly from slow grounded ice (where vertical shear dominates,
the SIA limit) to fast streaming ice and floating shelves (where membrane
stresses dominate, the SSA limit). It does so by retaining the depth-integrated
membrane stress balance of SSA while reintroducing vertical shear through a
closure on $\frac{\partial u}{\partial z}$ that is consistent with the Glen
flow law.

The implementation is in [`src/physics/velocity_diva.f90`][diva_src].

## Continuum equations

### Depth-integrated stress balance

The unknown is the depth-averaged horizontal velocity
$\bar{\mathbf u} = (\bar u, \bar v)$. The momentum balance is the
depth integral of the horizontal Stokes equations, with the assumption
that the longitudinal stresses do not vary in the vertical
($\frac{\partial u}{\partial z}$ contributes only to the strain-rate
invariant, not as an independent unknown):

$$
\frac{\partial}{\partial x}\!\left[\, 2\,\bar\mu\,H\,(2\,\bar\varepsilon_{xx} + \bar\varepsilon_{yy})\,\right]
+ \frac{\partial}{\partial y}\!\left[\, 2\,\bar\mu\,H\,\bar\varepsilon_{xy}\,\right]
- \beta_\mathrm{eff}\,\bar u
\;=\; \tau_{d,x},
$$

$$
\frac{\partial}{\partial y}\!\left[\, 2\,\bar\mu\,H\,(2\,\bar\varepsilon_{yy} + \bar\varepsilon_{xx})\,\right]
+ \frac{\partial}{\partial x}\!\left[\, 2\,\bar\mu\,H\,\bar\varepsilon_{xy}\,\right]
- \beta_\mathrm{eff}\,\bar v
\;=\; \tau_{d,y},
$$

with the depth-averaged effective viscosity

$$
\bar\mu \;=\; \frac{1}{H}\int_b^s \mu \,\mathrm dz.
$$

(The product $\bar\mu\,H$ is stored in the code as `visc_eff_int`.) The
horizontal strain rates are

$$
\bar\varepsilon_{xx} \;=\; \frac{\partial \bar u}{\partial x},
\qquad
\bar\varepsilon_{yy} \;=\; \frac{\partial \bar v}{\partial y},
\qquad
\bar\varepsilon_{xy} \;=\; \tfrac{1}{2}\!\left(\frac{\partial \bar u}{\partial y} + \frac{\partial \bar v}{\partial x}\right).
$$

The **driving stress** on the right-hand side is the depth-integrated
gravitational force,

$$
\tau_{d,x} \;=\; \rho_i\, g\, H\, \frac{\partial s}{\partial x},
\qquad
\tau_{d,y} \;=\; \rho_i\, g\, H\, \frac{\partial s}{\partial y},
$$

evaluated on ac-faces by [`calc_driving_stress`][diva_taud] in
`velocity_general.f90` and limited to `ydyn.taud_lim` (default 200 kPa). With
`ydyn.taud_gl_method` > 0, the surface slope at the grounding line is modified
(`calc_driving_stress_gl`): 1, slope weighted by the grounded fraction, with a
virtual slope for floating ice; 2, one-sided differences (Feldmann et al.,
2014); 3, linear interpolation (Gladstone et al., 2010). The default is 0 (no
modification).

At ice fronts, the depth-integrated imbalance between ice and water pressure
$\tau_{l,\mathrm{int}}$ enters as a Neumann condition on the front faces.
`ydyn.ssa_lat_bc` selects the fronts where it is applied: `"marine"`
(default; floating and marine-terminating fronts), `"floating"` (or
`"float"`; floating fronts), `"all"` or `"none"`. The other front faces are
solved as inner faces (see [Numerical solution](solvers.md)). A face between
ice and ice-free land is not an ocean front: next to floating ice it carries
no velocity, and next to marine ice it is treated as a grounded front.

### Effective viscosity from Glen's flow law

The effective viscosity follows Glen's law,

$$
\mu \;=\; \tfrac{1}{2}\, A^{-1/n}\, \dot\varepsilon_e^{\,(1-n)/n},
$$

computed in 3D with the rate factor $A(T'(z))$ and depth-integrated to
$\bar\mu H$. With `ydyn.visc_method = 1` (default) it is evaluated at Gaussian
quadrature points of each cell (`calc_visc_eff_3D_nodes`), with `visc_method = 2`
at the cell centres (`calc_visc_eff_3D_aa`), both in [`velocity_diva.f90`][diva_src];
`visc_method = 0` imposes the constant `visc_const`. The effective strain rate $\dot\varepsilon_e$ in DIVA uses the
horizontal strain rates plus a contribution from the (closure-defined)
vertical shear, regularised by $\varepsilon_0$ (`ydyn.eps_0`) to avoid singular
viscosity in stagnant ice:

$$
\dot\varepsilon_e^{\,2}
\;=\; \dot\varepsilon_{xx}^{\,2}
   + \dot\varepsilon_{yy}^{\,2}
   + \dot\varepsilon_{xx}\,\dot\varepsilon_{yy}
   + \dot\varepsilon_{xy}^{\,2}
   + \tfrac{1}{4}\!\left(\frac{\partial u}{\partial z}\right)^{\!2}
   + \tfrac{1}{4}\!\left(\frac{\partial v}{\partial z}\right)^{\!2}
   + \varepsilon_0^{\,2}.
$$

The two
$\left(\frac{\partial u}{\partial z}\right)^{\!2},
 \left(\frac{\partial v}{\partial z}\right)^{\!2}$
terms are what make DIVA different from SSA: they reintroduce the
vertical-shear contribution to the deformation rate, and hence soften
the ice in places where vertical shear is large (slow, grounded
interior).

### Vertical-shear closure

Because the depth-integrated equation does not solve for the full vertical
velocity profile, $\frac{\partial u}{\partial z}$ and
$\frac{\partial v}{\partial z}$ are obtained from a closure that is
consistent with the Stokes balance under DIVA's assumptions. With the
basal shear stress $\boldsymbol\tau_b$ acting at $z = b$ and vanishing
shear stress at the surface $z = s$,

$$
\frac{\partial u}{\partial z}(z) \;=\; \frac{\tau_{b,x}\,(s-z)}{\mu(z)\,H},
\qquad
\frac{\partial v}{\partial z}(z) \;=\; \frac{\tau_{b,y}\,(s-z)}{\mu(z)\,H}
$$

(Robinson et al. 2022 Eq. 21; routine `calc_vertical_shear_3D`). Define
the **F-integrals**

$$
F_m \;\equiv\; \int_b^s \frac{1}{\mu(z)}\!\left(\frac{s-z}{H}\right)^{\!m}\,\mathrm dz
$$

(Robinson et al. 2022 Eq. 15; routine `calc_F_integral`). These are
precisely the integrals needed to relate the depth-averaged velocity to
the basal velocity:

$$
\bar u \;=\; u_b \,+\, \tau_{b,x}\,F_2,
\qquad
\bar v \;=\; v_b \,+\, \tau_{b,y}\,F_2.
$$

### Basal stress closure and effective drag

With a basal friction law of the form $\boldsymbol\tau_b = \beta\,\mathbf u_b$
(see [Basal friction](../basal-friction.md) for the supported choices of
$\beta$), the relation
$\bar u = u_b + \tau_{b,x}\,F_2 = u_b\,(1 + \beta\,F_2)$ inverts to
express $\boldsymbol\tau_b$ directly in terms of the depth-averaged
velocity:

$$
\boldsymbol\tau_b \;=\; \beta_\mathrm{eff}\,\bar{\mathbf u},
\qquad
\beta_\mathrm{eff} \;\equiv\; \frac{\beta}{1 + \beta\,F_2}
$$

(Robinson et al. 2022 Eq. 19; routine `calc_beta_eff`). In the no-slip
limit ($\beta \to \infty$) this reduces to
$\beta_\mathrm{eff} = 1/F_2$ (Robinson et al. 2022 Eq. 20), which is the
form used with `ydyn.solver = "diva-noslip"`. The effective drag
$\beta_\mathrm{eff}$ — not the raw $\beta$ — is what appears in the
depth-integrated momentum equation above.

### Picard iteration

The viscosity $\mu(\bar{\mathbf u})$ and the effective drag
$\beta_\mathrm{eff}(\bar{\mathbf u})$ both depend on the solution, so the
momentum equation is nonlinear. Yelmo solves it by Picard iteration
(`calc_velocity_diva`). Each iteration:

1. computes the vertical shear $\partial u/\partial z$, $\partial v/\partial z$
   from the current basal stress and viscosity;
2. computes the 3D viscosity $\mu$ from $\bar{\mathbf u}$ and the vertical
   shear, relaxes it in log-space against the previous iterate (Sandip et al.,
   2023; `picard_relax_visc`), and integrates it to $\bar\mu H$;
3. computes $\beta$ from the current basal velocity (see
   [Basal friction](../basal-friction.md)), the F-integral $F_2$, and
   $\beta_\mathrm{eff}$ on the velocity faces. A non-zero $\beta$ on the
   faces is raised to at least `beta_min`, and $\beta_\mathrm{eff}$ is set to
   `beta_min` on grounded faces where it is zero (`set_beta_min_grounded`);
4. adds the [velocity-limit drag](solvers.md#velocity-limit), and assembles and
   solves the linear system for $\bar{\mathbf u}$ (see
   [Numerical solution](solvers.md));
5. relaxes $\bar{\mathbf u}$ against the previous iterate with the factor
   `ydyn.ssa_iter_rel` (default 0.7; `picard_relax_vel`);
6. updates the basal stress and the basal velocity;
7. with an effective pressure that depends on the basal velocity, re-evaluates
   $N$ from $u_b$, and with it the friction coefficient `c_bed` (`neff_hook`).
   This happens with the K24 hydrology (`yhyd.method_transport = 1`) and
   `yhyd.k24_N_ub_coupled = True` (default), or when a host has registered a
   callback with `yelmo_set_neff_callback` (C API). With
   `k24_N_ub_coupled = False`, $N$ is evaluated once per time step.

The iteration stops when the relative L2 change of the velocity,
$\lVert \bar{\mathbf u}^{k} - \bar{\mathbf u}^{k-1}\rVert / \lVert \bar{\mathbf u}^{k-1}\rVert$
(De Smedt et al., 2010), over the solved faces faster than $10^{-5}$ m a$^{-1}$,
is below `ydyn.ssa_iter_conv` (default $10^{-2}$), or after `ydyn.ssa_iter_max`
(default 20) iterations. With the hook of step 7, the relative L1 change of
`c_bed` must also be below `ssa_iter_conv`. The hook is used only with DIVA,
not with the SSA or hybrid solvers.

When the iteration converges, the basal stress
$\boldsymbol\tau_b = \beta_\mathrm{eff}\bar{\mathbf u}$ and basal velocity
$u_b = \bar u - \tau_{b,x}\,F_2$ are diagnosed, and the 3D horizontal
velocity field is reconstructed by integrating the vertical-shear closure
from the bed upward:

$$
u(z) \;=\; u_b \,+\, \tau_{b,x}\, \tilde F_1(z),
\qquad
\tilde F_1(z) \;\equiv\; \int_b^z \frac{1}{\mu(z')}\!\left(\frac{s-z'}{H}\right)\,\mathrm dz'
$$

(routine `calc_vel_horizontal_3D`).

## Limits

- **SSA limit**: as $\frac{\partial u}{\partial z},
  \frac{\partial v}{\partial z} \to 0$, the F-integrals enter only
  through $\beta_\mathrm{eff} \to \beta$ (no vertical-shear enhancement
  of $\dot\varepsilon_e$, no basal/depth-averaged offset). The
  depth-integrated stress balance reduces to the SSA equation — see
  [SSA](ssa.md).
- **SIA limit**: as membrane gradients become negligible relative to
  vertical-shear gradients, the balance collapses to
  $\boldsymbol\tau_b \approx \boldsymbol\tau_d$, and the closure on
  $\frac{\partial u}{\partial z}$ reproduces the SIA velocity profile —
  see [SIA](sia.md).

[diva_src]: https://github.com/fesmc/yelmo/blob/main/src/physics/velocity_diva.f90
[diva_taud]: https://github.com/fesmc/yelmo/blob/main/src/physics/velocity_general.f90
