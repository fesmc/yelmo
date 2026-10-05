# Thermodynamics

Yelmo solves for the energy of the ice in terms of enthalpy, which describes
cold and temperate ice with one variable (Aschwanden et al., 2012). The ice
column is coupled to a column of bedrock below. The thermodynamics module is
`calc_ytherm` in
[`src/yelmo_thermodynamics.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_thermodynamics.f90);
the column solvers are in
[`src/physics/ice_enthalpy.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/ice_enthalpy.f90)
and the material properties and boundary fluxes in
[`src/physics/thermodynamics.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/thermodynamics.f90).

## Methods

The method is selected with `ytherm.method`:

| `method` | Description |
|---|---|
| `"enth"` | Enthalpy solver (**default**): cold and temperate ice, water content $\omega$ |
| `"temp"` | Temperature solver (`calc_temp_column`): no water content, heat above the pressure melting point is drained as basal melt, centred vertical advection, `"wtil"` basal boundary condition |
| `"robin"` | Steady Robin solution in grounded columns with positive net mass balance; linear profiles elsewhere (base at $T_\mathrm{pmp}$ for grounded ice, 271.15 K for floating ice) |
| `"robin-cold"` | Robin solution averaged with a linear profile, to obtain a cold base |
| `"linear"` | Linear profile between the surface temperature and $T_\mathrm{pmp} - 10$ K (frozen bed) at the base |
| `"fixed"` | No update: the temperature/enthalpy fields stay as initialised |

`yelmo_init_state` initialises the thermodynamic state with `"linear"`,
`"robin"` or `"robin-cold"` (argument `thrm_method`). The Robin solution
assumes constant material properties and uses `const_kt` and `const_cp`
(whatever `use_const_kt` / `use_const_cp`), so that the basal gradient is
$-G/k$ with $k$ = `const_kt`.

## Enthalpy

The specific enthalpy $E$ [J kg$^{-1}$] is related to the temperature $T$ and
the water content $\omega$ by

$$
E = \begin{cases}
E_c(T) & T < T_\mathrm{pmp} \\
E_c(T_\mathrm{pmp}) + \omega L & T = T_\mathrm{pmp}
\end{cases}
$$

where $L$ is the latent heat of fusion (`L_ice`) and $T_\mathrm{pmp}$ is the
pressure melting point,

$$
T_\mathrm{pmp} = T_0 - \beta \rho_i g (s - z),
$$

with $\beta$ = `T_pmp_beta` (`input/yelmo_phys_const.nml`). The cold-ice
enthalpy $E_c$ depends on `ytherm.enth_cp_method`:

- `"integral"` (default): $E_c(T) = \int_0^T c(T')\,\mathrm dT'$ with
  $c(T) = 146.3 + 7.253\,T$ J kg$^{-1}$ K$^{-1}$ (Greve and Blatter, 2009,
  Eq. 4.39), so that $\partial E/\partial T = c(T)$ exactly;
- `"const"`: $E_c = c_\mathrm{ref}\,T$ with $c_\mathrm{ref}$ = 2009 J kg$^{-1}$ K$^{-1}$.

The thermal conductivity is $k(T) = 9.828\,e^{-0.0057\,T}$ W m$^{-1}$ K$^{-1}$
(Greve and Blatter, 2009, Eq. 4.37). With `use_const_kt`, the constant
`const_kt` is used instead. `use_const_cp` / `const_cp` apply to the other
methods only (the Robin methods always use `const_kt` and `const_cp`); in the enthalpy solver, the heat capacity follows
`enth_cp_method`. Yelmo expects
conductivities in J a$^{-1}$ m$^{-1}$ K$^{-1}$, i.e. W m$^{-1}$ K$^{-1}$
multiplied by `sec_year`.

## Ice column

In the terrain-following coordinate $\zeta$, the enthalpy equation solved in each
column is

$$
\frac{\partial E}{\partial t}
= \frac{1}{H^2}\frac{\partial}{\partial \zeta}\!\left(K \frac{\partial E}{\partial \zeta}\right)
- \frac{w^\star}{H}\frac{\partial E}{\partial \zeta}
- \mathbf u \cdot \nabla_\zeta E
+ \frac{\Phi}{\rho_i},
$$

where $w^\star$ is the sigma-relative vertical velocity `uz_star` (see
[Vertical velocity](vertical-velocity.md)) and $\Phi$ is the strain heating.
The diffusivity is

$$
K = \begin{cases}
k/(\rho_i c) & \text{cold ice} \\
\epsilon\, k/(\rho_i c) & \text{temperate ice}
\end{cases}
$$

with the conductivity ratio $\epsilon$ = `enth_cr` (default $10^{-3}$), which
represents the small diffusion of water in temperate ice. On the cold side of the
cold–temperate transition surface (CTS), the flux across the interface is
computed with the cold-ice diffusivity (Blatter and Greve, 2015, Eq. 25).

The water content is limited to `omega_max` (default 0.01). Water above this
limit is drained to the bed and added to the basal melt rate (`melt_int`). The
height of the CTS above the bed is diagnosed as `H_cts`.

### Discretisation

Each column is solved implicitly (backward Euler) as a tridiagonal system
(`calc_enth_column_internal`), on the vertical aa-nodes of the ice layers,
with the diffusivity at the layer faces from a weighted harmonic mean.

- **Vertical advection** is upwind with second-order accuracy: the matrix
  holds first-order upwind, and the right-hand side adds the difference to a
  minmod-limited second-order upwind gradient from the start-of-step enthalpy
  (deferred correction). Where the profile is smooth this removes the numerical
  diffusion of first-order upwind; at extrema, such as the CTS, the limiter
  reduces the scheme to first-order upwind.
- **Horizontal advection** is explicit and enters the column equation as a
  source term (`calc_advec_horizontal_3D`). It is computed in flux form,
  $\mathbf u\cdot\nabla E = \nabla\cdot(\mathbf u E) - E\,\nabla\cdot\mathbf u$,
  with face values from a van Leer MUSCL reconstruction
  (`advecxy_order = 2`, default) or plain upwind (`advecxy_order = 1`). The
  scheme is monotone and gives zero advection for a uniform field. When the
  maximum horizontal Courant number over the thermodynamic time step exceeds
  `advecxy_cfl` (default 0.5), the horizontal advection of the whole domain is
  sub-cycled, with $\lceil \mathrm{CFL}_\mathrm{max}/$`advecxy_cfl`$\rceil$
  sub-steps, at most `advecxy_nmax` (default 10).
- **Strain heating** is $\Phi = 4\mu\dot\varepsilon_e^2$ (Greve and Blatter, 2009,
  Eqs. 4.7, 5.65), from the 3D viscosity and effective strain rate of the
  material module (`strain_heating = "full"`, default). With
  `strain_heating = "sia"`, the SIA approximation is used instead, and with
  `"none"` strain heating is switched off.

The thermodynamics uses the column of the dynamics, with thickness `H_ice_dyn`
(the effective thickness `H_eff` in partial front cells). Columns are solved
where the cell is fully ice covered and thicker than `H_ice_thin` (default
10 m). Thinner columns get a linear temperature profile between the surface
and the basal temperature. Partially covered and ice-free cells next to the
ice take the mean of their fully covered neighbours, so that ice advected into
them starts with a realistic temperature.

The column solvers are tested against the Kleiner et al. (2015) benchmarks
with a standalone driver (see [Benchmarks](../benchmarks.md#enthalpy-column-tests)).

### Surface boundary condition

The enthalpy at the surface is prescribed from the surface temperature,
$E_s = E_c(\min(T_\mathrm{srf}, T_0))$.

### Basal boundary condition

The basal heat balance involves the geothermal (or bedrock) heat flux $G$
(`Q_rock`), the frictional heating $Q_b = \boldsymbol\tau_b \cdot \mathbf u_b$
(`Q_b`), the heat supplied by basal water $Q_w$ (dissipation and sensible heat,
from the hydrology model), and the conductive flux into the ice
$Q_{i,b} = -k\,\partial T/\partial z|_b$ (`Q_ice_b`).

**Grounded ice.** With `basal_bc_method = "capacity"` (default), the model
first computes the basal mass balance $\dot b^\star$ (`bmb_grnd_star`) that
would keep the base at the pressure melting point. It is compared with the
freeze-on capacity $C$, the rate at which the water at the bed can be frozen:

- if $\dot b^\star \le 0$ (melting) or $\dot b^\star \le C$, the base is held
  at the pressure melting point (Dirichlet condition);
- otherwise all available water freezes ($\dot b = C$), and the base cools with
  the flux condition
  $k\,\partial T/\partial z|_b = -(G + Q_b + Q_w + \rho_i L C)$.

With a temperate layer above the base (CTS above the first layer), the basal
enthalpy follows the layer above (zero gradient). The source of $C$ is set by
`cap_source`: `"hyd"` (computed by the hydrology model), `"till"` (all till
water above `cap_W_floor` frozen in one step), `"water"` (same, from the water
thickness), `"none"` ($C = 0$), or `"auto"` (default; `"hyd"` with water
transport, otherwise `"till"`). Capacities below `cap_eps` count as a dry bed.
The diagnostic `bc_b` records the condition used (1: pressure melting point,
2: flux). The older rule `basal_bc_method = "wtil"` decides from a predictor of
the till water thickness and is deprecated.

The grounded basal mass balance is then (Cuffey and Paterson, 2010, Eq. 9.38)

$$
\dot b_g = -\frac{G + Q_b + Q_w - Q_{i,b}}{\rho_i\,(L - \max(E_b - E_\mathrm{pmp}, 0))},
$$

where the latent heat is reduced by the water already stored in the basal ice.
Melt is only allowed for a near-temperate base ($T_b - T_\mathrm{pmp} > -1$ K),
and freeze-on is limited to $C$. The englacial drainage `melt_int` is added to
the basal melt.

**Floating ice.** The base is held at the freezing temperature of sea water at
the depth of the ice base (Jenkins, 1991, salinity 34.75 psu), limited to the
pressure melting point. Near the grounding line, the base is at the pressure
melting point for ice less than 100 m below flotation, and blends to the sea
water freezing point between 100 and 200 m below flotation. The
diffusivity at the base is the cold-ice value, so the ocean can cool the ice
above and form a cold basal boundary layer even when the shelf is temperate.
Partially grounded cells use the mean of the two temperatures, weighted by the
grounded fraction.

### Basal frictional heating

The frictional heating is formed on the velocity faces and averaged to the cell
centres, selected by `qb_method`: 1, face products; 2 (default), face products
at quadrature nodes; 3, staggering of $u_b$ and $\tau_b$ to the cell centre;
4, quadrature of the staggered fields. Methods 1 and 2 conserve the energy
dissipated by basal friction.

## Bedrock column

Below the ice, a bedrock column of thickness `H_rock` (default 2000 m) with
`nzr_aa` layers is treated according to `ytherm.rock_method`:

- `"equil"` (default): linear profile in equilibrium with the geothermal heat
  flux, $\partial T_r/\partial z = -Q_\mathrm{geo}/k_r$, and the basal ice
  temperature at the top;
- `"active"`: the heat equation
  $\rho_r c_r\, \partial T_r/\partial t = k_r\, \partial^2 T_r/\partial z^2$ is
  solved with $Q_\mathrm{geo}$ at the bottom and the basal ice temperature (the
  ocean freezing temperature under floating ice) at the top. The flux at the top
  of the bedrock, `Q_rock`, is the geothermal input to the ice base.
- `"fixed"`: no update.

The bedrock properties are `rhoc_rock` (volumetric heat capacity, default
$2\times10^6$ J m$^{-3}$ K$^{-1}$) and `kt_rock` (default
$6.3\times10^7$ J a$^{-1}$ m$^{-1}$ K$^{-1}$, i.e. 2 W m$^{-1}$ K$^{-1}$).
In SICOPOLIS and GRISLI, $c_r$ = 1000 J kg$^{-1}$ K$^{-1}$ is used (Rogozhina
et al., 2012; Greve, 2005). The conductivity follows Rogozhina et al. (2012)
for Greenland, and is consistent with Lösing et al. (2020) for Antarctica and
the 2–3 W m$^{-1}$ K$^{-1}$ of the upper crust (Cammarano and Guerri, 2017).
Older values are 3 W m$^{-1}$ K$^{-1}$ (Greve, 1997, 2005) and
3.3 W m$^{-1}$ K$^{-1}$ (GRISLI).

## Coupling to the flow

The temperature and water content set the Glen rate factor in the material
module (`calc_ymat`, `ymat.rf_method = 1`):

$$
A = E_f\,A_0\,\exp\!\left(-\frac{Q}{R\,T^\ast}\right),
$$

with the homologous temperature on the absolute scale
$T^\ast = T - T_\mathrm{pmp} + T_0$ (limited to $[220\,\mathrm K, T_0]$), and $A_0$ = 1.25671$\times10^{-5}$ Pa$^{-3}$ a$^{-1}$,
$Q$ = 60 kJ mol$^{-1}$ for $T^\ast \le 263.15$ K and
$A_0$ = 6.0422976$\times10^{10}$ Pa$^{-3}$ a$^{-1}$, $Q$ = 139 kJ mol$^{-1}$
above (Greve and Blatter, 2009). `rf_use_eismint2 = True` uses the EISMINT2
constants (Payne et al., 2000). With `rf_with_water = True`, $A$ is multiplied
by $1 + 181.25\,\omega$ (Lliboutry and Duval, 1985). $E_f$ is the enhancement
factor (`enh_method`, `enh_shear`, `enh_stream`, `enh_shlf`).

The basal temperature enters the basal friction through the frozen-bed
sliding factor (`ydyn.frz_scale`, see [Basal friction](basal-friction.md)), with
the basal homologous temperature $T'_b = T_b - T_\mathrm{pmp}$ (`T_prime_b`).
The temperate fraction of the base, `f_pmp` = $\exp(\max(T'_b, -20\,\mathrm K)/\gamma)$
with $\gamma$ = `ytherm.gamma` (0 and 1 below 0.01 and above 0.99; binary for
$\gamma = 0$; 1 for floating ice), is a diagnostic used in the bed mask.

## References

- Aschwanden, A., Bueler, E., Khroulev, C., and Blatter, H. (2012). An enthalpy
  formulation for glaciers and ice sheets. J. Glaciol., 58, 441–457.
- Blatter, H. and Greve, R. (2015). Comparison and verification of enthalpy
  schemes for polythermal glaciers and ice sheets with a one-dimensional model.
  Polar Sci., 9, 196–207.
- Cuffey, K. M. and Paterson, W. S. B. (2010). *The Physics of Glaciers*, 4th ed.
- Greve, R. and Blatter, H. (2009). *Dynamics of Ice Sheets and Glaciers.* Springer.
- Jenkins, A. (1991). A one-dimensional model of ice shelf–ocean interaction.
  J. Geophys. Res., 96, 20671–20677.
- Kleiner, T., Rückamp, M., Bondzio, J. H., and Humbert, A. (2015). Enthalpy
  benchmark experiments for numerical ice sheet models. The Cryosphere, 9, 217–228.
