# Vertical velocity — `uz` and `uz_star`

Yelmo carries **two** distinct vertical velocity fields, both defined on
vertical `ac`-nodes (cell edges, `nz_ac = nz_aa + 1`) and both owned by the
dynamics class (`dyn%now%uz`, `dyn%now%uz_star`):

| Field | Symbol | Meaning |
|---|---|---|
| `uz` | $w$ | the **physical** vertical velocity, in a fixed (Eulerian) frame |
| `uz_star` | $w^\star$ | the **sigma-relative** advective vertical velocity, $w^\star = H\,\mathrm{D}\zeta/\mathrm{D}t$ |

They are not interchangeable, and confusing them is easy because both have
units of $\mathrm{m\,a^{-1}}$ and both are "a vertical velocity". This page
derives each one, maps it onto the code, and states the rule for which to use
where.

Both are computed together in
[`src/physics/velocity_general.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/velocity_general.f90),
dispatched from
[`src/yelmo_dynamics.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_dynamics.f90)
by the `ydyn.uz_method` parameter.

## The vertical coordinate

Yelmo uses a terrain-following, time-dependent sigma coordinate

$$
\zeta \;=\; \frac{z - b(x,y,t)}{H(x,y,t)},
\qquad \zeta = 0 \ \text{at the base},\quad \zeta = 1 \ \text{at the surface},
$$

with $b$ the ice base, $s = b + H$ the surface, and $H$ the thickness. The
$\zeta$ surfaces both **tilt** (they follow $b$ and $s$) and **move in time**
(as $b$ and $s$ evolve). Every subtlety below follows from those two facts.

The metric terms we will need are

$$
\frac{\partial \zeta}{\partial z} = \frac{1}{H},
\qquad
\frac{\partial \zeta}{\partial x}\bigg|_z
= -\frac{1}{H}\Bigl[(1-\zeta)\,\frac{\partial b}{\partial x} + \zeta\,\frac{\partial s}{\partial x}\Bigr],
$$

and identically for $y$ and $t$. These are Greve and Blatter (2009),
Eqs. 5.131–5.132, and they appear verbatim in the code as the `c_x`, `c_y`,
`c_t` factors — **but with two different normalisations**, which is the single
most important detail on this page (see [below](#the-two-c_x-are-not-the-same)).

## `uz` — the physical vertical velocity

### Continuum form

Ice is treated as incompressible, so

$$
\frac{\partial u}{\partial x}\bigg|_z + \frac{\partial v}{\partial y}\bigg|_z + \frac{\partial w}{\partial z} \;=\; 0 .
$$

Integrating upward from the base, anchored by the basal kinematic boundary
condition (Greve and Blatter, 2009, Eq. 5.31),

$$
w_b \;=\; \frac{\partial b}{\partial t} \;+\; u_b\,\frac{\partial b}{\partial x} \;+\; v_b\,\frac{\partial b}{\partial y} \;+\; \dot b ,
$$

gives

$$
w(z) \;=\; w_b \;-\; \int_b^z \left( \frac{\partial u}{\partial x}\bigg|_{z'} + \frac{\partial v}{\partial y}\bigg|_{z'} \right) \mathrm{d}z' .
$$

Note the derivatives inside the integral are at **constant $z$** — true
partial derivatives. This is the crux of the discrete implementation.

### Discrete form

The model stores $u$ and $v$ on constant-$\zeta$ layers, so a naive horizontal
difference at fixed vertical index $k$ yields $\partial u/\partial x|_\zeta$,
*not* $\partial u/\partial x|_z$. The two are related by the chain rule:

$$
\frac{\partial u}{\partial x}\bigg|_z
\;=\; \frac{\partial u}{\partial x}\bigg|_\zeta
\;+\; \frac{\partial \zeta}{\partial x}\,\frac{\partial u}{\partial \zeta} .
$$

The code applies exactly this correction. In `calc_uz_3D` (`uz_method = 2`;
`calc_uz_3D_aa` is analogous):

```fortran
c_x = -H_inv * ( (1.0-zeta_now)*dzbdx_aa + zeta_now*dzsdx_aa )   ! = dzeta/dx
c_y = -H_inv * ( (1.0-zeta_now)*dzbdy_aa + zeta_now*dzsdy_aa )   ! = dzeta/dy

dudx_aa = <du/dx at constant zeta>  +  c_x*dudz_aa               ! -> du/dx at constant z
dvdy_aa = <dv/dy at constant zeta>  +  c_y*dvdz_aa               ! -> dv/dy at constant z
```

where `dudz_aa` $= \partial u/\partial\zeta$. Note the `H_inv`: here
`c_x` $= \partial\zeta/\partial x$.

The corrected divergence is then integrated upward (Greve and Blatter, 2009,
Eq. 5.95), with $H\,\Delta\zeta = \Delta z$:

```fortran
uz(i,j,1) = dzbdt_now + f_bmb*bmb(i,j) + ux_aa*dzbdx_aa + uy_aa*dzbdy_aa   ! basal KBC
uz(i,j,k) = uz(i,j,k-1) - H_now*(zeta_ac(k)-zeta_ac(k-1))*(dudx_aa+dvdy_aa)
```

with `dzbdt_now = tpo%now%dzbdt_kin`, the kinematic rate of the column base
(`dzsdt_kin` for the surface), computed by `calc_column_kinematic_rates` from
the tendencies of the topography step:

- vertical thickness change of the column, `dHidt_vert` = advection (including
  `mb_clip`) + smb + bmb (+ relaxation). Calving, frontal melt, discharge, front
  advance and removals are lateral and do not move the column surface or base;
- grounded ice: the base follows the bedrock, `dzbdt_kin = dz_bed/dt`;
- floating ice: the column floats, `dzbdt_kin = dz_sl/dt - (rho_ice/rho_sw)*dHidt_vert`;
- blended by the grounded fraction; `dzsdt_kin = dzbdt_kin + dHidt_vert`;
- zero in partial front cells (column from `H_eff`), cells on the `H_eff`
  floor, cells that became ice covered during the step, and ice-free cells.

The bedrock and sea-level rates are the changes since the previous call of
`yelmo_update` (`ybound_update_rates`). They are written to the restart file,
so a continued run uses the same rates; they are zero on the first call after
initialisation or after a restart that does not continue the run.

**The sigma-ness is removed, not embedded.** The coordinate corrections are
applied to the *horizontal derivatives* precisely so that what gets integrated
is the true constant-$z$ divergence. What comes out, `uz`, is honest physical
$w$ in a fixed frame.

### The four methods

`ydyn.uz_method` selects among four implementations. Methods 1–3 integrate the
horizontal divergence as above and differ in how it is discretized. They use
the column of the dynamics (`H_ice_dyn`), so partial front cells are full
columns with the effective thickness. Method 4 integrates the layer mass budget
instead (see [below](#uz-flux)):

| `uz_method` | Routine | Notes |
|---|---|---|
| 1 (`uz_aa`) | `calc_uz_3D_aa` | divergence from the two faces of each cell, $(u_{i+1/2}-u_{i-1/2})/\Delta x$, the same face differences as the thickness equation |
| 2 (`uz_nodes`) | `calc_uz_3D` | Gaussian-quadrature sub-node averaging (same effective stencil as method 3) |
| 3 (`uz_jac`) | `calc_uz_3D_jac` | **default**; uses the precomputed 3D velocity Jacobian `jvel` from `calc_jacobian_vel_3D_uxyterms` ([`deformation.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/deformation.f90)), with 3D quadrature. The constant-$z$ correction is applied in the Jacobian, so `uz` is integrated from the corrected divergence directly. |
| 4 (`uz_flux`) | `calc_uz_3D_flux` | `uz_star` from the layer mass budget, with the layer fluxes closed against the applied thickness step; `uz` recovered from `uz_star` |

Methods 2 and 3 average the divergence over neighbouring cells. Along the
flow, the centred face derivatives $(u_{i+3/2}-u_{i-1/2})/(2\Delta x)$ averaged
to the cell centre give a four-face stencil, which does not see a 2$\Delta x$
oscillation of the velocity; across the flow, the corner and quadrature
averaging weights the neighbouring rows by 1/4, 1/2, 1/4. This smooths grid-scale
noise in the divergence: in Greenland (16 km, 1 kyr), ISLAND4 (16 km, 1 kyr) and
TROUGH-F17 (4 km, 5 kyr), method 3 gives a 1.5–4 times smaller surface mismatch
`uz_srf_err` and a 2–5 times smoother `uz` than method 1. However, where the
divergence changes from one cell to the next, the vertical velocity of a column
is mixed with that of its neighbours: in the A4 benchmark (plug flow with
uniform thickness along the flow), method 3 gives a
surface `uz` off by up to a factor of 5 next to a jump in the divergence, while
method 1 is exact. At an ice margin, the Jacobian of method 3 uses second-order
one-sided differences on the ice side instead of the centred ones.

### Flux-consistent vertical velocity (`uz_method = 4`) {#uz-flux}

Method 4 computes `uz_star` directly from the mass budget of each layer $k$
between $\zeta_{k-1/2}$ and $\zeta_{k+1/2}$,

$$
w^\star_{k+1/2} \;=\; w^\star_{k-1/2} \;-\; \Delta\zeta_k\,\frac{\partial H}{\partial t} \;-\; D_k,
\qquad w^\star_{1/2} = \dot b ,
$$

where $D_k = \nabla\cdot(H\,\Delta\zeta_k\,\mathbf{u}_k)$ is the flux divergence of
layer $k$ and $\partial H/\partial t$ = `dzsdt_kin` − `dzbdt_kin` is the vertical
thickness change of the column. $D_k$ is formed with the stencil of the
thickness solver (face thickness of `H_ice` upwind by the sign of the
depth-averaged velocity). Since the trapezoid weights of `zeta_aa` are the layer
thicknesses $\Delta\zeta_k$, the sum of the layer fluxes is the flux divergence
$\nabla\cdot(H\bar{\mathbf{u}})$ of the current velocity. Where the column rate is
given by the applied thickness step (`tpo%now%mask_kin` = 1), the layer fluxes
are corrected additively,

$$
D_k \;\leftarrow\; D_k + \Delta\zeta_k\Bigl(-\dot H_\mathrm{dyn} - \sum_m D_m\Bigr),
$$

so that they sum to the applied transport $-\dot H_\mathrm{dyn}$ (`dHidt_dyn`).
The correction absorbs the differences between the thickness step and the
current velocity solution (time-filtered velocity, implicit face thickness,
predictor–corrector weighting) and is distributed over the column like the
layer thickness. The surface value is then

$$
w^\star_s \;=\; -\dot a \;-\; \dot H_\mathrm{clip} \;-\; \dot H_\mathrm{relax},
$$

the kinematic condition, with the clip of negative thickness (`mb_clip`) and the
relaxation (`mb_relax`) as additional surface terms. The physical `uz` is
recovered from `uz_star` with the coordinate terms of the next section, so the
two fields stay consistent.

`mask_kin` is 0 where the column is re-derived (partial front cells, cells on
the `H_eff` floor, cells that became fully ice covered during the step) and
wherever the thickness is not advanced (`topo_fixed`, initialization). There
the layer fluxes are those of the current velocity and $\partial H/\partial t = 0$,
so the surface condition holds only if the ice is in balance. For uniform
thickness along the flow (benchmark A4), the uncorrected layer fluxes are exact
and do not mix neighbouring columns.

The vertical velocity is computed in the dynamics, after the predictor step.
With `pc_use_H_pred = True` (default), the predictor is the applied thickness
step and the closure is exact. With `pc_use_H_pred = False`, the corrected
thickness is applied, and the closure differs from it by the
predictor–corrector truncation error.

### Practical caveats

- **The surface kinematic BC is not enforced (methods 1–3).** `uz` is anchored at the *base*
  and integrated upward, and the surface condition
  $w_s = \partial s/\partial t + u_s\,\partial s/\partial x + v_s\,\partial s/\partial y - \dot a$
  is not imposed. So `uz` *exactly satisfies the basal BC*, and any mismatch
  with the surface BC accumulates as a residual at the top rather than being
  spread through the column. The mismatch is diagnosed as `uz_srf_err` =
  $w^\star_s + \dot a$ in fully ice-covered cells. It is not zero, since the
  integrated divergence ($H\,\nabla\cdot\mathbf{u}$ plus the basal and
  coordinate terms) is not the discrete flux divergence of the thickness
  equation, which uses face thicknesses and the time-filtered velocity
  (`pc_filter_vel`). It is largest at grounding lines and outlet margins
  (several m/yr in Greenland at 16 km) for all three methods. Method 4 closes
  it by construction where `mask_kin` = 1.
- **No clamps.** `uz` and `uz_star` are not limited; in fast outlets $w$ reaches
  tens of m/yr. Values below `TOL_UNDERFLOW` are zeroed.
- **Ice-free points** get `uz = dzbdt - max(smb,0)` and `uz_star = uz`.

## `uz_star` — the sigma-relative advective velocity

### Why it exists

Consider advecting any scalar $X$ (enthalpy, temperature, age, a tracer). The
material derivative is a physical statement:

$$
\frac{\mathrm{D}X}{\mathrm{D}t}
= \frac{\partial X}{\partial t}\bigg|_z
+ u\,\frac{\partial X}{\partial x}\bigg|_z
+ v\,\frac{\partial X}{\partial y}\bigg|_z
+ w\,\frac{\partial X}{\partial z} .
$$

But a sigma-coordinate solver cannot evaluate $\partial X/\partial x|_z$
directly — differencing $X$ across neighbouring columns at fixed layer index
$k$ gives $\partial X/\partial x|_\zeta$. Substituting the chain rule
$\partial X/\partial f|_z = \partial X/\partial f|_\zeta + (\partial\zeta/\partial f)\,\partial X/\partial\zeta$
for $f \in \{t,x,y\}$, and $\partial X/\partial z = H^{-1}\partial X/\partial\zeta$,
all the leftover pieces collapse into a single vertical term:

$$
\frac{\mathrm{D}X}{\mathrm{D}t}
= \frac{\partial X}{\partial t}\bigg|_\zeta
+ u\,\frac{\partial X}{\partial x}\bigg|_\zeta
+ v\,\frac{\partial X}{\partial y}\bigg|_\zeta
+ \underbrace{\frac{\mathrm{D}\zeta}{\mathrm{D}t}}_{\textstyle \equiv\, w^\star/H}\,\frac{\partial X}{\partial \zeta} ,
$$

where

$$
\frac{\mathrm{D}\zeta}{\mathrm{D}t}
= \frac{\partial\zeta}{\partial t}
+ u\,\frac{\partial\zeta}{\partial x}
+ v\,\frac{\partial\zeta}{\partial y}
+ \frac{w}{H} .
$$

Multiplying through by $H$ defines

$$
\boxed{\;
w^\star \;\equiv\; H\,\frac{\mathrm{D}\zeta}{\mathrm{D}t}
\;=\; w
\;+\; u\,\underbrace{H\frac{\partial\zeta}{\partial x}}_{c_x}
\;+\; v\,\underbrace{H\frac{\partial\zeta}{\partial y}}_{c_y}
\;+\; \underbrace{H\frac{\partial\zeta}{\partial t}}_{c_t}
\;}
$$

which is Greve and Blatter (2009), Eq. 5.148, and exactly the code:

```fortran
c_x = -( (1.0-zeta_now)*dzbdx_aa  + zeta_now*dzsdx_aa )   ! = H * dzeta/dx
c_y = -( (1.0-zeta_now)*dzbdy_aa  + zeta_now*dzsdy_aa )   ! = H * dzeta/dy
c_t = -( (1.0-zeta_now)*dzbdt_now + zeta_now*dzsdt_now )  ! = H * dzeta/dt

uz_star(i,j,k) = uz(i,j,k) + ux_aa*c_x + uy_aa*c_y + c_t
```

The vertical advection term is then $w^\star\,\partial X/\partial z$, with
$\partial X/\partial z$ formed as $\Delta X / (H\,\Delta\zeta)$. This is why
the factor $H^{-1}$ is **deliberately not applied** when `uz_star` is built —
it is supplied by each consumer's advection step. `uz_star` therefore has
units of $\mathrm{m\,a^{-1}}$.

### What it represents

$w^\star$ is the parcel's velocity **relative to the moving sigma surfaces**,
expressed with the dimensions of a vertical velocity. Concretely:

- $w^\star = 0$ $\iff$ the parcel stays on its $\zeta$ layer (it may still be
  moving vertically in physical space, if the layer itself is moving).
- $w^\star$ is what tells you the rate at which ice *crosses* model layers.

This is the moving-mesh (ALE) mesh-relative velocity, and the direct analogue
of $\omega$ in atmospheric sigma coordinates or the diasurface velocity in
ocean models. It is a genuine geometric/kinematic quantity — **not** a
numerical artifact. It is, however, *coordinate-relative* rather than
*frame-relative*: it only means anything once you have chosen a
terrain-following vertical coordinate. The presence of $c_t$ makes this vivid —
that term exists purely because the coordinate surfaces move in time, which has
nothing to do with numerics.

## The two `c_x` are not the same {#the-two-c_x-are-not-the-same}

In `calc_uz_3D`, `c_x` is defined twice with different normalisations within
the *same subroutine*. This trips people up:

| Context | Code | Value | Multiplies |
|---|---|---|---|
| `uz` loop | `c_x = -H_inv * (...)` | $\partial\zeta/\partial x$ | $\partial u/\partial\zeta$ (a **velocity derivative**) |
| `uz_star` loop | `c_x = -(...)` | $H\,\partial\zeta/\partial x$ | $u$ (the **velocity itself**) |

In `calc_uz_3D_aa` and `calc_uz_3D_flux`, the `uz_star` loop includes the
velocity in `c_x` (`c_x = -ux_aa*(...)`), so there `c_x` is the product
$u\,H\,\partial\zeta/\partial x$ itself.

The corresponding physical distinction is worth stating outright, because it is
the most common misconception about these two fields:

- The corrections that build **`uz`** involve **derivatives of $u$ and $v$**
  ($\partial u/\partial\zeta$, $\partial v/\partial\zeta$). They arise from
  converting constant-$\zeta$ velocity gradients into constant-$z$ ones so that
  incompressibility can be integrated.
- The corrections that build **`uz_star`** involve **$u$ and $v$ themselves,
  undifferentiated**, multiplied by *geometry slopes* — plus $c_t$ from the
  coordinate's time dependence. They arise from evaluating the *scalar's*
  horizontal gradient at constant $\zeta$ instead of constant $z$.

Different corrections, different origins, same routine.

## Which one do I use?

The rule is determined entirely by **what coordinate your derivative or
integration is taken in** — not by whether the code is "numerical".

| Task | Use |
|---|---|
| Eulerian advection of a scalar with horizontal gradients at constant $\zeta$ | **`uz_star`** |
| Lagrangian particle trajectory integrating physical $z$: $\mathrm{d}z/\mathrm{d}t = w$ | **`uz`** |
| Lagrangian particle tracked by layer, integrating $\mathrm{d}\zeta/\mathrm{d}t$ | **`uz_star`**$/H$ |
| Reporting/diagnosing the physical vertical velocity, coupling to an external model | **`uz`** |
| Vertical CFL constraint for the 3D advection above | **`uz_star`** |

Put briefly: **integrate $z$ → `uz`; integrate $\zeta$ → `uz_star`.** Anything
reasoning in physical space wants `uz`; anything stepping through the model's
own layers wants `uz_star`.

A worked example of the second row: to hand $(u, v, w)$ to an offline
Lagrangian particle tracer that integrates
$\dot x = u,\ \dot y = v,\ \dot z = w$, pass `ux`, `uy`, `uz` — never
`uz_star`.

## Where each field is used in the code

| Consumer | Field | Location |
|---|---|---|
| Enthalpy/temperature vertical advection | `uz_star` | [`yelmo_thermodynamics.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_thermodynamics.f90) → `calc_ytherm_enthalpy_3D` |
| Eulerian deposition-time tracer (`t_dep`) | `uz_star` | [`yelmo_tracers.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_tracers.f90) → `calc_tracer_3D` ([`ice_tracer.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/ice_tracer.f90)) |
| Lagrangian particle tracer (`tracer` backend) | `uz` | [`yelmo_tracers.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_tracers.f90) |
| `enh_bnd` tracer advection (`*-tracer` enhancement methods) | `uz_star` | [`yelmo_material.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_material.f90) → `calc_tracer_3D` |
| Velocity Jacobian `jvel` | `uz` | [`deformation.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/deformation.f90) → `calc_jacobian_vel_3D_uzterms` |
| `uz_b`, `uz_s` diagnostics | `uz` | `yelmo_dynamics.f90` |
| Output / restart / C API | both | `yelmo_io.f90`, `yelmo_c_api.f90` |

Both fields are **diagnostic**: they are recomputed from the velocity solution
every dynamics step and carry no prognostic state, even though both are written
to the restart file.

::: {.callout-note}
## 3D CFL timestep

There is no 3D advective CFL limit: the adaptive time step
([Time stepping](timestepping.md)) only applies the 2D (depth-averaged)
advective CFL. The vertical advection of the thermodynamics is implicit, and
its horizontal advection is sub-cycled internally. A 3D constraint would have
to be driven by `uz_star`, since that is the velocity at which scalars actually
cross model layers.
:::

## References

- Greve, R. and Blatter, H. (2009). *Dynamics of Ice Sheets and Glaciers.*
  Springer. — Eq. 5.31 (basal kinematic BC), Eq. 5.95 (upward integration of
  $w$), Eqs. 5.131–5.132 (sigma metric terms), Eq. 5.148 ($w^\star$).
- [Robinson et al. (2020)](https://doi.org/10.5194/gmd-13-2805-2020) — Yelmo
  model description.
- The Glimmer ice-sheet model
  [documentation](https://www.geos.ed.ac.uk/~mhagdorn/glide/glide-doc/glimmer_htmlse9.html),
  whose algorithm the `uz` integration follows.
