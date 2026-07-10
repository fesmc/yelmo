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

The code applies exactly this correction. In `calc_uz_3D`:

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

with `dzbdt_now = dzsdt_now - dhdt_now` (the base rate is deduced from the
surface and thickness rates, so no explicit bedrock-uplift input is needed,
and the expression is valid for grounded *and* floating ice).

**The sigma-ness is removed, not embedded.** The coordinate corrections are
applied to the *horizontal derivatives* precisely so that what gets integrated
is the true constant-$z$ divergence. What comes out, `uz`, is honest physical
$w$ in a fixed frame.

### The three methods

`ydyn.uz_method` selects among three implementations that differ only in how
the horizontal derivatives are evaluated — all three produce the same
quantity:

| `uz_method` | Routine | Notes |
|---|---|---|
| 1 (`uz_aa`) | `calc_uz_3D_aa` | simplest; plain finite differences on `aa`-nodes |
| 2 (`uz_nodes`) | `calc_uz_3D` | Gaussian-quadrature sub-node averaging |
| 3 (`uz_jac`) | `calc_uz_3D_jac` | **default**; uses the precomputed 3D velocity Jacobian `jvel` from `calc_jacobian_vel_3D_uxyterms` ([`deformation.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/deformation.f90)). Most stable and most correct. |

### Practical caveats

- **The surface kinematic BC is not enforced.** `uz` is anchored at the *base*
  and integrated upward. The corresponding surface condition
  $w_s = \partial s/\partial t + u_s\,\partial s/\partial x + v_s\,\partial s/\partial y - \dot a$
  is computed but the redistribution correction is commented out in all three
  routines. So `uz` is *exactly divergence-consistent and exactly satisfies the
  basal BC*; any mismatch with the surface BC accumulates as a residual at the
  top rather than being spread through the column. "Consistent with mass
  conservation" should be read in that precise sense.
- **Stability clamps.** The basal value is floored at `uz_min = -10 m/yr`, and
  `calc_uz_3D_jac` additionally clamps $|w| \le$ `uz_lim = 10 m/yr` throughout.
  Values below `TOL_UNDERFLOW` are zeroed. These are numerical guards that only
  bite when something is already going wrong, but they do make `uz` mildly
  non-physical near the bed in those cases.
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

Within the *same subroutine*, `c_x` is defined twice with different
normalisations. This trips people up:

| Context | Code | Value | Multiplies |
|---|---|---|---|
| `uz` loop | `c_x = -H_inv * (...)` | $\partial\zeta/\partial x$ | $\partial u/\partial\zeta$ (a **velocity derivative**) |
| `uz_star` loop | `c_x = -(...)` | $H\,\partial\zeta/\partial x$ | $u$ (the **velocity itself**) |

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
`uz_star`. (Note the basal `uz_min` clamp above if trajectories approach the
bed.)

## Where each field is used in the code

| Consumer | Field | Location |
|---|---|---|
| Enthalpy/temperature vertical advection | `uz_star` | [`yelmo_thermodynamics.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_thermodynamics.f90) → `calc_ytherm_enthalpy_3D` |
| Age (`dep_time`) and `enh_bnd` tracer advection | `uz_star` | [`yelmo_material.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_material.f90) → `calc_tracer_3D` ([`ice_tracer.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/ice_tracer.f90)) |
| Velocity Jacobian `jvel` | `uz` | [`deformation.f90`](https://github.com/fesmc/yelmo/blob/main/src/physics/deformation.f90) → `calc_jacobian_vel_3D_*` |
| `uz_b`, `uz_s` diagnostics | `uz` | `yelmo_dynamics.f90` |
| Output / restart / C API | both | `yelmo_io.f90`, `yelmo_c_api.f90` |

Both fields are **diagnostic**: they are recomputed from the velocity solution
every dynamics step and carry no prognostic state, even though both are written
to the restart file.

::: {.callout-note}
## 3D CFL timestep

`set_adaptive_timestep` in
[`yelmo_timesteps.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_timesteps.f90)
currently takes `uz`, but the 3D advective CFL is **disabled**: `dt_adv3D` is
hardcoded to `1000.0` and both `calc_adv3D_timestep*` helpers are commented
out. Should that constraint ever be reactivated, it must be driven by
`uz_star`, since that is the velocity at which scalars actually cross model
layers.
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
