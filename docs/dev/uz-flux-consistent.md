# Flux-consistent vertical velocity (design)

*Status: design for discussion, 2026-10-05. Nothing here is implemented.*

## Problem

The vertical velocity in Yelmo is not consistent with the discrete thickness equation. The kinematic surface condition requires the sigma-relative vertical velocity at the surface to equal the negative surface mass balance, w*_s = −SMB. Yelmo diagnoses the mismatch as `uz_srf_err` = w*_s + SMB, and it is far from zero:

| Case | `uz_method = 1`, rms | `uz_method = 3`, rms |
|---|---|---|
| Greenland, 16 km, 1 kyr | 5.4 m a⁻¹ | 3.2 m a⁻¹ |
| ISLAND4, 16 km, 1 kyr | 0.69 m a⁻¹ | 0.17 m a⁻¹ |
| TROUGH-F17, 4 km, 5 kyr | 0.063 m a⁻¹ | 0.044 m a⁻¹ |

In all three cases the mismatch is concentrated at grounding lines and outlet margins, where it reaches several m a⁻¹ to tens of m a⁻¹. A mismatch of this size is comparable to the vertical velocity itself, so the vertical advection of heat (and of age and tracers) is wrong there. The choice of method also changes the thermal state substantially: in ISLAND4 the temperate basal fraction after 1 kyr is 0.17 with method 1 and 0.08 with method 3.

The two methods fail in different ways. Method 1 integrates H ∇·u with the face differences of the thickness equation and is exact when the thickness is uniform along the flow (benchmark A4), but it passes grid-scale noise of the divergence to w. Method 3 averages the divergence over neighbouring cells, which smooths w and lowers the mismatch in realistic flow, but it mixes the vertical velocity of a column with that of its neighbours, with errors of up to a factor of 5 in A4.

## Cause

Both methods integrate the continuum identity

$$w_s = w_b - \int_b^s \nabla\cdot\mathbf{u}\,dz, \qquad \int_b^s \nabla\cdot\mathbf{u}\,dz = \nabla\cdot\mathbf{q} - \mathbf{u}_s\cdot\nabla s + \mathbf{u}_b\cdot\nabla b,$$

with a discretization of their own. The thickness equation uses another one: first-order upwind face thicknesses, implicit in time, with the mean of the current and previous velocity solutions (`pc_filter_vel`), combined over the predictor and corrector steps. The discrete versions of H ∇·u + u·∇H and of ∇·(Hu) therefore differ by terms of order Δx ∇²(uH), which are largest where H, u or the bed change abruptly: at grounding lines, margins and outlets. No choice of derivative stencil in w removes this, since the thickness equation itself is not a centred scheme.

## Proposed formulation

The sigma-relative vertical velocity w* = H dζ/dt (`uz_star`) is the quantity that the thermodynamics and the tracer advection need. In sigma coordinates it follows from the mass budget of each layer k between the levels ζ_{k−1/2} and ζ_{k+1/2}:

$$\frac{\partial}{\partial t}\left(H\,\Delta\zeta_k\right) + \nabla\cdot\left(H\,\Delta\zeta_k\,\mathbf{u}_k\right) + w^*_{k+1/2} - w^*_{k-1/2} = 0,$$

with w*_{1/2} = −BMB at the base. This is the form used for the diasigma velocity in terrain-following ocean models. Summed over all layers, it gives w*_s = −(∂H/∂t + ∇·(H ū)) − BMB, which equals −SMB exactly when the same ∂H/∂t and the same flux divergence are used as in the thickness equation.

We propose to compute w* with the flux divergence of the thickness equation itself:

1. **Total divergence.** D = ∇·(H ū) is taken from the thickness update that was actually applied (the advective tendency, after the predictor–corrector weighting), together with the vertical thickness change ∂H/∂t = `dHidt_vert` that already defines `dzsdt_kin` and `dzbdt_kin`. These are the terms the surface condition contains, so the closure holds by construction.
2. **Vertical distribution.** The layer divergences D_k^raw = ∇·(H Δζ_k u_k) are computed with the same stencil as the thickness solver (upwind face thicknesses, face velocities of layer k with the same filtering). They are corrected additively so that their sum equals D,

   $$D_k = D_k^{\mathrm{raw}} + \Delta\zeta_k\left(D - \sum_m D_m^{\mathrm{raw}}\right),$$

   which keeps the vertical structure of the shear profile and needs no division by a possibly small total.
3. **Integration.** w*_{k+1/2} = w*_{k−1/2} − Δζ_k ∂H/∂t − D_k, upward from the base.
4. **Physical vertical velocity.** w (`uz`) is recovered from w* with the coordinate terms that are already used to form w* from w (u·∇z at constant ζ and ∂z/∂t), so that the two fields stay consistent. w enters the strain-rate tensor (∂w/∂z) and the diagnostics.

For plug flow with uniform thickness along the flow (A4), D_k = Δζ_k D and w* = −SMB ζ, the exact solution. Where the velocity carries grid-scale noise, w* carries the same noise as the thickness transport. This is the consistent result: any smoothing then belongs in the velocity or thickness solution, not in w*.

## Points to settle in the implementation

- **Time level.** w* must be formed from the tendencies of the thickness step that the thermodynamics follows. At present the dynamics computes `uz` before the corrector and advance steps, with `dzsdt_kin` from the previous advance. The order of the calls in `yelmo_update` and the stored tendencies (`dHidt_dyn_raw`, `dHidt_vert`) need to be checked, so that D and ∂H/∂t refer to the same step.
- **Depth-averaged velocity.** The correction in step 2 absorbs any difference between ū and the layer mean Σ Δζ_k u_k, but a large difference would distort the vertical profile. Both should be computed with the same vertical quadrature.
- **Partial front cells and the H_eff floor.** The kinematic rates are set to zero in partial front cells, cells on the H_eff floor and newly ice-covered cells, since their column is re-derived. The same cells need a defined w* (e.g., the current treatment), and `uz_srf_err` is evaluated only in fully ice-covered cells.
- **Floating ice.** The base of floating ice moves with the thickness (`dzbdt_kin`); the layer budget includes this through ∂H/∂t and w*_{1/2} = −BMB.
- **Other thickness solvers.** The explicit and second-order upwind solvers use other stencils. The layer divergence must use the stencil of the solver in use, or the method is restricted to `impl-lis` and `impl-upwind`.

## Implementation and tests

The formulation would be added as `ydyn.uz_method = 4` (`calc_uz_3D_flux` in `velocity_general.f90`), next to the existing methods, which allows a direct comparison before a change of the default.

Tests, in order:

1. A4 with one row per parameter set and all columns compared: exact to round-off at the surface, temperature as for the three-row design.
2. `uz_srf_err` in ISLAND4, TROUGH-F17 and Greenland: round-off in fully ice-covered cells.
3. Roughness of w and of the basal temperature against methods 1 and 3, and the ISLAND4 D4 symmetry in double precision.
4. EISMINT EXPF symmetry gate (thermomechanical feedback) and the regression set (CalvingMIP, MISMIP+, MISMIP3D, TROUGH-F17).
5. Effect on the Greenland and Antarctica thermal state and on the timing (the extra cost is one pass over the layers with the upwind stencil, which should be small).

## Alternatives considered

- **Method 3 with a compact divergence.** Replacing the averaged divergence of method 3 by the face differences of method 1 removes the mixing between neighbouring cells, but keeps the mismatch at grounding lines and margins, which comes from the different discretization of the thickness equation.
- **Correcting w to the surface condition.** Scaling or shifting w in each column so that it meets the surface condition (the commented-out line in `calc_uz_3D_aa`) enforces the closure, but distributes the error through the column without a physical basis and leaves the layer mass budget inconsistent.
