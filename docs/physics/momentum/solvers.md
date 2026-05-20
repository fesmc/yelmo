# Numerical solution of DIVA / SSA

The DIVA and SSA momentum balances share the same 2D depth-integrated
form: a nonlinear elliptic problem for the depth-averaged horizontal
velocity $\bar{\mathbf u} = (\bar u, \bar v)$, with coefficients
$(\mu, \beta_\mathrm{eff})$ that depend on the solution. Yelmo
discretises this problem on the Arakawa C-grid (the **ac-grid**) and
linearises it by Picard iteration, freezing
$(\mu, \beta_\mathrm{eff}, H)$ at the previous iterate, solving a linear
system for the new $\bar{\mathbf u}$, relaxing, and repeating until
convergence.

Two solvers are provided for the linear step. Both produce the same
solution at interior cells (the two systems are related by a sign and a
cell-area factor), but they differ in the structure of the assembled
matrix and in how boundary conditions are imposed. The choice is made
at runtime via the parameter `ydyn.ssa_solver = "residual" | "energy"`,
which the DIVA solver dispatches on.

## Discretisation conventions

Yelmo uses the Arakawa C / staggered AC grid:

- **aa-nodes** (cell centres): scalars — $H$, $s$, $\bar\mu$, $\bar\mu\,H$, $\beta$ (before staggering).
- **acx-nodes** (right face of each aa-cell): $\bar u$, $\tau_{d,x}$, $\beta_{\mathrm{eff},x}$.
- **acy-nodes** (top face): $\bar v$, $\tau_{d,y}$, $\beta_{\mathrm{eff},y}$.
- **ab-nodes** (cell corners): cross-coupling viscosity
  $(\bar\mu H)^{\mathrm{ab}}$, obtained by averaging
  $\bar\mu\,H$ from its four neighbouring aa-cells
  (`stagger_visc_aa_ab`).

The two unknowns per cell ($\bar u$ at the right face and $\bar v$ at the
top face) are interleaved into a single state vector
$\mathbf x = [\,\bar u_{1}, \bar v_{1}, \bar u_{2}, \bar v_{2}, \dots\,]$,
so the assembled linear system has dimension $2 N_\mathrm{cells}$. The
matrix is stored in CSR format and solved with
[LIS](http://www.ssisc.org/lis/) (Library of Iterative Solvers for
Linear Systems); the iterative method and preconditioner are configured
at runtime via `ydyn.ssa_lis_opt`.

Per-row solver masks (`ssa_mask_acx`, `ssa_mask_acy`) classify each
ac-node as one of:

- `0`: Dirichlet zero velocity (frozen bed, edge of domain, etc.);
- `-1`: Dirichlet prescribed velocity (e.g. observed input);
- `1`: solve (interior);
- `3`: lateral / calving-front Neumann row driven by the
  depth-integrated lateral stress $\tau_{l,\mathrm{int}}$.

The two assemblers share the same argument list so they are drop-in
interchangeable from the Picard loop in `calc_velocity_diva`.

## Solver A — residual form (Yelmo v1)

The residual assembler ([`solver_ssa_ac.f90`][resid_src]) builds a
non-symmetric matrix $A_\mathrm{res}$ and right-hand side
$\mathbf b_\mathrm{res}$ directly from the strong form of the SSA PDE.
Writing $N \equiv \bar\mu\,H$ for brevity, the interior $\bar u$ row
is a 5-point stencil in $\bar u$ plus cross-coupling terms in $\bar v$
on the neighbouring acy-faces, balanced against the basal drag and the
driving stress. Schematically (cell $(i,j)$, with $N^{\mathrm{aa}}$ on
aa-nodes and $N^{\mathrm{ab}}$ on ab-nodes):

$$
\sum_{(i',j')}\alpha^{x}_{i',j'}\,\bar u_{i',j'}
\;+\; \sum_{(i',j')}\gamma^{xy}_{i',j'}\,\bar v_{i',j'}
\;-\; \beta_\mathrm{eff}\,\bar u_{i,j}
\;=\; \tau_{d,x}(i,j).
$$

The $\bar u$ stencil coefficients carry the membrane stretching and
shearing,

$$\alpha^{x}_{i-1,j} \;=\; \tfrac{4}{\Delta x^{2}}\,N^{\mathrm{aa}}_{i,j},$$

$$\alpha^{x}_{i+1,j} \;=\; \tfrac{4}{\Delta x^{2}}\,N^{\mathrm{aa}}_{i+1,j},$$

$$\alpha^{x}_{i,j-1} \;=\; \tfrac{1}{\Delta y^{2}}\,N^{\mathrm{ab}}_{i,j-1},$$

$$\alpha^{x}_{i,j+1} \;=\; \tfrac{1}{\Delta y^{2}}\,N^{\mathrm{ab}}_{i,j},$$

with the centre coefficient $\alpha^{x}_{i,j}$ being minus the sum of
its four neighbours. The $\gamma^{xy}$ coefficients couple the row to
four neighbouring $\bar v$-faces via the $\partial \bar v/\partial x$
and $\partial \bar v/\partial y$ cross-terms in the membrane stress.
The $\bar v$ row has the analogous structure.

The right-hand side is the driving stress itself (no cell-area factor).
The system
$A_\mathrm{res}\,\mathbf x = \mathbf b_\mathrm{res}$ is non-symmetric,
and is solved with a Krylov method such as BiCGStab plus an algebraic
preconditioner.

- Lateral / calving-front rows substitute the membrane-stress balance
  with the prescribed depth-integrated lateral stress
  $\tau_{l,\mathrm{int}}$ in a row-specific stencil.
- Dirichlet rows replace the equation by $\bar u_{i,j} = u^*$ with the
  corresponding column kept in place (non-symmetric).
- Free-slip, no-slip and periodic conditions are applied per side via
  the boundary-code helper `get_neighbor_indices_bc_codes`.

This is the formulation inherited from Yelmo v1 and is the legacy
default.

## Solver B — energy form (new)

The energy assembler ([`solver_ssa_ac_energy.f90`][energy_src]) builds
the Hessian of a discrete energy functional and solves
$K\,\mathbf x = \mathbf b$ for the velocity that minimises that energy.
With $(\mu, \beta, H)$ frozen during each Picard step the energy is
quadratic, so $K = \frac{\partial^{2} W}{\partial \mathbf x^{\,2}}$ is
symmetric positive (semi-)definite and the linear step can use a
symmetric Krylov method — CG with an AMG preconditioner — in place of
BiCGStab.

### Energy density

The continuum energy density underlying the SSA momentum balance is the
sum of a membrane (deformation) term, a basal-drag term, and a
gravitational potential-energy term:

$$
\begin{aligned}
W \;=\;
\bar\mu\,H\,&\biggl(
       2\!\left(\frac{\partial \bar u}{\partial x}\right)^{\!2}
     + 2\!\left(\frac{\partial \bar v}{\partial y}\right)^{\!2}
     + 2\,\frac{\partial \bar u}{\partial x}\,\frac{\partial \bar v}{\partial y}
     + \tfrac{1}{2}\!\left(\frac{\partial \bar u}{\partial y} + \frac{\partial \bar v}{\partial x}\right)^{\!2}
\,\biggr) \\[4pt]
&+\; \tfrac{1}{2}\,\beta\,(\bar u^{\,2} + \bar v^{\,2})
\;+\; \rho_i\,g\,H\,\!\left(\bar u\,\frac{\partial s}{\partial x} + \bar v\,\frac{\partial s}{\partial y}\right).
\end{aligned}
$$

The first line is $2\,\bar\mu\,H\,\dot{\bar\varepsilon}_{ij}\,\dot{\bar\varepsilon}_{ij}$
written out for the depth-averaged horizontal strain rates. Stationarity
of the integral $\mathcal W = \int W \,\mathrm dx\,\mathrm dy$ with
respect to $(\bar u, \bar v)$ reproduces exactly the SSA / DIVA strong
form, so any critical point of $\mathcal W$ is a solution of the
momentum balance.

### Discrete assembly

Yelmo's discrete energy is the cell-by-cell evaluation of $W$ on the
C-grid, with the derivatives
$\frac{\partial \bar u}{\partial x}, \frac{\partial \bar v}{\partial y},
\frac{\partial \bar u}{\partial y}, \frac{\partial \bar v}{\partial x}$
expressed as the natural finite differences between adjacent ac-nodes.
The Hessian $K$ then has the same stencil graph as the residual matrix
but is symmetric in $(\bar u, \bar v)$. At inner cells the two
formulations are related by an exact algebraic identity (documented in
the header of the energy assembler):

$$
K_\mathrm{inner} \;=\; -\,A_\mathrm{res, inner}\cdot \Delta x\,\Delta y,
\qquad
\mathbf b_\mathrm{inner} \;=\; -\,\boldsymbol\tau_d \cdot \Delta x\,\Delta y.
$$

So the energy formulation is, at interior cells, the residual
formulation rescaled by the cell area and a sign — the physical
solution at interior cells is identical to machine precision. Where the
two solvers differ is at boundaries:

- **Lateral / calving-front BC**: the front stress enters the energy as
  a boundary-work term $\pm\,\tau_{l,\mathrm{int}}\,\Delta y$ on the
  RHS, with the sign determined by the outward normal. This is the
  variational form of the Neumann condition and is symmetric by
  construction.
- **Dirichlet rows**: prescribed values are imposed by **static
  condensation** — the prescribed column is multiplied by the known
  velocity and moved to the RHS — instead of by row replacement. This
  preserves the symmetry of $K$ and so allows CG / AMG to be used
  without spoiling the SPD structure.

The viscosity staggering aa $\to$ ab is the same routine
(`stagger_visc_aa_ab`) used by the residual assembler.

### Why bother?

Two practical advantages flow from the SPD structure:

1. **Linear solver choice**: CG with AMG is typically faster and more
   robust than BiCGStab + ILU for large, well-conditioned SPD systems,
   and converges monotonically in the energy norm.
2. **Physical interpretability and discrete consistency**: every term
   in $K$ and $\mathbf b$ corresponds to a contribution to a discrete
   energy. Boundary conditions that are natural for the continuum
   functional (e.g. Neumann front stress) become natural for the
   discrete one. This makes it easier to add new physics — e.g.
   alternative friction laws, additional body forces — in a way that
   provably preserves the variational structure.

The Picard loop, the viscosity update, the F-integral closure and the
basal-stress and 3D velocity diagnostics are identical between the two
solvers: only the linear-system assembly differs.

[resid_src]: https://github.com/fesmc/yelmo/blob/main/src/physics/solver_ssa_ac.f90
[energy_src]: https://github.com/fesmc/yelmo/blob/main/src/physics/solver_ssa_ac_energy.f90
