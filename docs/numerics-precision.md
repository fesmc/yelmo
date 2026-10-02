# Numerical precision

## Current state

Yelmo and its Fortran dependencies are compiled in single precision. The
working precision `wp` is hard-coded in each package:

| Package | File | Setting |
|---|---|---|
| yelmo | `src/yelmo_defs.f90` | `wp = sp` |
| fesm-utils | `src/precision.f90` | `wp = sp` |
| FastHydrology | `src/fast_hydrology.f90`, `src/bucket.f90`, `src/closures.f90` | `wp = sp` |
| tracer | `src/tracer_precision.f90` | `wp = sp`, `prec_time = sp` |
| elsa | `src/elsa_precision.f90` | `wp = dp` |

There is no build switch to change this. Interfaces between packages pass
`wp` arrays, so all packages above must agree on `wp`.

Some solvers already work internally in double precision:

- Lis (`solver_linear.F90`) solves in double precision; the matrix is
  assembled in `wp` and converted.
- FastHydrology's K24 transport copies its `wp` inputs to `dp`, solves, and
  converts the result back (`hydro_update`, `TRANSPORT_K24` branch). This is
  the template for the plan below.

## Why it matters: symmetry regression check

The EISMINT-moving symmetry check (`symmetry_check` in the `&ctrl` group,
`src/yelmo_symmetry.f90`) compares the solution with its mirror images,
normalised by the maximum ice thickness, against the tolerance `symmetry_tol`.
EISMINT-moving and EXPA check the final state at 1e-3. EXPF's symmetric state
is thermomechanically unstable and shows transient symmetry-breaking bursts
with either thermal solver, so it checks the mean over the last 20 kyr
(`symmetry_avg_time`) at 2e-2; persistent asymmetry ("spokes", ~0.1) still
fails.

With the SIA solver the check passes in single precision. With
`ydyn.solver = "diva"` it fails (Linf/Hmax ≈ 1.7e-3), although every DIVA
operation is symmetric: with exactly symmetric inputs the solver output is
symmetric to 4e-7 (sp) and 2e-13 (dp). The cause is dynamical. Thin grounded
frontal cells switch on and off every 100-300 yr (grounded cells have
`f_ice = 1` whenever `H_ice > 0`), and the eight mirror-image cells slowly
drift out of phase with an e-folding time of ~1.7 kyr. Single-precision
round-off (~1e-7) grows to a phase slip after ~18 kyr. In double precision
the asymmetry at 25 kyr is ~4e-12.

The symmetry check is therefore meaningful for DIVA only in double
precision.

## Building in double precision (for regression checks)

Build a separate copy; do not modify the shared checkouts.

1. Copy fesm-utils, FastHydrology, tracer and yelmo into a scratch tree and
   link them to each other there (`fesm-utils`, `FastHydrology`, `tracer`
   symlinks inside each package, as in a normal setup). elsa can be linked
   as is.
2. Set `wp = dp` in every file listed in the table above, and also
   `prec_time = dp` in tracer.
3. Add `-fdefault-real-8 -fdefault-double-8` to `FFLAGS_BASE` in the yelmo
   `Makefile`. This is currently required because some calls pass
   default-real literals to `wp` dummy arguments (e.g. in
   `solver_ssa_ac.f90`), which fail to compile with `wp = dp` otherwise.
4. Build in dependency order: `make fesmutils-static`, `make
   fasthydro-static`, `make tracer-static`, then `make yelmo-static
   benchmarks` in yelmo.

Run the check with, for example:

```bash
runme -r -e benchmarks -o output/sym-diva-dp -n par/yelmo_EISMINT_moving.nml -p ctrl.time_end=25e3 ydyn.solver=diva
```

## Plan: double-precision internals, `wp` interfaces

Rather than switching the whole model to double precision, the plan is to
run the precision-sensitive internal calculations in `dp` while all
interfaces (derived types, module arguments, I/O) stay in `wp`, following the
K24 pattern:

```fortran
! inputs wp -> dp
H_dp = real(H_ice, dp)
...
call solve_something_dp(x_dp, H_dp, ...)
! result dp -> wp
x = real(x_dp, wp)
```

Candidates, in order of expected benefit:

- Picard iteration of the SSA/DIVA momentum balance, including viscosity and
  effective-drag updates and the convergence measure.
- DIVA vertical integrals (F-integrals) and basal velocity.
- Vertical column solvers in thermodynamics (tridiagonal enthalpy and
  temperature solves).
- Implicit advection matrix assembly (the Lis solve is already `dp`).

Prerequisites:

- Make yelmo `wp`-clean, so that it compiles with `wp = dp` without
  `-fdefault-real-8`. Literals passed to `wp` arguments need a `_wp` kind.
- Add a double-precision build of the symmetry check to the regression
  tests.
