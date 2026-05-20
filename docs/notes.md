# Notes

## `timeout` module

Example parameters for using the `timeout` module. The name of the section
should be specified when calling `timeout_init`.

```fortran
&tm_1D
    method          = "file"            ! "const", "file", "times"
    dt              = 1.0
    file            = "input/timeout_ramp_100kyr.txt"
    times           = -10, -5, 0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55
/

```

```fortran
! Get output times
call timeout_init(tm_1D,path_par,"tm_1D","small", time_init,time_end)
```

If we are loading the desired output times from a file, the format is one
time per line, or a range of times using the format `t0:dt:t1`:

```bash
0:10:200
200:20:300
300:50:500
500:100:1000
1000:200:5000
5000:500:10000
10e3:1e3:20e3
20e3:2e3:200e3
200e3:5e3:1e6
```

Duplicate times will be removed, as well as times outside of the range of
`time_init` and `time_end`. In this way, once `timeout_init` is called,
we know how many timesteps of output will be generated. This can help
confirm that we designed the experiment well, and how much data to expect.

Then during the timeloop, simply use the function `timeout_check` to
determine if the current time should be written to output:

```fortran
if (timeout_check(tm_1D,time)) then  
    !call write_step_1D_combined(yelmo1,hyst1,file1D_hyst,time=time)
end if 
```

## `timing` module

All calls to the intrinsic routine `cpu_time()` have been replaced by timing
calculations performed in the new `timing` module. This has the benefit
of ensuring timing will work properly for parallel and serial programs,
and allows us to keep track of multiple timing objectives with one
simple object and subroutine.

The control of a timing object is handled via `timer_step`:

First, initialize and reset the `timer` object:

```fortran
call timer_step(tmrs,comp=-1)
```

Then, e.g., within the timeloop, get the timing for isostasy calls and for Yelmo calls:

```fortran
call timer_step(tmrs,comp=0) 

! == ISOSTASY ==========================================================
call isos_update(isos1,yelmo1%tpo%now%H_ice,yelmo1%bnd%z_sl,time,yelmo1%bnd%dzbdt_corr) 
yelmo1%bnd%z_bed = isos1%now%z_bed

call timer_step(tmrs,comp=1,time_mod=[time-dtt_now,time]*1e-3,label="isostasy") 

! Update ice sheet to current time 
call yelmo_update(yelmo1,time)
            
call timer_step(tmrs,comp=2,time_mod=[time-dtt_now,time]*1e-3,label="yelmo")      
```

The option `comp` tells us which component the timing is being calculated for, and
we can additionally provide a label to associate with this component. This is useful for
printing a table later.

After all components have been calculated, we can print to a summary file:

```fortran
if (mod(time_elapsed,10.0)==0) then
    ! Print timestep timing info and write log table
    call timer_write_table(tmrs,[time,dtt_now]*1e-3,"m",tmr_file,init=time_elapsed .eq. 0.0)
end if 
```

The resulting file will look something like this, here for 4 components measured during
during the time loop:

```bash
     time      dt    yelmo isostasy  climate       io    total     rate  
    0.000   0.010    0.000    0.000    0.016    0.051    0.067    6.694
    0.010   0.010    0.000    0.000    0.025    0.000    0.025    2.533
    0.020   0.010    0.000    0.000    0.024    0.000    0.025    2.458
```

Based on the options supplied, the time units are in `[m]` and the model time in `[kyr]`. The
rate is then calculated as `[m/kyr]` - this is the inverse of what we used to measure `[kyr/hr]`.
The rate as defined now is easier to manage in terms of summing the contribution of different
components, and so is preferred moving forward. To recover `[kyr/hr]`, simply take 60/rate.

## `master` to `main`

Following updated conventions, the default branch is now called `main` and the branch `master` has been deleted.

To update a working copy locally that already contains a `master` branch and therefore points to it as the default branch, the following steps should be applied:

1. Get the branch `main`.
2. Delete the local branch `master`.
3. Make sure your local repository sees `main` as the default branch.

```bash
# Get all branch information from the origin (github):
git fetch --all

# Get onto the new default branch:
git checkout main 

# Delete the branch master:
git branch -d master

# Clean up any branches that no longer exist at origin:
git fetch --prune origin

# Set the local 'head' to whatever is specified at the origin (which will be main):
git remote set-head origin -a
```

Done! Now your local copy should work like normal, with `main` instead of `master`.

## Thermodynamics equations

### Ice column

Prognostic equation:

$$
\frac{\partial T}{\partial t} = \frac{k}{\rho c} \frac{\partial^2 T}{\partial z^2} - u \frac{\partial T}{\partial x} - v \frac{\partial T}{\partial y} - w \frac{\partial T}{\partial z} + \frac{\Phi}{\rho c}
$$

Ice surface boundary condition:

$$
T(z=z_{\rm srf}) = {\rm min}(T_{\rm 2m},T_0)
$$

Ice base (temperate) boundary condition:

$$
T(z=z_{\rm bed}) = T_{\rm pmp}
$$

Ice base (frozen) boundary condition:

$$
k \frac{\partial T}{\partial z} = k_r \frac{\partial T_r}{\partial z}
$$

Note, the following internal Yelmo variables are defined for convenience:

$$
Q_{\rm ice,b} = -k \frac{\partial T}{\partial z}; \quad
Q_{\rm rock} = -k_r \frac{\partial T_r}{\partial z}
$$

### Bedrock column

Prognostic equation:

$$
\frac{\partial T_r}{\partial t} = \frac{k_r}{\rho_r c_r} \frac{\partial^2 T_r}{\partial z^2}
$$

Bedrock surface boundary condition:

$$
T_r(z=z_{\rm bed}) = T(z=z_{\rm bed})
$$

Bedrock base boundary condition:

$$
\frac{\partial T_r}{\partial z} = -\frac{Q_{\rm geo}}{k_r}
$$

## Equilibrium bedrock

In this case, the bedrock temperature profile is prescribed to the equilibrium linear temperature profile. The slope follows:

$$
\frac{\partial T_r}{\partial z} = -\frac{Q_{\rm geo}}{k_r}
$$

and the bedrock surface temperature is given by the ice temperature at its base:

$$
T_r(z=z_{\rm bed}) = T(z=z_{\rm bed})
$$

## Active bedrock

Yelmo calculates the temperature in the lithosphere along with the ice temperature. This can be achieved by assuming equilibrium conditions in the bedrock, i.e., that the temperature profile in the bedrock is always linear with `T_lith_s = T_ice_b` and the slope equal to `dT/dz = -Q_geo / k_lith`. Or, the temperature equation can be solved in the lithosphere together with the temperature in the ice column.

The parameter block `ytherm_lith` controls how the lithosphere is calculated with `ytherm_lith.method=['equil','active']` deciding the two cases above.

### Density of the upper lithosphere

### Heat capacity of the upper lithosphere

In both SICOPOLIS and GRISLI, a value of `cp = 1000.0 [J kg-1 K-1]` is used (referenced in Rogozhina et al., 2012; Greve, 2005; Greve, 1997). This value is adopted in Yelmo as well.

```fortran
cp = 1000.0    ! [J kg-1 K-1] 
```

### Heat conductivity of the upper lithosphere

Note, Yelmo expects input parameter values in units of `[J a-1 m-1 K-1]`, while much literature uses `[W m-1 K-1]`. Given the number of seconds in a year `sec_year = 31536000.0`, `kt [W m-1 K-1] * sec_year = kt [J a-1 m-1 K-1]`.

Rogozhina et al. (2012) use `kt = 2 [W m-1 K-1]` for Greenland:

```fortran
kt = 6.3e7     ! [J a-1 m-1 K-1]
```

This value is supported by Lösing et al. (2020), who perform a Bayesian inversion for GHF in Antarctica. Assuming exponentially decreasing heat production with depth, lower values of `kt` are supported (see Fig. 7b). In a study on the global thermal characteristics of the lithosphere, Cammarano and Guerri (2017) adopt an upper crust thermal conductivity of `kt = 2.5 [W m-1 K-1]`.

To do: This study is also potentially relevant:
[https://link.springer.com/article/10.1186/s40517-020-0159-y](https://link.springer.com/article/10.1186/s40517-020-0159-y). They show ranges of on the order of `kt = 2-3 [W m-1 K-1]` for the Canadian shield.

The above value of `kt = 2 [W m-1 K-1] = 6.3e7 [J a-1 m-1 K-1]` is adopted as the default thermal conductivity of the upper crust in Yelmo. For historical context, see other estimates below.

From Greve (1997) and Greve (2005):

```fortran
kt = 9.46e7    ! [J a-1 m-1 K-1]
```

which is equivalent to `kt = 3 [W m-1 K-1]`. The source of this value is not known.

From GRISLI:

```fortran
kt = 1.04e8    ! [J a-1 m-1 K-1]
```

which is equivalent to `kt = 3.3 [W m-1 K-1]`. The source of this value is not known.

## How to read `yelmo_check_kill` output

The subroutine `yelmo_check_kill` is used to see if any instability is arising in the model. If so, then a restart file is written at that moment (the earlier in the instability, the better), and the model is stopped with diagnostic output to the log file.

Note that `pc_eps` is the parameter that defines our target error tolerance in the time stepping of ice thickness evolution. At each time step, the diagnosed model error `pc_eta` is compared with `pc_eps`. If `pc_eta >> pc_eps`, this is interpreted as instability and the model is stopped.

## Margin-front mass balance

Following Pollard and DeConto (2012,2016), an ice-margin front melting scheme has been implemented that accounts for the melt rate along the vertical face of ice submerged by seawater.

The frontal mass balance ($\dot{f}$, m yr$^{-1}$) is calculated as:

$$
\dot{f} = \dot{b}_{\rm eff} \frac{A_f}{A_{\rm tot}} \theta_f
$$

where $\dot{b}_{\rm eff}$ is the effective basal mass balance (the mean of the basal mass balance calculated for the ice-free neighbors), $A_{\rm tot}=\Delta x \Delta x$ is the horizontal grid area and $A_f$ is the area of the submerged faces (i.e., the sum of the depth of submerged ice for each face of the grid cell adjacent to an ice-free cell -- potentially four faces in total). $\theta_f=10$ is a scaling coefficient that implies the face mass balance should be ~10 times higher than the basal mass balance (Pollard and DeConto, 2016, appendix).

## Calving schemes

Here is a summary of calving schemes.

### Lipscomb et al. (2019)

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

## Vertical velocity

$$
w = u_b \frac{\partial b}{\partial x} + v_b \frac{\partial b}{\partial y} - \int_b^z \left( \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} \right) dz'
$$

## Ice margin, calving rates, mass conservation

ajr, 2021-06-22

Through `v1.42`, ice margins were not fully consistently treated in Yelmo. This has been thoroughly revised. Now the following should be true:

1. The variable `f_ice` contains information of the ice area fraction of a grid cell. If `f_ice=0`, no ice is present, if `f_ice=1`, the cell is fully ice covered, and for a fractional value, this cell is designated an ice margin point with partial ice cover. To determine `f_ice`, we need to calculate the "effective ice thickness" of a grid point. For floating cells at the margin, the effective ice thickness is equivalent to either the ice thickness or the minimum ice thickness of the neighboring cell, whichever is larger. For grounded cells, the effective ice thickness must at least be that of 1/2 of the minimum ice thickness of a neighboring cell. With effective ice thickness known, `f_ice = H_ice / H_eff`.
2. Any grid cell with fractional ice cover `0<f_ice<1` is designated dynamically inactive on its outer borders (borders with ice-free points). Thus before a new cell can be populated with ice, first the fractional cell must be filled to reach `f_ice=1`.
3. Velocity and the ssa solver now explicitly only treat cells with `f_ice=1`. The dynamic-inactive borders are also enforced explicitly in the mass_conservation routine for safety.
4. The calving rate is diagnosed as a lateral flux assuming a thickness of `H_eff`, but then is converted to a horizontal mass balance component applied to the whole cell.
5. Additionally, it is possible to ensure 'residual calving' is applied to the upstream grid point of a margin cell, if the calving rate is large enough. This is done via the new routine `calc_calving_residual`. In the case of the `simple` and `flux` methods, applying full residual calving rates causes far too much calving. Rather a small amount of residual calving is applied to conservatively diminish the dynamic status of the upstream cell. The latter was implemented following CISM.

Other important changes

1. Major bug fix with thermodynamics and vertical velocity. The vertical velocity correction that was being applied to account for the sigma coordinates was incorrect. A correction must be applied to get the vertical velocity itself. Then another correction must be applied when calculating the vertical advection in the thermodynamics routine. Now this is hopefully done correctly. EISMINT EXPA and EXPF appear to work well. This may hopefully prove crucial for Javi's advance/retreat issues in Antarctica.
2. Modified staggering of 3D viscosity to be done on horizontal ab-nodes, then averaged to the center. This follows a 'quadrature' approach and appears to be more rigorously close to actually integrating over the area. It's not clear how this might affect stability of ice streams (i.e., Daniel's results).
3. Two-step mass balance with updates to `f_ice` in between, with more explicit tracking of all mass changes. Now a new variable `mb_resid` holds any threshold changes applied at the end of the mass balance update. I haven't checked in detail, but I hope all mass changes are now accounted for explicitly in the model.
