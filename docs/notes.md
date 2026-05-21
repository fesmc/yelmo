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
