# Notes

## `timeout` module

The `timeout` module of fesm-utils defines the output times of a driver
program. Each output stream has its own parameter group, whose name is given to
`timeout_init`. For example, initmip uses `&t1D`, `&t2Dsm`, `&t2D` and `&trst`
(restarts):

```fortran
&t2D
    method          = "const"           ! "none", "const", "file", "times"
    dt              = 1000.0
    file            = "none"
    times           = 0
/
```

```fortran
! Get output times
call timeout_init(t2D,path_par,"t2D","heavy",time_init,time_end)
```

With `method = "file"`, the times are loaded from a file with one time per
line, or a range of times using the format `t0:dt:t1`:

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
`time_init` and `time_end`, and `time_end` is always included. In this way, once `timeout_init` is called,
we know how many timesteps of output will be generated. This can help
confirm that we designed the experiment well, and how much data to expect.

Then during the timeloop, simply use the function `timeout_check` to
determine if the current time should be written to output:

```fortran
if (timeout_check(t2D,time)) then
    call yelmo_write_step(yelmo1,file2D,time)
end if
```

## `timer` module

The `timer` module of fesm-utils measures the wall time of several components
of a driver program (e.g. isostasy, climate and Yelmo in yelmox) with one
object. Yelmo itself measures its timing with `yelmo_cpu_time`.

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
rate is then calculated as `[m/kyr]`, which can be summed over components. To
obtain `[kyr/hr]`, take 60/rate.

## How to read `yelmo_check_kill` output

The subroutine `yelmo_check_kill` is used to see if any instability is arising in the model. If so, then a snapshot of the model state is written to `yelmo_killed.nc` at that moment (the earlier in the instability, the better), and the model is stopped with diagnostic output to the log file.

The criteria are listed in [Time stepping](physics/timestepping.md#instability-checks). The error measure is `pc_eta`, the norm of the predictor–corrector truncation error [1/yr]: `pc_eps` is its target value for the adaptive time step, and `pc_tol` the value above which a time step is redone with a smaller dt. If the mean of the stored `pc_eta` values (the last three steps) exceeds `10*pc_tol`, this is interpreted as instability and the model is stopped. The checks are switched off with `yelmo.disable_kill = True`.
