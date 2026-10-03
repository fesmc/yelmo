# Yelmo

Yelmo is a 3D ice-sheet-shelf model solving
for the coupled dynamics and thermodynamics of the ice sheet system. Yelmo
can be used for idealized simulations, stand-alone ice sheet simulations
and fully coupled ice-sheet and climate simulations.

The physics and design of the model are described in the following article:

> Robinson, A., Alvarez-Solas, J., Montoya, M., Goelzer, H., Greve, R., and Ritz, C.: Description and validation of the ice-sheet model Yelmo (version 1.0), Geosci. Model Dev., 13, 2805–2823, [https://doi.org/10.5194/gmd-13-2805-2020](https://doi.org/10.5194/gmd-13-2805-2020), 2020.

The model documentation is provided to help with proper use of the model,
and can be found at:

 [https://fesmc.github.io/yelmo](https://fesmc.github.io/yelmo)
 
While the model has been designed to be easy to use, there
are many parameters that require knowledge of ice-sheet
physics and numerous parameterizations. It is not recommended to use the ice
sheet model as a black box without understanding of the key parameters that
affect its performance.

The test cases shown by Robinson et al. (2020) can be run with the current code
following the instructions below in the section "Test cases".

To get started with compiling and running the model, see the quick-start
instructions below, or the documentation: [https://fesmc.github.io/yelmo/getting-started.html](https://fesmc.github.io/yelmo/getting-started.html).

# Getting started

Here you can find the basic information and steps needed to get **Yelmo** running.

## Super-quick start

A summary of commands to get started is given below. For more detailed information see subsequent sections.

```bash
# Install the configme and runme commands (one time, system-wide)
pip install git+https://github.com/fesmc/configme
pip install git+https://github.com/fesmc/runme

# Clone Yelmo with the packages it needs (fesm-utils, FastHydrology, elsa, tracer),
# generate the Makefiles for this machine and compiler, and build fesm-utils
configme install yelmo
cd yelmo

# Create the local runme config (.runme/config.toml)
runme config init

# Compile the benchmarks program
make clean
make benchmarks

# Run a test simulation of the EISMINT1-moving experiment
runme -r -e benchmarks -o output/eismint1-moving -n par/yelmo_EISMINT_moving.nml

# Compile the initmip program and run a simulation of Antarctica
make initmip
runme -r -e initmip -o output/ant-pd -n par/yelmo_initmip.nml -p ctrl.set_nm="set_ant_pd" yelmo.domain="Antarctica" yelmo.grid_name="ANT-32KM"
```

## Dependencies

See: [Installation](https://fesmc.github.io/yelmo/getting-started.html) for installation tips.

- NetCDF library (preferably version 4.0 or higher), with the Fortran interface. This is the only library to install yourself.
- LIS: [Library of Iterative Solvers for Linear Systems](http://www.ssisc.org/lis/) and FFTW, built inside `fesm-utils` by `configme install`.
- The Fortran packages [fesm-utils](https://github.com/fesmc/fesm-utils), [FastHydrology](https://github.com/fesmc/FastHydrology), [elsa](https://github.com/fesmc/elsa) and [tracer](https://github.com/fesmc/tracer), cloned and linked into the checkout by `configme install yelmo`.
- Python 3 with the [`configme`](https://github.com/fesmc/configme) package (Makefile generation and installation) and the [`runme`](https://github.com/fesmc/runme) package (running single simulations and ensembles, job submission). Install both with `pip` as shown above.

## Parameter configuration tool (`yelmo-config`)

Yelmo reads a canonical default parameter set from `input/yelmo_defaults.nml`; a
user parameter file only needs to list the parameters it wants to override.
`yelmo-config` is a command-line tool to **discover, manage, compare and
validate** Yelmo parameter files against those defaults and the model's built-in
consistency checks.

Install it system-wide straight from this repository:

```bash
pip install -U "git+https://github.com/fesmc/yelmo#subdirectory=tools/yelmo-config"
```

A snapshot of the defaults and constraints is bundled, so it works without a
local checkout; run it inside a checkout to use that copy's live files. Then,
e.g.:

```bash
yelmo-config list ydyn                   # browse the dynamics parameters
yelmo-config check par/yelmo_initmip.nml # validate a parameter file
yelmo-config diff runA.nml runB.nml      # compare two runs
yelmo-config update                      # self-update to the latest version
```

See [tools/yelmo-config/README.md](tools/yelmo-config/README.md) and the
[documentation](https://fesmc.github.io/yelmo/yelmo-config.html) for the full
command set.

## Directory structure

```fortran
    config/
        Makefile templates (config/Makefile, used by configme) and legacy
        host configurations (config/legacy/).
    docs/
        Documentation (Quarto site).
    input/
        Default parameters (yelmo_defaults.nml), physical constants and the
        variable tables used for output.
    libs/
        Auxiliary libraries nesecessary for running the model.
    libyelmo/
        Compiled files: include/ (object and module files, libyelmo.a) and bin/ (executables).
    output/
        Default location for model output.
    par/
        Parameter files of the test programs.
    src/
        Source code for Yelmo.
    tests/
        Source code and analysis scripts for specific model benchmarks and tests.
    tools/
        Auxiliary tools, including yelmo-config for parameter management.
```

## Usage

Follow the steps below to (1) obtain the code, (2) configure the Makefile for your system,
(3) compile the Yelmo static library and an executable program and (4) run a test simulation.

### 1. Get the code.

The recommended way is `configme install yelmo` (see above), which clones
[https://github.com/fesmc/yelmo](https://github.com/fesmc/yelmo) together with the
packages it needs and configures them. To pick the machine and compiler, the clone
protocol or the directory explicitly:

```bash
configme install yelmo -m dkrz_levante -c ifx   # machine + compiler
configme install yelmo -d https                 # clone over HTTPS
configme install yelmo --dir ~/models/yelmo     # checkout location
```

Run `configme list` for the supported machines and compilers.

If you plan to make changes to the code, it is wise to check out a new branch:

```bash
git checkout -b user-dev
```

You should now be working on the branch `user-dev`.

### 2. Create the system-specific Makefile.

`configme install` already generated the Makefiles. To regenerate them (e.g. after
changing machine or compiler, or after `config/Makefile` changed), run from the
checkout:

```bash
configme config yelmo -m macbook -c gfortran
```

This writes `Makefile` in the Yelmo root (from the template `config/Makefile`) and
the Makefiles of `fesm-utils`, `FastHydrology`, `elsa` and `tracer`, with the netCDF
paths detected automatically. A new machine is added with `configme new`. The old
`python config.py config/legacy/<host>` route is kept for reference only.

### 3. Compile the code.

Now you are ready to compile Yelmo as a static library:

```bash
make clean    # This step is very important to avoid errors!!
make yelmo-static [debug=1] [openmp=1]
```
This compiles the static libraries of FastHydrology, elsa and tracer, then all of the
Yelmo modules (as defined in `config/Makefile_yelmo.mk`), and links them in the static
library `libyelmo/include/libyelmo.a`. All compiled files can be found in the folder `libyelmo/`.

Once the static library has been compiled, it can be used inside of external Fortran programs and modules
via the statement `use yelmo`.
To include/link yelmo-static during compilation of another program, its location must be defined:

```
INC_YELMO = -I${YELMOROOT}/libyelmo/include
LIB_YELMO = -L${YELMOROOT}/libyelmo/include -lyelmo
```

Alternatively, several test programs exist in the folder `tests/` to run Yelmo
as a stand-alone ice sheet, for example:

```bash
make benchmarks    # libyelmo/bin/yelmo_benchmarks.x: EISMINT, HALFAR, BUELER
make mismip        # libyelmo/bin/yelmo_mismip.x:     MISMIP3D
make trough        # libyelmo/bin/yelmo_trough.x:     TROUGH-F17, MISMIP+, SLAB-S06
make ismiphom      # libyelmo/bin/yelmo_ismiphom.x:   ISMIP-HOM
make calving       # libyelmo/bin/yelmo_calving.x:    CalvingMIP
make initmip       # libyelmo/bin/yelmo_initmip.x:    realistic domains (initMIP Greenland/Antarctica)
```

`make usage` lists all targets. The option `debug=1` compiles with debugging flags
(e.g., `make benchmarks debug=1`); the code then runs much slower, so this option
is not recommended unless necessary. The option `openmp=1` compiles with OpenMP.

### 4. Run the model.

Once an executable has been created, you can run the model. This can be
achieved via the `runme` command, installed via `pip` (see above). Before the first
run, create the local runme config with `runme config init` (this writes
`.runme/config.toml` from `.runme/config.default.toml`; set `hpc` and `account` there
to submit jobs). The following steps are carried out by `runme`:

1. The output directory is created.
2. The executable is copied to the output directory
3. The relevant parameter files are copied to the output directory.
4. Links to the input data paths (`input`, `ice_data` and `maps`) are created in the output directory. Note that many simulations, such as benchmark experiments, do not depend on these external data sources, but the links are made anyway.
5. The executable is run from the output directory, either as a background process or it is submitted to the queue via `sbatch` (the SLURM workload manager).

To run a benchmark simulation, for example, use the following command:

```bash
runme -r -e benchmarks -o output/test -n par/yelmo_EISMINT_moving.nml
```

where the option `-r` implies that the model should be run as a background process. If this is omitted, then the output directory will be populated, but no executable will be run, while `-s` instead prepares a job script for the cluster queue system and `-rs` also submits it (`-q` selects the queue alias, see `runme queues`). The option `-e` lets you specify the executable. For the standard programs, shortcuts are defined in `.runme/info.json`:

```
benchmarks = libyelmo/bin/yelmo_benchmarks.x
calving    = libyelmo/bin/yelmo_calving.x
mismip     = libyelmo/bin/yelmo_mismip.x
initmip    = libyelmo/bin/yelmo_initmip.x
trough     = libyelmo/bin/yelmo_trough.x
ismiphom   = libyelmo/bin/yelmo_ismiphom.x
mask_ice   = libyelmo/bin/yelmo_mask_ice.x
regridding = libyelmo/bin/yelmo_test_regridding.x
```
The last two mandatory arguments `-o OUTDIR` and `-n PAR_PATH` are the output/run directory and the parameter file to be used for this simulation, respectively. In the case of the above simulation, the output directory is defined as `output/test`, where all model parameters (loaded from the file `par/yelmo_EISMINT_moving.nml`) and model output can be found.

It is also possible to modify parameters inline via the option `-p KEY=VAL [KEY=VAL ...]`. The parameter should be specified with its namelist group and its name. E.g., to change the resolution of the EISMINT benchmark experiment to 10km, use:

```bash
runme -r -e benchmarks -o output/test -n par/yelmo_EISMINT_moving.nml -p ctrl.dx=10
```

To run an ensemble, pass comma-separated values to `-p` (e.g. `-p ctrl.dx=10,20,40`); `runme` creates one run directory per combination under `-o`.

See `runme -h` for more details, or the [runme README](https://github.com/fesmc/runme). 

## Test cases

The published model description includes several test simulations for validation
of the model's performance. The following section describes how to perform these
tests with the current code; results differ from the model version documented in
the article (see CHANGELOG.md). From this point, it is assumed that the model is
configured for your system (see above) and ready to compile. The script
`run_batch_gmd.sh` runs these and further benchmarks as a batch, and
[docs/benchmarks.md](docs/benchmarks.md) describes them.

### 1. EISMINT1 moving margin experiment
To perform the moving margin experiment, compile the benchmarks
executable and call it with the EISMINT parameter file:

```bash
make benchmarks
runme -r -e benchmarks -o output/eismint-moving -n par/yelmo_EISMINT_moving.nml
```

### 2. EISMINT2 EXPA
To perform Experiment A from the EISMINT2 benchmarks, compile the benchmarks
executable and call it with the EXPA parameter file:

```bash
make benchmarks
runme -r -e benchmarks -o output/eismint-expa -n par/yelmo_EISMINT_expa.nml
```

### 3. EISMINT2 EXPF
To perform Experiment F from the EISMINT2 benchmarks, compile the benchmarks
executable and call it with the EXPF parameter file:

```bash
make benchmarks
runme -r -e benchmarks -o output/eismint-expf -n par/yelmo_EISMINT_expf.nml
```

### 4. MISMIP RF
To perform the MISMIP3D rate factor experiment, compile the mismip executable
and call it with the MISMIP3D parameter file and the three parameter permutations of interest (default, subgrid and subgrid+gl-scaling):

```bash
make mismip
runme -r -e mismip -o output/mismip-rf-0 -n par/yelmo_MISMIP3D.nml -p ctrl.experiment="RF" ydyn.beta_gl_stag=0 ydyn.beta_gl_scale=0
runme -r -e mismip -o output/mismip-rf-1 -n par/yelmo_MISMIP3D.nml -p ctrl.experiment="RF" ydyn.beta_gl_stag=3 ydyn.beta_gl_scale=0
runme -r -e mismip -o output/mismip-rf-2 -n par/yelmo_MISMIP3D.nml -p ctrl.experiment="RF" ydyn.beta_gl_stag=3 ydyn.beta_gl_scale=2
```
To additionally change the resolution of the simulations change the parameter `ctrl.dx` [km], e.g. for the default simulation with 10km resolution, call:

```bash
runme -r -e mismip -o output/mismip-rf-0-10km -n par/yelmo_MISMIP3D.nml -p ctrl.experiment="RF" ydyn.beta_gl_stag=0 ydyn.beta_gl_scale=0 ctrl.dx=10
```

### 5. Age profile experiments
To perform the age profile experiments, compile the Fortran program `tests/test_icetemp.f90`
and run it with the number of vertical points, the conductivity ratio, the solver and
the experiment as arguments (output is written to `output/`):

```bash
make icetemp
./libyelmo/bin/test_icetemp.x 12 1e-3 temp eismint
./libyelmo/bin/test_icetemp.x 32 1e-3 temp eismint
./libyelmo/bin/test_icetemp.x 52 1e-3 temp eismint
```

To compare single and double precision, recompile after changing the working
precision `wp` (`sp` or `dp`) in the file `src/yelmo_defs.f90`.

### 6. Antarctica present-day and glacial simulations
To perform the Antarctica simulations as presented in the paper, it is necessary
to compile the `initmip` executable and run `par/yelmo_initmip.nml` with the
present-day (`set_ant_pd`) and glacial (`set_ant_lgm`) settings:


```bash
make initmip
runme -r -e initmip -o output/ant-pd  -n par/yelmo_initmip.nml -p ctrl.set_nm="set_ant_pd"  yelmo.domain="Antarctica" yelmo.grid_name="ANT-32KM"
runme -r -e initmip -o output/ant-lgm -n par/yelmo_initmip.nml -p ctrl.set_nm="set_ant_lgm" yelmo.domain="Antarctica" yelmo.grid_name="ANT-32KM"
```
