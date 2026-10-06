# Installation

Here you can find the basic information and steps needed to get **Yelmo** running.

## Quick start

Yelmo is configured and built with [`configme`](https://github.com/fesmc/configme),
a Python tool that detects your netCDF installation, configures every package of
the stack for your machine and compiler, and clones, links and builds them with
one command. Simulations are run with [`runme`](https://github.com/fesmc/runme).
Install both once, system-wide:

```bash
pip install git+https://github.com/fesmc/configme
pip install git+https://github.com/fesmc/runme
```

Then, from the directory where the checkout should live:

```bash
configme install yelmo
```

This clones Yelmo together with the packages it needs (`fesm-utils`,
`FastHydrology`, `elsa` and `tracer`), generates the Makefiles for your machine
and compiler, links the packages, and builds `fesm-utils` (LIS, FFTW, SHTns and
the utilities, which can take 10–30 min; `configme` asks before building). If
`configme` can detect your machine from the hostname it does so, otherwise it
prompts you.

Finally, create the local runme configuration in the checkout:

```bash
cd yelmo
runme config init     # writes .runme/config.toml from .runme/config.default.toml
```

The only system dependency you must install yourself is **netCDF** (see
[Dependencies](#dependencies)). Once the install finishes you are ready to
compile and run; see [Usage](#usage) below.

### configme install options

Common options for `configme install yelmo`:

```bash
configme install yelmo -m dkrz_levante -c ifx   # pick the machine + compiler explicitly
configme install yelmo -d https                 # clone over HTTPS (no GitHub SSH key needed)
configme install yelmo --dir ~/models/yelmo     # put the checkout here instead of ./yelmo
configme install yelmo --overwrite              # re-clone over an existing checkout
configme install yelmo --build-deps             # build the dependency packages without prompting
```

Run `configme list` for the supported machines and compilers, and
`configme --help` for the full command surface. `configme install` also writes
the commands it ran to a `.install.sh` script in the checkout; see
[configme install details](configme-install-details.md).

If the `configme` or `runme` command is not found after installation, your
Python user bin directory is probably not on your `PATH`; add it in your
`~/.bashrc` / `~/.zshrc`:

```bash
export PATH="${PATH}:${HOME}/.local/bin"
```

To update `configme` later, run `configme update`.

### Existing checkout

If you cloned Yelmo by hand
([https://github.com/fesmc/yelmo](https://github.com/fesmc/yelmo)), run
`configme install yelmo` from inside the checkout: it uses the existing
checkout and clones, configures and links the packages it needs. To regenerate
the Makefiles later (e.g. for another compiler, or after `config/Makefile`
changed), run from the checkout:

```bash
configme config yelmo -m macbook -c gfortran
```

This writes the Makefiles of Yelmo, `fesm-utils`, `FastHydrology`, `elsa` and
`tracer`, with the netCDF paths detected automatically.

### Input data

Simulations of realistic domains (e.g. initmip) read their input from an
`ice_data` directory in the checkout, which is not part of the repository.
Link it to the data on your system:

```bash
ln -s <path-to>/ice_data ice_data
```

See [HPC notes](hpc-notes.md) for the locations on the supported clusters.

## Dependencies

The Yelmo stack depends on:

- [NetCDF](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)
  (preferably version 4.0 or higher), with the Fortran interface;
- [LIS](http://www.ssisc.org/lis/) (Library of Iterative Solvers for Linear
  Systems), FFTW and SHTns, built inside `fesm-utils`;
- the Fortran packages [fesm-utils](https://github.com/fesmc/fesm-utils),
  [FastHydrology](https://github.com/fesmc/FastHydrology),
  [elsa](https://github.com/fesmc/elsa) and [tracer](https://github.com/fesmc/tracer);
- Python 3 with the `configme` and `runme` packages.

Only **netCDF** must be installed on your system beforehand; everything else is
managed by `configme`.

### Installing NetCDF

The NetCDF library is typically available with different distributions (Linux, Mac, etc).
Along with installing `libnetcdf`, it will be necessary to install the package `libnetcdf-dev`.
Installing the NetCDF viewing program `ncview` is also recommended.

If you want to install NetCDF from source, then you must install both the
`netcdf-c` and subsequently `netcdf-fortran` libraries. The source code and
installation instructions are available from the Unidata website:

[https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)

`configme` finds netCDF through `nf-config` and `nc-config`, so these must be on
your `PATH` (on clusters, load the netCDF modules first).

## Directory structure

```fortran
    analysis/
        Analysis scripts and data of model tests.
    config/
        Makefile templates (config/Makefile, used by configme) and legacy
        host configurations (config/legacy/).
    docs/
        Documentation (Quarto site).
    input/
        Default parameters (yelmo_defaults.nml), physical constants and the
        variable tables used for output.
    libs/
        Auxiliary libraries necessary for running the model.
    libyelmo/
        Compiled files: include/ (object and module files, libyelmo.a) and
        bin/ (executables).
    maps/
        Grid description files and scripts to generate remapping weights.
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

### 1. Compile the code

Compile Yelmo as a static library:

```bash
make clean    # This step is very important to avoid errors!!
make yelmo-static [debug=1] [openmp=0]
```

This compiles the static libraries of FastHydrology, elsa and tracer, then all
of the Yelmo modules (as defined in `config/Makefile_yelmo.mk`), and links them
in the static library `libyelmo/include/libyelmo.a`. `make clean` also cleans
FastHydrology, elsa and tracer. `fesm-utils` is not built by Yelmo's Makefile.

Once the static library has been compiled, it can be used inside of external
Fortran programs and modules via the statement `use yelmo`. To include/link
yelmo-static during compilation of another program, its location must be
defined:

```bash
INC_YELMO = -I${YELMOROOT}/libyelmo/include
LIB_YELMO = -L${YELMOROOT}/libyelmo/include -lyelmo
```

The program must also link the libraries of the dependencies (elsa, tracer,
FastHydrology, fesm-utils, LIS, FFTW and netCDF); see `LFLAGS` in the
generated `Makefile`.

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

`make usage` lists all targets. The option `debug=1` compiles with debugging
flags (e.g., `make benchmarks debug=1`); the code then runs much slower, so this
option is not recommended unless necessary. Yelmo is compiled with OpenMP by
default; `openmp=0` compiles without it.

### 2. Run the model

Once an executable has been created, you can run the model with `runme`. The
following steps are carried out by `runme`:

1. The output directory is created.
2. The executable is copied to the output directory.
3. The relevant parameter files are copied to the output directory.
4. Links to the input data paths (`input`, `ice_data` and `maps`) are created
   in the output directory. Many simulations, such as benchmark experiments, do
   not depend on these external data sources, but the links are made anyway.
5. The executable is run from the output directory, either as a background
   process or as a job submitted to the queue via `sbatch` (the SLURM workload
   manager).

To run a benchmark simulation, for example, use the following command:

```bash
runme -r -e benchmarks -o output/test -n par/yelmo_EISMINT_moving.nml
```

The option `-r` runs the model as a background process. Without it, the output
directory is populated, but no executable is run. The option `-s` prepares a
job script for the cluster queue system instead, and `-rs` also submits it
(`-q` selects the queue, `-w` the wall time and `--omp` the number of OpenMP
threads). To submit jobs, set `hpc` and `account` in `.runme/config.toml`.

The option `-e` sets the executable. For the standard programs, shortcuts are
defined in `.runme/info.json`:

```bash
benchmarks = libyelmo/bin/yelmo_benchmarks.x
calving    = libyelmo/bin/yelmo_calving.x
mismip     = libyelmo/bin/yelmo_mismip.x
initmip    = libyelmo/bin/yelmo_initmip.x
trough     = libyelmo/bin/yelmo_trough.x
ismiphom   = libyelmo/bin/yelmo_ismiphom.x
mask_ice   = libyelmo/bin/yelmo_mask_ice.x
regridding = libyelmo/bin/yelmo_test_regridding.x
```

The arguments `-o OUTDIR` and `-n PAR_PATH` are the output/run directory and
the parameter file to be used for this simulation, respectively. In the case of
the above simulation, the output directory is `output/test`, where all model
parameters (loaded from the file `par/yelmo_EISMINT_moving.nml`) and model
output can be found.

Parameters can be modified inline with `-p KEY=VAL [KEY=VAL ...]`, with the
parameter given by its namelist group and its name. E.g., to change the
resolution of the EISMINT benchmark experiment to 10 km, use:

```bash
runme -r -e benchmarks -o output/test -n par/yelmo_EISMINT_moving.nml -p ctrl.dx=10
```

For ensembles, pass comma-separated values to `-p` (e.g. `-p ctrl.dx=10,20,40`);
`runme` creates one run directory per combination under `-o`. See `runme -h`
for more details, or the [runme README](https://github.com/fesmc/runme).
