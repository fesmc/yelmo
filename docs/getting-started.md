# Installation

Here you can find the basic information and steps needed to get **Yelmo** running.

## How to get Yelmo

The Yelmo code repository is here: [https://github.com/fesmc/yelmo](https://github.com/fesmc/yelmo). To get a local copy:

```bash
git clone git@github.com:fesmc/yelmo.git
cd yelmo
```

## Install configme (one time)

Yelmo is configured and built with [`configme`](https://github.com/fesmc/configme), a small Python tool that detects your netCDF installation, configures every package in the stack for your machine and compiler, and clones/links/builds the whole thing with one command. It is installed once, globally, and provides the `configme` command on your `PATH`:

```bash
pip install git+https://github.com/fesmc/configme
```

To upgrade it later, add `--upgrade` to the same command. If the `configme` command is not found afterwards, your Python user bin directory is probably not on your `PATH`; add it in your `~/.bashrc` / `~/.zshrc`:

```bash
export PATH="${PATH}:${HOME}/.local/bin"
```

The only system dependency you must install yourself is **netCDF** (see [Dependencies](#dependencies)). Everything else — LIS, FFTW, the `fesm-utils` libraries, and `runme` — is managed by `configme`.

## Quick start

With `configme` installed, build Yelmo and the packages it needs with a single command from the directory where you want the checkout to live:

```bash
configme install yelmo
```

This clones Yelmo together with `fesm-utils`, configures each for your machine and compiler, links them, and builds `fesm-utils` (LIS + FFTW + utils, which can take 10-30 min). If `configme` can detect your machine from the hostname it does so, otherwise it prompts you.

Common options:

```bash
configme install yelmo -m dkrz_levante -c ifx   # pick the machine + compiler explicitly
configme install yelmo -d https                 # clone over HTTPS (no GitHub SSH key needed)
configme install yelmo --dir ~/models/yelmo     # put the checkout here instead of ./yelmo
configme install yelmo --overwrite              # re-clone over an existing checkout
configme install yelmo --build-deps             # rebuild dependency packages without prompting
```

Run `configme list` for the supported machines and compilers, and `configme --help` for the full command surface. The exact clone/configure/link/build commands `configme install yelmo` runs for you are recorded in a `.install.sh` script in the checkout, and are also shown for context on the [configme install details](configme-install-details.md) page — these are shown for reference only; `configme install` is the recommended path and you do not need to run them by hand.

Once the install finishes you are ready to compile and run; see [Usage](#usage) below.

## Dependencies

The Yelmo stack depends on the following libraries:

- [NetCDF](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)
- [Library of Iterative Solvers for Linear Systems](http://www.ssisc.org/lis/)
- FFTW (ver. 3.9+)
- ['runme' Python package (fesmc)](https://github.com/fesmc/runme)

Of these, only **NetCDF** must be installed on your system beforehand. LIS, FFTW (built via `fesm-utils`) and `runme` are all managed for you by `configme`. Installation tips for netCDF and `runme` can be found below.

### Installing NetCDF (preferably version 4.0 or higher)

The NetCDF library is typically available with different distributions (Linux, Mac, etc).
Along with installing `libnetcdf`, it will be necessary to install the package `libnetcdf-dev`.
Installing the NetCDF viewing program `ncview` is also recommended.

If you want to install NetCDF from source, then you must install both the
`netcdf-c` and subsequently `netcdf-fortran` libraries. The source code and
installation instructions are available from the Unidata website:

[https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)

### Installing runme

`runme` prepares and runs Yelmo simulations (single runs and ensembles). `configme install` installs it for you if it is missing, so you normally do not need to install it by hand. To install it manually into your system's Python installation via `pip`:

```bash
pip install git+https://github.com/fesmc/runme
```

That's it! Ensemble support is built in, so no separate `runner` package is needed. Now check that the system command `runme` is available by running `runme -h`. If the command is not found, it means that the Python bin directory is not available in your PATH. To add it, typically something like this is needed in your .profile or .bashrc file:

```bash
PATH=${PATH}:${HOME}/.local/bin
export PATH
```

## Directory structure

```fortran
    config/
        Configuration files for compilation on different systems.
    input/
        Location of any input data needed by the model.
    libs/
        Auxiliary libraries nesecessary for running the model.
    libyelmo/
        Folder containing all compiled files in a standard way with
        lib/, include/ and bin/ folders.
    output/
        Default location for model output.
    par/
        Default parameter files that manage the model configuration.
    src/
        Source code for Yelmo.
    tests/
        Source code and analysis scripts for specific model benchmarks and tests.
```

## Usage

After `configme install yelmo` (see [Quick start](#quick-start)) you have a configured Yelmo checkout with its `Makefile`, linked libraries and `.runme_config` already in place. The steps below cover (1) compiling the code and (2) running a simulation.

### 1. Compile the code

Now you are ready to compile Yelmo as a static library:

```bash
make clean    # This step is very important to avoid errors!!
make yelmo-static [debug=1]
```

This will compile all of the Yelmo modules and libraries (as defined in `config/Makefile_yelmo.mk`),
and link them in a static library. All compiled files can be found in the folder `libyelmo/`.

Once the static library has been compiled, it can be used inside of external Fortran programs and modules
via the statement `use yelmo`.
To include/link yelmo-static during compilation of another program, its location must be defined:

```bash
INC_YELMO = -I${YELMOROOT}/include
LIB_YELMO = -L${YELMOROOT}/include -lyelmo
```

Alternatively, several test programs exist in the folder `tests/` to run Yelmo
as a stand-alone ice sheet.
For example, it's possible to run different EISMINT benchmarks, MISMIP benchmarks and the
ISIMIP6 INITMIP simulation for Greenland, respectively:

```bash
make benchmarks    # compiles the program `libyelmo/bin/yelmo_benchmarks.x`
make mismip        # compiles the program `libyelmo/bin/yelmo_mismip.x`
make initmip       # compiles the program `libyelmo/bin/yelmo_initmip.x`
```

The Makefile additionally allows you to specify debugging compiler flags with the option `debug=1`, in case you need to debug the code (e.g., `make benchmarks debug=1`). Using this option, the code will run much slower, so this option is not recommended unless necessary.

### 2. Run the model

Once an executable has been created, you can run the model. This can be
achieved via the `runme` command (installed via `pip`, see [Installing runme](#installing-runme)). The following steps
are carried out by `runme`:

1. The output directory is created.
2. The executable is copied to the output directory
3. The relevant parameter files are copied to the output directory.
4. Links to the input data paths (`input` and `ice_data`) are created in the output directory. Note that many simulations, such as benchmark experiments, do not depend on these external data sources, but the links are made anyway.
5. The executable is run from the output directory, either as a background process or it is submitted to the queue via `sbatch` (the SLURM workload manager).

To run a benchmark simulation, for example, use the following command:

```bash
runme -r -e benchmarks -o output/test -n par/yelmo_EISMINT.nml
```

where the option `-r` implies that the model should be run as a background process. If this is omitted, then the output directory will be populated, but no executable will be run, while `-s` instead will submit the simulation to cluster queue system instead of running in the background. The option `-e` lets you specify the executable. For some standard cases, shortcuts have been created:

```bash
benchmarks = libyelmo/bin/yelmo_benchmarks.x
mismip     = libyelmo/bin/yemo_mismip.x
initmip    = libyelmo/bin/yelmo_initmip.x
```

The last two mandatory arguments `-o OUTDIR` and `-n PAR_PATH` are the output/run directory and the parameter file to be used for this simulation, respectively. In the case of the above simulation, the output directory is defined as `output/test`, where all model parameters (loaded from the file `par/yelmo_EISMINT.nml`) and model output can be found.

It is also possible to modify parameters inline via the option `-p KEY=VAL [KEY=VAL ...]`. The parameter should be specified with its namelist group and its name. E.g., to change the resolution of the EISMINT benchmark experiment to 10km, use:

```bash
runme -r -e benchmarks -o output/test -n par/yelmo_EISMINT.nml -p ctrl.dx=10
```

For ensembles, pass comma-separated values to `-p` (e.g. `-p ctrl.dx=10,20,40`); `runme` creates one run directory per combination under `-o`. See `runme -h` for more details, or the [runme README](https://github.com/fesmc/runme).
