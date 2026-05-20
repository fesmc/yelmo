# Dependencies

Yelmo is dependent on the following libraries:

- [NetCDF](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)
- [Library of Iterative Solvers for Linear Systems](http://www.ssisc.org/lis/)
- ['runner' Python library (fesmc version)](https://github.com/fesmc/runner)

YelmoX is additionally dependent on the following library:

- FFTW (ver. 3.9+)

Installation tips for each dependency can be found below.

## Installing NetCDF (preferably version 4.0 or higher)

The NetCDF library is typically available with different distributions (Linux, Mac, etc).
Along with installing `libnetcdf`, it will be necessary to install the package `libnetcdf-dev`.
Installing the NetCDF viewing program `ncview` is also recommended.

If you want to install NetCDF from source, then you must install both the
`netcdf-c` and subsequently `netcdf-fortran` libraries. The source code and
installation instructions are available from the Unidata website:

[https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)

## Install LIS, FFTW and utils

These packages could be installed individually and linked into the main directories of Yelmox and Yelmo. However, to ensure the right versions are used, etc., we have now made a separate repository for managing the installation of LIS and FFTW from the versions available in that repository, as well as some custom modules. This repository is managed as part of the Fast Earth System Model Community (FESMC).

Please download the code from this repository and see the README for installation instructions:
[https://github.com/fesmc/fesm-utils](https://github.com/fesmc/fesm-utils)

## Installing runner

1. Install `runner` to your system's Python installation via `pip`, along with dependency `tabulate`.

```bash
pip install https://github.com/fesmc/runner/archive/refs/heads/master.zip
```

That's it! Now check that system command `job` is available by running `job -h`. If the command is not found, it means that the Python bin directory is not available in your PATH. To add it, typically something like this is needed in your .profile or .bashrc file:

```bash
PATH=${PATH}:${HOME}/.local/bin
export PATH
```
