# Shared build configuration for yelmo (dependency wiring).
#
# Loaded *after* the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> netCDF -> common). This file references
# variables those provide: FFLAGS / FFLAGS_OPENMP (compiler) and LIB_NC
# (machine or auto-detected netCDF).

# Dependency paths (serial build by default).
FESMUTILSROOT = fesm-utils/utils
INC_FESMUTILS = -I${FESMUTILSROOT}/include-serial
LIB_FESMUTILS = -L${FESMUTILSROOT}/include-serial -lfesmutils

LISROOT = fesm-utils/lis-serial
INC_LIS = -I${LISROOT}/include
LIB_LIS = -L${LISROOT}/lib/ -llis

# PETSc is an optional linear solver (enabled with `make petsc=1`). It is not
# managed by configme: PETSCROOT / PETSC_DIR are site-specific — set PETSC_DIR
# in the environment, or override PETSCROOT here for your machine.
PETSCROOT = /opt/local/lib/petsc
INC_PETSC = -I $(PETSC_DIR)/include
LIB_PETSC = -L${PETSC_DIR}/lib -lpetsc

# OpenMP build (make openmp=1): swap the serial dependency builds for their
# OpenMP variants and append the compiler's OpenMP flag (FFLAGS_OPENMP, set in
# the compiler fragment).
ifeq ($(openmp), 1)
    INC_FESMUTILS = -I${FESMUTILSROOT}/include-omp
    LIB_FESMUTILS = -L${FESMUTILSROOT}/include-omp -lfesmutils

    LISROOT = fesm-utils/lis-omp
    INC_LIS = -I${LISROOT}/include
    LIB_LIS = -L${LISROOT}/lib/ -llis

    FFLAGS += $(FFLAGS_OPENMP)
endif

# Linear solvers to include: LIS is always required, PETSc is optional.
INC_LINEAR = $(INC_LIS)
LIB_LINEAR = $(LIB_LIS)
ifeq ($(petsc), 1)
    INC_LINEAR = $(INC_LIS) $(INC_PETSC)
    LIB_LINEAR = $(LIB_LIS) $(LIB_PETSC)
endif

# Extra link flags. -Wl,-zmuldefs works around duplicate symbols in the static
# deps (the default on Linux). A machine fragment can disable it by setting
# `LFLAGS_EXTRA =` (macOS ld rejects -zmuldefs, so the macbook fragment does).
LFLAGS_EXTRA ?= -Wl,-zmuldefs

LFLAGS = $(LIB_NC) $(LIB_FESMUTILS) $(LIB_LINEAR) $(LFLAGS_EXTRA)
