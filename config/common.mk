# Shared build configuration for yelmo (dependency wiring).
#
# Loaded *after* the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> netCDF -> common). This file references
# variables those provide: FFLAGS / FFLAGS_OPENMP (compiler) and LIB_NC
# (machine or auto-detected netCDF).

# Dependency paths (serial build by default).
FESMUTILSROOT = fesm-utils
INC_FESMUTILS = -I${FESMUTILSROOT}/include-serial
LIB_FESMUTILS = -L${FESMUTILSROOT}/include-serial -lfesmutils

LISROOT = fesm-utils/lis/lis-serial
INC_LIS = -I${LISROOT}/include
LIB_LIS = -L${LISROOT}/lib/ -llis

# FastHydrology: subglacial hydrology library. Built in-tree as a sibling
# checkout (yelmo/FastHydrology), exposing module/lib under include/. The
# directory casing matches the upstream repo and what `configme install`
# creates as the link target.
FASTHYDROROOT = FastHydrology
INC_FASTHYDRO = -I${FASTHYDROROOT}/include
LIB_FASTHYDRO = -L${FASTHYDROROOT}/include -lfasthydro

# elsa: englacial layer-tracing model, a passive-tracer backend for %trc. Built
# in-tree as a sibling checkout (yelmo/elsa), exposing module/lib under
# libelsa/include/. Directory casing matches the upstream repo and the link
# target `configme install` creates.
ELSAROOT = elsa
INC_ELSA = -I${ELSAROOT}/libelsa/include
LIB_ELSA = -L${ELSAROOT}/libelsa/include -lelsa

# tracer: Lagrangian particle-tracing model, a passive-tracer backend for %trc.
# Built in-tree as a sibling checkout (yelmo/tracer), exposing module/lib under
# libtracer/include/.
TRACERROOT = tracer
INC_TRACER = -I${TRACERROOT}/libtracer/include
LIB_TRACER = -L${TRACERROOT}/libtracer/include -ltracer

# FFTW: required transitively by FastHydrology. Built by fesm-utils into a
# sibling tree; swapped to the OpenMP variant in the openmp block below.
FFTWROOT = fesm-utils/fftw/fftw-serial
INC_FFTW = -I${FFTWROOT}/include
LIB_FFTW = -L${FFTWROOT}/lib -lfftw3 -lm

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

    LISROOT = fesm-utils/lis/lis-omp
    INC_LIS = -I${LISROOT}/include
    LIB_LIS = -L${LISROOT}/lib/ -llis

    FFTWROOT = fesm-utils/fftw/fftw-omp
    INC_FFTW = -I${FFTWROOT}/include
    LIB_FFTW = -L${FFTWROOT}/lib -lfftw3_omp -lfftw3 -lm

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

LFLAGS = $(LIB_NC) $(LIB_ELSA) $(LIB_TRACER) $(LIB_FESMUTILS) $(LIB_FASTHYDRO) $(LIB_FFTW) $(LIB_LINEAR) $(LFLAGS_EXTRA)
