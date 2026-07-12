.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = src
objdir = libyelmo/include
bindir = libyelmo/bin
libdir = libs
testdir = tests

# Command-line options at make call
debug  ?= 0
openmp ?= 0 
petsc  ?= 0 

## COMPILER CONFIGURATION ##
# Loaded from a config/<host> fragment by `python config.py config/<host>`,
# or assembled by `configme` in the modern path. The fragment defines the
# compiler (FC, FFLAGS, FFLAGS_OPENMP, DFLAGS_*), netCDF wiring (INC_NC,
# LIB_NC), and any host-specific overrides (e.g. LFLAGS_EXTRA).

<COMPILER_CONFIGURATION>

## SHARED DEPENDENCY WIRING ##
# Pulled in after the compiler/machine fragment so it can reference FFLAGS,
# FFLAGS_OPENMP, LIB_NC, and LFLAGS_EXTRA set there. Defines the FESMUTILS,
# LIS, FASTHYDRO, FFTW, and (optional) PETSc paths, the openmp= / petsc=
# build switches, and the final LFLAGS.

include config/common.mk

# Determine whether to use normal flags or debugging flags
DFLAGS = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS = $(DFLAGS_PROFILE)
endif

###############################################
##							
## List of yelmo rules and source files
##
###############################################

include config/Makefile_yelmo.mk

###############################################
##
## Compilation of complete programs
##
###############################################

check-hash:
	./check_githash > git_yelmo.txt

# FastHydrology's static archive provides fast_hydrology.mod, which yelmo's
# sources `use` (transitively via yelmo_defs, which declares a hydro_class
# member of yelmo_class), so this must build before yelmo-static.
fasthydro-static:
	$(MAKE) -C $(FASTHYDROROOT) fasthydro-static openmp=$(openmp)

# elsa/tracer static archives provide the .mod files yelmo `use`s (transitively
# via yelmo_defs, which declares elsa_class/tracer_class members of ytrc_class),
# so they must build before yelmo-static.
elsa-static:
	$(MAKE) -C $(ELSAROOT) elsa-static openmp=$(openmp)

tracer-static:
	$(MAKE) -C $(TRACERROOT) tracer-static openmp=$(openmp)

# Static library compilation
yelmo-static: fasthydro-static elsa-static tracer-static $(yelmo_libs) $(yelmo_physics) $(yelmo_base) $(yelmo_tests)
	ar rc $(objdir)/libyelmo.a $(yelmo_libs) $(yelmo_physics) $(yelmo_base) $(yelmo_tests)
	ranlib $(objdir)/libyelmo.a
	./check_githash > git_yelmo.txt
	@echo " "
	@echo "    $(objdir)/libyelmo.a is ready."
	@echo " "

yelmo-c: yelmo-static
	gfortran -shared -fPIC $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(srcdir)/yelmo_c_api.f90 -o $(objdir)/libyelmo_c_api.so \
	-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
	@echo " "
	@echo "    $(objdir)/libyelmo_c_api.so is ready."
	@echo " "

# Test programs that use yelmo-static
initmip : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/yelmo_initmip.x tests/yelmo_initmip.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    yelmo_initmip.x is ready."
		@echo " "

benchmarks : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/yelmo_benchmarks.x tests/yelmo_benchmarks.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    yelmo_benchmarks.x is ready."
		@echo " "

calving : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/yelmo_calving.x tests/yelmo_calving.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    yelmo_calving.x is ready."
		@echo " "

symmetry : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/test_symmetry.x tests/test_symmetry.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    test_symmetry.x is ready."
		@echo " "

mismip : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/yelmo_mismip.x tests/yelmo_mismip.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    yelmo_mismip.x is ready."
		@echo " "

ismiphom : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/yelmo_ismiphom.x tests/yelmo_ismiphom.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    yelmo_ismiphom.x is ready."
		@echo " "

regridding : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/yelmo_ismiphom.x tests/yelmo_test_regridding.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    yelmo_test_regridding.x is ready."
		@echo " "

icetemp : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/test_icetemp.x tests/test_icetemp.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    test_icetemp.x is ready."
		@echo " "

enthalpy : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/test_enthalpy.x tests/test_enthalpy.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    test_enthalpy.x is ready."
		@echo " "

variables : tests/test_variables.f90
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_variables.x tests/test_variables.f90 \
			$(LFLAGS)
		@echo " "
		@echo "    test_variables.x is ready."
		@echo " "

ssa_energy : yelmo-static tests/test_ssa_energy.f90
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/test_ssa_energy.x tests/test_ssa_energy.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    test_ssa_energy.x is ready."
		@echo " "

ssa_energy_sym : yelmo-static tests/test_ssa_energy_sym.f90
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/test_ssa_energy_sym.x tests/test_ssa_energy_sym.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    test_ssa_energy_sym.x is ready."
		@echo " "

ssa_energy_lr_sym : yelmo-static tests/test_ssa_energy_lr_sym.f90
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/test_ssa_energy_lr_sym.x tests/test_ssa_energy_lr_sym.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    test_ssa_energy_lr_sym.x is ready."
		@echo " "

opt : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/yelmo_opt.x tests/yelmo_opt.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    yelmo_opt.x is ready."
		@echo " "

trough : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/yelmo_trough.x tests/yelmo_trough.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    yelmo_trough.x is ready."
		@echo " "

slab : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/yelmo_slab.x tests/yelmo_slab.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    yelmo_slab.x is ready."
		@echo " "

mask_ice : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) -o $(bindir)/yelmo_mask_ice.x tests/yelmo_mask_ice.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS)
		@echo " "
		@echo "    yelmo_mask_ice.x is ready."
		@echo " "

check_sim :
		$(FC) $(DFLAGS) $(INC_FESMUTILS) $(FFLAGS) -o ./check_sim.x check_sim.f90 \
			$(LFLAGS)
		@echo " "
		@echo "    check_sim.x is ready."
		@echo " "

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make initmip    : compiles yelmo_initmip.x, for running yelmo on a given domain defined in param file."
	@echo " make benchmarks : compiles yelmo_benchmarks.x, for running various EISMINT domain experiments.""
	@echo " make mismip     : compiles yelmo_mismip.x, for running various MISMIP domain experiments.""
	@echo " make opt        : compiles yelmo_opt.x, for running various MISMIP domain experiments.""
	@echo " make trough     : compiles yelmo_trough.x, for running trough domain experiments.""
	@echo " make clean      : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
	rm -f ./git_yelmo.txt
	$(MAKE) -C $(FASTHYDROROOT) clean


# Experimental
yelmo-shared: $(yelmo_libs) $(yelmo_physics) $(yelmo_base) $(yelmo_tests)
	$(FC) -shared -o $(objdir)/libyelmo.so -fPIC $(yelmo_libs) $(yelmo_physics) $(yelmo_base) $(yelmo_tests) \
		$(LFLAGS) $(INC_LIS) 
	./check_githash > git_yelmo.txt
	@echo " "
	@echo "    $(objdir)/libyelmo.so is ready."
	@echo " "
