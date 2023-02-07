SRC_DIR = source
MOD_DIR = modules
OBJ_DIR = objects
BIN_DIR = bin

MKDIR_P = mkdir -p

GMCTDH = $(BIN_DIR)/macaw

all: directories $(GMCTDH)
directories : $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR)

# Compiler
FC = gfortran

# Compiler options
## Optimal performance
FCFLAGS =  -O3 -J$(MOD_DIR) -fdefault-integer-8 -freal-4-real-8 -fimplicit-none 
# For debug:
#FCFLAGS = -g -Wall -fcheck=all -fbacktrace -J$(MOD_DIR) -fdefault-integer-8 -fimplicit-none -m64 

# Standard Lapack (change the path)
LIBS= /lib/lapack/cygblas-0.dll /lib/lapack/cyglapack-0.dll

#
# Dependencies
#

GMCTDH_OBJS = $(addprefix $(OBJ_DIR)/, \
	sysparam.o globals.o timingmod.o input.o psidef.o dvr.o init_calc.o hamilton.o \
	storage.o gmctdh_objects.o propagator.o derivatives.o integrator.o \
       	cmf2_base.o cmf2_integration.o cmf2_derivatives.o \
       	output.o maths.o)

$(OBJ_DIR)/sysparam.o:   $(SRC_DIR)/sysparam.f90
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@ 

$(OBJ_DIR)/timingmod.o:   $(SRC_DIR)/timingmod.f90 $(OBJ_DIR)/sysparam.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/globals.o:   $(SRC_DIR)/globals.f90  $(OBJ_DIR)/sysparam.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/psidef.o:   $(SRC_DIR)/psidef.f90  
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/dvr.o:   $(SRC_DIR)/dvr.f90  $(OBJ_DIR)/sysparam.o  $(OBJ_DIR)/psidef.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/init_calc.o:   $(SRC_DIR)/init_calc.f90  \
	$(OBJ_DIR)/sysparam.o  $(OBJ_DIR)/globals.o  $(OBJ_DIR)/psidef.o  \
	$(OBJ_DIR)/storage.o   $(OBJ_DIR)/hamilton.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/hamilton.o:   $(SRC_DIR)/hamilton.f90  \
	$(OBJ_DIR)/sysparam.o   $(OBJ_DIR)/psidef.o     $(OBJ_DIR)/globals.o \
	$(OBJ_DIR)/dvr.o        $(OBJ_DIR)/storage.o    $(OBJ_DIR)/timingmod.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/storage.o:   $(SRC_DIR)/storage.f90  \
	$(OBJ_DIR)/sysparam.o  $(OBJ_DIR)/psidef.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/gmctdh_objects.o:   $(SRC_DIR)/gmctdh_objects.f90  \
	$(OBJ_DIR)/sysparam.o $(OBJ_DIR)/storage.o  $(OBJ_DIR)/psidef.o \
	$(OBJ_DIR)/hamilton.o $(OBJ_DIR)/globals.o   $(OBJ_DIR)/timingmod.o \
	$(OBJ_DIR)/dvr.o      $(OBJ_DIR)/maths.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/output.o:   $(SRC_DIR)/output.f90   \
	$(OBJ_DIR)/sysparam.o  $(OBJ_DIR)/globals.o    $(OBJ_DIR)/psidef.o   \
	$(OBJ_DIR)/storage.o   $(OBJ_DIR)/hamilton.o   $(OBJ_DIR)/gmctdh_objects.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/propagator.o:   $(SRC_DIR)/propagator.f90  \
	$(OBJ_DIR)/globals.o      $(OBJ_DIR)/psidef.o  \
	$(OBJ_DIR)/storage.o      $(OBJ_DIR)/hamilton.o    $(OBJ_DIR)/derivatives.o  \
	$(OBJ_DIR)/integrator.o   $(OBJ_DIR)/cmf2_base.o \
	$(OBJ_DIR)/gmctdh_objects.o   $(OBJ_DIR)/output.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/integrator.o:   $(SRC_DIR)/integrator.f90  \
	$(OBJ_DIR)/sysparam.o  $(OBJ_DIR)/globals.o      $(OBJ_DIR)/psidef.o  \
	$(OBJ_DIR)/derivatives.o   $(OBJ_DIR)/dvr.o     $(OBJ_DIR)/timingmod.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/cmf2_base.o:   $(SRC_DIR)/cmf2_base.f90  \
	$(OBJ_DIR)/sysparam.o  $(OBJ_DIR)/globals.o      $(OBJ_DIR)/psidef.o  \
	$(OBJ_DIR)/cmf2_integration.o   
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/cmf2_integration.o:  $(SRC_DIR)/cmf2_integration.f90 \
	$(OBJ_DIR)/sysparam.o  $(OBJ_DIR)/globals.o  $(OBJ_DIR)/psidef.o \
	$(OBJ_DIR)/storage.o   $(OBJ_DIR)/dvr.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/cmf2_derivatives.o:  $(SRC_DIR)/cmf2_derivatives.f90 \
	$(OBJ_DIR)/sysparam.o  $(OBJ_DIR)/globals.o  $(OBJ_DIR)/psidef.o \
	$(OBJ_DIR)/hamilton.o  $(OBJ_DIR)/dvr.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

$(OBJ_DIR)/derivatives.o:   $(SRC_DIR)/derivatives.f90  \
	$(OBJ_DIR)/sysparam.o  $(OBJ_DIR)/psidef.o   $(OBJ_DIR)/storage.o  \
	$(OBJ_DIR)/hamilton.o  $(OBJ_DIR)/dvr.o      $(OBJ_DIR)/globals.o   \
	$(OBJ_DIR)/timingmod.o $(OBJ_DIR)/maths.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

# Algebra
$(OBJ_DIR)/maths.o:    $(SRC_DIR)/maths.f90 \
	$(OBJ_DIR)/sysparam.o
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@

# Input/Output
$(OBJ_DIR)/input.o:   $(SRC_DIR)/input.f90  \
	$(OBJ_DIR)/sysparam.o   $(OBJ_DIR)/globals.o   $(OBJ_DIR)/psidef.o  \
	$(OBJ_DIR)/dvr.o   $(OBJ_DIR)/hamilton.o 
	$(FC) $(FCFLAGS) $(PROFILE) -c $< -o $@


# Programs
$(GMCTDH): $(SRC_DIR)/main.f90   $(GMCTDH_OBJS)
	$(FC) $(FCFLAGS) $^ $(LIBS) -o $@

#
$(OBJ_DIR):
	$(MKDIR_P) $(OBJ_DIR)

$(BIN_DIR):
	$(MKDIR_P) $(BIN_DIR)

$(MOD_DIR):
	$(MKDIR_P) $(MOD_DIR)

clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR)

