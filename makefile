#------------------------------------------------------------------
#| Makefile to create the BCS-model executable			  | 
#------------------------------------------------------------------
#| Requires the Linear algebra libraries    Lapack and Blas       |
#|                                                                |
#------------------------------------------------------------------
#
# ----------------------------------------------------------------------
#              Default:   Intel Fortran  compiler 
#              		LINUX
# ----------------------------------------------------------------------
FOR=ifort
FLAGS= -g -check all -fp-stack-check -heap-arrays 
#         Intel's math kernel library, for LINUX
LIBS= -L/opt/intel/Compiler/11.0/069/mkl/lib/em64t/ -lmkl_lapack -lmkl_core -lmkl_em64t -lmkl_intel_thread -lmkl_intel_lp64 -liomp5
# 
# ----------------------------------------------------------------------
#              Default:   Intel Fortran  compiler 
#              		MACITOSH	
# ----------------------------------------------------------------------
#FOR=ifort
#FLAGS= -g -check all -fp-stack-check -heap-arrays 
#         Intel's math kernel library, for MAC
#MKLPATH=/opt/intel/compilers_and_libraries_2019.1.144/mac/mkl/lib
#LIBS= -mkl -Wl,-rpath,$(MKLPATH) -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -liomp5 -I$(MOD_DIR) 
# ----------------------------------------------------------------------
#            Object files and Modules                                   
# ----------------------------------------------------------------------
MOD_DIR=modules
MODS=-I$(MOD_DIR)
#  
OBJ_DIR=$(MOD_DIR)
OBJECTS_SRCS:=$(wildcard $(OBJ_DIR)/*.f90)
OBJECTS:=$(OBJECTS_SRCS:%.f90=%.o)
#
MAIN= BCSmodel.f90
# ----------------------------------------------------------------------
#            Targets
# ----------------------------------------------------------------------
# "make" will build all
all: $(OBJECTS) $(MAIN)
	$(FC) $(FLAGS) $(MODS) -o BCS $(OBJECTS) $(MAIN)
%.o: %.f90	
	$(FC) $(FLAGS) $(LIBS) $(MODS) -c $<  -o $@
# Utility targets
clean: 
	rm -f $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod *.o *.mod BCS

