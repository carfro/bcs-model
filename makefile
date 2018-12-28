#------------------------------------------------------------------
#| Makefile to create the BCS-model executable			  | 
#------------------------------------------------------------------
#| Uses intels MKL lib (uses lapack and blas)			  | 
#------------------------------------------------------------------
#-----------------------------------------------------------------------
#              Default:   Intel Fortran  compiler 
#              		LINUX
# ----------------------------------------------------------------------
FC=ifort
FLAGS= -g -check all -fp-stack-check -heap-arrays 
##         Intel's math kernel library, for LINUX
#LIBS= -L/opt/intel/Compiler/11.0/069/mkl/lib/em64t/ 
#LIBS= -llapack -lblas 
LIBS= -mkl -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -liomp5 
PFAPATH=/nfs/users3/carlfrost/Documents/Code/Pfaffian/pfapack/fortran
# ----------------------------------------------------------------------
#              Default:   Intel Fortran  compiler 
#              		MACINTOSH	
# ----------------------------------------------------------------------
#FC=ifort
#FLAGS= -g -check all -fp-stack-check -heap-arrays 
##         Intel's math kernel library, for MAC
#MKLPATH=/opt/intel/compilers_and_libraries_2019.1.144/mac/mkl/lib 
#PFAPATH=/Users/carlfrostenson/Documents/1_UNI/1_MasterThesis/Fortran/Routines/Pfaffian/pfapack/fortran
##
#LIBS= -mkl -Wl,-rpath,$(MKLPATH) -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -liomp5 
# ----------------------------------------------------------------------
#            Object files and Modules                                   
# ----------------------------------------------------------------------
MOD_DIR=modules
PFAFLAGS=-I$(PFAPATH) -L$(PFAPATH)
MFLAG=-module $(MOD_DIR)
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
	$(FC) $(FLAGS) $(PFAFLAGS) $(MFLAG) -o BCS $(OBJECTS) $(MAIN) $(LIBS) -lpfapack
%.o: %.f90	
	$(FC) $(FLAGS) $(LIBS) $(MFLAG) -c $<  -o $@
# Utility targets
clean: 
	rm -f $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod *.o *.mod BCS

