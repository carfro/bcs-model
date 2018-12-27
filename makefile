#------------------------------------------------------------------
#| Makefile to create the executables fro the .f90 files in the main directory | 
#------------------------------------------------------------------
#| Uses intels MKL lib (uses lapack and blas)			  | 
#------------------------------------------------------------------
#-----------------------------------------------------------------------
#              Default:   Intel Fortran  compiler 
#              		LINUX
# ----------------------------------------------------------------------
<<<<<<< ed64cb91f3bfcf01ea739c71ab7d2bf93469f696
FC=ifort
FFLAGS= -g -check all -fp-stack-check -heap-arrays 
##         Intel's math kernel library, for LINUX
#LIBS= -L/opt/intel/Compiler/11.0/069/mkl/lib/em64t/ 
#LIBS= -llapack -lblas 
LIBS= -mkl -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -liomp5 
#PFAPATH=/nfs/users3/carlfrost/Documents/Code/Pfaffian/pfapack/fortran
PFAPATH=lib/fortran
=======
#FC=ifort
#FLAGS= -g -check all -fp-stack-check -heap-arrays $(IFLAGS)
##         Intel's math kernel library, for LINUX
#LIBS= -L/opt/intel/Compiler/11.0/069/mkl/lib/em64t/ 
>>>>>>> Last commit was uncomplete. Added libs and scrapbook containing old cod that might be useful
# ----------------------------------------------------------------------
#              Default:   Intel Fortran  compiler 
#              		MACINTOSH	
# ----------------------------------------------------------------------
<<<<<<< ed64cb91f3bfcf01ea739c71ab7d2bf93469f696
#FC=ifort
#FFLAGS= -g -check all -fp-stack-check -heap-arrays 
##         Intel's math kernel library, for MAC
#MKLPATH=/opt/intel/compilers_and_libraries_2019.1.144/mac/mkl/lib 
#PFAPATH=/Users/carlfrostenson/Documents/1_UNI/1_MasterThesis/Fortran/Routines/Pfaffian/pfapack/fortran
##
#LIBS= -mkl -Wl,-rpath,$(MKLPATH) -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -liomp5 
=======
FC=ifort
FLAGS= -g -check all -fp-stack-check -heap-arrays $(IFLAGS)
#         Intel's math kernel library, for MAC
MKLPATH=/opt/intel/compilers_and_libraries_2019.1.144/mac/mkl/lib 
PFAPATH=/Users/carlfrostenson/Documents/1_UNI/1_MasterThesis/Fortran/Routines/Pfaffian/pfapack/fortran
#
LIBS= -mkl -Wl,-rpath,$(MKLPATH) -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -liomp5 
>>>>>>> Last commit was uncomplete. Added libs and scrapbook containing old cod that might be useful
# ----------------------------------------------------------------------
#            Object files and Modules                                   
# ----------------------------------------------------------------------
MOD_DIR=modules
<<<<<<< ed64cb91f3bfcf01ea739c71ab7d2bf93469f696
PFAFLAGS=-I$(PFAPATH) -L$(PFAPATH)
=======
IFLAGS=-I$(MOD_DIR) -I$(PFAPATH) $(PFAPATH)/libpfapack.a
>>>>>>> Last commit was uncomplete. Added libs and scrapbook containing old cod that might be useful
MFLAG=-module $(MOD_DIR)
#  
OBJ_DIR=$(MOD_DIR)
OBJECTS:=$(patsubst %.f90,%.o,$(wildcard $(OBJ_DIR)/*.f90))
#
MAIN=$(wildcard *.f90)
#MAIN=$(patsubst %.f90,%,$(wildcard *.f90))
# ----------------------------------------------------------------------
#            Targets
# ----------------------------------------------------------------------
# "make" will build all
<<<<<<< ed64cb91f3bfcf01ea739c71ab7d2bf93469f696
all: $(MAIN)
	
$(MAIN):$(OBJECTS) 
	$(FC) $(FFLAGS) $(PFAFLAGS) $(MFLAG) -o $(patsubst %.f90,%,$@) $? $@ $(LIBS) -lpfapack
%.o: %.f90
	$(FC) $(FFLAGS) $(LIBS) $(MFLAG) -c $<  -o $@
=======
all: $(OBJECTS) $(MAIN)
	$(FC) $(FLAGS) $(LIBS) $(MFLAG) -o BCS $(OBJECTS) $(MAIN)
%.o: %.f90	
	$(FC) $(FLAGS) $(LIBS) $(MFLAG) -c $<  -o $@
>>>>>>> Last commit was uncomplete. Added libs and scrapbook containing old cod that might be useful
# Utility targets
clean: 
	rm -f $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod *.o *.mod $(MAIN_TARGET)

