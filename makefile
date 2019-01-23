#------------------------------------------------------------------
#| Makefile to create the executables fro the .f90 files in the main directory | 
#------------------------------------------------------------------
#| Uses intels MKL lib (uses lapack and blas)			  | 
#------------------------------------------------------------------
#-----------------------------------------------------------------------
#              Default:   Intel Fortran  compiler 
#              		LINUX
# ----------------------------------------------------------------------
<<<<<<< f8ecfa0a313073c6540362884910bffeee3ffa04
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
=======
FC=ifort
FFLAGS= -g -check all -fp-stack-check -heap-arrays 
##         Intel's math kernel library, for LINUX
#LIBS= -L/opt/intel/Compiler/11.0/069/mkl/lib/em64t/ 
#LIBS= -llapack -lblas 
LIBS= -mkl -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -liomp5 
<<<<<<< 39081309af0c5314a057e32da0b2ba81764cc1d9
PFAPATH=/nfs/users3/carlfrost/Documents/Code/Pfaffian/pfapack/fortran
>>>>>>> Finally managed to get the linking of the pfapack-library to workgit status!
=======
#PFAPATH=/nfs/users3/carlfrost/Documents/Code/Pfaffian/pfapack/fortran
PFAPATH=lib/fortran
>>>>>>> The pfaffian performance test is now finished. Currently in the process of creating particle nbr projector
# ----------------------------------------------------------------------
#              Default:   Intel Fortran  compiler 
#              		MACINTOSH	
# ----------------------------------------------------------------------
<<<<<<< f8ecfa0a313073c6540362884910bffeee3ffa04
<<<<<<< ed64cb91f3bfcf01ea739c71ab7d2bf93469f696
#FC=ifort
#FFLAGS= -g -check all -fp-stack-check -heap-arrays 
=======
#FC=ifort
<<<<<<< 39081309af0c5314a057e32da0b2ba81764cc1d9
#FLAGS= -g -check all -fp-stack-check -heap-arrays 
>>>>>>> Finally managed to get the linking of the pfapack-library to workgit status!
=======
#FFLAGS= -g -check all -fp-stack-check -heap-arrays 
>>>>>>> The pfaffian performance test is now finished. Currently in the process of creating particle nbr projector
##         Intel's math kernel library, for MAC
#MKLPATH=/opt/intel/compilers_and_libraries_2019.1.144/mac/mkl/lib 
#PFAPATH=/Users/carlfrostenson/Documents/1_UNI/1_MasterThesis/Fortran/Routines/Pfaffian/pfapack/fortran
##
#LIBS= -mkl -Wl,-rpath,$(MKLPATH) -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -liomp5 
<<<<<<< f8ecfa0a313073c6540362884910bffeee3ffa04
=======
FC=ifort
FLAGS= -g -check all -fp-stack-check -heap-arrays $(IFLAGS)
#         Intel's math kernel library, for MAC
MKLPATH=/opt/intel/compilers_and_libraries_2019.1.144/mac/mkl/lib 
PFAPATH=/Users/carlfrostenson/Documents/1_UNI/1_MasterThesis/Fortran/Routines/Pfaffian/pfapack/fortran
#
LIBS= -mkl -Wl,-rpath,$(MKLPATH) -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -liomp5 
>>>>>>> Last commit was uncomplete. Added libs and scrapbook containing old cod that might be useful
=======
>>>>>>> Finally managed to get the linking of the pfapack-library to workgit status!
# ----------------------------------------------------------------------
#            Object files and Modules                                   
# ----------------------------------------------------------------------
MOD_DIR=modules
<<<<<<< f8ecfa0a313073c6540362884910bffeee3ffa04
<<<<<<< ed64cb91f3bfcf01ea739c71ab7d2bf93469f696
PFAFLAGS=-I$(PFAPATH) -L$(PFAPATH)
=======
IFLAGS=-I$(MOD_DIR) -I$(PFAPATH) $(PFAPATH)/libpfapack.a
>>>>>>> Last commit was uncomplete. Added libs and scrapbook containing old cod that might be useful
=======
PFAFLAGS=-I$(PFAPATH) -L$(PFAPATH)
>>>>>>> Finally managed to get the linking of the pfapack-library to workgit status!
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
<<<<<<< 39081309af0c5314a057e32da0b2ba81764cc1d9
<<<<<<< ed64cb91f3bfcf01ea739c71ab7d2bf93469f696
=======
>>>>>>> The pfaffian performance test is now finished. Currently in the process of creating particle nbr projector
all: $(MAIN)
	
$(MAIN):$(OBJECTS) 
	$(FC) $(FFLAGS) $(PFAFLAGS) $(MFLAG) -o $(patsubst %.f90,%,$@) $? $@ $(LIBS) -lpfapack
%.o: %.f90
	$(FC) $(FFLAGS) $(LIBS) $(MFLAG) -c $<  -o $@
<<<<<<< 39081309af0c5314a057e32da0b2ba81764cc1d9
=======
all: $(OBJECTS) $(MAIN)
	$(FC) $(FLAGS) $(PFAFLAGS) $(MFLAG) -o BCS $(OBJECTS) $(MAIN) $(LIBS) -lpfapack
%.o: %.f90	
	$(FC) $(FLAGS) $(LIBS) $(MFLAG) -c $<  -o $@
>>>>>>> Last commit was uncomplete. Added libs and scrapbook containing old cod that might be useful
=======
>>>>>>> The pfaffian performance test is now finished. Currently in the process of creating particle nbr projector
# Utility targets
clean: 
	rm -f $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod *.o *.mod $(MAIN_TARGET)

