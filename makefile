# makefile for the BCS-model
#
# Defining variables -------------------------------------------------
# Compiler
#FC= gfortran 
FC= ifort
# Compiler flags
MKLINCLUDE=/opt/intel/compilers_and_libraries_2019.1.144/mac/mkl/include
MKLPATH=/opt/intel/compilers_and_libraries_2019.1.144/mac/mkl/lib
# 
# gfortran compiler flags
#FLAGS=-fcheck=all -framework accelerate -J $(MOD_DIR) 
# iFort compiler flags
FLAGS= -g -check all -fp-stack-check -heap-arrays 
LIBS= -mkl -Wl,-rpath,$(MKLPATH) -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -liomp5 -I$(MOD_DIR) 
# Directories for the modules and object files
MOD_DIR=modules
OBJ_DIR=$(MOD_DIR)
MODS=-I$(MOD_DIR)
#  
MAIN= BCSmodel.f90
OBJECTS_SRCS:=$(wildcard $(OBJ_DIR)/*.f90)

OBJECTS:=$(OBJECTS_SRCS:%.f90=%.o)
#MAIN:=$(MAIN:%.f90=%.o)
#----------Actual makefile commands----------
# "make" will build all
all: $(OBJECTS) $(MAIN)
	$(FC) $(FLAGS) $(MODS) -o BCS $(OBJECTS) $(MAIN)
%.o: %.f90	
	$(FC) $(FLAGS) $(LIBS) $(MODS) -c $<  -o $@
# Utility targets
clean: 
	rm -f $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod *.o *.mod BCS
#$(OBJ_DIR)/%.o: $(OBJ_DIR)/%.f03
#	$(FC) $(FLAGS) -c $<  -o $@
