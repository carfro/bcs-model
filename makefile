# makefile for the BCS-model
#
# Defining variables -------------------------------------------------
#Compiler
FC= gfortran 
# Compiler flags
FLAGS=-fcheck=all -framework accelerate -J $(MOD_DIR) 
# Directories for the modules and object files
MOD_DIR=modules
OBJ_DIR=modules
#  
MAIN= BCSmodel.f03
OBJECTS_SRCS:=$(wildcard $(OBJ_DIR)/*.f03)

OBJECTS:=$(OBJECTS_SRCS:%.f03=%.o)
MAIN:=$(MAIN:%.f03=%.o)
#----------Actual makefile commands----------
# "make" will build all
all: $(OBJECTS) $(MAIN)
	$(FC) $(FLAGS) -o BCS $(OBJECTS) $(MAIN)
%.o: %.f03	
	$(FC) $(FLAGS) -c $<  -o $@
# Utility targets
clean: 
	rm -f $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod *.o *.mod
#$(OBJ_DIR)/%.o: $(OBJ_DIR)/%.f03
#	$(FC) $(FLAGS) -c $<  -o $@
