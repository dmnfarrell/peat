CPP = g++  -g -Wall
FLAGS = -fPIC -I../ -I../Energies -I../Modelling -I../Utilities -I../Structure
PLATFORM ?= linux

all:	Modelling  

Modelling:	Model.o ModelIII.o Model_global.o parse_aadef.o quatfit.o quatfit_class.o Rotamers.o assign.o parse_hydrogens.o

.SUFFIXES:	.o .cpp
.cpp.o:	
	$(CPP) $(FLAGS) -c -O $< 


