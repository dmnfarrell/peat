CPP = g++  -g -Wall
FLAGS = -fPIC -I../ -I../Energies -I../Modelling -I../Utilities -I../Structure
PLATFORM ?= linux

all:	Structure

Structure:	atom_class.o chain_class.o fff.o residue_class.o 

.SUFFIXES:	.o .cpp
.cpp.o:	
	$(CPP) $(FLAGS) -c -O $< 


