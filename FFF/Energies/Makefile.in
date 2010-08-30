CPP = g++  -g -Wall
FLAGS = -fPIC -I../ -I../Energies -I../Modelling -I../Utilities -I../Structure
PLATFORM ?= linux

all:	Energies  

Energies:	energy.o Boxes.o Lennard_Jones.o Coulomb.o Titration.o

.SUFFIXES:	.o .cpp
.cpp.o:	
	$(CPP) $(FLAGS) -c -O $< 


