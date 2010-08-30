CPP = g++  -g -Wall
FLAGS = -fPIC -I../ -I../Energies -I../Modelling -I../Utilities -I../Structure
PLATFORM ?= linux

all:	Utils  

Utils:	ostools.o string_tools.o                              

.SUFFIXES:	.o .cpp
.cpp.o:	
	$(CPP) $(FLAGS) -c -O $< 


