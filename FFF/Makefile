CPP = g++  -g -Wall
FLAGS = -fPIC  -IEnergies/ -IModelling/ -IUtilities -IStructure -I./
PLATFORM ?= linux

SUBDIRS = Energies Structure Utilities Modelling
# Default Target

all:	FFF pyFFF

compile_subdirs:
	rm -f fff _FFFcontrol.so FFFcontrol.py	
	@echo	
	@echo "#######################################"
	@echo "### BUILDING ALL TARGETS ###"
	@echo "#######################################"
	@echo
	for i in $(SUBDIRS) ; do \
	( echo "Now compiling: " $$i ;cd $$i ; make ; echo "----------------") ; \
	done

clean:	
	rm -f *.o **/*.o fff _FFFcontrol.so FFFcontrol.py

FFF:	compile_subdirs main.o pKa.o
	$(CPP) $(FLAGS) **/*.o main.o pKa.o   -o FFF

pyFFF:	compile_subdirs pKa.o FFFcontrol.o
	swig -python -c++ FFFcontrol.i
	$(CPP) $(FLAGS) -c FFFcontrol_wrap.cxx -I/usr/include/python2.6/
	
ifeq ($(OSTYPE),linux)
	$(CPP)  -export-dynamic -shared **/*.o pKa.o FFFcontrol.o FFFcontrol_wrap.o -o _FFFcontrol.so
else
	$(CPP) -u _PyMac_Error -framework Python  -flat_namespace -shared pKa.o FFFcontrol.o **/*.o FFFcontrol_wrap.o -o _FFFcontrol.so
endif

.SUFFIXES:	.o .cpp
.cpp.o:	
	$(CPP) $(FLAGS) -c -O $< 

	
