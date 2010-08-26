
SHELL=/bin/sh	# hvilken shell-type
srcdir=.

INCDIR=  -I/usr/local/boost_1_35_0 -I/sw/include/boost-1_33_1/ -I/usr/include -I$(srcdir)  # include mapper
LIBDIR = -L/usr/lib -lm #MPI mappe

#CC = g++ -pg # for profiling
CC = g++ 


#CFLAGS = ${INCDIR} -fast -xipo -xarch=v8plusb -xchip=ultra4
#CFLAGS = ${INCDIR} -g0 -Wall -pipe -ansi
#CFLAGS = ${INCDIR} 




ifeq ($(MAKECMDGOALS),isomerase)
	CFLAGS = ${INCDIR} -O3 -ffast-math -funroll-loops -fomit-frame-pointer  #standard optimization
	TARGET_DIR = /home/people/chresten/run/.
else
	ifeq ($(MAKECMDGOALS),peat)
		CFLAGS = ${INCDIR} -O3 -ffast-math -funroll-loops -fomit-frame-pointer 
		TARGET_DIR = /usr/local/bin/.
	else
		CFLAGS = ${INCDIR} -O3 -ffast-math -mtune=i686 -fomit-frame-pointer -funroll-loops  #optimized for desktops   -ggdb
		TARGET_DIR = /usr/local/bin/.
	endif
endif

methods = simplelj.o distances.o peoe.o bondmatrix.o method.o complementarity.o nnpot.o targets.o nnmodel.o elementsnn.o prepinreader.o distancenn.o energy.o vdw.o es.o pdbwriter.o vdw_flat.o vdw_line.o hbonds.o protein_entropy.o generalparameters.o ligand_entropy.o correct_atomnames.o nn_es.o nn_vdw.o nn_hbond.o elementsnn2.o topology.o nn_ls.o metal.o hp.o nn_fitter.o bfactor.o water_entropy.o up_complementarity.o genericnn.o nn_discrimination.o mutant_stability_change.o write_file.o matrix_io.o desolvation.o fragmental_volume.o simple_parameters.o chemical_shift_changes.o reverse_chemical_shift_changes.o backbone_entropy.o water_bridges.o crystallized_water_bridges.o Boxes.o 
tools = NN.o


all: compile install

isomerase: compile install

peat: compile install

compile: algebra.o soupmanip.o stringtools.o atom.o soup.o soupobject.o water.o ion.o molecule.o residue.o proteinchain.o control.o taskmaster.o main.o $(methods) $(tools)
	${CC} -o Uffbaps algebra.o soupmanip.o stringtools.o atom.o soup.o soupobject.o water.o ion.o molecule.o residue.o proteinchain.o control.o taskmaster.o main.o $(methods) $(tools) $(LIBDIR)


atom.o: atom.cpp atom.h
	${CC} ${CFLAGS} -c atom.cpp

stringtools.o: stringtools.cpp stringtools.h
	${CC} ${CFLAGS} -c stringtools.cpp

soupmanip.o: soupmanip.cpp soupmanip.h
	${CC} ${CFLAGS} -c soupmanip.cpp

soup.o: soup.cpp soup.h
	${CC} ${CFLAGS} -c soup.cpp

soupobject.o: soupobject.cpp soupobject.h
	${CC} ${CFLAGS} -c soupobject.cpp

water.o: water.cpp water.h
	${CC} ${CFLAGS} -c water.cpp

molecule.o: molecule.cpp molecule.h
	${CC} ${CFLAGS} -c molecule.cpp

ion.o: ion.cpp ion.h
	${CC} ${CFLAGS} -c ion.cpp

residue.o: residue.cpp residue.h
	${CC} ${CFLAGS} -c residue.cpp

proteinchain.o: proteinchain.cpp proteinchain.h
	${CC} ${CFLAGS} -c proteinchain.cpp

control.o: control.cpp control.h
	${CC} ${CFLAGS} -c control.cpp

taskmaster.o: taskmaster.cpp taskmaster.h
	${CC} ${CFLAGS} -c taskmaster.cpp

main.o: main.cpp main.h
	${CC} ${CFLAGS} -c main.cpp


simplelj.o: methods/simplelj.cpp methods/simplelj.h
	${CC} ${CFLAGS} -c methods/simplelj.cpp


peoe.o: methods/peoe.cpp methods/peoe.h
	${CC} ${CFLAGS} -c methods/peoe.cpp

bondmatrix.o: methods/bondmatrix.cpp methods/bondmatrix.h
	${CC} ${CFLAGS} -c methods/bondmatrix.cpp

distances.o: methods/distances.cpp methods/distances.h
	${CC} ${CFLAGS} -c methods/distances.cpp

method.o: methods/method.cpp methods/method.h
	${CC} ${CFLAGS} -c methods/method.cpp

NN.o:   tools/NN/NN.cpp tools/NN//NN.h
	${CC} ${CFLAGS} -c tools/NN/NN.cpp

complementarity.o: methods/complementarity.cpp methods/complementarity.h
	${CC} ${CFLAGS} -c methods/complementarity.cpp

nnpot.o: methods/nnpot.cpp methods/nnpot.h
	${CC} ${CFLAGS} -c methods/nnpot.cpp

targets.o: methods/targets.cpp methods/targets.h
	${CC} ${CFLAGS} -c methods/targets.cpp

nnmodel.o: methods/nnmodel.cpp methods/nnmodel.h
	${CC} ${CFLAGS} -c methods/nnmodel.cpp

elementsnn.o: methods/elementsnn.cpp methods/elementsnn.h
	${CC} ${CFLAGS} -c methods/elementsnn.cpp

distancenn.o: methods/distancenn.cpp methods/distancenn.h
	${CC} ${CFLAGS} -c methods/distancenn.cpp

prepinreader.o: methods/prepinreader.cpp methods/prepinreader.h
	${CC} ${CFLAGS} -c methods/prepinreader.cpp

energy.o: methods/energy.cpp methods/energy.h
	${CC} ${CFLAGS} -c methods/energy.cpp

vdw.o: methods/vdw.cpp methods/vdw.h
	${CC} ${CFLAGS} -c methods/vdw.cpp

es.o: methods/es.cpp methods/es.h
	${CC} ${CFLAGS} -c methods/es.cpp

pdbwriter.o: methods/pdbwriter.cpp methods/pdbwriter.h
	${CC} ${CFLAGS} -c methods/pdbwriter.cpp

vdw_flat.o: methods/vdw_flat.cpp methods/vdw_flat.h
	${CC} ${CFLAGS} -c methods/vdw_flat.cpp

vdw_line.o: methods/vdw_line.cpp methods/vdw_line.h
	${CC} ${CFLAGS} -c methods/vdw_line.cpp

hbonds.o: methods/hbonds.cpp methods/hbonds.h
	${CC} ${CFLAGS} -c methods/hbonds.cpp

protein_entropy.o: methods/protein_entropy.cpp methods/protein_entropy.h
	${CC} ${CFLAGS} -c methods/protein_entropy.cpp

generalparameters.o: methods/generalparameters.cpp methods/generalparameters.h
	${CC} ${CFLAGS} -c methods/generalparameters.cpp

ligand_entropy.o: methods/ligand_entropy.cpp methods/ligand_entropy.h
	${CC} ${CFLAGS} -c methods/ligand_entropy.cpp

correct_atomnames.o: methods/correct_atomnames.cpp methods/correct_atomnames.h
	${CC} ${CFLAGS} -c methods/correct_atomnames.cpp

nn_es.o: methods/nn_es.cpp methods/nn_es.h
	${CC} ${CFLAGS} -c methods/nn_es.cpp

nn_vdw.o: methods/nn_vdw.cpp methods/nn_vdw.h
	${CC} ${CFLAGS} -c methods/nn_vdw.cpp

nn_hbond.o: methods/nn_hbond.cpp methods/nn_hbond.h
	${CC} ${CFLAGS} -c methods/nn_hbond.cpp

elementsnn2.o: methods/elementsnn2.cpp methods/elementsnn2.h
	${CC} ${CFLAGS} -c methods/elementsnn2.cpp

topology.o: methods/topology.cpp methods/topology.h
	${CC} ${CFLAGS} -c methods/topology.cpp

nn_ls.o: methods/nn_ls.cpp methods/nn_ls.h
	${CC} ${CFLAGS} -c methods/nn_ls.cpp

metal.o: methods/metal.cpp methods/metal.h
	${CC} ${CFLAGS} -c methods/metal.cpp

hp.o: methods/hp.cpp methods/hp.h
	${CC} ${CFLAGS} -c methods/hp.cpp

nn_fitter.o: methods/nn_fitter.cpp methods/nn_fitter.h
	${CC} ${CFLAGS} -c methods/nn_fitter.cpp

bfactor.o: methods/bfactor.cpp methods/bfactor.h
	${CC} ${CFLAGS} -c methods/bfactor.cpp

water_entropy.o: methods/water_entropy.cpp methods/water_entropy.h
	${CC} ${CFLAGS} -c methods/water_entropy.cpp

up_complementarity.o: methods/up_complementarity.cpp methods/up_complementarity.h
	${CC} ${CFLAGS} -c methods/up_complementarity.cpp

genericnn.o: methods/genericnn.cpp methods/genericnn.h
	${CC} ${CFLAGS} -c methods/genericnn.cpp

nn_discrimination.o: methods/nn_discrimination.cpp methods/nn_discrimination.h
	${CC} ${CFLAGS} -c methods/nn_discrimination.cpp

mutant_stability_change.o: methods/mutant_stability_change.cpp methods/mutant_stability_change.h
	${CC} ${CFLAGS} -c methods/mutant_stability_change.cpp

write_file.o: methods/write_file.cpp methods/write_file.h
	${CC} ${CFLAGS} -c methods/write_file.cpp

matrix_io.o: methods/matrix_io.cpp methods/matrix_io.h
	${CC} ${CFLAGS} -c methods/matrix_io.cpp

desolvation.o: methods/desolvation.cpp methods/desolvation.h
	${CC} ${CFLAGS} -c methods/desolvation.cpp

fragmental_volume.o: methods/fragmental_volume.cpp methods/fragmental_volume.h
	${CC} ${CFLAGS} -c methods/fragmental_volume.cpp

simple_parameters.o: methods/simple_parameters.cpp methods/simple_parameters.h
	${CC} ${CFLAGS} -c methods/simple_parameters.cpp

chemical_shift_changes.o: methods/chemical_shift_changes.cpp methods/chemical_shift_changes.h
	${CC} ${CFLAGS} -c methods/chemical_shift_changes.cpp

reverse_chemical_shift_changes.o: methods/reverse_chemical_shift_changes.cpp methods/reverse_chemical_shift_changes.h
	${CC} ${CFLAGS} -c methods/reverse_chemical_shift_changes.cpp

backbone_entropy.o: methods/backbone_entropy.cpp methods/backbone_entropy.h
	${CC} ${CFLAGS} -c methods/backbone_entropy.cpp

crystallized_water_bridges.o: methods/crystallized_water_bridges.cpp methods/crystallized_water_bridges.h
	${CC} ${CFLAGS} -c methods/crystallized_water_bridges.cpp

water_bridges.o: methods/water_bridges.cpp methods/water_bridges.h
	${CC} ${CFLAGS} -c methods/water_bridges.cpp

algebra.o: algebra.cpp algebra.h
	${CC} ${CFLAGS} -c algebra.cpp

Boxes.o:   methods/Boxes.cpp methods/Boxes.h
	${CC} ${CFLAGS} -c methods/Boxes.cpp


clean: 	
	rm -f *.o *.bak *~ \#*\#  methods/*~ 

clean-bck:
	rm -f *~ *.bak \#*\#

install:
	mv -f ./Uffbaps ${TARGET_DIR}








