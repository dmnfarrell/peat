#! /usr/bin/env python
'''Calculates the electrostatic solvation energy of a ligand in water, in complex with a protein, and the difference.

The protein and ligand are supplied in separate files. 
Ensure that the ligand has partial charges assigned and is in a bound conformation
with respect to the protein structure'''
import os
import PEATSA.Core as Core
import optparse, sys

usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-p", "--protein", dest="protein",
                  help="A pdb file containing the protein structure", metavar="PROTEIN")
parser.add_option("-l", "--ligand",
                  dest="ligand",
                  help="The ligand bound to the protein in mol2 format",
                  metavar='LIGAND')		  		  		  		  		  		  
		  
(options, args) = parser.parse_args()

if options.protein is None:
        print 'PDB file must be provided'
        sys.exit(1)
		  
if options.ligand is None:
	print 'A ligand must be provided'
        sys.exit(1)

options.protein = os.path.abspath(options.protein)
options.ligand = os.path.abspath(options.ligand)

solver = Core.Programs.PBSolver()
result = solver.solvationEnergy(pdbFile=options.protein, ligandFile=options.ligand)
print '\nEnergies are in kj.mol-1\n'
print 'Components are ', zip(['Water-Ligand(Free) Energy', 'Water-Ligand(Complex) Energy','Protein-Ligand Energy', 'Desolvation'], result)
print 'Water-Ligand Solvation Energy %lf' % result[0]
print 'Protein-Ligand Solvation Energy (Pre-Org) %lf' % (result[1] + result[2])
print 'Protein-Ligand Solvation Energy (No-Org) %lf' % (result[1] + 0.5*result[2])
print 'Delta Delta Solvation (Pre-Org) %lf' % (result[2] + result[3])
print 'Delta Delta Solvation (No-Org) %lf' % (0.5*result[2] + result[3])
