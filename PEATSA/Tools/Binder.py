#! /bin/env python
#
# Protein Engineering Analysis Tool Structure Analysis (PEATSA)
# Copyright (C) 2010 Michael Johnston & Jens Erik Nielsen
#
# Author: Michael Johnston
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact information:
# Email: Jens.Nielsen_at_gmail.com
# Normal mail:
# Jens Nielsen
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
#

#
#  Binder.py
#  ProteinDesignTool
#
#  Created by Michael Johnston on 15/07/2009.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
#
'''Calculates the effect of mutations on binding free-energies'''
import PEATSA.Core as Core
import os, operator
import optparse, sys

usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-b", "--bound-conf", dest="bound",
                  help="The bound configuration of the protein in pdb format.", metavar="BOUND")
parser.add_option("-u", "--unbound-conf",
                   dest="unbound",
                  help="The unbound configuration of the protein - optional",
                  metavar='UNBOUND')
parser.add_option("-l", "--ligand",
                  dest="ligand",
                  help="Mol2 file of the ligand",
                  metavar='LIGAND')
parser.add_option("-m", "--mutationList",
			dest="mutationList",
			help="Mutations to create",
			metavar='MUTATION-LIST')

parser.add_option("-p", "--prefix",
			dest="prefix", default="",
			help="Prefix for the results directories",
			metavar='PREFIX')
		  		  		  		  
(options, args) = parser.parse_args()

if options.bound is None:
        print 'Bound configuration must be provided'
        sys.exit(1)
		  
if options.ligand is None:
        print 'Ligand must be provided'
        sys.exit(1)
	
if options.mutationList is None:
        print 'Mutation list must be provided'
        sys.exit(1)	
	
configuration = Core.Environment.Configuration(writeFile=False)
workingDirectory = Core.Environment.WorkingDirectory(location='.')
environment = Core.Environment.Environment()

#Clean the pdbs

boundStructure = os.path.abspath(options.bound)
ligand = os.path.abspath(options.ligand)
	
#Parse mutation list
mutationListFile = Core.Data.MutationListFile(options.mutationList)

#Create output dir
output = Core.Data.DataSet(name="%sBinderResults" % options.prefix, location=".", overwrite=True)

#Create bound mutants
environment.output("[MONGREL] Creating bound mutants")
boundMutants = Core.Data.MutantCollection(pdbFile=boundStructure,
					name="%sBoundMutants" % options.prefix,
					mutationList=mutationListFile.mutantList(),
					ligandFiles=[ligand],
					location=".",
					maxOverlap=configuration.getfloat('PKA SCAN PARAMETERS', 'mutation_quality'),
					clean=True)
environment.output("[MONGREL] Complete - %s" % boundMutants)

#Calculate est binding component
#First the wild-type
pbSolver = Core.Programs.PBSolver()
environment.output("\n[MONGREL] Calculating wild-type electrostatic binding component")
wildTypeSolvationEnergy = pbSolver.solvationEnergy(pdbFile=boundStructure, ligandFile=ligand)
environment.output("[MONGREL] Result %s" % str(wildTypeSolvationEnergy))

#Now the mutants
environment.output("\n[MONGREL] Calculating mutant electrostatic binding components")
mutantSolvationEnergies = pbSolver.mutantSolvationEnergies(mutantCollection=boundMutants, ligandFile=ligand)
newRow = ['WT']
newRow.extend(wildTypeSolvationEnergy)
mutantSolvationEnergies.addRow(newRow, index=0)
output.addMatrix(mutantSolvationEnergies, name="AbsoluteSolvationEnergies")
environment.output("[MONGREL] Done")

#Calculate nonbonded binding components
#Use boundMutants.pdb() since that is cleaned.
environment.output("\n[MONGREL] Calculating non-polar binding component")
uffbaps = Core.Programs.UFFBAPS(sourcePDB=boundMutants.pdbFile,
					executablePath=configuration.uffbapsLocation(),
					workingDirectory=workingDirectory)
uffbaps.binding(comparePDBs=boundMutants.mutantFiles(), ligand=ligand, verbose=True)	
nonPolarResults = uffbaps.deltaBindingResults()
output.addMatrix(matrix=nonPolarResults, name="DeltaNonPolar")
environment.output("\n[MONGREL] Done")

if options.unbound is not None:
	unboundStructure = os.path.abspath(options.unbound)
	#Calculate reorganisation
	environment.output("\n[MONGREL] Calculating changes in reorganisation energy")
	environment.output("[MONGREL] Creating unbound mutants")
	unboundMutants = Core.Data.MutantCollection(pdbFile=unboundStructure,
					name="%sUnboundMutants" % options.prefix,
					mutationList=mutationListFile.mutantList(),
					ligandFiles=None,
					location=".",
					maxOverlap=configuration.getfloat('PKA SCAN PARAMETERS', 'mutation_quality'),
					clean=True)
	environment.output("[MONGREL] Complete - %s" % unboundMutants)

	environment.output("[MONGREL] Calculating change in unbound stability")
	uffbaps = Core.Programs.UFFBAPS(sourcePDB=unboundMutants.pdbFile,
						executablePath=configuration.uffbapsLocation(),
						workingDirectory=workingDirectory)
						
	uffbaps.deltaStability(comparePDBs=unboundMutants.mutantFiles(), verbose=True)	
	unboundStabilityResults = uffbaps.deltaStabilityResults()
	output.addMatrix(matrix=unboundStabilityResults, name="UnboundStabilityResults")
	
	environment.output("\n[MONGREL] Calculating change in bound stability")
	uffbaps = Core.Programs.UFFBAPS(sourcePDB=boundMutants.pdbFile,
						executablePath=configuration.uffbapsLocation(),
						workingDirectory=workingDirectory)
						
	uffbaps.deltaStability(comparePDBs=boundMutants.mutantFiles(), verbose=True)	
	boundStabilityResults = uffbaps.deltaStabilityResults()
	output.addMatrix(matrix=boundStabilityResults, name="BoundStabilityResults")
	
	environment.output("\n[MONGREL] Finished reorganisation calculation")

else:
	environment.output('[MONGREL] No unbound structure provided - not calculating reorganisation energy')

preorgResult = None
noorgResult = None
if environment.isRoot():

	print '[MONGREL] Combining results ...\n'
	
	#Ligand solvation energy - Use value returned by wt calc
	ligandSolvation = wildTypeSolvationEnergy[0]
	
	#Wild-Type solvation Data
	wildTypeInteractionEnergy = wildTypeSolvationEnergy[2]
	wildTypeProteinWaterSolvation = wildTypeSolvationEnergy[1]
	wildTypePreorgCharging = wildTypeInteractionEnergy + wildTypeProteinWaterSolvation - ligandSolvation
	wildTypeNoorgCharging = 0.5*wildTypeInteractionEnergy + wildTypeProteinWaterSolvation - ligandSolvation
	print 'Widl Typ Pre org', wildTypePreorgCharging
	print 'Widl Typ No org', wildTypeNoorgCharging
	
	preorg = []
	noorg = []
	vanDerWaals = nonPolarResults.indexOfColumnWithHeader('Van der Waals')
	ligandEntropy = nonPolarResults.indexOfColumnWithHeader('Ligand Entropy')
	for mutationSet in mutationListFile.mutantList():
		#First calculate change in protein-ligand est interaction energy
		try:
			print '[MONGREL] %s - Calculating change in binding due to change in electrostatics' % mutationSet
			mutantData = mutantSolvationEnergies.dataForMutationSet(mutationSet)
			interactionEnergy = mutantData[3]
			proteinWaterSolvation = mutantData[2]
			proteinPreorgCharging = interactionEnergy + proteinWaterSolvation - ligandSolvation
			proteinNoorgCharging = 0.5*interactionEnergy + proteinWaterSolvation - ligandSolvation
			print proteinPreorgCharging
			print proteinNoorgCharging
			#Calculate absolute and delta values
			#Skip absolute for now
			diffPreorg =  proteinPreorgCharging - wildTypePreorgCharging
			diffNoorg = proteinNoorgCharging - wildTypeNoorgCharging
			print 'Mutant diff Pre org %lf (%lf %lf)' % (diffPreorg, wildTypePreorgCharging, proteinPreorgCharging)
			print 'Mutant diff No org %lf (%lf %lf)' % (diffNoorg, wildTypeNoorgCharging, proteinNoorgCharging)
		except IndexError, data:
			print '[MONGREL] No Solvation Data available for mutant %s - Skipping' % mutationSet
			continue
	
		#Now non-polar components
		try:
			print '[MONGREL] %s - Calculating change in binding due to change in non-polar interactions' % mutationSet
			mutantData = nonPolarResults.dataForMutationSet(mutationSet)
			deltaVanDerWaals = mutantData[vanDerWaals]
			deltaLigandEntropy = mutantData[ligandEntropy]
			deltaNonPolar = deltaVanDerWaals + deltaLigandEntropy
			result = [diffPreorg, deltaVanDerWaals, deltaLigandEntropy]
		except IndexError, data:
			print '[MONGREL] No non-polar data available for mutant %s - Skipping' % mutationSet
			continue	
			
		#The reorganisation component if present
		if options.unbound is not None:
			try:
				print '[MONGREL] %s - Calculating change in binding due to protein stability' % mutationSet
				boundData = boundStabilityResults.dataForMutationSet(mutationSet)
				unboundData = unboundStabilityResults.dataForMutationSet(mutationSet)
				reorganisation = boundData[-1] - unboundData[-1]
				result.append(reorganisation)
			except IndexError, data:
				print '[MONGREL] No reorganisation data available for mutant %s - Skipping' % mutationSet
				continue
				
		#Combine the results		
		total = reduce(operator.add, result)
		result.append(total)
		result.insert(0, mutationSet.codeString(pdb=boundMutants.pdb))
		print '[MONGREL] %s Pre-org ' % mutationSet, result
		preorg.append(result)
		
		result = [diffNoorg, deltaVanDerWaals, deltaLigandEntropy]
		if options.unbound is not None:
			result.append(reorganisation)

		total = reduce(operator.add, result)
		result.append(total)
		result.insert(0, mutantData[0])
		print '[MONGERL] %s No-org ' % mutationSet, result
		noorg.append(result)

	headers = ["Mutations", "EST(Pre-Org)", "VanDerWaals", "LigandEntropy", "Total"]
	if options.unbound is not None:
		headers.insert(4, 'Reorganisation')
		
	preorgResult = Core.Matrix.Matrix(rows=preorg, headers=headers)
	headers[1] = "EST(No-Org)"
	noorgResult = Core.Matrix.Matrix(rows=noorg, headers=headers)

output.addMatrix(preorgResult, name="PreorgResults")
output.addMatrix(noorgResult, name="NoorgResults")


