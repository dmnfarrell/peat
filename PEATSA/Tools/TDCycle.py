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
# Five different types of free-energy differences can be specified.
# A. Stability - ΔGstability(A) = G(Fold A) - G(Unfold A)
# B. Binding - ΔGbinding(AB) = G(AB) - G(Bound A + Ligand)
# C. Interaction - Protein-Protein interaction:  ΔGpp(A + B -> AB) = ΔGstability(AB) - ΔGstability(A) - ΔGstability(B)  
# D. BindingWithChange -  ΔGbwc(A + B -> AB) = ΔGbinding(AB) + ΔGstability(A Bound Conf) - ΔGstability(A Unbound Conf)
# F. Mutation - ΔG(WT -> M) = ΔGstability(M) - ΔGstability(WT) + G(WT Unfold -> M unfold)
#	However the last term is assumed to be zero. Also this does not work with bound ligand states.

# Each difference can be calculated for the WT and any number of mutants.
# This allows WT and Mutant cycles to be compared.
#! /bin/env python
import Environment, Utilities
from ProteinDesignTool import ProteinDesignTool

def InteractionFreeEnergy(structureA, structureB, complex, mutantListA, mutantListB);

	'''Calculates the free-energy difference between two structures and their complex.
	
	ΔGstability(AB) - ΔGstability(A) - ΔGstability(B) '''
	
	proteinDesignTool.dataDirectory = interactionData
	
	#ΔGstability(AB)
	proteinDesignTool.setPDB(complex)
	proteinDesignTool.runStabilityCalculation(mutantFiles=mutantFiles, resultsName='StabilityAB')
	
	proteinDesignTool.setPDB(structureA)
	proteinDesignTool.runStabilityCalculation(mutantFiles=mutantFiles, resultsName='StabilityA')
	
	proteinDesignTool.setPDB(structureB)
	proteinDesignTool.runStabilityCalculation(mutantFiles=mutantFiles, resultsName='StabilityB')
	
	#Totals for WT and each mutant
	for mutation in mutations:
		ab = interactionData.stabilityAB.dataForMutation[mutation] 
		a = interactionData.stabilityA.dataForMutation[mutation] 
		b = interactionData.stabilityB.dataForMutation[mutation]
		#Compute the resulting totals
		result = ab - a - b
		

#BindingWithChange -  ΔGbwc(A + B -> AB) = (ΔΔGstability(AU->AB) + ( ΔGstability(H) - ΔGstability(AU))
def BindingWithChange(holoStructure, apoStructure, ligand, mutantions)
	
	#
	#ΔGstability(Unbound Conf)
	#
	
	#Pass the HOLO pdb file to the ProteinDesignTool instance.
	#This cleans the pdb and copies it to the working directory
	proteinDesignTool.setPDB(pdb=holoStructure, dataDirectoryName='BindingWithChange')
	resultsData = proteinDesignTool.dataDirectory
	
	#Create the mutants - NOTE: Must use cleaned pdb file
	mutantCollection = Data.MutantCollection(pdbFile=proteinDesignTool.pdbFile,
					mutationList=mutations,
					location=proteinDesignTool.outputDirectory,
					maxOverlap=proteinDesignTool.configuration.getfloat('PKA SCAN PARAMETERS', 'mutation_quality'),
					clean=False)
		
	#Run the calculation							
	proteinDesignTool.runStabilityCalculation(mutantFiles=mutantCollection.mutantFiles(), 
				resultsName='StabilityHolo')
	
	#
	#ΔGstability(Bound Conf - No Ligand)
	#
	
	#Pass the APO pdb file to the tool. NOTE: The ligand is removed here during the cleaning process
	#Then create mutants with the ligand in place
	proteinDesignTool.setPDB(pdb=apoStructure)
	mutantCollection = Data.MutantCollection(pdbFile=proteinDesignTool.pdbFile,
					mutationList=mutations,
					ligandFiles=[ligand],
					location=proteinDesignTool.outputDirectory,
					maxOverlap=proteinDesignTool.configuration.getfloat('PKA SCAN PARAMETERS', 'mutation_quality'),
					clean=False)

	proteinDesignTool.runStabilityCalculation(mutantFiles=mutantCollection.mutantFiles(), resultsName='StabilityApoUnbound')
	
	#ΔGstability(Bound Conf - With Ligand)
	#Pass the APO pdb to the tool this time keeping the ligand.
	#NOTE: The mutant files must also contain the ligand coordinates
	proteinDesignTool.setPDB(apoStructure, removeLigand=False)
	proteinDesignTool.runStabilityCalculation(mutantFiles=mutantCollection.mutationFiles(), resultsName='StabilityApoBound')
	
	#Do the sum FreeEnergy = ΔGStability(AB) + ΔGstability(A Bound Conf) - ΔGstability(A Unbound Conf)
	
	results = []
	headers = ['Mutation', 'DeltaBind', 'DeltaDeltaApoBoundStability', 'DeltaDeltaApoUnboundStability', 'DeltaDeltaHoloStability'] 
	for i in range(len(resultsData.stabilityHolo.numberOfRows())):
		data = [resultsData.stabilityApoBound.total[i], resultsData.stabilityApoUnbound.total[i], resultsData.stabilityHolo.total[i]]
		freeEnergy = data[0] + data[1] - data[2]
		staticEnergy = 
		data.insert(0, freeEnergy)
		data.insert(0, resultsData.stabilityApoBound.mutation[i])
		results.append(data)
			
	matrix = Data.PEATSAMatrix(rows=results, headers=headers, name='DeltaBinding')
	resultsData.addMatrix(matrix, name='DeltaBinding')
	
	
def main():
	
	'''Main function for the TDCycle program.
	
	Parses the command line and starts the jobs requested'''
	
	environment = Environment.Environment()
	
	#Read arguments
	parser = Utilities.CommandLineParser()
	parser.parseCommandLine()
	
	if parser.helpRequested() is True:
		print __doc__
		sys.exit(0)

	if parser.createConfigurationRequested() is True:				
		print 'Creating default configuration file, proteinDesignTool.conf, in the current directory'
		print 'See the programs README file for information on the options'
		Environment.Configuration(searchDefaultLocations=False, writeFile=True)
		sys.exit(0)
	
	#Create the ProteinDesignTool instance
	tool = ProteinDesignTool(parser.configurationFile(), 
			workingDirectory=parser.workingDirectory(),
			pdbFile=parser.pdbFile(), 
			outputDirectory=parser.outputDirectory())
								








	
