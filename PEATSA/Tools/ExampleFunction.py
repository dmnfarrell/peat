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
'''Contains an example function which calculates stability differences between a protein and arbitrary mutants use PEATSA.Core classes'''
import os
import PEATSA.Core as Core

def DeltaStability(inputFile, mutationList, configurationFile=None, workingDirectory=os.getcwd(), outputDirectory=os.getcwd()):

	'''Calculates the stability difference between a protein and set of mutants
	
	Parameters:
		inputFile: A PDB file of the protein
		mutationList: A list of Data.MutationSet instances. Each represents a mutant of the protein.
		configurationFile: The location of a proteinDesignTool.conf file - defaults to home directory.
		workingDirectory: Where the calculation will be run. 
		outputDirectory: Where the results will be written.
		
	Returns
		A Data.DataSet instance containing one matrix, stabilityResults.
		Each row of this matrix corresponds to a mutant defined in the mutationList argument.'''

	#Create the ProteinDesignTool instance
	tool = Core.ProteinDesignTool.ProteinDesignTool(configurationFile, 
			workingDirectory=workingDirectory,
			pdbFile=inputFile, 
			outputDirectory=outputDirectory)
	
	#The above cleans the pdb file and copies it to the working directory.
	#Use this pdb from now on.
	inputFile = tool.pdbFile
	
	#Create the mutants
	mutantCollection = Core.Data.MutantCollection(pdbFile=inputFile,
					mutationList=mutationList,
					location=outputDirectory,
					temporary=True)
					
	#Run stability calculation
	#The results are added to the ProteinDesignTool instance's dataDirectory attribute
	#This is an instance of Data.DataSet class
	tool.runStabilityCalculation(mutantFiles=mutantCollection.mutantFiles())

	#Clean up - Deletes files copied to the working directory for Uffbaps
	tool.cleanUp()
	
	return tool.dataDirectory

if __name__ == "__main__":

	#Create MutationSet instances for the mutations specfied
	#This class provides an easy way to represent mutants with multiple mutations
	#We'll create the mutation A:0001:LYS+A:0042:GLY
	set = Core.Data.MutationSet()
	set.addMutation(chain="A", residueIndex=1, mutation='LYS')
	set.addMutation(chain="A", residueIndex=42, mutation='GLY')
	
	print 'The mutation is %s' % set.codeString() 

	#If the configuration file is in the directory where you run the script 
	#you don't have to pass it - it will be read automaticalyl
	results = DeltaStability("2AQU.pdb", [set])
	
	#Just to print out the results.
	print results.stabilityResults.csvRepresentation()
	results.delete()
