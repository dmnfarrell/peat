'''Contains an example function which calculates stability differences between a protein and arbitrary mutants use PEAT_SA.Core classes'''
import os
import PEAT_SA.Core as Core

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
