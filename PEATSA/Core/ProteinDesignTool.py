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

'''Protein Design Tool - Creates, and runs a number of calculations on, a set of mutants of given protein

	Written by Michael Johnston
	This programs copyright - Michael Johnston and Jens Nielsen 2008-2010.
	For copyright of consituent programs check each

	NOTE: For detailed instructions see the programs README file.

	Usage: ProteinDesignTool -p [pdbname] -w [directory]  [optional arguments]
	       ProteinDesignTool --create-configuration

	-p - The path to the pdb the tool is to be run on - can be relative
	
	-w - The directory where the program will be run. It is created if it doesn't exist
	     
	--create-configuration 
		This outputs a default configuration file in the current directory.
		The configuration file should be modified to reflect your setup.
	
	The mutants to create can be specified using the --mutation, --mutationList 
	or --mutants arguments (see below).
	
	Calculation Arguments
	
	The following are the specific calculation arguments. 	
	
	--pKa	- Performs a pKa calculation on the protein. 
		      Note: Cannot be run in parallel.

	--delta-pKa
		Calculate the change in pKa of a set of ionisable residues due to the specified mutations. 
		The set to use can be specified using the --ionisablGroups option. 	
		If the --ionisableGroups option is not present all ionisable residue are used.
		
		Note: If the tool detects an intial pKa calculation hasn't been done 
		it does one (if running in serial).

	--stability - 
		Calculates the change in stability due to the specified mutations.
		Note: If a ligand is supplied the muations will be modelled with it present
	
	--ligand=[LigandFile]   
		Calculates the change in binding-free energy for the given ligand for
		the specified of mutants.
	
	Mutation Arguments
	
	If no mutation argument is specified --mutation=ALA is assumed.
	Note: Mutatants can be created independant of calculations if desired.
	
	--mutation=[MutationType]	            
		MutationType is a three letter amino acid code.
		The program creates a set of mutants where each residue in the protein
		is mutated to MutationType.
		Overridden by --mutationList if that argument is also specified.

	--mutationList=[MutationListFile]
		The program creates the set of mutants specified in the mutationList.	
		See the program guide for details on the file format.	
		
	--mutants=[MutantCollection]
		The program uses the set of mutants in the specified mutant collection.
		This will have been created by a previous run of the program
		      
	Optional Arguments
	
	-o	- Directory where the results of the various runs will be placed.
		  This defaults to the directory where the program was executed.
		  NOTE: This directory must already exist
		  
	-n	- The name used for the programs outputs.
		  This defaults to the name of the pdb file used.
		  
	-v	- Print more verbose output
	
	--configurationFile	
		 The location of the configuration file to use.   
		 
	--ionisableGroups=[ResidueCodes]
		ResidueCodes is a list of the residue codes of the ionisable groups to use in a
		delta-pKa calculation. 		 
	'''
 
import sys, os.path, datetime
import Environment, Exceptions, Programs, Utilities, Data
			
class ProteinDesignTool:
	
	def __init__(self, configurationFile, workingDirectory, pdbFile, outputDirectory, dataName=None):

		'''Creates a new PDT instance for running calculations on pdbName.
		
		Parameters
			configurationFile - File specifying information on the type of calculations to run.
			
			workingDirectory - Where to run the jobs. If you want to use the output of previous 
			pKa calculation as the input for a scan, the workingDirectory must contain it
			
			pdbFile - The file to run the jobs on.
			
			outputDirectory - Where to output information on the calculations and the results
			
		The results of any calculation are accesible through the returned objects
		dataDirectory attribute which is an instance of the Data.DataSet class.'''
		
		self.outputDirectory = outputDirectory
		self.pdbFile = pdbFile	
		self.environment = Environment.Environment()
		
		self.environment.output("[PEAT-SA] Beginning - %s" % datetime.datetime.now().strftime('%H:%M:%S'))
		#Readin the configuration file if specified
		#Otherwise search for one in the default locations (user home & current directory)
		self.configuration = Environment.Configuration(configurationFile, writeFile=False)
					
		#Create the working directory object
		self.workingDirectory = Environment.WorkingDirectory(location=workingDirectory, configuration=self.configuration)
		
		#Create a name for the data directory if none was provided.
		if dataName is None:
			dataName = os.path.split(self.pdbFile)[1]
			dataName = os.path.splitext(dataName)[0]
		
		#Set the pdb - This copies the pdb to the working directory and cleans it.
		#If also created the dataDirectory if requested.
		self.setPDB(pdbFile=pdbFile, dataDirectoryName=dataName)	
					
		#These will be created when they're needed. 
		#If some calculations are run on self.pdbFile before these are used
		#the file may change. Instantiating them with the unmodified file and using them
		#later may cause problems.
		self.scanner = None
		self.uffbaps = None

	def runPKACalculation(self):
	
		'''Runs a pKa calculation on the current pdb'''
	
		if not self.environment.isParallel:
			#Unlike the other tools PKA uses the WorkingDirectory class
			pKaRun = Programs.PKARun(workingDirectory=self.workingDirectory, 
						outputDirectory=self.outputDirectory,
						configuration=self.configuration,
						filename=self.pdbFile)
			pKaRun.run(cleanPDB=False)
		else:
			raise Exceptions.EnvironmentError, "Cannot run a pKa calculation in a parallel environment"

	def runScan(self, mutantCollection, ionisableGroups=None, resultsName='ScanResults', verbose=False):
	
		'''Calculates the effect of the mutations in mutantCollection on the pKa of ionisable groups.
		
		Parameters:
			mutantCollection - A MutantCollection instance representing the mutations to be scanned.
			ionisableGroups - A list of the residue codes of the groups that will be checked
			resultsName - The results matrix will be added to the dataDirectory with this name.
			verbose - If True the scan will output more verbose information on the
				progress of the calculation. Default - False.'''
			
		scanner = self.scanTool()
		if ionisableGroups is not None:
			scanner.setIonisableResidues(ionisableGroups)
	
		runpka = False
		try:
			scanner.scan(mutantCollection=mutantCollection, verbose=verbose)
			self.environment.output("[PEAT-SA] Finished scan - %s" % datetime.datetime.now().strftime('%H:%M:%S'))
		except Exceptions.MissingPKADataError, data:
			self.environment.output('\n[PEAT-SA] Detected error setting up scan\n%s.' % data)
			self.environment.output('\n[PEAT-SA] Initial pKa calculation required - running now.')
			runpka = True
			raise
		
		if runpka is True:
			try:
				#The initial design run wasn't done - Do it now
				pKaRun = Programs.PKARun(self.workingDirectory, 
						self.outputDirectory,
						configuration=self.configuration,
						filename=self.pdbFile)
				pKaRun.run(cleanPDB=False)
				#Retry
				scanner.scan(mutantCollection=mutantCollection, verbose=verbose)
			except BaseException, data:
				self.environment.output('\n[PEAT-SA] Caught exception while running scan (%d)' % data)
				raise
			
		self.dataDirectory.addMatrix(matrix=scanner.scanResults(), name=resultsName)

	def runStabilityCalculation(self, mutantFiles=[], resultsName='StabilityResults', verbose=False):
	
		'''Runs a stability calculation
			
		Parameters:
			mutantFiles=A list of path to pdb files containing the mutants.
			
		Returns: Nothing
			Access the results of the scan through the receivers dataDirectory attribute.
			i.e dataDirectory.scanResults'''
	
		if len(mutantFiles) is 0:
			raise Exceptions.MutantCollectionError, 'No mutant files available for %s' % os.path.split(self.pdbFile)[1]
		
		uffbaps = self.uffbapsTool()		
		uffbaps.deltaStability(comparePDBs=mutantFiles, verbose=verbose)
		self.environment.output("[PEAT-SA] Finished stability - %s" % datetime.datetime.now().strftime('%H:%M:%S'))

		
		self.dataDirectory.addMatrix(matrix=uffbaps.deltaStabilityResults(), name=resultsName)
		self.environment.output('\n[STABILITY] Stability run complete.')
		
	def runBindingCalculation(self, ligandFile, mutantFiles=[], 
			bindingResults='BindingResults', 
			deltaBindingResults='DeltaBindingResults', 
			verbose=False):
	
		'''Runs a binding calculation
		
		Parameters:
			ligandFile=Path to  a mol2 file containing the ligand
			mutantFiles=A list of pdb files containing the mutants.'''

		if len(mutantFiles) is 0:
			raise Exceptions.MutantCollectionError, 'No mutant files available for %s' % os.path.split(self.pdbFile)[1]
		
		uffbaps = self.uffbapsTool()
		uffbaps.binding(comparePDBs=mutantFiles, ligand=ligandFile, verbose=verbose)
		self.environment.output("[PEAT-SA] Finished binding - %s" % datetime.datetime.now().strftime('%H:%M:%S'))

		
		self.dataDirectory.addMatrix(matrix=uffbaps.bindingResults(), name=bindingResults)
		self.dataDirectory.addMatrix(matrix=uffbaps.deltaBindingResults(), name=deltaBindingResults)
			
		self.environment.output('\n[STABILITY] Binding run complete.')

	def scanTool(self):
	
		'''Returns the scanner instance used by the receiver'''
	
		#Create a scanner object
		#The scanner is needed regardles of wheather or not a scan is performed.
		#This is because we may need to access the data from a previous scan for
		#e.g. a stability calculation.
		
		if self.scanner is None:
			self.scanner = Programs.Scanner(pdbName=self.pdbFile,
					workingDirectory=self.workingDirectory,
					configuration=self.configuration, 
					outputDirectory=self.outputDirectory,
					ionisableResidues=None)

		return self.scanner
		
	def uffbapsTool(self):
	
		'''Returns the UFFBAPS instance used by the receiver'''
	
		if self.uffbaps is None:
			
			self.uffbaps = Programs.UFFBAPS(sourcePDB=self.pdbFile,
						executablePath=self.configuration.uffbapsLocation(),
						workingDirectory=self.workingDirectory)	

		return self.uffbaps	
			
	def cleanUp(self):
	
		self.workingDirectory.clean()
		
	def setPDB(self, pdbFile, dataDirectoryName=None, removeLigand=True):	
	
		'''Set the receiver to work on the given pdb file.
		
		The pdbFile is copied to the recievers workingDirectory and cleaned (remove waters, add hydrogens, correct structure).
		
		Parameters:
			dataDirectoryName: If None the current data directory is used.
			removeLigand: If True all ligands are removed from the pdb as well as standard cleaning. Defaults to True'''
		
		#Copy the pdbfile to the working directory - this is the version of the file
		#we will use with all subsequent processes.
		self.pdbFile = self.workingDirectory.copyFileToDirectory(pdbFile)
		self.environment.output('[PEAT-SA] Using pdb file at %s as working copy' % self.pdbFile)	
					
		#Create a data dir in the output directory if a directory name is specified.
		if dataDirectoryName is not None:	
			self.dataDirectory = Data.DataSet(name=dataDirectoryName+'.peatsa', location=self.outputDirectory)
						
		#Clean the pdb
		if removeLigand is True:
			self.environment.output('[PEAT-SA] Cleaning pdb - Removing water and ligands, correcting structure, adding hydrogens.')
			if self.environment.isRoot():
				Utilities.CleanPDB(inputFile=self.pdbFile, outputFile=self.pdbFile)
		else:
			self.environment.output('[PEAT-SA] Cleaning pdb - Removing water, correcting structure, adding hydrogens.')
			if self.environment.isRoot():
				Utilities.CleanPDB(inputFile=self.pdbFile, outputFile=self.pdbFile, removeLigand=False)	
	
		self.environment.output("[PEAT-SA] Finished clean - %s" % datetime.datetime.now().strftime('%H:%M:%S'))

	
		self.environment.wait()	
		self.environment.output('[PEAT-SA] Done')


def CreateMutants(pdbFile, configuration, parser, environment, cleanPDB=False):

	'''Creates a set of mutant pdbs based on the command line options'''

	environment.output('[PEAT-SA] Generating mutants now ...')
	ligandFile = None
	if parser.ligandFile() is not None:
		ligandFile = [parser.ligandFile()]
		
	#Check if a mutation list was provided
	#If not create a scan list
	if parser.mutationListFile() is not None:
		filename = parser.mutationListFile()
		object = Data.MutationListFile(filename)
		#Validate - Will raise an exception if it fails
		errorData = {}
		if not object.validateOnPDB(pdbFile, errorData):
			raise Exceptions.MutationListValidationError, errorData["errorString"] 
		else:	
			mutationList = object.mutantList()
	else:
		skipResidueTypes = []
		if parser.scanMutation() is 'ALA' or parser.scanMutation() is 'GLY':
			skipResidueTypes.append('ALA')
			skipResidueTypes.append('GLY')
		else:
			skipResidueTypes.append(parser.scanMutation())	

		mutationList = Data.CreateScanList(pdbFile=pdbFile, mutation=parser.scanMutation(), skipResidueTypes=skipResidueTypes)	
		
	#Create the collection		
	mutantCollection = Data.MutantCollection(pdbFile=pdbFile,
					name=parser.outputName(),
					mutationList=mutationList,
					ligandFiles=ligandFile,
					location=parser.outputDirectory(),
					maxOverlap=configuration.getfloat('PKA SCAN PARAMETERS', 'mutation_quality'),
					clean=cleanPDB)
				
				
	environment.output("[PEAT-SA] Finished mutant creation - %s" % datetime.datetime.now().strftime('%H:%M:%S'))
				
	environment.output('[PEAT-SA] Complete - %s' % mutantCollection)
	return mutantCollection
		
def main():
	
	'''Main function for the ProteinDesignTool program.
	
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
			outputDirectory=parser.outputDirectory(),
			dataName=parser.outputName())
	
	#Check if an initial pKa run requested
	if parser.pKa() is True:	
		tool.runPKACalculation()
		
	#Check if a previously created MutantCollection was specified
	mutantCollection = parser.mutants()
	if mutantCollection is None:	
		environment.output('[PEAT-SA] No mutants available from previous run.\n')
		#Note we use the ProteinDesignTool instances pdbFile attribute.
		#This points to the cleaned pdbFile the tool operates on.
		#Therefore we must use this to generate the mutants.
		mutantCollection = CreateMutants(tool.pdbFile, tool.configuration, parser, environment)	
		
	#Scan 	
	if parser.scan() is True:
		tool.runScan(mutantCollection=mutantCollection, 
				ionisableGroups=parser.ionisableGroups(),
				verbose=parser.verbose())

	#Stability	
	if parser.stability() is True:
		tool.runStabilityCalculation(mutantFiles=mutantCollection.mutantFiles(), 
				verbose=parser.verbose())

	#Binding
	if parser.binding() is True:
		tool.runBindingCalculation(mutantFiles=mutantCollection.mutantFiles(), 
			ligandFile=parser.ligandFile(), 
			verbose=parser.verbose())
			
	if parser.modes():
		doCalculation = True
		if parser.mutationListFile() is not None:
			object = Data.MutationListFile(parser.mutationListFile())
			doCalculation = object.isSinglePointMutationFormat()
				
		if doCalculation:
			mutations = mutantCollection.mutations()
			residueIndexes = []
			for mutation in mutations:
				index = Utilities.ParseResidueCode(mutation.residueCodes()[0])[1]
				residueIndexes.append(index)
			
			residueIndexes.sort()
			import Goodvibes.Main
			environment.output('[MODES] Calculating the effect of %d mutations on modes.\n', len(residueIndexes))
			modeCalculator = Goodvibes.Main.vibration()	
			results = modeCalculator.main('PEAT-SA',tool.pdbFile,residueIndexes=residueIndexes)
			print residueIndexes
		else:
			environment.output('[MODES] Can only perform mode calcualtion for SPMs.\n')

		
	if parser.transitionPath():
		print '[PATH] Transition path calculation not implemented yet - skipping'
	
	tool.cleanUp()
	
if __name__ == "__main__":

	environment = Environment.Environment()
	
	try: 
		main()
	except Exceptions.ProteinDesignToolException, data:
		print "\nError - %s" % data.errorString
		print "Description - %s\n" % data
		#If its an instance of Exception.ArgumentError print out the help docs.
		if data.__class__ == Exceptions.ArgumentError:
			print __doc__
		
		#In serial this does nothing
		#In parallel if all processes are here they exit
		#In addition if only one process raised the exception and is here, 
		#the others will be waiting elsewhere in the code
		#In this case, passing True here causes all the processes to exit.
		environment.wait(error=True)
		environment.exit(1)
	except SystemExit, data:
		#This is raised by sys.exit - check was the arg 0
		#If so then its a normal exit - don't print any error message
		if data.message != 0:
			print "\nEncountered an unexpected error - This is likely a bug"
			print "Please contact the developers with the following traceback information",
			print "along with details on the calculation being run (arguments, configuration file etc.)\n"
			print "Program exited with exit code %d" % data.message
			print "Traceback %s" % Utilities.GetTraceback()
			environment.wait(error=True)
			environment.exit(data.message)
		else:
			environment.wait(error=True)
			environment.exit(data.message)
			
	except BaseException, data:
		print "\nEncountered an unexpected error - This is likely a bug"
		print "Please contact the developers with the following traceback information",
		print "along with details on the calculation being run (arguments, configuration file etc.)\n"
		print "Error info %s" % data
		print "Traceback %s" % Utilities.GetTraceback()
		environment.wait(error=True)
		environment.exit(1)
			
	
