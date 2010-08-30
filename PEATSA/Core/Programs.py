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


'''Contains classes representing the programs used by PEAT-SA to perform various calculations'''
import os, sys, shutil, subprocess, copy, datetime, operator
import Environment, Exceptions, Matrix, Utilities

stabilityCommentStrings = ["",	"Values are in KJ/mol",
				"For details on the various columns see the programs documentation",
				"This file is best viewed in a spreadsheet program e.g. Excel",
				""]

def ResidueCodeCompare(residue1, residue2):
	
	'''Compares two residue codes which have the form ChainID:ResidueNumber
	
	Note: Actually can compare any code where the components are separated by colons
	and the second element is the residue number
	Also note this does not take into accound possible different chain ids.
		
	Returns an integer - Meaning is the same as the inbuilt cmp function'''	

	residue1 = int(residue1.split(':')[1])
	residue2 = int(residue2.split(':')[1])
	
	if residue1 > residue2:
		return 1
	elif residue1 < residue2:
		return -1
	else:
		return 0

def scannerUsingStoredConfiguration(filename):

	'''Creates a new Scanner object using configuration information stored in a file
	
	Parameters: 
		filename: The name of a file containing a scan options
				
	Exceptions:
		ConfigurationError: If the file isn't a valid configurationFile'''

	configuration = Environment.Configuration(filename)
	return Scanner(configuration)	

class PKARun:

	'''Class representing the PKARun tool.
	
	This class runs an initial pKa calculation on a pdb'''
	
	def __init__(self, workingDirectory, outputDirectory, configuration=None, parallel=False, filename=None):
	
		'''Returns a new PKARun instance
		
		Note: Instances of this class used WHATIF. There must be a .WHATIF file in the
		users home directory or the environment variable WIF must be set with its location.
		
		Parameters:
			workingDirectory - A ProteinDesignTool.Core.WorkingDirectory instance representing
				the directory where the calculation will be run.
			outputDirectory - The directory where the class will store the log of the pKa calculation.	
			configuration - A ProteinDesignToo.Core.Configuration instance from which the instance
				will obtain parameters for the design run. If this is none no options are 
				passed to the run.
			parallel - Set to True if the job is to be run in parallel. Default is False
			filename - The path to the pdb file the job is to be run on.
			
		Exceptions:
		
			Raises an ArgumentError if the specified pdbfile doesn't exist. 	
		'''		
	
		self.outputDirectory = outputDirectory
		self.workingDirectory = workingDirectory
		self.configuration = configuration
		self.pdbFile = filename
		self.parallel=parallel
		self.environment = Environment.Environment()

	def run(self, cleanPDB=True):
	
		'''Runs the pKa calculation. This method does not return until the run is finished which may take some time.
		
		Parameters:
			cleanPDB - If True the all waters and ligands are removed from the pdbfile before the run.
				Default is True.'''
	
		#Import these here rather than at the top of the module
		#to avoid the program immediately crahsing if they are not present
		#FIXME: Possibly work around this in top-level ...
		import pKarun.pKarun_base
	
		#Reroute the WHATIF output to the log file
		#The file is called PDBNAME.log
		standardOut = sys.stdout
		scanLogFile = os.path.join(self.outputDirectory, '%s%d.log' % (os.path.split(self.pdbFile)[1], self.environment.rank()))
		#Stream is line buffered
		scanLogStream = open(scanLogFile, 'w', 1)
		sys.stdout = scanLogStream
	
		standardOut.write('[PKA CALC] Setting up pKa calculation for %s\n' % self.pdbFile)
		standardOut.write('[PKA CALC] Calculation log %s\n' % scanLogFile)
	
		#Setup the working directory for a pKa calculation
		self.workingDirectory.setupPKARun()
		#Copy the pdbfile to the directory if its not already there
		if not self.workingDirectory.containsFile(self.pdbFile) and self.environment.isRoot():
			self.pdbFile = self.workingDirectory.copyFileToDirectory(self.pdbFile)
		
		#Clean the pdb file if specified
		if cleanPDB is True:
			standardOut.write('[PKA CALC] Cleaning %s - Removing ligands and water, correcting pdb\n' % self.pdbFile)
			try:
				if self.environment.isRoot():
					Utilities.CleanPDB(self.pdbFile, self.pdbFile)
					
				self.environment.wait()	
				standardOut.write('[PKA CALC] Cleaning done.\n')
			except:
				#Reset streams
				scanLogStream.close()
				sys.stdout = standardOut
				raise
		
		#Get the options for the calculation from the configuration object (if it exists)
		if self.configuration is not None:
			options = []
			#See if any options were specified
			if self.configuration.has_section('PKA RUN'):
				standardOut.write('[PKA CALC] Reading options from configuration\n' % self.pdbFile)
				for option in self.configuration.options('PKA RUN'):
					options.append(option)
					options.append(self.configuration.get('PKA RUN', option))
				standardOut.write('[PKA CALC] Options are %s\n' % " ".join(options))	
			else:
				standardOut.write('[PKA CALC] No options provided. Defaults will be used\n')
			
			options.insert(0, 'pKarun.py')
			options.insert(1, self.pdbFile) 
		
		#Do the calculation
		try: 
			#I think this is just a group of miscellaneous functions
			#wrapped in a class
			functionObject = pKarun.pKarun_base.pKamisc()
			parameters = functionObject.parse_parameterlist(options)
			pKarunInstance = pKarun.pKarun_base.pKarun(os.getcwd(), self.pdbFile, parameters)
							
			standardOut.write('[PKA CALC] Running job in parallel\n')
			pKarunInstance.runseq()
		except:
			#Reset the output streams and reraise
			scanLogStream.close()
			sys.stdout = standardOut
			raise
			
		standardOut.write('[PKA CALC] Done\n')	
		
	def setPDBFile(self, filename):
	
		'''Sets the pdb file the job will be run on
		
		No checking is done on filename until the run() is called'''
		
		self.pdbFile = filename
	
	def setConfiguration(self, configuration):
	
		'''Sets the configuration object the instance obtains the run options from'''
		
		self.configuration = configuration	
	
	def setWorkingDirectory(self, workingDirectory):
	
		'''Sets the directory to run the job in'''
		
		self.workingDirectory = workingDirectory
	
	def pdbFile(self):
	
		'''Returns the path to the currently set pdb.
		
		Exceptions -
			Raises an ArgumentError if the specified pdbfile doesn't exist. '''
			
		return self.pdbFile	
	
	def configuration(self):
	
		'''Returns a copy of the configuration object the instance is using'''
	
		return copy.deepcopy(self.configuration)
		
	def workingDirectory(self):
	
		'''Returns the working directory object the instance is using'''
	
		return self.workingDirectory

class UFFBAPS():

	'''Class wrapping runs of the UFFBAPS/genialtNav tool'''
	
	#The contributions calculated by a delta stability run.
	#This is also corresponds to the order of the values in the list returned by _chopDeltaStabilityResultLine
	stabilityResultsHeaders = ['Van der Waals', 'Electrostatic', 'H-Bond', ' Desolvation', 'Backbone Entropy', 'Water-Bridges', 'Total']

	#The contributions calculated by a ligand binding run.
	#This is also corresponds to the order of the values in the list returned by __chopBindingResultLine__
	bindingResultsHeaders = ['Van der Waals',  'Electrostatic', 'H-Bond', 'Desolvation', 'Backbone Conformation',
		'Protein Entropy', 'Ligand Entropy', 'Water Bridges', 'Total']

	def __init__(self, sourcePDB, executablePath, workingDirectory):
	
		'''Creates a new UFFBAPS object
		
		Parameters:
			sourcePDB - The full path to a PDB file.
				The pdb in this file will be what all other pdbs will be compared against
				
			executablePath - The absolute path to the directory containing Uffbaps executable
			
			workingDirectory - A Environment.WorkingDirectory instance
			
		Exceptions:
					
			Raises a Exceptions.StabilityError if either the value supplied for the workingDirectory
			is not valid (i.e. not a directory) or the executablePath is invalid.
			The existance of the pdbs and necessary files is only checked when run() is called'''	
				
		self.sourcePDB = sourcePDB
		self.executablePath = executablePath
		self.resultsMatrix = None
		self.bindingMatrix = None
		self.deltaBindingMatrix = None
		self.workingDirectory = workingDirectory
		
		#Set up the working directory
		#This creates a subdir in the working dir containing the necessary files for running UFFBAPS
		self.workingDirectory.setupUFFBAPS()
		self.directoryPath = workingDirectory.uffbapsRunDirectory()

		#Check sourcePDB exists
		if not os.path.isfile(self.sourcePDB):
			raise Exceptions.StabilityError, "Supplied source pdb (%s) invalid - Cannot run stability calculation", self.sourcePDB
		
		if not os.path.isdir(self.executablePath):
			raise Exceptions.StabilityError, "Supplied location of Uffbaps invalid - Not a directory"

		self.environment = Environment.Environment()
	
	def absoluteBinding(self, ligand):
	
		#Copy the wildType and ligand to the working dir
		pdbName = os.path.split(self.sourcePDB)[1]
		ligandName = os.path.split(ligand)[1]
		if self.environment.isRoot():
			shutil.copyfile(self.sourcePDB, os.path.join(self.directoryPath, pdbName))
			shutil.copyfile(ligand, os.path.join(self.directoryPath, ligandName))
			
		self.environment.wait()	
	
		#Create the output stream
		outputFile = os.path.join(self.directoryPath, "output%d.txt" % self.environment.rank())
		outputStream = open(outputFile, "a+")
		
		#The executable
		executable = os.path.join(self.executablePath, "Uffbaps")
		
		self.environment.output('\n[BINDING] Running initial wild-type calculation for %s' % ligandName)
		self.environment.output('[BINDING] Using Uffbaps at %s' % self.executablePath)
		
		results = self._runBindingCalculation(pdbName, ligandName, executable, outputStream)

		self.environment.output('[BINDING] Results - %s' % zip(UFFBAPS.bindingResultsHeaders, results))
	
		outputStream.close()
		
		matrix = Matrix.PEATSAMatrix(rows=[results], headers=UFFBAPS.bindingResultsHeaders, name="Binding")
		
		return matrix
	
	def absoluteStability(self, pdbFiles=[], verbose=False):
	
		'''Calculates the absolute stability of a number of structures
		
		FIXME: Incomplete'''
				
		#If there are no compare pdbs supplied raise an exception
		if len(pdbFiles) == 0:
			raise Exceptions.StabilityError, "Unable to perform stability calculation - No files supplied"
			
		#The run file name
		runFile = os.path.join(self.directoryPath, "run%d.txt" % self.environment.rank())
		outputFile = os.path.join(self.directoryPath, "output%d.txt" % self.environment.rank())
		outputStream = open(outputFile, "a+")
		
		#The executable
		executable = os.path.join(self.executablePath, "Uffbaps")
		
		self.environment.output('\n[STABILITY] Running absolute stability calculation for %d structures' % len(pdbFiles))
		self.environment.output('[STABILITY] Using Uffbaps at %s' % self.executablePath)
		self.environment.output('[STABILITY] Beginning stability calculation ')
		
		#Force output of last print line - required when using ","
		sys.stdout.flush()
	
		#Variables for informing user of progress
		ten = int(len(pdbFiles)/10.0)
		#Avoid module zero if their are less than 10 pbds
		if ten == 0:
			ten = 1
		count = 0
	
		#Accumulate the results for each file into this list
		#It will be used at the end to create a matrix
		rows = []
		
		#Ask the environment to split the pdbFiles array if this is a parallel run
		pdbFiles = self.environment.splitArray(pdbFiles)
		for pdbFile in pdbFiles:
			if verbose:
				self.environment.output('[STABILITY] Structure %s' % os.path.split(pdbFile)[1], rootOnly=False)
			
			#Copy the mutant to the working dir
			pdbName = os.path.split(pdbFile)[1]
			copiedFile = self.workingDirectory.copyFileToUffbapsDirectory(pdbFile)
			
			#Create the run file 
			self.__writeAbsoluteStabilityRunFile__(pdbName, filename=runFile)
			
			#Run
			process = subprocess.Popen([executable, runFile], cwd=self.directoryPath, stdout=outputStream)
			process.wait()
		
			components = self.__processAbsoluteStabilityResults__()
						
			if verbose:
				self.environment.output(
					'[STABILITY] %s' % zip(UFFBAPS.stabilityResultsHeaders, components), 
					rootOnly=False)
			
			#Add the name of the pdbFile and then the stability values 
			#to the 'rows' variable
			row = [os.path.splitext(pdbName)[0]]
			row.extend(components)	
			rows.append(row)
			
			#Reset the stream
			outputStream.seek(0)
			outputStream.truncate()
			
			#Get rid of the mutant pdb file we copied
			os.remove(os.path.join(self.directoryPath, copiedFile))
			
			if not verbose:
				count = count + 1
				if count%ten == 0:
					#A ugly space is left if print is used
					sys.stdout.write(".")
					sys.stdout.flush()
								
		outputStream.close()
		#Combine all the rows
		rows = self.environment.combineArray(rows)
		
		if self.environment.isRoot():
			headers = ['Mutations']
			headers.extend(UFFBAPS.stabilityResultsHeaders)
			self.resultsMatrix = Matrix.PEATSAMatrix(rows=rows, headers=headers, name="Stability")
	
			#Add some comments
			comments = stabilityCommentStrings[:]
			date = datetime.datetime.today()
			comments.insert(1, "Generated  by PEAT_SA on %s" % date.strftime('%A %d/%m/%Y %H:%M:%S %Z'))	

											
	def _chopDeltaStabilityResultLine(self, resultLine):
	
		'''Chops the result line from a delta stability run and returns the essential information
		
		Note: Private method for internal use only
		
		Returns:
			A list containing the stability information
			'''
		
		#Chop up the result
		components = resultLine.split(" ")
		#Get rid of the pdb name which is the first string
		#Then make floats from all the values
		components.pop(0)
		components = [element.strip() for element in components]
		components = [element for element in components if element != "|" and element != ""]
		components = [float(element) for element in components]
	
		#Now just take the last 7 float values which are the 5 components of 
		#the overall change in stability, plus the target and total stability. 
		#Get rid of the target because its meaningless here.
		
		components = components[-7:]
		components.pop(5)
		
		#The water-bridges term is not output by the msc task - calculate it here
		subtotal = reduce(operator.add, components[:-1])
		bridges = components[-1] - subtotal
		components.insert(-1, bridges)
		
		return components
		
	def _writeDeltaStabilityRunFile(self, pdbOne, pdbTwo, filename=None):
	
		'''Writes a input file for Uffbaps for comparing the stability of two proteins
		
		 Note: Private method for internal use only'''
		 
		#Create the run file 
		run = "load %s %s\n" % (pdbOne, pdbOne)
		run = run + "load %s %s\n" % (pdbTwo, pdbTwo)
		run = run + "task msc {%s}{%s}\n" % (pdbOne, pdbTwo)
		
		runStream = open(filename, 'w')
		runStream.write(run)
		runStream.close()
		
	def deltaStability(self, comparePDBs=[], verbose=False):
	
		'''Calculate the difference in stability between the source pdb and the compare pbds
		
		Parameters
			comparePDBs - A list of paths to PDB files.
				All the pdbs in this list will be compared to the sourcePDB.
				The names of the compare pdbs are expected to conform to the following
				.....
				(FIXME: Possibly only have to be of a certain format to do a 'mutantStability' run?)
		
		Exceptions:
			Raises an Exceptions.UFFBAPSError if the source pdb does not exist.
			Raises an Exceptions.UFFBAPSError if any of the comparePDBS does not exist.
			However the calculation is performed for all existing pdbs first.'''
	
				
		#If there are no compare pdbs supplied raise an exception
		if len(comparePDBs) == 0:
			raise Exceptions.StabilityError, "Unable to perform stability calculation - No mutants supplied"
	
		#Copy the wildType to the working dir
		pdbName = os.path.split(self.sourcePDB)[1]
		#Uffbaps requires the pdb to have a .pdb extension - check this
		components = os.path.splitext(pdbName)
		if components[1] != 'pdb':
			pdbName = '%s.pdb' % components[0]
			
		if self.environment.isRoot():
			shutil.copyfile(self.sourcePDB, os.path.join(self.directoryPath, pdbName))
		self.environment.wait()	
		
		#The run file name
		runFile = os.path.join(self.directoryPath, "run%d.txt" % self.environment.rank())
		outputFile = os.path.join(self.directoryPath, "output%d.txt" % self.environment.rank())
		outputStream = open(outputFile, "a+")
		
		#The executable
		executable = os.path.join(self.executablePath, "Uffbaps")
		
		self.environment.output(('\n[STABILITY] Running stability comparison between %s and %d mutants' 
			% (os.path.split(self.sourcePDB)[1], len(comparePDBs))))
		self.environment.output('[STABILITY] Using Uffbaps at %s' % self.executablePath)
		self.environment.output('[STABILITY] Beginning stability calculation in %s' % self.directoryPath)
		
		#Force output of last print line - required when using ","
		sys.stdout.flush()
	
		#Variables for informing user of progress
		ten = int(len(comparePDBs)/10.0)
		#Avoid module zero if their are less than 10 pbds
		if ten == 0:
			ten = 1
		count = 0
	
		#Accumulate the results for each mutation into this list
		#It will be used at the end to create a matrix
		rows = []
		
		#Ask the environment to split the comparePDB array if necessary 
		#That is if this is a parallel run
		comparePDBs = self.environment.splitArray(comparePDBs)
		for pdbFile in comparePDBs:
			if verbose:
				self.environment.output('[STABILITY] Mutant %s' % os.path.split(pdbFile)[1], rootOnly=False)
			
			#Copy the mutant to the working dir
			compareName = os.path.split(pdbFile)[1]
			shutil.copyfile(pdbFile, os.path.join(self.directoryPath, compareName))
			
			#Create the run file 
			self._writeDeltaStabilityRunFile(pdbName, compareName, filename=runFile)
			
			#Run
			process = subprocess.Popen([executable, runFile], cwd=self.directoryPath, stdout=outputStream)
			process.wait()
		
			#Get the results
			outputStream.seek(0)
			lines = outputStream.readlines()
			result = lines[-1:]
			
			try:
				#If an exception is raised here it means something unexpected
				#was returned by UFFBAPS
				components = self._chopDeltaStabilityResultLine(result[0])
				
				#Add the name of the mutation and then the stability values 
				#to the 'rows' variable
				row = [os.path.splitext(compareName)[0]]
				row.extend(components)	
				rows.append(row)
				
				#Get rid of the mutant pdb file we copied
				os.remove(os.path.join(self.directoryPath, compareName))
				
				if not verbose:
					count = count + 1
					if count%ten == 0:
						#A ugly space is left if print is used
						sys.stdout.write(".")
						sys.stdout.flush()
				else:
					self.environment.output(
						'[STABILITY] %s' % zip(UFFBAPS.stabilityResultsHeaders, components), 
						rootOnly=False)		

			except BaseException:
				self.environment.output(
					'[STABILITY] Unable to parse unexpected UFFBAPS output\n%s\n' % lines, 
					rootOnly=False,
					stream=sys.stderr)
				self.environment.output(
					'[STABILITY] Skipping mutant %s\n' % os.path.split(pdbFile)[1], 
					rootOnly=False,
					stream=sys.stderr)
	
			#Reset the stream
			outputStream.seek(0)
			outputStream.truncate()
			
		outputStream.close()
		#This is required in order to catch if any the processes exited the 
		#above loop due to an exception	
		self.environment.wait()	
		#Combine all the rows
		rows = self.environment.combineArray(rows)
		
		if self.environment.isRoot():
			headers = ['Mutations']
			headers.extend(UFFBAPS.stabilityResultsHeaders)
			self.resultsMatrix = Matrix.PEATSAMatrix(rows=rows, headers=headers, name="Stability")
	
			#Add some comments
			comments = stabilityCommentStrings[:]
			date = datetime.datetime.today()
			comments.insert(1, "Generated  by PEAT_SA on %s" % date.strftime('%A %d/%m/%Y %H:%M:%S %Z'))
			
	def _chopLigandResultLine(self, resultLine):
	
		'''Chops the result line from a  binding run and returns the essential information
		
		Note: Private method for internal use only
		
		Returns:
			A list containing the binding information'''
	
		#Chop up the result
		components = resultLine.split(" ")
	
		#Get rid of the pdb name which is the first string
		#Then make floats from all the values
		components.pop(0)
		components = [element.strip() for element in components]
		components = [element for element in components if element != "|" and element != ""]
		components = [float(element) for element in components]
	
		#There are 8 components for the binding energy plus the prediction and total stability. 
		#Get rid of the prediction because its meaningless here.
		
		components.pop(8)
		
		return components		
		
	def _writeBindingRunFile(self, pdbOne, pdbTwo, filename=None):
	
		'''Writes a file that acts as input for Uffbaps fo comparing the stability of two proteins
		
		 Note: Private method for internal use only'''
		 
		#Create the run file 
		run = "load %s %s\n" % (pdbOne, pdbOne)
		run = run + "load %s %s\n" % (pdbTwo, pdbTwo)
		run = run + "task energy {%s}{%s}\n" % (pdbOne, pdbTwo)
		
		runStream = open(filename, 'w')
		runStream.write(run)
		runStream.close()
	
	def _runBindingCalculation(self, pdbName, ligandName, executable, outputStream):
	
		'''Runs a binding calculation for ligandName to pdbName 
		
		Note: This is a private method - runBindingCalculation() should be used by external callers
		
		The method writes a Uffbaps run file for calcualting the binding energy of ligandName to pdbName.
		It then runs Uffbaps and waits for it to finish.
		Finally it reads the last line output, extracts the information from it, and returns the results.
		
		Parameters:
			pdbName - Name of a pdb file in the directory where Uffbaps is being run
			
			ligandName - Name of a ligand file in the directory where Uffbaps is being run
			
			exectutable - Path to Uffbaps executable
			
			outputStream - Where the output should be written.
			Note this is empty on return.
			
		Return:
			Returns a list of the components of the binding energy and the total value.
			The order of the elements in the list correpsonds to the order of elements
			in bindingResultsHeaders'''
		
		runFile = os.path.join(self.directoryPath, "run%d.txt" % self.environment.rank())	
					
		#Create the run file 
		self._writeBindingRunFile(pdbName, ligandName, runFile)

		#Run
		process = subprocess.Popen([executable, runFile], cwd=self.directoryPath, stdout=outputStream)
		process.wait()
		
		#Get the results
		outputStream.seek(0)
		
		lines = outputStream.readlines()
		result = lines[-1:]
		try:
			components = self._chopLigandResultLine(result[0])
		except BaseException:
			self.environment.output(
				'[STABILITY] Unable to parse unexpected UFFBAPS output\n%s\n' % lines, 
				rootOnly=False,
				stream=sys.stderr)
			self.environment.output(
				'[STABILITY] Skipping mutant %s\n' % os.path.split(pdbName)[1], 
				rootOnly=False,
				stream=sys.stderr)
			components = None	
			
		#Reset the stream
		outputStream.seek(0)
		outputStream.truncate()									
		
		return components
	
	def binding(self, comparePDBs=[], ligand=None, verbose=False):
	
		'''Calculates the binding free-energy of a ligand to the source pdb and the compare pdbs.
		
		Also calculates the change in binding free-energy for each to the source pdb
		
		Parameters
			comparePDBs - A list of paths to PDB files.
				All the pdbs in this list will be compared to the sourcePDB.
				The names of the compare pdbs are expected to conform to the following
				.....
				(FIXME: Possibly only have to be of a certain format to do a 'mutantStability' run?)
				
			ligand - The path to a pdb file containing the ligand to check'''
	
		#If there are no mutants supplied raise an exception
		if len(comparePDBs) == 0:
			raise Exceptions.StabilityError, "Unable to perform stability calculation - No mutants supplied"

		#Copy the wildType and ligand to the working dir
		pdbName = os.path.split(self.sourcePDB)[1]
		ligandName = os.path.split(ligand)[1]
		if self.environment.isRoot():
			shutil.copyfile(self.sourcePDB, os.path.join(self.directoryPath, pdbName))
			shutil.copyfile(ligand, os.path.join(self.directoryPath, ligandName))
			
		self.environment.wait()	
	
		#Create the output stream
		outputFile = os.path.join(self.directoryPath, "output%d.txt" % self.environment.rank())
		outputStream = open(outputFile, "a+")
		
		#The executable
		executable = os.path.join(self.executablePath, "Uffbaps")
		
		self.environment.output('\n[BINDING] Running initial wild-type calculation for %s' % ligandName)
		self.environment.output('[BINDING] Using Uffbaps at %s' % self.executablePath)
		
		wildTypeResults = self._runBindingCalculation(pdbName, ligandName, executable, outputStream)

		self.environment.output('[BINDING] Wild-Type results - %s' % zip(UFFBAPS.bindingResultsHeaders, wildTypeResults))
		self.environment.output('\n[BINDING] Beginning mutant binding calculations \n')
		
		#Force output of last print line - required when using ","
		sys.stdout.flush()
	
		#Variables for informing user of progress
		ten = int(len(comparePDBs)/10.0)
		#Avoid module zero if their are less than 10 pbds
		if ten == 0:
			ten = 1
		count = 0
	
		#Accumulate the results for each mutation into this list
		#It will be used at the end to create a matrix
		rows = []
		comparePDBs = self.environment.splitArray(comparePDBs)
		for pdbFile in comparePDBs:
			if verbose:
				self.environment.output('[BINDING] Mutant %s' % os.path.split(pdbFile)[1], rootOnly=False)
			
			#Copy the mutant to the working dir
			compareName = os.path.split(pdbFile)[1]
			shutil.copyfile(pdbFile, os.path.join(self.directoryPath, compareName))
			
			#Do the calculation
			components = self._runBindingCalculation(compareName, ligandName, executable, outputStream)
			
			if components is not None:
				if verbose:
					self.environment.output(
						'[BINDING] Results - %s' % zip(UFFBAPS.bindingResultsHeaders, components),
						rootOnly=False)
				else:
					count = count + 1
					if count%ten == 0:
						#A ugly space is left if print is used
						sys.stdout.write(".")
						sys.stdout.flush()		
			
				#Add the name of the mutation and then the ligand values 
				#to the 'rows' variable
				row = [os.path.splitext(compareName)[0]]
				row.extend(components)	
				rows.append(row)
			
				#Get rid of the mutant pdb file we copied
				os.remove(os.path.join(self.directoryPath, compareName))
								
		outputStream.close()
		
		#This is required in order to catch if any the processes exited the 
		#above loop due to an exception	
		self.environment.wait()	
		#Combine all the rows
		rows = self.environment.combineArray(rows)
		
		if self.environment.isRoot():
			
			headers = ['Mutations']
			headers.extend(UFFBAPS.bindingResultsHeaders)
			self.bindingMatrix = Matrix.PEATSAMatrix(rows=rows, headers=headers, name="Binding")
			
			#Create a delta stability matrix Gmut - Gwt
			deltaRows = []
			for row in rows:
				#Get just the numbers
				values = row[1:]
				#Add the residue information to a new row entry
				newRow = [row[0]]
				#Extend it with  (values - wildTypeResults)
				try:
					newRow.extend(map(operator.sub, values, wildTypeResults))
					deltaRows.append(newRow)
				except:
					print 'Unable to subtract %s from %s' % (values, wildTypeResults)

			self.deltaBindingMatrix = Matrix.PEATSAMatrix(rows=deltaRows, headers=headers, name="DeltaBinding")
		
	def deltaStabilityResults(self):
	
		'''Returns the results of the last stability comparison run.
		
		Return:
			Returns None if no successful stability comparison run has been performed.'''
		
		#Return the results of the last run
		return copy.deepcopy(self.resultsMatrix)
		
	def bindingResults(self):
	
		'''Returns the results of the last binding  run.
		
		Return:
			Returns None if no successful binding comparison run has been performed.'''
	
		return copy.deepcopy(self.bindingMatrix)
		
	def deltaBindingResults(self):
	
		'''Returns the delta binding results of the last binding  run.
		
		Return:
			Returns None if no successful binding comparison run has been performed.'''
	
		return copy.deepcopy(self.deltaBindingMatrix)

class Scanner:

	'''Instances calculate the change in pKa of ionisable residues due to mutating other residues.
	
	Terminology - 
		scan - The calculation of the change in pKa of all ionisable residues due to a single mutation.
		i.e. Means 'scan all ionisable residues after making a given change'
	
	Each objects main attribute is a set of parameters which define how the scan is to be performed
	
	Note a scan requires the output of a previous pKA run to be present in the workingDirectory.'''

	isInitialised = False

	def __init__(self, workingDirectory, configuration=None, pdbName=None, outputDirectory=None, ionisableResidues=None,):
	
		'''Creates a new Scanner instance.
		
		Parameters
			pdbName: The name of the pdb file containing the protein the scan is to be performed on.
			The pdb is assumed to be already in the working directory along with the other required
			output of a pKa calculation.
			
			outputDirectory: The object will write the log of the scan to this directory.
				Defaults to the working directory.
				
			workingDirectory: A Environment.WorkingDirectory instance representing the directory where the
				scan will be run.
			
			ionisableResidues: A list containing the names of residues to be scanned.
				The names must be in the format A:xxxx:YYY where A is the chain number
				xxxx is the residue number e.g. 0040, and YYY is its three letter code.
				The chain number may be omitted
				If this is None or is not specified all ionisable residues are selected.
			
			configurationObject: An Environment.Configuration instance which represents
				the parameters to be used. If this is None the default configuration is used.
		
		Note
			The scan is always run in the working directory.
			This will be expected to contain the necessary auxillary input files for the pdb.
			If the current directory is not the working directory it is changed to it.
								
		Exceptions
			TypeError: Raised if the format of the entries in residues is not valid'''
			
		#Initialise these ivars to None in case something goes wrong in setPDB
		self.pdbName = None
		self.pdb = None						
		self.setPDB(pdbName)
		self.designer = None
		self.mutantCollection=None
		self.environment = Environment.Environment()
		
		if workingDirectory is None:
			raise ValueError, 'WorkingDirecotry argument cannot be None'
		else:
			self.workingDirectory = workingDirectory
							
							
		#if no ionisable residue are explicitly passed use them all.
		if ionisableResidues == None:
			#Contrary to its name get_titratable_groups does not return anything.
			#Have to access the data directly through the ivar titratable_groups
			self.pdb.get_titratable_groups()
			#Keys are the chain + residue numbers e.g. T:0003. 
			#The corresponding value is a list containing a dictionary (strange) with three keys (currently)
			#These are charge, name and atoms. name gives the second part of the ionisable residue id.
			ionisableResidues = []
			string = ""
			keys = self.pdb.titratable_groups.keys()
			keys.sort(ResidueCodeCompare)
			for key in keys:
				name = self.pdb.titratable_groups[key][0]['name']
				#Skip NTERM and CTERM
				#FIXME - the extra colon in NTERM is due to a bug in Protool.
				if name != 'NTERM' and name != 'CTERM':
					residueCode = "%s:%s" % (key, self.pdb.titratable_groups[key][0]['name'])
					ionisableResidues.append(residueCode)
					string = string + "%s " % residueCode

		self.setIonisableResidues(ionisableResidues)	
		
		#Set the output directory - default to current if none is supplied
		if outputDirectory == None:
			outputDirectory = os.getcwd()
	
		self.outputDirectory = os.path.abspath(outputDirectory)
		
		#If no configuration object was supplied use default
		if configuration == None:
			self.configuration = Environment.Configuration()
		else:
			self.configuration = configuration
		
	def _singleScan(self, mutationSet, defaults, mutantCollection):
	
		'''Does a scan of the pdb defined by mutationSet
		
		This is a private method. Do not call externally
		
		Parameters
			mutationSet is a Data.MutationSet instance
			defaults are the parameters to pass to calc_dpKa
			mutantCollection is the Data.MutantCollection instance the mutationSet
			is part of.'''
		
		#Import this module here so it won't cause the program to immediatey fail
		#on start if it missing. If it is missing then Environment.Configuration instance
		#will have already noticed and dealt with the situtation
		import pKa
		
		# Set what the residue is being mutated to - usually ALA
		defaults['mutations'][0] = mutationSet.codeString(self.pdb, reduced=False)
		self.environment.output('Mutation %s' % defaults['mutations'][0], rootOnly=False)
		
		#This was a conversion done internally in run_opt.
		#Now that its here the conversion to the defaults ivar format looks
		#redundant. Will be changed at some later stage
		params={}
		for key in defaults.keys():
			params[key]=defaults[key][0]

		print defaults
		#Create a new designer if none was created so far
		if self.designer is None:
			self.designer = pKa.Design_pKa.Design_pKa(params)
		else:
			# If we have an instance then just copy parameters
			self.designer.params=params.copy()
			self.designer.parse_parameters()
 
		filename = mutantCollection.fileForMutant(mutationSet)
		self.designer.name = pKa.Design_pKa.print_setup(params)
		self.designer.mutfile_names={mutationSet.codeString(self.pdb, reduced=False):filename}
		self.environment.output('Mutation File %s' % self.designer.mutfile_names, rootOnly=False)
		solutions = self.designer.calc_dpKa()
		
		#Convert the information from the (annoying ;-) format its in
		#Note: Key is a string rep of a list with the code in it e.g. "['A:0005:GLU:ALA']"
		key = repr(mutationSet.mutationCodes(self.pdb, reduced=False))
		aList = [mutationSet.codeString(self.pdb, reduced=True)]
		for residue in self.ionisableResidues:
			print 'The results for %s are %s' % (residue, solutions)
			try:
				aList.append(solutions[key][residue])
			except KeyError:
				#This reisdue is the mutation residue - just add 0 for this
				aList.append(0)
			except TypeError:
				#The result (solutions[key]) is None
				aList.append('Error')	

		print 'The converted list %s' % aList
				
		return aList

	def _processTasks(self, taskArray, results, standardOut, scanLogStream, verbose, checkInterval, skewThreshold):
	
		'''Process the tasks in taskArray and accumulates the answers in results.
		
		Every checkInterval tasks this method checks if the loadSkew
		(the balance of tasks across nodes) is greater than skewThreshold.
		If it is, this method exits, returning an array containing the uncompleted tasks'''
	
		#Get the scan options from the configuration object
		#in the format used by Design_pKa.run_opt()
		#Add whats missing
		defaults = self.configuration.designPKAOptions()
		defaults['pdb'][0] = self.pdbName
		defaults['pKas'][0] = self.dummyPKAArg
	
		#Count how many mutations there are
		numberMutations = len(taskArray)
	
		#Variables for informing user of progress
		ten = int(numberMutations/10)
		if ten == 0:
			ten = 1
		completedTasks = 0

		start = datetime.datetime.now()
		checkTime = 600
	
		for element in taskArray:
		
			standardOut.flush()
			#On the first run there is a long wait while various files are read. 
			#The first mutation is then performed. 
			if completedTasks == 0:
				self.environment.output(
					'[SCAN] (%d on %s) Initialising ...\n' % (self.environment.rank(), self.environment.processorName()),
					stream=standardOut, 
					rootOnly=False)
		
			try:
				if verbose:	
					data = (self.environment.rank(), self.environment.processorName(), element.codeString(self.pdb))
					self.environment.output(
						'[SCAN] (%d on %s) Calculating for Mutation %s\n' % data,
						stream=standardOut,
						rootOnly=False)

				aList = self._singleScan(element, defaults, self.mutantCollection)
				if verbose: 
					remain = len(taskArray) - completedTasks
					self.environment.output(
						'[SCAN] (%d on %s) Done %d remaining\n' %  (self.environment.rank(), self.environment.processorName(), remain), 
						stream=standardOut, rootOnly=False)
			except BaseException, reason:
				self.environment.output('[SCAN] Caught exception - Data %s. Reraising' % reason, stream=standardOut)
				scanLogStream.close()
				sys.stdout = standardOut
				raise
				
			#Cheating slightly
			#This line simply informs the user the initialisation is over is starting.
			if completedTasks == 0:
				self.environment.output('[SCAN] Scanning ')
				
			results.append(aList)
			
			#Print some dots to show progress
			completedTasks = completedTasks + 1
			if not verbose:
				if completedTasks%ten == 0:
					standardOut.write(".")

			#Check if we require load balancing - we really need this to be asynchronous
			#As it is every node has to reach this point to calculate the skew
			#We wait to check if we have been running longer than the checkInterval
			#OR if the time to wait now is less than the time we estimate we will have
			#exceeded the checkInterval if we go ahead.
			runningTime = datetime.datetime.now() - start
			seconds = 86400*runningTime.days + runningTime.seconds
			timePerStep = seconds/float(completedTasks)
			timeToCheck = checkTime - seconds
			estimatedNextTime = seconds + timePerStep

			self.environment.output('[SCAN] (%d) Running for %s - %lf to check interval\n' % (self.environment.rank(), runningTime, timeToCheck), 
				stream=standardOut,
				rootOnly=False)

			self.environment.output('[SCAN] (%d) Time per step %lf - Estimated next check %lf.\n' % (self.environment.rank(), timePerStep, estimatedNextTime), 
				stream=standardOut,
				rootOnly=False)

			self.environment.output('[SCAN] (%d) Estimated next time to check interval %lf.\n' % (self.environment.rank(), checkTime - estimatedNextTime), 
				stream=standardOut,
				rootOnly=False)
				
			if seconds >= checkTime or (abs(timeToCheck) < abs(checkTime - estimatedNextTime)):
				startWait = datetime.datime.now()
				self.environment.output('[SCAN] Checking skew\n', stream=standardOut, rootOnly=False)
				skew = self.environment.loadSkew(len(taskArray) - completedTasks)
				self.environment.output('[SCAN] Skew is %lf\n' % skew, stream=standardOut, rootOnly=False)
				endWait = datetime.datime.now()
				self.environment.output('[SCAN] Waited for %s\n' % (endWait - startWait), stream=standardOut, rootOnly=False)
				checkTime += 600
				if  skew > skewThreshold:
					break

		return taskArray[completedTasks:]

	def scan(self, mutantCollection, verbose=False, printList=True):
	
		'''Performs a scan for each mutation specified in data
		
		Parameters: 
			mutantCollection - MutantCollection containing the mutants to be scanned
			verbose - True if more verbose output is required. False by default.
			printList - Outputs a list of all the mutations being scanned if True.	
				
		Returns:
			Nothing.
			
			To access the scan results use the scanResults method
			
		Exceptions:
			AttributeError is raised if no pdb or ionisable residues have been provided.
			Exceptions.MissingPKADataError is raised if the working directory is missing required data.
			Exceptions.ScanError is raised if the scan can't be performed for any other reason.'''
	
		if len(mutantCollection.mutations()) == 0:
			raise Exceptions.ArgumentError, "No mutations supplied"
		
		self.environment.output('[SCAN] Ionisable residues - %s' % " ".join(self.ionisableResidues))
		self.mutantCollection = mutantCollection
		data = self.mutantCollection.mutations()
		if verbose is True:
			self.environment.output('[SCAN] Scanning the following mutations - ')	
			for element in data:
				self.environment.output('%s' % (element.__str__()))
				
		#Setup the workding directory
		#This will check to see if the necessary files are present
		#and subsequently set up the directory so the scan process uses the files in the mutantCollection.
		#This method will raise an Exceptions.MissingPKADataError if the first requirement is not satisfied.
		#First make sure we are in the working directory
		os.chdir(self.workingDirectory.path())
		self.workingDirectory.setupScan(pdbName=self.pdbName, mutantCollection=self.mutantCollection)
		
		#Redirect to scan output to file
		self.scanLogFile = os.path.join(self.outputDirectory, '%s_Scan%d.log' % (os.path.splitext(self.pdbName)[0], self.environment.rank()))
		self.environment.output('\n[SCAN] Scan output redirected to file Scan.log ...')
		self.environment.output('[SCAN] Reading necessary scan information')
						
		standardOut = sys.stdout
		scanLogStream = open(self.scanLogFile, 'w')
		sys.stdout = scanLogStream
		
		#Each element of data is a (ResidueCode, MutationList) pair.
		#Make each mutation in mutationList and scan the effect of the mutation on each ionisable residue
		results = []

		#Split the tasks across the nodes
		taskArray = self.environment.splitArray(data)
		codes = [el.codeString(self.pdb) for el in data]
		info = (self.environment.rank(), self.environment.processorName(), " ".join(codes))
		self.environment.log(
			'\n[SCAN] (%d on %s) Mutations %s\n' % info)

		while len(taskArray) != 0:
			#Process the tasks - If the loadSkew exceeds the threshold this 
			#method will return an array of the uncompleted tasks
			taskArray = self._processTasks(taskArray, results, standardOut, scanLogStream, verbose, checkInterval=5, skewThreshold=1.0)
			self.environment.output('[SCAN] (%d) Rebalancing load ...' % self.environment.rank())
			#Re-distributed the remaining tasks across the nodes - this should be asynchronous somehow
			taskArray = self.environment.balanceArrays(taskArray)
			self.environment.output('[SCAN] (%d) Continuing with %d tasks'% (self.environment.rank(), len(taskArray)))

		#This is required in order to catch if any the processes exited the 
		#above loop due to an exception	
		self.environment.wait()	
		results = self.environment.combineArray(results)
	
		#Store the results in a matrix
		headers = ['Mutations']
		headers.extend(self.ionisableResidues)
		self.resultsMatrix = Matrix.Matrix(rows=results, 
					headers=headers)
		
		#Redirect log back to stdout
		scanLogStream.close()
		sys.stdout = standardOut
		sys.stdout.flush()
		
		self.workingDirectory.cleanUpScan()
		self.environment.output('\n[SCAN] Scan completed.')	
	
	def scanLogFile(self):
	
		'''Returns the path to the log file for the last scan'''
		
		return self.scanLogFile

	def setIonisableResidues(self, residues):
	
		'''Sets the ionisable residues used.
		
		Parameters
			residues: A list of residue names. Can be in any format.'''
		
		self.pdb.get_titratable_groups()
		#This is a list of Chain:Number pairs
		titratableGroups = self.pdb.titratable_groups.keys()	
		filteredList = []	
		errorList = []
		for code in residues:
			try:
				#Parse the code - it may be reduced, extended etc.
				components = list(Utilities.ParseResidueCode(code))
			except Exceptions.ResidueCodeFormatError, data:
				error = {}
				error['code'] = code
				error['reason'] = data	
				errorList.append(error)	
				continue
			
			#Check that the code refers to a titratable group
			#Note: The code in the titratableGroups list don't contain the residue name		
			keyCode = Utilities.CreateResidueCode(chain=components[0], number=components[1])
			if keyCode in titratableGroups:
				#If the name was supplied werify it is correct
				#If not its an error
				name = self.pdb.titratable_groups[keyCode][0]['name']
				if (len(components) == 3) and (components[2] != name):
					error = {}
					error['code'] = code
					error['reason'] = "Incorrect name (%s) specified for residue %s. Correct to %s" % (components[2], keyCode, name)	
					errorList.append(error)		
				else:
					code = Utilities.CreateResidueCode(chain=components[0], 
							number=components[1], 
							name=name)
					filteredList.append(code)
			else:
				error = {}
				error['code'] = code
				error['reason'] = "Supplied residue code does not refer to an ionsiable group"	
				errorList.append(error)	
			
		if len(errorList) != 0:
			data = ""
			for el in errorList:
				data = data + "\tCode %s. Problem: %s\n" % (el['code'], el['reason'])
				
			raise Exceptions.ArgumentError, "Some of the ionisable residue codes you supplied are invalid\n%s" % data
					
		self.ionisableResidues = filteredList
			
		#Have to make a dummy value for the 'pKas' argument when doing the scan
		#This is a list of the ionisable residues plus how much their pKa
		#should be changed (although no change is done)
	
		self.dummyPKAArg =''
		for target in self.ionisableResidues:
			self.dummyPKAArg = self.dummyPKAArg + target + '=0.0,' # Dummy pKa value
			
		#Remove the last comma
		self.dummyPKAArg = self.dummyPKAArg[:-1]
				
	def setPDB(self, pdbName):
	
		'''Sets the pdb to be scanned.
		
		Parameters
			pdbName - The name of the pdb file containing the protein structure the scan 
				will be performed on
			
		Exceptions:
			Raises an IOError when the pdb can't be read
			Raises an Exceptions.FileFormatError when the pdb format is invalid.
		
		NOTE: The ionisable residues are not changed and must be done so separately if required'''	

		#In case the path was passed
		pdbName = os.path.split(pdbName)[1]
		
		#Create pdb structure object
		import Protool
		pdb = Protool.structureIO()
		pdb.readpdb(pdbName)
		try:
			pdb = Protool.structureIO()
			pdb.readpdb(pdbName)
		except Exception, data:	
			raise Exceptions.FileFormatError, 'Format of specified PDB file %s not valid.\nUnderlying error - %s' % (pdbName, data)
			
		self.pdb = pdb
		self.pdbName = pdbName
		
	def PDB(self):
	
		return self.pdbName
								
	def scanResults(self):
	
		'''The results of the last scan performed.
		
		Returns
			A Matrix instance.
			The first four columns are:
			'Chain, Residue Code, Residue Number, Mutation'
			
			There is one remaining colomn for each ionisable residue the scan was defined for.
			Each entry in these columns gives the change in pKa due to the mutation defined
			by the entries in the first four columns in the given row.
			
			This method return None if no scan has been performed'''
		
		return self.resultsMatrix
				

class PBSolver:

	'''Calculates electrostatic solvation energies for protein ligand complexes using PB equation'''
	
	
	def __init__(self):
	
		self.environment = Environment.Environment()
	
	def _getNewSID(self, tag):

		'''Returns a unique string for use as a Session ID'''
	
		import time, md5, random
		
		t1 = time.time()
		t2 = t1 + random.random()
		base = md5.new( tag + str(t1 +t2) )
		sid = tag + '_' + base.hexdigest()
		return sid
		
	def _convertProtonationData(self, protonationStates):
	
		'''Replaces PEAT type residue codes with pdb2pqr residue codes'''
		
		#Convert protonation states
		converted = {}
		for residueCode in protonationStates.keys():
			data = Utilities.ParseResidueCode()
			pdb2pqrCode = '%s %s %d' % (data[2], data[0], data[1])
			converted[pdb2pqrCode] = protonationStates[residueCode]
			
		return converted																					
				
	def solvationEnergy(self, pdbFile, ligandFile, protonationStates={}, verbose=True):
	
		'''Returns the electrostatic solvation energies for the protein ligand complex
		
		Params:
			ligandFile: A ligand in mol2 format
			pdbFile: A pdb the ligand binds to
		
		Returns:
			A list of energies. The elements of the list are
			(Desolvation, Interaction, Protein Solvation, Water Solvation)
			The Protein and Water solvation energies have been scaled by 0.5 ala Born'''
		
		try:	
			import pdb2pqr.pdb2pka.pka
		except ImportError, data:
			raise Exceptions.EnvironmentError, 'Unable to compute solvation energies - %s' % data
		
		dir = self._getNewSID('PB')
		dir = os.path.abspath(dir)
		try:
			os.mkdir(dir)
			os.chdir(dir)
		except OSError, data:
			print data
			dir = self._getNewSID('PB')
			dir = os.path.abspath(dir)
			os.mkdir(dir)
			os.chdir(dir)
		
		#Replace PEAT mutation codes with pdb2pqr mutation codes
		if protonationStates is not None and len(protonationStates.keys()) != 0:
			protonationStates = self._convertProtonationData(protonationStates)

		#Copy the ligand and pdb files to the working dir since
		#the pdb2pka function will modify them
		try:
			if verbose is False:
				stdout = sys.stdout
				newout = open('PB.out', 'a+')
				sys.stdout = newout
				stderr = sys.stderr
				newerr = open('PB.err', 'a+')
				sys.stderr = stderr
		
			pdbCopy = os.path.join(dir, os.path.basename(pdbFile))
			shutil.copyfile(pdbFile, pdbCopy)
			ligandCopy = os.path.join(dir, os.path.basename(ligandFile))
			shutil.copyfile(ligandFile, ligandCopy)
						
			result = pdb2pqr.pdb2pka.pka.get_res_energies(pdbCopy, 
					ligandCopy ,
					'L:0001:LIG',
					 fix_states=protonationStates)
		except Exception, data:
			print 'Encountered an exception with file %s' % pdbFile 
			print data	
			raise		 
		finally:
			if verbose is False:
				newout.close()
				sys.stdout = stdout
				newerr.close()
				sys.stderr = stderr

			os.chdir('..')
			shutil.rmtree(dir)
			#if os.path.isdir('phidir'):
			#	shutil.rmtree('phidir')

		#Convert the result to kjoule per mol  (its in kT) - assume 300K - 2.494
		result = map(lambda x: x*2.494, result)
	
		result = list(result)
		result.reverse()	
		return result
		
	def mutantSolvationEnergies(self, mutantCollection, ligandFile, protonationStates={}):
	
		'''Calculates electrostatic solvation energies for all the mutant/ligand complexes
		
		Params:
			ligandFile: A ligand in mol2 format
			mutantCollection: A Core.Data.MutantCollection instance
		
		Returns:
			A matrix of results. The matrix columns are
			(MutationCode, Desolvation, Interaction, Protein Solvation, Water Solvation)'''
		
		rows = []
		mutations = mutantCollection.mutations()
		mutations = self.environment.splitArray(mutations) 
		for set in mutations:
			mutant = mutantCollection.fileForMutant(set)
			result = self.solvationEnergy(mutant, ligandFile, protonationStates)
			result.insert(0, "+".join(set.reducedMutationCodes(pdb=mutantCollection.pdb)))
			rows.append(result)
	
		rows = self.environment.combineArray(rows)
	
		return Matrix.PEATSAMatrix(rows=rows, headers=["Mutations", "Solvation W", "Solvation P+W", "P-L Interaction", "Delta Solvation"])

