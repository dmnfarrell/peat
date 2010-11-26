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

'''Contains functions and classes related to handling PEAT-SA's environment'''

import Exceptions, ConfigParser, os, sys, shutil, math, time, tempfile, datetime, operator

#Assign the location of the module to a variable
moduleLocation = os.path.dirname(__file__)
moduleLocation = os.path.abspath(moduleLocation)

#Add it to the path so we can always access the modules
sys.path.insert(0, moduleLocation)

#Create the resource path
resources = os.path.join(moduleLocation, 'Resources')

#Required files for any alanine scan 
requiredPKADesignFiles = ['DELCRG.DAT', 'DELRAD.DAT', 'TOPOLOGY.H']
#Required files for the specific pdb
proteinScanFilesExtensions = ['BACKGR.DAT', 'DESOLV.DAT', 'MATRIX.DAT', 'PKA.DAT', 'TITCURV.DAT']

#The global configuration instance to be used by default by all other objects.
appConfiguration = None

def resourceDirectory():
	
	'''Returns the directory containing the resources used by the ProteinDesignTool'''

	return resources

def UffbapsResourceDirectory():

	'''Returns the directory containing the files necessary for running genialtNav from UFFBAPS
	
	genialtNav must be run in this directory'''
	
	return os.path.join(resourceDirectory(), 'UFFBAPSRun')
	
def RequiredScanFilesDirectory():
	
	'''Returns the paths to the directory containing the required scan files'''
	
	return os.path.join(resourceDirectory(), 'RequiredScanFiles')

class Configuration:

	'''Creates instances representing an options for a design run.
	
	The returned instance is polymorphic to python ConfigParser.Config parser class.
	That is all that classes methods can also be used with this class.'''
	
	scanSections = ['PKA SCAN PARAMETERS', 'PKA SCAN METHOD', 'PKA SCAN OUTPUT']
	requiredOptions = {'PKA SCAN PARAMETERS':['pKMCsteps', 'recalc_intpka', 'recalc_intpka_dist',
							'use_titration_curves', 'calc_dpka', 'mutation_quality',
							'generate_mutations'],
				'PKA SCAN METHOD':['dpKa_method', 'tabulated', 'MCsteps', 'PBEsolver'], 
				'PKA SCAN OUTPUT':['verbose','save_solutions'], 
				'PATHS':['pKa Tools','uffbaps']}
	
	def __validateConfiguration(self, object, errorList=None):
	
		'''Checks if object is a valid configuration object
		
		Parameters:
			object - A ConfigParser instance
			errorList - An empty list
			
		Returns:
			True if object contains all the required options.
			False otherwise.
			If False and errorList was set if contains a list of what was missing'''
		
		if errorList == None:
			errorList = []
			
		requiredSections = Configuration.requiredOptions.keys()
		
		#Check all sections are present
		#and all required options in that section are present
		for section in requiredSections:
			if object.has_section(section):
				options = object.options(section)
				required = Configuration.requiredOptions[section]
				for option in required:
					if options.count(option) == 0:
						errorList.append('Missing option %s' % option)
			else:
				errorList.append('Missing section %s' % section)
					
		if len(errorList) == 0:
			retval = True
		else:
			retval = False
		
		return retval	

	def __readConfigurationsFromDefaultLocations(self):
	
		'''Checks if a configuration file exists at set of default locations and reads it
		
		The locations are (in order)
			- the current directory
			- the value of PEATSA_CONFIGURATION_PATH environment variable
			- the users home directory
		
		The first found configuration file is used.
		
		Returns:
			The configuration object if one was successfully created.
			None otherwise.
		
		Exceptions:
			Raises an Exceptions.ConfigurationError if the configuration is
			missing any of the required sections.'''
			
		default = ConfigParser.SafeConfigParser()
		default.optionxform = str
	
		#First check if any configuration files exist in default locations
		current = os.path.abspath('proteinDesignTool.conf')
		home = os.path.join(os.path.expanduser('~'), 'proteinDesignTool.conf')
		array = [home, current]
		
		if os.environ.has_key('PEATSA_CONFIGURATION_PATH'):
			env = os.environ['PEATSA_CONFIGURATION_PATH']
			env = os.path.join(env, 'proteinDesignTool.conf')
			array.insert(1, env)
		
		#Note this first read $HOME if it exists and then overrides it with current
		result = default.read(array)
		if len(result) == 0:
			default = None
		else:
			#Check if all necessary sections are present
			list = []
			if not self.__validateConfiguration(default, list):
				raise Exceptions.ConfigurationError, "Configuration file %s not valid.\n%s" % (result, list)	
	
		return default

	def __defaultConfiguration(self, writeFile=True):

		'''Creates a default ConfigParser object instances of this class can use
		
		Parameters -
			writeFile - If True writes a default configuration file in the current directory'''
	
		default = ConfigParser.SafeConfigParser()
		default.optionxform = str

		#Add default configuration values
		default.add_section('PKA SCAN PARAMETERS')
		default.add_section('PKA SCAN METHOD')
		default.add_section('PKA SCAN OUTPUT')
		default.add_section('PATHS')
		default.add_section('WORKING DIRECTORY')
	    
		# pKa calculation parameters
		default.set('PKA SCAN PARAMETERS', 'pKMCsteps', str(200000))
		default.set('PKA SCAN PARAMETERS', 'pKMCsteps', str(200000))
		default.set('PKA SCAN PARAMETERS', 'recalc_intpka', str(1))
		default.set('PKA SCAN PARAMETERS', 'recalc_intpka_dist', str(20))
		default.set('PKA SCAN PARAMETERS', 'use_titration_curves', str(1))
		default.set('PKA SCAN PARAMETERS', 'calc_dpka', str(1))
		default.set('PKA SCAN PARAMETERS', 'generate_mutations', str(False))
		default.set('PKA SCAN PARAMETERS', 'save_temp_files', str(False))
		default.set('PKA SCAN PARAMETERS', 'mutation_quality', str(0.5))

		# Method
		default.set('PKA SCAN METHOD', 'dpKa_method', 'MC')
		default.set('PKA SCAN METHOD', 'tabulated', str(1))
		default.set('PKA SCAN METHOD', 'MCsteps', str(0))
		default.set('PKA SCAN METHOD', 'PBEsolver', 'DelPhi')

		# Be not-so-noisy
		default.set('PKA SCAN OUTPUT', 'verbose', str(1))
		# Do not save the solutions
		default.set('PKA SCAN OUTPUT', 'save_solutions', str(None))
			
		#Set a default location for the pKa Tools
		default.set('PATHS', 'pKa Tools', os.path.join(os.path.expanduser('~'), 'PEAT'))
		default.set('PATHS', 'uffbaps', os.path.join(os.path.expanduser('~'), 'PEAT/UFFBAPS'))
		
		#Working directory configuration (optional)
		#Users may not want to overwrite files in an existing working dir
		#These options configure this behaviour
		default.set('WORKING DIRECTORY', 'copyIfExists', str(0))
		default.set('WORKING DIRECTORY', 'copyLocation', os.getcwd())
		default.set('WORKING DIRECTORY', 'overwriteExistingCopies', str(0))
		default.set('WORKING DIRECTORY', 'useUniqueID', str(0))
	
		#If requested write out this default configuration 
		#wherever the program is being run from
		if writeFile:
			file = open('proteinDesignTool.conf', 'w')
			default.write(file)
			file.close()
		
		return default	
		
	def __init__(self, filename=None, searchDefaultLocations=True, writeFile=False):
	
		'''Initialises a new configuration object.
		
		Parameters
			filename: Optional name of a file containing valid options.
			If a filename is not specified default values are used.
			
			searchDefaultLocations: If filename is None and this is True
			then a number of default places are searched for a configuration file
			
			writeFile: If True a default configuration file will be written
			in the current directory if no filename is passed AND either no configuration file
			is found in a default location OR searchDefaultLocations = False
			
		Exceptions:
			Raises an Exceptions.EnvironmentError if a filename is given and it 
			contains no options'''
	
		self.configuration = None
		self.environment = Environment()
		#If a filename was provided read from it.
		#Otherwise search in default locations or create default options
		if filename == None:
			if searchDefaultLocations:
				self.configuration = self.__readConfigurationsFromDefaultLocations()
			
			#If the above failed or searchDefaultLocations was false create a default configuration		
			if self.configuration == None:	
				print 'No configuration file found - creating default'
				self.configuration = self.__defaultConfiguration(writeFile=writeFile)
		else:
			self.configuration = ConfigParser.SafeConfigParser()
			self.configuration.optionxform = str
			self.configuration.read([filename])
			#Validate the file
			list = []
			if not self.__validateConfiguration(self.configuration, list):
				raise Exceptions.ConfigurationError, "Configuration file %s not valid.\n%s" % (filename, list)
		
		#Add specified pKa Tool path to sys.path
		modulePath = self.configuration.get('PATHS', 'pKa Tools')
		modulePath = os.path.abspath(modulePath)
		self.pKaToolAvailable = True
		if modulePath != None:
			#Check is a directory
			if not os.path.isdir(modulePath):
				self.pKaToolAvailable = False
				self.environment.output("Invalid value given for location of pKa Tools (%s) - file %s. pKa function disabled" % (modulePath, filename))
			#Append it to the path if its not already present
			elif sys.path.count(modulePath) == 0: 
				sys.path.insert(0, modulePath)	
				
		#If this is the first Configuration object created set it as the scripts global configuration
		global appConfiguration
		if appConfiguration is None:
			appConfiguration = self
			
	def __cmp__(self, object):
	
		'''Compares by checking if every section and option in the configuration objects are the same.'''
		
		equal = 0
		
		#If at any stage the objects are detected to be inequal a KeyError is raised.
		try: 
			#Check sections are the same
			sections = object.sections()
			if sections != self.configuration.sections():
				raise KeyError
		
			#Check all the option names are the same in each section
			for section in sections:
				#Order of elements could be different - use sets to compare
				if set(self.configuration.options(section)) == set(object.options(section)):
					#If all the option names are the same check all the option values are the same
					#Use same list of options for both so we get the same order of values
					#in the list comprehensions
					options = self.configuration.options(section)
					valuesOne = [self.configuration.get(section, option) for option in options]
					valuesTwo = [object.get(section, option) for option in options]
					if valuesOne != valuesTwo:
						raise KeyError
				else:
					raise KeyError
		except KeyError:
			#Specifically happens when the object is polymorphic to a ConfigParser
			#But hasn't the same options
			equal = -1
		except BaseException:
			#Happens if the object is not polymorphic
			equal = -1
			
		return equal
	
	def __getattr__(self, aname):
		
		'''Forward all attribute access we don't know about to the configuration object'''
		
		return getattr(self.configuration, aname)
	
	def writeToFile(self, filename):
	
		'''Writes the configuration object to a file
		
		Exceptions
			Raises an IOError if file cannot be created or written to'''		
		
		file = open(filename, 'w')
		self.configuration.write(file)
		file.close
		
	def designPKAOptions(self):
	
		'''Returns the options for the scan calculation in the format used by Design_pKa.run_opt
		
		Exceptions
			Raises an Exceptions.EnvironmentError if the pKa modules cannot be found'''
		
		try:
			import pKa.Design_pKa as Design_pKa
		except ImportError:
			raise Exceptions.ConfigurationError, "Cannot located Design_pKa module. Check correct path specified in configuration"
		
		#Get the default options
		defaults = Design_pKa.get_defaults()
		
		#Go through sections and set values
		for section in Configuration.scanSections:
			for option in self.configuration.options(section):
				#print option
				#print defaults[option]
				#Do some translating
				if self.configuration.get(section, option) == 'None':
					defaults[option][0] = None
					#print 'None conversion - Changed value of %s to %s' % (option, defaults[option])
				elif type(defaults[option][0]) == float:
					defaults[option][0] = float(self.configuration.get(section, option))
					#print 'Float conversion - Changed value of %s to %s' % (option, defaults[option])
				elif type(defaults[option][0]) == int:
					defaults[option][0] = int(self.configuration.get(section, option))
					#print 'Int conversion - Chaged value of %s to %s' % (option, defaults[option])
				elif type(defaults[option][0]) == bool:
					confValue = self.configuration.get(section, option)
					#Accepts 1,0, True and False
					if confValue == 'True':
						confValue = True
					elif confValue == 'False':
						confValue = False
					else:
						confValue = bool(int(confValue))	

					defaults[option][0] = confValue
					#print 'Bool - Changed value of ', option, ' to ', defaults[option]
				else:
					defaults[option][0] = self.configuration.get(section, option)
					#print 'No conversion - Changed value of %s to %s' % (option, defaults[option])
		
		return defaults	

	def uffbapsLocation(self):
	
		'''Returns the location of Uffbaps as given by the configuration file
		
		Defaults to /usr/local/bin if nothing is supplied'''
		
		return self.configuration.get('PATHS', 'uffbaps')
		
	def pkaModuleLocation(self):
	
		'''Returns the location of pKaTool.
		
		Note: This is the path provided by the configuration file. It may not be valid.
		User pKaCalculationsAvailable() to check'''
			
		return self.configuration.get('PATHS', 'pKa Tools')
		
	def pKaCalculationsAvailable(self):
		
		'''Returns True if pKaCalculations are available'''
		
		return self.pKaToolAvailable
			
class WorkingDirectory:

	'''Handles setup of working directory for run of the protein design tool
	
	This involves checking that the necessary files for a given protein are present
	along with generic run files e.g. TOPOLOGY.H etc.
	If generic files are not present they are copied there'''
	
	def __init__(self, location, configuration=None):
	
		'''Creates an instance for performing a run in a specified directory.
		
		The various classes representing the programs to be run do not depend on this class.
		It merely sets up the necessary conditions.
		
		Note that the current working directory is changed to the specified directory on instantiation.
		Also note that internally the class always uses absolute paths so if the current working directory
		is changed it will still function except for the fact that path() will not return the working
		directory.
		
		Many of this classes methods only work in the root process in a parallel environment.
		
		Parameters
		
			location: The path to the working directory
			configuration: An optional Environment.Configuration instance'''
	
		#Create environment instance
		self.environment = Environment()

		#set default behaviour if working directory already exists
		copyIfExists = False
		ignoreCopyError = False
		overwriteExistingCopies = False
		useUniqueID = False
		copyLocation = os.getcwd()

		#override defaults with configuration values if available
		if configuration is not None:
			#Note: Casting a string directly to a bool always gives True
			copyIfExists = int(configuration.get('WORKING DIRECTORY', 'copyIfExists'))
			copyIfExists = bool(copyIfExists)
			useUniqueID = int(configuration.get('WORKING DIRECTORY', 'useUniqueID'))
			useUniqueID = bool(useUniqueID)
			overwriteExistingCopies = bool(configuration.get('WORKING DIRECTORY', 'overwriteExistingCopies'))
			copyLocation = configuration.get('WORKING DIRECTORY', 'copyLocation')
			copyLocation = os.path.abspath(copyLocation)

			try:
				ignoreCopyErrors = int(configuration.get('WORKING DIRECTORY', 'ignoreCopyErrors'))
				ignoreCopyErrors = bool(ignoreCopyErrors)
			except ConfigParser.NoOptionError, data:
				self.environment.output(
					'[PEAT-SA] Configuration file does not contain ignoreCopyErrors option - defaulting to false',
					rootOnly=True)
				ignoreCopyErrors = False
		
		self.directory = os.path.abspath(location)
			
		#If the working dir doesn't exist create it
		#If it does and copyIfExists is True, copy it
		#Otherwise use it as is
		creatingDir = False 
		if not os.path.exists(self.directory):
			creatingDir = True
			if self.environment.isRoot():
				#Create it
				try:
					os.mkdir(self.directory)
				except BaseException:
					raise Exceptions.WorkingDirectoryError, "Working directory - %s - does not exist and cannot be created" % self.path()
		elif copyIfExists is True:
			creatingDir = True
			#Get name for copied dir - all processes enter this code
			destinationName = os.path.split(self.directory)[1]
			if useUniqueID is True:
				if self.environment.isRoot():
					#Cumbersome way of getting a unique file necessary since mktmp was deprecated
					#We use the original dir name so we can id the tmp dir later
					temp = tempfile.NamedTemporaryFile(prefix='%s_' % destinationName, dir=copyLocation)
					destinationName = temp.name 
					temp.close() 
				destinationName = self.environment.broadcast(data=destinationName, process=0)

			originalLocation = self.directory
			#Reassign self.directory to new value
			self.directory = os.path.join(copyLocation, destinationName)
			
			self.environment.output(
				'[PEAT-SA] Specified working directory exists - creating copy at %s' % self.directory,
				rootOnly=True)

			#only root copies the directory
			if self.environment.isRoot():				
				if os.path.exists(self.directory) and not overwriteExistingCopies:
					raise Exceptions.WorkingDirectoryError, "Directory %s exists and cannot create copy at %s" % (originalLocation, self.directory())
				else:
					#Delete existing dir
					shutil.rmtree(self.directory, ignore_errors=True)
				
				#Copy the directory to the new destination
				try:
					shutil.copytree(originalLocation, self.directory , symlinks=True)
				except shutil.Error, data:
					self.environment.output(
						'[PEAT-SA] Encountered error while copying',
						rootOnly=True)
					if ignoreCopyErrors is True:
						self.environment.output( '[PEAT-SA] Errors were:', rootOnly=True)
						for copyError in data:
							self.environment.output( '\t%s' % data, rootOnly=True)
	
						self.environment.output(
							'[PEAT-SA] Ignoring errors and trying to continue',
							rootOnly=True)
					else:
						raise
						 
		else:
			self.environment.output(
				'[PEAT-SA] Specified working directory exists - using it',
				rootOnly=True)
					
		#Wait for roots copy/creation to finish
		self.environment.wait()				
			
		#Never executed if the directory creation fails
		if creatingDir:
			self.environment.log('Root finished creation of directory %s - checking if its appeared' % self.directory)
			while not os.path.isdir(self.directory):
				self.environment.log('Directory hasnt appeared yet')
				time.sleep(0.5)			

		#Switch current working directory to the supplied directory
		os.chdir(self.path())	
	
	def _removeErrorHandler(self, function, path, execInfo):
	
		'''Convience method for creating an remove directory exception
		
		Note: Do not use as handler to shutil.rmtree as the exception raised here
		cannot be caught then'''
		
		data = 'Error removing directory %s\n' % path
		data = data + 'Exception information - %s\n' % str(execInfo)
		
		raise Exceptions.WorkingDirectoryError, data
	
	def _checkFile(self, pdbFile):
	
		'''Checks if the path pdbFile refers to an existing file.
		FIXME - Make changing of current working directory understandable.
		
		Parameters:
			pdbFile - A path to a file - Can be absolute or not.
			
		Returns:
			The absolute path to the file.
			
		Exceptions:
			Raises an exception if the file does not exist or if its not a file'''
			
		#Make the filename absolute
		if not os.path.isabs(pdbFile):
			pdbFile = os.path.abspath(pdbFile)
		
		#Check the specified file exists.
		if not os.path.exists(pdbFile):
			raise Exceptions.WorkingDirectoryError, 'File %s does not exist' % pdbFile
			
		#Check its a file
		if not os.path.isfile(pdbFile):
			raise Exceptions.WorkingDirectoryError, 'Object at %s is not a file' % pdbFile	
		
		return pdbFile
		
	def _copyFileToDirectory(self, filename, directory, overwrite=False):
	
		'''For internal use only. Copies filename to directory.
		
		Note that since this is for internal use error checking is less vigourous.
		
		Parameters
			filename - A path to a file. Can be relative or absolute.
			If the path is relative its is resolved w.r.t the specified directory.
			directory - The directory to copy the file to.
			overwrite - If True, if a file exists in the directory with filename
				it is overwritten if possible. Otherwise that file is used.
				Default - False
		Returns:
			The full path to the copied file
			
		Exceptions:
			Raises an Exceptions.ArgumentError if the file does not exist.
			Raises an Excpetions.WorkingDirectoryError if the file could not be copied.'''
	
		#Only the root process does the copying
		destination = os.path.join(directory, os.path.split(filename)[1])
		self.environment.log('Copying %s to %s' % (filename, destination))
		if self.environment.isRoot():
			filename = self._checkFile(filename)
			
			#Check if the specified file to is in the working directory
			#If it is we don't do anything
			containingDirectory = os.path.split(filename)[0]
			if containingDirectory != directory:
				#Check if a copy of the file already is present
				#Again do nothing if this is so
				directoryContents = os.listdir(directory)
				if directoryContents.count(filename) is 1 and overwrite is True:
					#This automatically overwrites destination if it exists 
					#and the operation is possible
					shutil.copyfile(filename, destination)
				elif directoryContents.count(filename) is 0:
					shutil.copyfile(filename, destination)
		
		self.environment.wait()
		return destination
		
	def copyFileToDirectory(self, filename, overwrite=False):
	
		'''Copies filename to the working directory.
		
		Parameters
			filename - A path to a file. Can be relative or absolute.
			If the path is relative its is resolved w.r.t the working directory.
			overwrite - If True, if a file called filename exists in the working directory
				it is overwritten if possible. Otherwise that file is used.
				Default - False
		Returns:
			The full path to the copied file
			
		Exceptions:
			Raises an Exceptions.ArgumentError if the file does not exist.
			Raises an Excpetions.WorkingDirectoryError if the file could not be copied.'''
	
		return self._copyFileToDirectory(filename=filename, directory=self.path(), overwrite=overwrite)
	
	def copyFileToUffbapsDirectory(self, filename, overwrite=False):
	
		'''Copies filename to the uffbaps run directory.
		
		Parameters
			filename - A path to a file. Can be relative or absolute.
			If the path is relative its is resolved w.r.t the uffbaps run directory.
			overwrite - If True, if a file called filename exists in the uffbaps run directory
				it is overwritten if possible. Otherwise that file is used.
				Default - False
		Returns:
			The full path to the copied file
			
		Exceptions:
			Raises an Exceptions.ArgumentError if the file does not exist.
			Raises an Excpetions.WorkingDirectoryError if the file could not be copied.'''
	
		return self._copyFileToDirectory(filename=filename, directory=self.uffbapsRunDirectory(), overwrite=overwrite)	
		
	def containsFile(self, filename):
	
		'''Returns True if the directory contains filename - False otherwise
		
		Note: Only check in the top-level directory'''
		
		filename = os.path.abspath(filename)
		path = os.path.split(filename)[0]
		isPresent = False
		
		if path == self.directory:
			isPresent = True
		
		return isPresent
	
	def setupPKARun(self):
	
		'''Checks that the directory contains all the necessary files for a pKa calculation
		
		Exceptions:
			Raises Exceptions.WorkingDirectoryError if there was a problem with the setup.'''

		#Only the root process does the copying
		if self.environment.isRoot():
			directoryContents = os.listdir(self.path())
		
			#Check that all the pka non-protein specific files are present
			#If any is copy the resource dir version
			for requiredFile in requiredPKADesignFiles:
				if directoryContents.count(requiredFile) == 0:
					copiedRequiredFiles = True
					resourceVersion = os.path.join(RequiredScanFilesDirectory(), requiredFile)
					shutil.copyfile(resourceVersion, os.path.join(self.path(), requiredFile))
					
		self.environment.wait()			
	
	def setupScan(self, pdbName, mutantCollection):
					
		'''Checks that the directory contains all the necessary files for a scan run on pdbName
		
		Note that is the mutantCollection is None does not cause this method to fail.
		
		Parameters:
			pdbName - The name of the pdbFile the scan will be run on e.g. 5P21.pdb
			A pdb with the name must be in the working directory.
			mutantCollection - A MutantCollection instance containing the mutant files the scan will be run on
		
		Exceptions:
			Raises an Exceptions.MissingPKADataError if the required pKa data is missing.
			Raises Exceptions.WorkingDirectoryError if there was any other problem with the setup.'''

		#Only the root process does the copying
		if self.environment.isRoot():

			#Accumulate all errors into one string so the user finds out 
			#everything that is wrong in one go
			errorString = ""
		
			#Get the directory contents
			directoryContents = os.listdir(self.path())
			
			#Flag tracking if we had to copy the required protein design files to the working directory
			copiedRequiredFiles = False
		
			#Check that all the pka non-protein specific files are present
			#If any is copy the resource dir version
			for requiredFile in requiredPKADesignFiles:
				if directoryContents.count(requiredFile) == 0:
					copiedRequiredFiles = True
					resourceVersion = os.path.join(RequiredScanFilesDirectory(), requiredFile)
					shutil.copyfile(resourceVersion, os.path.join(self.path(), requiredFile))
			
			#Check that all the pka protein specific files are present
			for extension in proteinScanFilesExtensions:
				filename = pdbName + "." + extension
				if directoryContents.count(filename) == 0:
					errorString = errorString + "Required scan file %s not present\n" % filename
						
			if len(errorString) != 0:
				#If we copied the required protein design files remove them
				if copiedRequiredFiles:
					for file in requiredPKADesignFiles:
						os.remove(os.path.join(self.path(), file))
				self.clean()
				raise Exceptions.MissingPKADataError, errorString
				
			#Check the pdbFile is present
			path = os.path.join(self.directory, pdbName)
			if not self._checkFile(path):
				raise Exceptions.WorkingDirectoryError, "No pdb file called %s in the working directory (%s)" % (pdbName, self.directory)
		
			#Link the mutants to the required scan directory
			if mutantCollection is not None:
				mutantPath = os.path.join(self.directory, "%s.pdbs" % pdbName)
				self.environment.output('Linking created mutants to %s' % mutantPath)
					
				#Check if a link exists - if it does delete it
				#Move previous scan dirs to a new path
				if os.path.islink(mutantPath):
					os.remove(mutantPath)
				elif os.path.isdir(mutantPath):
					renamedPath = os.path.join(self.directory, "%s.pdbs_old" % pdbName)
					if not os.path.isdir(renamedPath):
						os.rename(mutantPath, renamedPath)	
						self.environment.output('Moved existing mutants at %s to %s' % (mutantPath, renamedPath))
					else:
						self.environment.output('Leaving  mutants present at %s' % renamedPath)
						self.environment.output('Deleting mutants present at %s' % mutantPath)
						shutil.rmtree(mutantPath, ignore_errors=True)
				
				#Create a link to the mutant collection dir
				os.symlink(os.path.join(mutantCollection.location, "Mutants"), mutantPath)
				self.environment.output('Link created')
		
		self.environment.wait()		
		
		#If we are in parallel every process copies the necessary files
		#to its own subdir of this directory and changes to work in it
		if self.environment.isParallel:
			self.environment.output('Current path %s. Current dir %s' % (self.path(), os.getcwd()), rootOnly=True)

			destination = os.path.join(self.path(), "Process%d" % self.environment.rank())
			#Check is this already exists - could have been a problem with a previous clean
			self.environment.log('Creating process specific directory at %s' % destination)
			if os.path.exists(destination):
				self.environment.output('Deleting previous process specific directory at %s' % destination, rootOnly=True)
				shutil.rmtree(destination, ignore_errors=True)

			template = os.path.join(self.path(), "template")
				
			if self.environment.isRoot():
				if os.path.exists(template):
					self.environment.output('Deleting previous template directory at %s' % template, rootOnly=True)
					shutil.rmtree(template, ignore_errors=True)
					
				#Create process temp copy without the sim links
				shutil.copytree(self.path(), template, symlinks=True)

			self.environment.wait()			
			shutil.copytree(template, destination, symlinks=True)

			self.environment.wait()
			if self.environment.isRoot():
				shutil.rmtree(template, ignore_errors=True)

			os.chdir(destination)
			self.environment.output('After move: Current path %s. Current dir %s' % (self.path(), os.getcwd()), rootOnly=True)
	
	def cleanUpScan(self):
	
		'''Performs clean-up of the WorkingDirectory after a scan
		
		This basically involves removing any per-parallel-process copies
		of the pKa data'''
	
		if self.environment.isParallel:
			os.chdir(self.path())
			directory = os.path.join(self.path(), "Process%d" % self.environment.rank())
			try:
				shutil.rmtree(directory, ignore_errors=True)
			except Exception, data:
				print Exception, data
			
	def setupUFFBAPS(self):
		
		'''Checks that the directory contains all the necessary files for a UFFBAPS run.
		
		Exceptions:
			Raises Exceptions.WorkingDirectoryError if there was a problem with the setup.'''
	
		if self.environment.isRoot():
			#There is a chance the directory will already be present - delete it
			if os.path.isdir(self.uffbapsRunDirectory()):
				shutil.rmtree(self.uffbapsRunDirectory(), ignore_errors=True)
			
			#Copy all the neccessary stuff for running Chrestens stability tool
			try:
				shutil.copytree(UffbapsResourceDirectory(), self.uffbapsRunDirectory())
				errorString = ''
			except BaseException, data:
				print 'Encountered an error when copying UFFBAPS data - ', data
				errorString = "Unable to move neccessary files for stability calculation to working directory"
				raise Exceptions.WorkingDirectoryError, errorString
				
		self.environment.wait()		
			
	def setup(self, pdbFile, mutantCollection=None, scan=True, stability=True, modes=False):
	
		'''Checks that the directory contains all the necessary files for the calculations specified on initialisation
		
		Note: Unlike the other setup functions overwrite defaults to True in here.
				
		Parameters:
			pdbFile - The absoloute path to the pdbFile the scan will be run on
			scan - True if a scan is to be performed. Default True
			
			stability - True if a stability calculation is to be performed
				Default True. May not be possible if a scan is not performed.
			
			modes - True of a modes calculation is to be performed. Default False
				Note: No effect currently	

		Exceptions:
			Raises Exceptions.WorkingDirectoryError if there is a problem with the setup.'''
	
		if scan:
			self.setupScan(pdbFile, mutantCollection=mutantCollection)
				
		if stability:
			self.setupUFFBAPS()
		
	def clean(self):
	
		'''Cleans the working directory of unnecessary files.'''
		
		if self.environment.isRoot():
			shutil.rmtree(self.uffbapsRunDirectory(), ignore_errors=True)
			
		#Remove scan files?
		
		
	def path(self):
	
		'''Returns the path of the working directory'''
		
		return self.directory	
	
	def uffbapsRunDirectory(self):
	
		'''Returns the path of the subdir where genialtNav will be run'''

		return os.path.join(self.path(), 'UFFBAPSRun')
					

class Environment:

	'''Class representing the application environment.
	
	Environment objects hide details about the number of processes in the environment from other Core objects.
	That is Core objects use Environment objects in the same way if the program is running in serial or parallel.
	This means all MPI code is contained in this class.
	The other Core objects use Environment objects to output text, split arrays etc.'''

	isInitialised = False
	isParallel = False
	
	def __init__(self, verbose=True):
	
		'''Initialises new Environment objects'''
	
		self.outputDirectory = os.getcwd()
		self.verbose = verbose
	
		#Do the subsequent things once only
		if not Environment.isInitialised:
			#See if a parallel environment is available
			try:
				import pypar
				#If its imported we might have a parallel environment
				Environment.isParallel = True
				self.output('[PEAT-SA] Parallel environment available')
				#Check environment size. 
				#If there is more than one processor there must be a parallel environment
				#If there's only one then its not parallel.
				if pypar.size() == 1:
					self.output('[PEAT-SA] Only one processor - parallel environment disabled')
					Environment.isParallel = False
				else:
					self.output('[PEAT-SA] Parallel environment enabled with %d processors' % pypar.size())
			except BaseException:
				#Importing pypar caused an exception - No parallel environment
				Environment.isParallel = False
				self.output('[PEAT-SA] Parallel environment disabled.\n')

			Environment.isInitialised = True	 

	def rank(self):
	
		'''Returns the rank of the process in the environment
		
		This is 0 if there is only one process and for the root processor'''
	
		if Environment.isParallel:
			import pypar
			return pypar.rank()
		else:
			return 0
			
	def isRoot(self):
	
		'''Returns True if the processes rank is 0'''
			
		if self.rank() == 0:
			return True
		else:
			return False
			
	def _divideArray(self, array):
	
		'''Divides an array roughly equally depending on the environment size
		
		Returns: A list with one entry for each node.
		The entry is a tuple giving the start and end elements in array 
		that should be assigned to that node.'''
	
		import pypar
	
		#Divide evenly then add remainder elements to processors
		maxElements = int(math.floor(len(array)/pypar.size()))
		remainder = len(array) - maxElements*pypar.size()
		start = 0
		end = 0
		divisions = []
		for i in range(pypar.size()):
			start = end
			end = end + maxElements
			if remainder != 0:
				end = end + 1
				remainder = remainder - 1
			divisions.append((start, end))
			
		return divisions									
			
	def splitArray(self, array):
	
		'''Splits array between all the processes in the environment.
		
		Each process will be returned a different section of the array to work on'''
	
		if Environment.isParallel:
			import pypar
			#Split the array into sections and return the section for this processor
			divisions = []
			if self.isRoot():
				#Root does the splitting - we send each processor the start and end index
				#NOTE: pypar broadcast won't work even when setting vanilla
				#It always returns message trucated error.
				divisions = self._divideArray(array)	
				for i in range(1,pypar.size()):
					pypar.send(divisions[i], i)
				start = divisions[0][0]
				end = divisions[0][1]
			else:	
				indexes = pypar.receive(0)
				start = indexes[0]
				end = indexes[1]
				
			return array[start:end]
		else:
			return array
			
	def combineArray(self, arrayFragment):
	
		'''Combines a set of arrayFragments from each processor into one array'''

		if Environment.isParallel:
			import pypar
			if self.isRoot():
				completeArray = arrayFragment
				for i in range(1, pypar.size()):
					fragment = pypar.receive(i)
					completeArray.extend(fragment)
				
				#Send the array
				for i in range(1, pypar.size()):
					pypar.send(completeArray, i)
			else:
				#Send the fragment
				pypar.send(arrayFragment, 0)
				#Retrieve the array
				completeArray = pypar.receive(0)
		else:
			completeArray = arrayFragment
		
		return completeArray	
		
																				
	def output(self, string, stream=None, rootOnly=True):
	
		'''Prints a string to a stream if the calling process is the root process
		
		Parameters:
			string: The string to be written
			stream: The stream to write to. If None or not specified defaults to stdout.
			rootOnly: If True only the root process writes the string. Default True.'''
	
		if rootOnly is True and self.isRoot() is False:
			return
		
		if stream is None:
			print string
			#With verbose on flush everything immediately
			if self.verbose:
				sys.stdout.flush()	
		else:
			stream.write(string)
			if self.verbose:
				stream.flush()

	def processorName(self):

		if Environment.isParallel:
			import pypar
			return pypar.get_processor_name()
		else:
			return "localhost"
	
	def log(self, string):
	
		'''Logs a string to a specific file for the calling process
		
		The file is called ProcessorX.log where X is the rank of the process.
		The string is appended to this file'''

		from inspect import stack
	
		filename = os.path.join(self.outputDirectory, "Processor%d.log" % self.rank())
		stream = open(filename, 'a+')
		stream.write(string)	
		stream.write(' (Logged at line %d of %s at %s)' % (stack()[1][0].f_lineno, 
				stack()[1][0].f_code.co_filename, datetime.datetime.now().strftime('%H:%M:%S')))
		stream.write("\n")	
		stream.close()		
		
	def logError(self, string):
	
		'''Logs a string to a specific error file for the calling process
		
		The file is called ProcessorX.error where X is the rank of the process.
		The string is appended to this file'''
		
		filename = os.path.join(self.outputDirectory, "Processor%d.error", self.rank())
		stream = open(filename, 'a+')
		stream.write(string)		
		stream.close()
		
	def wait(self, error=False):
	
		'''This method will not return until all process in the environment have called it.
		
		This is a wrapper around MPI_Barrier which handles the case where MPI is not available'''

		from inspect import stack
	
		if self.verbose is True:
			string = '(%s) Waiting at line %d of %s' % (datetime.datetime.now().strftime('%H:%M:%S'),
					 stack()[1][0].f_lineno, stack()[1][0].f_code.co_filename)
			self.log(string)
			
		if Environment.isParallel:
			import pypar
			pypar.barrier()	
			#Because MPI_ABORT doesn't work in pypar if called from one process
			#we need a way for process to communicate to each other if an error occurred
			#during the code they executed before this barrier. We do a scatter/gather of
			#the error parameter - This isn't very efficient but it's all we can do now
			errors = self.combineArray([error])
			if True in errors:
				self.exit(1)
				
		if self.verbose is True:
			string = '(%s) Finished waiting' % (datetime.datetime.now().strftime('%H:%M:%S'))
			self.log(string)		
			
	def exit(self, code):
	
		'''This method exits the simulation.
		
		In a parallel environment calls MPI_ABORT.
		In serial, calls sys.exit().
		
		Code is the exit code. Only used in serial processes.'''
		
		if Environment.isParallel:
			import pypar
			return pypar.abort(code)
		else:
			return sys.exit(code)
			
	def broadcast(self, data, process):
	
		'''Broadcasts data from process to all other nodes'''
		
		if Environment.isParallel:
			import pypar
			if self.rank() == process:
				#NOTE: pypar broadcast won't work even when setting vanilla
				#It always returns message trucated error.
				for i in range(pypar.size()):
					if i != self.rank():
						pypar.send(data, i)
			else:	
				data = pypar.receive(process)
		
		return data
		
	def balanceArrays(self, arrayFragment):
	
		'''Redistributes the elements in a set of arrays equally across the nodes'''
		
		if Environment.isParallel:
			import pypar
			if self.isRoot():
				completeArray = arrayFragment
				for i in range(1, pypar.size()):
					fragment = pypar.receive(i)
					completeArray.extend(fragment)
				
				#Divide it up
				divisions = self._divideArray(completeArray)
				
				#Send the fragments
				for i in range(1, pypar.size()):
					start, end = divisions[i]
					pypar.send(completeArray[start:end], i)
					
				self.output('[ENV] Rebalanced array divisions %s' % divisions)	
					
				#Assign root fragment	
				start, end = divisions[0]	
				arrayFragment = completeArray[start:end]
			else:
				#Send the fragment
				pypar.send(arrayFragment, 0)
				#Retrieve the array
				arrayFragment = pypar.receive(0)
		else:
			completeArray = arrayFragment
		
		return arrayFragment
		
	def loadSkew(self, numberTasks):
	
		'''Computes the skew in the number of tasks to be processed by each node
		
		The skew is the standard deviation of the task number across the nodes'''				
										
		if Environment.isParallel:
			import pypar
			if self.isRoot():
				taskDistribution = [numberTasks]
				for i in range(1, pypar.size()):
					numberTasks = pypar.receive(i)
					taskDistribution.append(numberTasks)
				
				mean = reduce(operator.add, taskDistribution)
				mean = mean/float(len(taskDistribution))
				
				#Std. dev
				stdev = 0
				for el in taskDistribution:
					stdev = stdev + math.pow((el - mean), 2)
					
				skew = stdev/float(len(taskDistribution))
				skew = math.sqrt(skew)	
				
				for i in range(1, pypar.size()):
					pypar.send(skew, i)
				
			else:
				#Send the fragment
				pypar.send(numberTasks, 0)
				#Retrieve the array
				skew = pypar.receive(0)
		else:
			skew = 0								
										
		return skew
											
		
