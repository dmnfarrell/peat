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

'''Module containing classes representing the PEAT_SA's model data'''
import sys, os, copy, shutil, pickle, tempfile, StringIO
import Matrix, Exceptions, Environment, Utilities

peatMatrices = ['BindingResults', 'DeltaBindingResults', 'ScanResults', 'StabilityResults']

def scanCollection(pdbFile, mutation, skipResidueTypes=['ALA', 'GLY'],
			ligandFiles=[], location=None, name=None,
			clean=True, temporary=False, overwrite=True):	
			
	'''Convience function for creating a set of scan mutants.'''
	
	mutationList = CreateScanList(pdbFile, mutation=mutation, skipResidueTypes=skipResidueTypes)
	return MutantCollection(pdbFile, mutationList=mutationList, ligandFiles=ligandFiles, location=location, 
				name=name, clean=clean, temporary=temporary, overwrite=overwrite)						

def ParseMutationList(filename):
	
	'''Parses the mutation lists in filename.
	
	Note: Currently no checking is performed to see if the format or content is valid'''
	
	stream = open(filename)
	#Simple parse
	lines = [line.strip('\n') for line in stream.readlines()]
	stream.close()
	lines = [line.split(",") for line in lines]
	mutationList = []
	for line in lines:
		code = line[0].strip()
		chain, residueIndex = Utilities.ParseResidueCode(code)
		mutations = [element.strip() for element in line[1:]]
		for mutation in mutations:
			listElement = MutationSet()
			listElement.addMutation(chain, residueIndex, mutation)
			mutationList.append(listElement)
		
	return mutationList
		
def CreateScanList(pdbFile, mutation='ALA', skipResidueTypes=['ALA', 'GLY']):

	'''Creates a list representing the mutation of each residue to the specified residue
	
	Errors:
		Raises an FileFormatError if the pdb file is not valid'''

	if mutation == "ALL":
		return CreateCompleteScanList(pdbFile)

	#Create a protool instance of the pdb
	import Protool
	try:
		pdb = Protool.structureIO()
		pdb.readpdb(pdbFile)
	except Exception, data:	
		raise Exceptions.FileFormatError, 'Format of specified PDB file %s not valid.\nUnderlying error - %s' % (pdbFile, data)

	#Create a list where each element is a residue, mutation pair.
	#The mutation in each case is the same and is given by the mutation parameter	
	#The keys returned in the next line have the format ChainID:ResidueNumber	
	residues=pdb.residues.keys()
	residues.sort()
			
	#Convert the passed in list or tuple to a set
	skipResidueTypes = set(skipResidueTypes)
	#Remove from residues
	residues = filter(lambda x: pdb.resname(x) not in skipResidueTypes, residues)
	
	mutations = [[mutation]]*len(residues)
	data = zip(residues, mutations)
	mutationList = []
	for residueCode, list in data:
		chain, residueIndex = Utilities.ParseResidueCode(residueCode)
		for element in list:
			listElement = MutationSet()
			listElement.addMutation(chain=chain, residueIndex=residueIndex, mutation=element)
			mutationList.append(listElement)
	
	return mutationList
	
def CreateCompleteScanList(pdbFile):

	#Create a protool instance of the pdb
	import Protool
	try:
		pdb = Protool.structureIO()
		pdb.readpdb(pdbFile)
	except Exception, data:	
		raise Exceptions.FileFormatError, 'Format of specified PDB file %s not valid.\nUnderlying error - %s' % (pdbFile, data)

	#Create a list where each element is a residue, mutation pair.
	#The mutation in each case is the same and is given by the mutation parameter	
	#The keys returned in the next line have the format ChainID:ResidueNumber	
	residues=pdb.residues.keys()
	residues.sort()
			
	mutations = [copy.copy(Utilities.allResidues) for element in residues]	
	data = zip(residues, mutations)
	
	for entry in data:
		name = pdb.resname(entry[0])
		index = entry[1].index(name)
		entry[1].pop(index)
	
	mutationList = []
	for residueCode, mutations in data:
		chain, residueIndex = Utilities.ParseResidueCode(residueCode)
		for mutation in mutations:
			listElement = MutationSet()
			listElement.addMutation(chain=chain, residueIndex=residueIndex, mutation=mutation)
			mutationList.append(listElement)
	
	return mutationList

def mutationListFileFromStream(stream):

	'''Returns a mutation list instance initialised from stream
	
	The returned MutationListFile instance is associated with a temporary filename.
	This file is not created unless writeToFile() is called on the returned object'''
	
	return MutationListFile(filename=None, mutationData=stream.read(), update=False)

class MutationListFile:

	'''Represents a set of mutants specified in a specially formatted file
	
	There are two possible formats: SinglePointMutation and Standard.
	This class can be used to read an existing, or to create a new, mutation list file.'''
	
	def _isSinglePointMutationFormat(self, stream):
	
		'''Returns True if filename is in single point mutation format
		
		Note: This is a very simple check.
		It only looks to see if the first line from stream is in a valid SPM format'''
		
		isSPM = True
		
		stream.seek(0)		
		line = stream.readline()
	
		components = line.split(",")
		#Check the first component
		try:
			if len(components) != 0:
				code = components[0].strip()
				parts = Utilities.ParseResidueCode(code)
				#should only have two components Chain:Number for SPM
				if len(parts) == 3:
					isSPM = False
			else:
				isSPM = False
		except BaseException:
			isSPM = False
			
		#Check the mutation codes	
		if isSPM is not False:
			mutations = [element.strip() for element in components[1:]]	
			for mutation in mutations:
				if Utilities.allResidues.count(mutation.upper()) is not 1:
					isSPM = False
					break
					
		return isSPM				
	
	def _parseSinglePointMutationList(self, stream):
	
		'''Parses a mutation list in simple format
		
		Assumes the data read from stream is in single point mutation format.
		
		Exceptions:
			Raises an Exceptions.MutationListFileFormatError if the file format is invalid'''
		
		stream.seek(0)
		lines = [line.strip('\n') for line in stream.readlines()]
		lines = [line.split(",") for line in lines]
		lines = [line for line in lines if line != ""]
		for line in lines:
			code = line[0].strip()
			try:
				#Use ParseStandardResidueCode as we know that this is the format it should be in
				#A standard code with a missing : could be misidentified as reduced e.g. A0053
				chain, residueIndex = Utilities.ParseStandardResidueCode(code)
			except Exceptions.ResidueCodeFormatError, data:
				raise Exceptions.MutationListFileFormatError, "Detected an error with a residue code in single point mutation list.\n%s" % data
			except ValueError, data:
				pass
				
			mutations = [element.strip() for element in line[1:]]
			for mutation in mutations:
				listElement = MutationSet()
				try:
					listElement.addMutation(chain, residueIndex, mutation)
					self.addMutant(listElement, autoUpdate=False)
				except Exceptions.MutationCodeFormatError, data: 
					raise Exceptions.MutationListFileFormatError, "Detected an error with a mutation code in single point mutation list.\n%s" % data
		
	def _parseStandardMutationList(self, stream):
	
		'''Parses a mutation list in standard format'''
		
		stream.seek(0)
		lines = [line.strip('\n') for line in stream.readlines()]
		lines = [line.strip() for line in lines]
		lines = [line for line in lines if line != ""]

		for line in lines:
			try:
				set = MutationSet(code=line)
				self.addMutant(set, autoUpdate=False)
			except Exceptions.MutationCodeFormatError, data:
				raise Exceptions.MutationListFileFormatError, "Detected an error with a mutation code in standard mutation list.\n%s" % data
	
	def _parseMutationList(self, stream):
	
		'''Parses the mutation list in stream. 
		
		The instance variable self.mutants is filled during the parse
		
		Exceptions
			Raises an Exceptions.MutationListFileFormatError if there is an error with the files format.
			Raises an Exceptions.MutationListDuplicateEntryError if the file contains duplicates'''
			
		if self._isSinglePointMutationFormat(stream):
			self._parseSinglePointMutationList(stream)
		else:
			self._parseStandardMutationList(stream)
	
	def __init__(self, filename, create=False,  mutationData=None, update=False):
	
		'''Initialises a new MutationListFile instance 
		
		Parameters: 
			filename -
				The name of the file to read/write the data to/from
				If None a temporary filename is generated
			
			create -		
				Indicates if the supplied filename already exists or not
				If create is false this file must exist and contain data in a valid format
				If create is true filename is used as the default argument for writeToFile().
				Note in this second case the file is not created until the first write is requested.
				This parameter is automatically set to True if filename is None
										
			mutationData - 
				A mutation list as a string. The codes in this string are parsed
				after the file (it it exists). 
				Note: If you have mutation set you wish to add to the list use addMutant
				
			update -
				If True the current mutation list is written to the file.
				The file is created if it doesn't exist.
				Note: If the file exists and mutationData is None, the files contents won't change
			
		Exceptions:
			Raises an Exceptions.MutationListError if the file doesn't exists and create is false'''
	
		self.mutants = []
		
		if filename is None:
			#Override create
			create = True
			temp = tempfile.NamedTemporaryFile()
			filename = temp.name
			temp.close()
	
		self.filename = filename

		#First parse the file if it is present
		#The try block will raise an IOError if the file does not exist
		try: 
			Utilities.CheckFile(filename)
			self.fileExists = True
			#self.mutants is filled during the parsing
			stream = open(filename)
			self._parseMutationList(stream)
			stream.close()
		except IOError:
			if create is True:
				self.fileExists = False
			else:
				raise Exceptions.MutationListError, "Specified mutation list file, %s, does not exist" % filename
		
		#Now parse any string data
		if mutationData != None:
			mutationDataStream = StringIO.StringIO(mutationData)
			self._parseMutationList(mutationDataStream)
			mutationDataStream.close()
		
		#If update is True write everything to the file
		if update is True:
			self.writeToFile()	
					
	def __str__(self):
	
		'''Returns a descriptive string containing information on the receiver'''
		
		if self.fileExists:
			string = 'Mutation List File at %s\n' % self.filename
		else:
			string = 'Mutation List File at %s (creation pending)\n' % self.filename
			
		string = string + 'Format %s\n' % self.format()
		string = string + 'Contains %d entries' % len(self.mutants)	
		
		return string
	
	def __cmp__(self, object):
	
		'''Returns true if object describes the same set of mutants as the receiver'''
	
		retval = 0
		if hasattr(object, 'mutantList'):
			for mutant in self.mutantList():
				if object.mutantList().count(mutant) != 1:
					retval = -1
					break
		else:
			retval = -1
			
		return retval	

	def addMutant(self, mutant, autoUpdate=True, ignoreDuplicates=False):
	
		'''Adds an entry for mutant to the mutation list file.
		
		Parameters: 
			mutant - A MutationSet instance.
			autoUpdate - The underlying file is updated automatically.
			
		Errors: 
			If the file is not writable an IOError is raised.
			If mutant is not a MutationSet instance an Exceptions.MutationListEntryError is raised'
			If the receiver already contains the mutant an Exceptions.MutationListDuplicatEntryError is raised
			'''
		if not ignoreDuplicates:	
			if self.mutants.count(mutant) is not 0:
				raise Exceptions.MutationListDuplicatEntryError, "Mutant %s is already present in the list" % mutant
	
		self.mutants.append(mutant)

		if autoUpdate:
			try:
				self.writeToFile()
			except Exceptions.MutationListEntryError:
				self.mutants.remove(mutant)
				raise
		
	def removeMutant(self, mutant):
	
		'''Removes the entry for mutant from the file	
		
		Does nothing if no entry for mutant exists'''
		
		if self.mutants.count(mutant) == 1:
			self.mutants.remove(mutant)
			self.writeToFile()
	
	def mutantList(self):
	
		'''Returns a list of MutationSet instances corresponding to the mutants specified in the file'''
		
		return copy.copy(self.mutants)
		
	def stringRepresentation(self):
	
		'''Returns the string that would be written if writeToFile() was called'''
		
		stream = StringIO.StringIO()
		self.writeToStream(stream)
		stream.flush()
		data = stream.getvalue()
		stream.close()	
		
		return data
		
	def removeDuplicates(self, autoUpdate=True):
	
		'''Removes any duplicates in the current in memory list and optionally the underlying file'''
		
		seen = {} 
		result = [] 
		for set in self.mutants: 
			#Check if set is in the seen dictionary				
			if set.filename(reduced=True) in seen:
				continue 
			#First time seeing this set - add it to seen	
			seen[set.filename(reduced=True)] = 1 
			result.append(set) 
		
		self.mutants = result
		if autoUpdate is True:
			self.writeToFile()		

	def _writeSinglePointMutationListToStream(self, stream):
	
		'''Write the mutants to stream in single point mutation format'''
	
		residues = {}
		#Find which residues have SPM's
		for set in self.mutantList():
			residueCode = set.residueCodes(reduced=False)[0]
			chain, index, mutation = Utilities.ParseMutationCode(set.mutationCodes(reduced=False)[0])
			if not residues.has_key(residueCode):
				residues[residueCode] = []
			residues[residueCode].append(mutation)
		#Write out the mutants of each residue	
		for key in residues.keys():
			string = "%s" % key
			for mutation in residues[key]:
				string = string + " ,%s" % mutation
			string = string + '\n'
			stream.write(string)
			
	def _writeStandardMutationListToStream(self, stream, reduced=True):
	
		'''Writes the mutations to stream in standard format'''
		
		for set in self.mutantList():
			stream.write("+".join(set.mutationCodes(reduced=True)))
			stream.write("\n")

	def writeToFile(self, filename=None, reduced=True):
	
		'''Writes the contents of the receiver to filename as a mutation list.
		
		Parameters
			filename: The file to write to.
			If None the information is written to the default file provided on initialisation
			reduced: If True list is written in reduced format (if it is a Standard Format list)
			
		Exceptions:
			Raises an IOError if the file can't be written to for any reason'''
	
		if filename is None:
			filename = self.filename
	
		stream = open(filename, "w+")
		self.writeToStream(stream, reduced=reduced)
		stream.close()
		
		#If we wrote to the file provided on initialisation
		#set fileExists to True
		if filename == self.filename:
			self.fileExists = True
		
	def writeToStream(self, stream, reduced=True):
	
		'''Writes the contents of the receiver to stream as a mutation list.
		
		Parameters
			filename: The stream to write to.
			
		Exceptions:
			Raises an IOError if the file can't be written to for any reason'''
	
		if self.isSinglePointMutationFormat():
			self._writeSinglePointMutationListToStream(stream)
		else:
			self._writeStandardMutationListToStream(stream, reduced=reduced)
		
	def isSinglePointMutationFormat(self):
	
		'''Returns True if the list is in SPM format.
		
		A list is in SPM format if every entry refers to a SPM.'''
		
		result = False
		#If each MutationSet instance ony contains one residue code then they 
		#are single point mutations
		mutationsPerSet = [len(set.residueCodes()) for set in self.mutantList()]
		if mutationsPerSet.count(1) == len(self.mutants):
			result = True
			
		return result	
		
	def isStandardMutationFormat(self):
	
		'''Returns True if the list is in full format.
		
		A list is in full format if any entry does not refer to a SPM'''
		
		return not self.isSinglePointMutationFormat()

	def format(self):

		'''Returns the format of the receivers list

		Returns:
			A string. Standard or SinglePointMutation'''

		if self.isStandardMutationFormat():
			return 'Standard'
		
		return 'SinglePointMutation'
		
	def validateOnPDB(self, pdbFile, errors={}):
	
		'''Validates the mutants defined by the list on a pdb file
		
		Returns NO if there are problems and the dictionary errors contains the following keys
		errorString - A string describing the validation problems.
		errorList - A list continaing the mutation codes that failed validation
		
		Raises Exceptions.FileFormatError if there is a problem with the pdb'''
		
		import Protool
		
		try:
			pdb = Protool.structureIO()
			pdb.readpdb(pdbFile)
		except Exception, data:	
			raise Exceptions.FileFormatError, 'Format of specified PDB file %s not valid.\nUnderlying error - %s' % (pdbFile, data)
		
		errorList = []
		invalid = 0
		reason = ""
		retval = True
		for set in self.mutantList():
			try: 
				set.codeString(pdb=pdb)
			except Exceptions.MutationCodeApplicationError, data:
				if invalid == 0:
					first = data.mutationCode
					retval = False
					reason = data.errorString
				invalid = invalid + 1
				errorList.append(data.mutationCode)
		
		if len(errorList) != 0:
			errors["errorList"] = errorList
			errors["errorString"] = "There were %d failures in total - the first failure was with %s.\n\n%s" % (len(errorList), first, reason)
			
		return retval	
										
class DataSet:

	'''Class representing a collection of csv files stored in a directory.

	Instances of this class can create, read and manipulate these directories and the data in them.
	This class is used to create the output data of the PEAT_SA command line tool.
	
	Attributes:
		name - The name of the data directory
		location - The path to the directory containing the data
		path - The full path to the data'''

	def _writeMatrixNames(self):
	
		'''Writes the names of the current matrices to the on-disk archive'''
		
		matrixNames = self.data.keys()
		
		if self.environment.isRoot():
			stream = open(os.path.join(self.path, 'matrixNames'), 'w')
			pickle.dump(matrixNames, stream)
			stream.close()
		
		return  matrixNames		

	def _removeErrorHandler(self):
	
		'''Handler for shutil.rmtree - Does nothing at the moment except raise an exception'''
		
		raise Exceptions.UnableToRemoveDataError, "Unable to remove .peatsa data directory"		

	def _loadMatrix(self, name):

		'''Unarchives the matrix called name from the data dir and returns it.
		
		Return None if no matrix called name exists or if the data could not be read.'''

		filename = os.path.join(self.path, name)
		try:
			matrix = Matrix.matrixFromCSVFile(filename)
		except IOError:
			print 'Expected matrix %s - missing ...' % filename
			matrix = None
		except TypeError:
			print 'Unable to read matrix %s - Possible file corruption. Ignoring ...' % filename
			matrix = None
			
		return matrix	

	def _readData(self):
	
		'''Reads the available data'''
		
		try:
			stream = open(os.path.join(self.path, 'matrixNames'), 'r+')
			matrixNames = pickle.load(stream)
			stream.close()
			for name in matrixNames:
				matrix = self._loadMatrix(name + '.csv')
				if matrix is not None:
					self.data[name] = matrix
		except IOError:
			#If the file doesn't exist we have old PEATSA data.
			#See which of peatMatrices are present and then write the 
			#missing matrixNames files with their names.
			for name in peatMatrices:
				matrix = self._loadMatrix(name + '.csv')
				if matrix is not None:
					self.data[name] = matrix
			
			#Root process output the names for next time
			matrixNames = self._writeMatrixNames()
		
		#Check if the expected contents are the same as the actual contents
		if matrixNames != self.data.keys():
			print 'Descrepency detected between expected and actual contents'
			print 'Updating to reflect current contents'
			self._writeMatrixNames()
		
	def __init__(self, name, location=None, overwrite=True):
	
		'''Initialise a new DataSet instance. The instance is created if it doesn't exist.
		
		Parameters:
			name - The name of the directory corresponding to the data set
			location - The location of the results directory.
			If None assumes current directory.'''
			
		self.environment = Environment.Environment()	
			
		if location is None:
			location = os.getcwd()
			
		if name is None:
			raise TypeError, 'Name must be a string'
			
		self.name = name
		self.location = location
		self.path = os.path.join(self.location, self.name)
		self.overwrite = overwrite
		self.data = {}
	
		#Check if the directory exists
		#If not create it
		if os.path.exists(self.path):
			if not os.path.isdir(self.path):	
				raise ValueError, 'Incorrect file type already exists at %s' % self.path
			else:
				self._readData()	
		else:
			#Make sure every process gets here before the directory is created.
			#Otherwise some might enter the above loop and attempt to read the matrixNames
			#pickle as its being written
			self.environment.wait()
			#If we have to create the dir append a ".peatsa' extension			
			#to the name if its not present
			components = os.path.splitext(self.name)	
			if components[1] == "":
				self.name = name + '.peatsa'
			elif components[1] is not '.peatsa':
				self.name = components[0] + '.peatsa'
			else:
				self.name = name
			
			self.path = os.path.join(self.location, self.name)
			
			#Create the directory - but don't do it until every process has reached this point
			self.environment.wait()
			if self.environment.isRoot():
				try: 
					os.mkdir(self.path)
				except BaseException, data:
					raise ValueError, 'Invalid location given (%s). Unable to create data directory. Reason %s' % (self.path, data)
					
		#Write out empty matrix names list
		#This distinguishes this directory from previous versions which don't have matrixNames file at all
		self._writeMatrixNames()			
				
	def __getattr__(self, name):
	
		'''Overrides normal __getattr__ to allow matrices to be accessed as properties'''
	
		#Avoid recursively trying to access self.data if its not defined.
		if name is 'data':
			raise AttributeError, name
			
		capName = name.capitalize()[0] + name[1:]	
		
		if capName in self.data:
			return self.data[capName]
		else:
			raise AttributeError, name
				
	def __str__(self):
	
		'''Returns a description string'''
	
		description = 'PEAT-SA data at %s. Contents :\n' % self.path
		contents = ""
		for name in self.data.keys():
			if self.data.has_key(name):
				contents = contents + '\t%s\n' % name	
		
		if contents == "":
			contents = 'None'
		
		return description + contents
	
	def addMatrix(self, matrix, name):
	
		'''Adds a matrix called name to the data dir.
		
		Note: Name is capitalized if its not already
		
		If data for name is already present, and overwriting is disabled,
		an Exceptions.DataAlreadyExistsError is raised.
		Otherwise the data is added to the self.data dictionary.'''
	
		if self.data.has_key(name):
			if not self.overwrite:
				raise Exceptions.DataAlreadyExistsError, 'A matrix called %s already present' % name
		
		name = name.capitalize()[0] + name[1:]
		
		if self.environment.isRoot():
			filename = os.path.join(self.path, name + '.csv')
			file = open(filename, 'w')
			file.write(matrix.csvRepresentation())
			file.close()
			
			#Update the matrixNames list with the new name
			matrixNames = self.data.keys()
			matrixNames.append(name)
			matrixNames.sort(lambda a,b: cmp(a.lower(), b.lower()))
			stream = open(os.path.join(self.path, 'matrixNames'), 'w+')
			pickle.dump(matrixNames, stream)
			stream.close()
		
		self.environment.wait()
		self.data[name] = matrix
		
	def matrixNames(self):
	
		'''Returns a list containing the names of all the matrices in the receiver.
		
		These can be used as attributes of the receiver to access the matrices'''
		
		names = self.data.keys()
		names.sort(lambda a, b: cmp(a.lower(), b.lower()))
		
		return names
	
	def allMatrices(self):
	
		'''Returns a dictionary whose keys are the stored matrix names and values are matrices'''
		
		return copy.copy(self.data)
		
	def delete(self):
	
		shutil.rmtree(self.path, ignore_errors=1, onerror=self._removeErrorHandler);
		
	def synchronize(self):
	
		'''Synchronizes the in memory contents of the data sets matrices with the disk versions'''
		
		if self.environment.isRoot():
			for name in self.data.keys():
				matrix = self.data[name]
				filename = os.path.join(self.path, name + '.csv')
				file = open(filename, 'w')
				file.write(matrix.csvRepresentation())
				file.close()
		
		self.environment.wait()


class MutantCollection:

	'''Represents a collection of mutants of a specified pdb file which are stored on disk.

	Instances of this class act like a list containing the paths to the mutant pdbs.
	
	MutantCollection objects can create the mutants if they do not exist already or add new mutants to the collection
	They are also parallelised - if running in a parralel environment the mutations will be created in parallel.
	
	Accessible ivars
		pdbFile - Location of the pdb file the mutants were modelled on'''

	def _createDirectories(self, pdbFile, ligandFiles):
	
		'''Does all the necessary directory creation for the instance'''
		
		#Check if the output directory already exits
		#If it does and self.overwrite is True delete it
		if os.path.isdir(self.location):
			if self.overwrite is True and self.environment.isRoot():
				self.environment.output('[MUTATE] Warning - Removing previously existing directory at %s' % self.location)
				shutil.rmtree(self.path(), ignore_errors=1, onerror=self._removeErrorHandler)
			else:
				raise  Exceptions.MutantCollectionError, 'Unable to create output directory for mutant collection - Directory already exists'
		
		#Create the collection directory and sub-directories
		try:
			os.mkdir(self.location)
			os.mkdir(os.path.join(self.location, 'Mutants'))
			os.mkdir(os.path.join(self.location, 'Ligands'))
		except BaseException, data:
			self.isTemporary = True
			raise  Exceptions.MutantCollectionError, 'Unable to create output directory for mutant collection\nReason - %s' % data
				
		#Copy original file to the collection directory
		try:
			basename = os.path.basename(pdbFile)
			destination = os.path.join(self.location, basename)
			shutil.copyfile(pdbFile, destination)
		except BaseException, data:
			#Make self temporary so on dealloc the directory will be deleted
			self.isTemporary = True
			raise Exceptions.MutantCollectionError, 'Unable to copy pdb file to collection location'
		
		if ligandFiles is not None:
			for file in ligandFiles:
				try:
					basename = os.path.basename(file)
					destination = os.path.join(self.location, 'Ligands')
					destination = os.path.join(destination, basename)
					shutil.copyfile(file, destination)
				except BaseException, data:
					self.isTemporary = True
					string = 'Unable to copy ligand file at %s to collection location' % file
					raise Exceptions.MutantCollectionError, string

				
	def _removeErrorHandler(self):
	
		'''Handler for shutil.rmtree - Does nothing at the moment except raise an exception'''
		
		raise Exceptions.MutantCollectionError, "Unable to remove mutant collection directory"		

	def __initNewCollection__(self,  pdbFile, ligandFiles=[], location=None, name=None, clean=True):
	
		'''Called on instantiation to create a new MutantCollection'''
		
		self.collection = []
		
		#Check pdbFile exists.
		#CheckFile will raise an exception if something is wrong
		try:
			pdbFile = Utilities.CheckFile(pdbFile)
		except IOError, data:
			raise Exceptions.MutantCollectionError, 'Problem with pdbfile - %s' % data
		
		#Check the pdb file has the correct format using protool.
		import Protool
		try:
			pdb = Protool.structureIO()
			pdb.readpdb(pdbFile)
		except Exception, data:	
			raise Exceptions.FileFormatError, 'Format of specified PDB file %s not valid.\nUnderlying error - %s' % (pdbFile, data)
		
		#Create the collection name
		if name is not None:
			name = "%s.mutants" % name
			self.location = os.path.join(location, name)
		else:
			name = os.path.basename(pdbFile)
			name = os.path.splitext(name)[0]
			name = "%s.mutants" % name
			self.location = os.path.join(location, name)
		
		#Create all the necessary directories and copy the required files there
		#Only the root process does the directory creation
		if self.environment.isRoot():
			self._createDirectories(pdbFile, ligandFiles)
		
		self.environment.wait()
		
		#Redirect any whatif output to a log file
		filename = os.path.join(self.location, 'Mutate.log')
		mutateLog = open(filename, 'w+')
		standardOut = sys.stdout
		sys.stdout = mutateLog		
		self.environment.output('[MUTATE] Mutant collection log at %s\n' % filename, stream=standardOut)
		standardOut.flush()

		#Set up ivars for the new locations of the pdb and ligand files
		self.pdbFile = os.path.join(self.location, os.path.basename(pdbFile))
		self.ligands = []
		if ligandFiles is not None:
			#Add the paths to the moved ligand files to self.ligands
			#Its necessary to do this here since only the root process
			#does the copying step.
			for file in ligandFiles:
				file = os.path.abspath(file)
				basename = os.path.basename(file)
				destination = os.path.join(self.location, 'Ligands')
				destination = os.path.join(destination, basename)
				self.ligands.append(destination)
		
				
		#Clean the pdb if requested 
		try:
			if clean is True and self.environment.isRoot():
				self.environment.output('[MUTATE] Cleaning pdb %s\n' % self.pdbFile, stream=standardOut)
				standardOut.flush()
				Utilities.CleanPDB(inputFile=self.pdbFile, outputFile=self.pdbFile)
			
			self.environment.wait()	
		except	BaseException, data:
			self.isTemporary = True
			raise Exceptions.MutantCollectionError, 'Failed to clean pdb file - Underlying error - %s' % data	
		
		mutateLog.close()
		sys.stdout = standardOut
		
		#Create a protool instance of the pdb
		self.pdb = Protool.structureIO()
		self.pdb.readpdb(self.pdbFile)	
					
	def __initFromLocation__(self, location, name):
	
		'''Called on instantiation when the location of a existing mutant collection is given'''
	
		self.location = os.path.join(location, name)

		#Find the name of the pdbfile - 
		contents = os.listdir(self.location)
		contents = [os.path.split(file)[1] for file in contents]
		pdbs = []
		for file in contents:
			if os.path.splitext(file)[1] == '.pdb':
				pdbs.append(file)
				
		if len(pdbs) != 1:
			raise Exceptions.MutantCollectionError, 'Cannot determine pdb file from collection %s' % self.location
		else:
			self.pdbFile = os.path.join(self.location, pdbs[0])
		
		#Create a protool instance of the pdb and check if it has a vaid format
		import Protool
		try:
			self.pdb = Protool.structureIO()
			self.pdb.readpdb(self.pdbFile)
		except Exception, data:	
			raise Exceptions.FileFormatError, 'Format of specified PDB file %s not valid.\nUnderlying error - %s' % (self.pdbFile, data)
					
		#Read in the mutation list
		stream = open(os.path.join(self.location, 'mutationList'), 'r')
		archiver = pickle.Unpickler(stream)
		self.mutationList = archiver.load()
		stream.close()
		
		#The default naming scheme is given by any element of the mutation list
		self.useReducedNaming = self.mutationList[0].reducedDefault
		
		#Create the list of mutant pdb names from the mutation list
		mutantDir = os.path.join(self.location, 'Mutants')
		self.collection = [os.path.join(mutantDir, element.defaultFilename(self.pdb)) for element in self.mutationList]
	
		#Read in the ligand files
		ligandDir = os.path.join(self.location, 'Ligands')
		ligandFiles = os.listdir(ligandDir)
		self.ligands = [os.path.join(ligandDir, file) for file in ligandFiles]
		
		#Read in the scores of the mutants if presnet
		filename = os.path.join(self.location, 'Scores.csv')
		if os.path.isfile(filename):
			self.mutationScores = Matrix.matrixFromCSVFile(filename)
		
	def __init__(self, location=None, name=None, mutationList=[],  pdbFile=None, ligandFiles=[], 
			maxOverlap=0.5, clean=True, useReducedNaming=True, temporary=False, overwrite=True):
	
		'''Initialises a new MutantCollection instance.
		
		The new instance can be initialised with a previous mutant collection that exists on disk, 
		or a new one can be created.
		
		In the former case the name and location of the file are required options.
		In the latter case a pdbfile is required.
		To create mutants in the latter case either pass a mutationList to __init__ or use the methods
		createMutant() or addMutant() at a later time
		
		Parameters:			
			location - Path to the directory where the collection will be stored or exists.
				If you are opening a previous collection this cannot be None.
				
				However it can be None if you are creating a new collection.
				In this case or if its not specified the current working directory is used.

			pdbFile - Path to the pdb file that will be mutated.
			
			ligandFiles - A list of ligands that will be added to the pdb file when the
				mutations are made.
			
			mutationList - A list MutationSet objects representing the mutants to create.
			
			 - If True mutants will be stored using the reduced mutation code.
				This code uses 1 letter amino-acid codes and no padding. 
				See Utilities.ConvertResidueCodeToReducedFormat() for more.
				Default is True.
				Note this can not be changed once set.
						
			clean - If True all ligands and water will be removed from the pdb.
				Default is True.
			
			temporary - If true the directory will be deleted on restart.
			
			overwrite - If true an existing collection instance at the same location is deleted'''
		
		self.environment = Environment.Environment()
		self.isTemporary = temporary
		self.overwrite = overwrite
		self.mutationList = []
		self.maxOverlap = maxOverlap
		self.useReducedNaming = useReducedNaming 
		self.mutationScores = None
		
		#Check the location of the output directory if one is supplied
		#CheckDirectory will raise an exception if something is wrong
		if location is not None:
			try:
				location = Utilities.CheckDirectory(location)
			except IOError, data:
				raise Exceptions.MutantCollectionError, 'Problem with mutant collection location - %s' % data
		else:
			location = os.getcwd()

		#If no pdbfile is specified then the location and name of a previous collection must be supplied	
		#Otherwise we are creating a new collection
		if pdbFile is None and location is not None and name is not None:
			self.isNew = False
			self.__initFromLocation__(location, name)		
		elif pdbFile is not None:
			self.isNew = True
			self.__initNewCollection__(pdbFile=pdbFile, 
					ligandFiles=ligandFiles, 
					location=location, 
					name=name, 
					clean=clean)
		else:
			raise Exceptions.MutantCollectionError, "Neither a PDB file or the location of previous collection were provided"

		#Create the requested mutants if this is a new instance
		if self.isNew and len(mutationList) is not 0:
			self.createMutants(mutationList=mutationList)
		
	def __del__(self):
	
		'''Deletes the directory containing the created pdbs if requested on instantiation'''
	
		if self.environment.isRoot() and self.isTemporary:
			shutil.rmtree(self.path(), ignore_errors=1, onerror=self._removeErrorHandler)

	def __cmp__(self, object):
	
		'''Compares two MutantCollection instances.
		
		Return
			Returns 0 if the data in the matrices ns identical and they both have identical headers.
			Otherwise returns -1'''
			
		if isinstance(object, MutantCollection):
			#Cmp returns 0 if two objects are equal
			return cmp(self.collection, object)
		else:
			return -1
			
	def __getitem__(self, index):
	
		'''collection[i] returns the i'th mutant file name'''
	
		return self.collection(index)
		
	def __str__(self):
	
		'''Returns an information string'''
		
		string = "Mutant collection with %d mutants stored at %s" % (len(self.collection), self.path())
		return string
		
	def _updateMutationList(self):
	
		'''Private method - Updates the on-disk pickle of the mutationList ivar'''
		
		if self.environment.isRoot():
			stream = open(os.path.join(self.location, 'mutationList'), 'w+')
			archiver = pickle.Pickler(stream)
			archiver.dump(self.mutationList)
			stream.close()
		
		self.environment.wait()	
		
	def _createMutant(self, element, location, successes, failures, qualities, standardOut):
	
		'''Private method - creates the mutant represented by mutationCode.
		
		Note: Other objects should never call this method
		
		Params:
			element - A MutationSet instance representing the mutant
			location - The mutant output directory
			successes - A list: If successful the filename of the created mutant is appended to this list
			failures - A list: If modelling failed the mutants name is appended to this list
			qualities - A list: The modelling quality is appended to this list
			standardOut - Where to write output'''	
		
		import Protool.mutate
		
		mutationCodes = element.mutationCodes(self.pdb, reduced=False)
		name = element.codeString(pdb=self.pdb, reduced=False)		
		self.environment.output('[MUTATE] (%d) Creating mutant %s\n' % (self.environment.rank(), name), 
			stream=standardOut, rootOnly=False)
			
		#Skip mutations that can't be modelled
		totalBump = len(mutationCodes)*self.maxOverlap
		mutant, score = Protool.mutate.Model_Mutations(self.pdbFile, 
					self.ligands, 
					mutationCodes, 
					max_overlap=self.maxOverlap, 
					max_totalbump=totalBump,
					return_score=True)
		self.environment.output('[MUTATE] (%d) Done - Score %lf. Total-bump %lf\n' % (self.environment.rank(), score, totalBump), 
			stream=standardOut, rootOnly=False)
		
		if mutant is False:
			qualities.append([name, score, 'Failed'])			
			self.environment.output('[MUTATE] Unable to create mutant %s - Could not be modelled.\n' %
					name,
					stream=standardOut,
					rootOnly=False)	
			failures.append(name)
		else:
			qualities.append([name, score, 'Created'])			
			mutant = mutant.PI
			#Remove the ligand
			mutant.remove_atoms_with_tag('LIGAND')
			name = element.defaultFilename(pdb=self.pdb)
			name = os.path.join(location, name)
			mutant.writepdb(name)
			#Add hydrogens to new residue
			self.environment.output('[MUTATE] (%d) Starting clean\n' % self.environment.rank(), 
				stream=standardOut, rootOnly=False)
			self.environment.output('[MUTATE] Cleaning mutant %s\n' % name, 
				stream=standardOut, rootOnly=False)
			Utilities.CleanPDB(inputFile=name, outputFile=name, 
					removeWater=False, 
					removeLigand=False, 
					removeAltLoc=False,
					addHydrogens=True,
					correct=True)
			successes.append(name)
			self.environment.output('[MUTATE] (%d) Done\n' % self.environment.rank(), 
				stream=standardOut, rootOnly=False)


	def createMutants(self, mutationList=[]):
	
		'''Creates all the mutants represented by the MutationSet instances in mutationList and adds them to the collection.
		
		If modelling a mutant fails it is skipped - see mutationList() for more.
		Note: The default naming of the MutationSets is overridden by the naming scheme of the receiver which is
		set on initialisation'''
	
		print 'Creating mutants'
		#Redirect any whatif output to a log file
		filename = os.path.join(self.location, 'Mutate.log')
		mutateLog = open(filename, 'a+')
		standardOut = sys.stdout
		standardError = sys.stderr
		sys.stdout = mutateLog		
		standardOut.flush()
		
		if len(mutationList) == 0:
			raise Exceptions.MutantCollectionError, 'Cannot create mutants - No mutations provided'
		
		filenames = []
		failedMutants = []
		qualities = []
		
		location = os.path.join(self.location, 'Mutants')
		listFragment = self.environment.splitArray(mutationList)
	
		for element in listFragment:
			#Set the element naming to the default
			element.reducedDefault = self.useReducedNaming
			self.environment.output('[MUTATE] (%d) Next element %s\n' % (self.environment.rank(), element), 
				stream=standardOut, rootOnly=False)
			try:
				self._createMutant(element, location, filenames, failedMutants, qualities, standardOut)
			except KeyboardInterrupt:
				raise
			except Exception, data:
				name = element.codeString(pdb=self.pdb, reduced=False)		
				self.environment.output('[MUTATE] Unable to create mutant %s. Encountered exception - Reason %s\n' % (name, data), 
					stream=standardOut,
					rootOnly=False)	
				failedMutants.append(name)	
				self.environment.output('[MUTATE] Traceback %s' % Utilities.GetTraceback(), 
					stream=standardOut,
					rootOnly=False)
			except BaseException, data:
				self.environment.output('Caught base exception %s %s' % (BaseException, data), 
					stream=standardOut, rootOnly=False)
				standardOut.flush()
				raise
			except:
				self.environment.output('Caught unknown exception (%d)' % self.envrionment.rank(), 
					stream=standardOut, rootOnly=False)
				raise	
				
			self.environment.output('[MUTATE] Process rank %d has processed %d mutants: %d failures, %d left\n' % (self.environment.rank(), len(filenames), len(failedMutants), len(listFragment) - len(filenames) - len(failedMutants)), stream=standardOut, rootOnly=False) 
			standardOut.flush()
			
		#This is required in order to catch if any the processes exited the 
		#above loop due to an exception	
		self.environment.wait()			
		
		#Combine arrays
		filenames = self.environment.combineArray(filenames)
		failedMutants = self.environment.combineArray(failedMutants)	
		qualities = self.environment.combineArray(qualities)
		
		noFailures = len(failedMutants)
		noSuccesses = len(filenames)
		if noFailures is not 0:
			self.environment.output('[MUTATE] Failed to model %d of %d mutations\n' % (noFailures, noSuccesses + noFailures),
						stream=standardError)	
			self.environment.output('[MUTATE] Failures = %s\n' % " ".join(failedMutants), stream=standardError)			
			self.environment.output('[MUTATE] Check Mutation.log for details\n', stream=standardError)
		
		mutateLog.close()
		sys.stdout = standardOut	
		self.collection.extend(filenames)
		
		#Update the mutationList ivar
		#Filter out the failed mutants
		mutationList = filter(lambda x:  x.codeString(self.pdb, reduced=False) not in failedMutants, mutationList)
		
		self.mutationList.extend(copy.deepcopy(mutationList))
		self._updateMutationList()
		self.environment.output('[MUTATE] Done: %s' % self.__str__()) 
		
		#Update qualities matrix
		if self.mutationScores is None:
			self.mutationScores = Matrix.PEATSAMatrix(rows=qualities,
						 headers=['Mutations', 'Score', 'Result'])
		else:				 
			for row in qualities:
				self.mutationScores.addRow(row)
		
		f = open(os.path.join(self.location, 'Scores.csv'), 'w+')		
		f.write(self.mutationScores.csvRepresentation())
		f.close()	
		
		
	def addMutant(self, mutant, mutationSet):
	
		'''Adds a mutant (represented by a Protool instance) to the collection. 
		
		Note: Its is assumed that mutant is an mutant of the pdb provided on initialisation
		
		Parameters:
			mutant - A protool instance representing the mutant pdb.
			chain - The chain the mutated residue is in
			residueIndex - The index of the residue that was mutated
			residueName - The three letter code of the residue
			mutation - A list of three letter code of the mutant (what the residue was mutated to)'''
			
		mutant.remove_atoms_with_tag('LIGAND')
		self.mutationList.append(copy.deepcopy(mutationSet))	
		self._updateMutationList()	
		
		#Get the mutationCode used by pKa modules to refer to this entry
		mutationSet.reducedDefault = self.useReducedNaming
		filename = mutationSet.defaultFilename(pdb=self.pdb)
		
		#Write out the mutant pdb to the mutants directory
		mutantDirectory = os.path.join(self.location, 'Mutants')
		filename = os.path.join(mutantDirectory, filename)
		mutant.writepdb(filename)	
		
		#Add the filename to the collection ivar.
		self.collection.append(filename)
		
	def path(self):
	
		'''Returns the path to the collection'''
		
		return self.location
	
	def mutantFiles(self):
	
		'''Returns a list containing the paths to all the mutant pdbs
		
		This the contents of the Mutants subdirectory'''
	
		return copy.copy(self.collection)
		
	def ligandFiles(self):
	
		'''Returns a list containing the paths to ligand files used
		
		This is the contents of the Ligand subdirectory'''
	
		return copy.copy(self.ligands)
		
	def mutant(self, mutationSet):
	
		'''Returns a Protool instance for the given mutation'''
		
		import Protool
		filename = self.fileForMutation(mutationSet)
		pdb = Protool.structureIO()
		pdb.readpdb(filename)
		
		return pdb
		
	def fileForMutant(self, mutationSet):
	
		'''Returns the path to the pdb file for mutationSet.
		
		If no pdb corresponding to muationSet exists this method returns None'''
	
		try:
			self.mutationList.index(mutationSet)
			filename = mutationSet.filename(reduced=self.useReducedNaming, pdb=self.pdb)
			mutantDirectory = os.path.join(self.location, 'Mutants')
			filename = os.path.join(mutantDirectory, filename)
		except ValueError:
			filename = None
		
		return filename		
		
	def mutations(self):
	
		'''Returns a list of the mutations in the collection.
		
		Each element in the list is a Data.MutationSet instance.
		
		Note that any mutatants that the receiver failed to model will not'
		be in this list.'''
	
		return copy.deepcopy(self.mutationList)	
		
	def mutationListFile(self, filename):
	
		'''Returns a MutationListFile instance detailing the mutants in the receiver
		
		MutationListFile objects represent files containing a list of mutation codes.
		
		Params:
			filename - The name of the file associated with the list
			
		Note: This method does not create the file.
		To do so writeToFile() must be called on the returned object'''
		
		mutations = self.mutations()
		list = MutationListFile(filename=filename, create=True)
		for set in mutations:
			list.addMutant(set, autoUpdate=False)
		
		return list
		
def mutationSetFromFilename(filename):

	'''Factory function for MutationSet class. Creates a mutation set from filename.
	
	Filename should have been created by a call to one of the MutationSet filename methods.'''

	mutationSet = MutationSet()
	codes = os.path.splitext(filename)[0]
	codes = codes.split(MutationSet.mutationSeparator)	
	for code in codes:
		data = Utilities.ParseResidueCode(code)		
		#If the residue code was in 
		if len(data) == 4:
			mutationSet.addMutation(chain=data[0], residueIndex=data[1], mutation=data[3])
		elif len(data) == 3:
			mutationSet.addMutation(chain=data[0], residueIndex=data[1], mutation=data[2])
			
	return mutationSet
	
def mutationSetFromSequences(initialSequence, targetSequence, offset=0, chain='A', name=None):

	'''Returns a MutationSet representing the mutations necessary to change initialSequence to targetSequence
	
	Parameters:
		initialSequence - "WildType"
		targetSequence - "Mutant"
		offset - The difference between the index of the first aa in the initial sequence and 1
			e.g. if the first amino-acid in initialSequence string is actually the 5th
			in the real sequence then the offset is 4.
		chain - The chain the sequence corresponds to - defaults to A
		name - A name to associate with the created mutation set	  
	
	Note: If this requires insertions or deletions this method returns None'''
	
	if len(initialSequence) != len(targetSequence):
		raise Exceptions.ProteinDesignToolException, 'Sequences are not the same length'

	if initialSequence.find('-') != -1:
		return None
	elif targetSequence.find('-') != -1:
		return None	

	set = MutationSet(name=name)
	for count in range(len(initialSequence)):
		initialAA=initialSequence[count]
		targetAA=targetSequence[count]
		if initialAA != targetAA:
			set.addMutation(chain, count + offset + 1, targetAA) 
	
	return set
	
class MutationSet:

	'''Instances of this class represent a set of mutations
	
	Attributes:
		mutations - A list of mutation codes the receiver represents'''

	#The separator used in pdb filenames
	mutationSeparator = '+'

	def _extendedMutationCodes(self, pdb):
		
		'''Private method. Returns the mutations in extended format
		
		Parameters:
			pdb - A protool instance representing the pdb that the receiver is a mutant of.
			
		Errors: 
			Raises an MutationCodeApplicationError if 
			- the residue index specified by a mutation code is not in the pdb.
			- the residue chain specified by a mutation code is not in the pdb'''
		
		mutations = []
		for code in self.mutationCodeList:
			chain, residueIndex, mutation = Utilities.ParseMutationCode(code)
			residueCode = Utilities.CreateResidueCode(chain, residueIndex)
			try:
				extendedResidueCode = Utilities.CreateResidueCode(chain=chain,
								number=residueIndex,
								name=pdb.resname(residueCode))
			except Exception:
				#Check if the error was because of chain or index
				if pdb.chains.keys().count(chain) == 0:
					message = "There is no chain '%s' in the pdb file" % chain
				else: 
					#Its a residue index problem
					message = "There is no residue %d in chain '%s' in the pdb file" % (residueIndex, chain)
			
				raise Exceptions.MutationCodeApplicationError(mutationCode=code, message=message)
				
			extendedMutationCode = extendedResidueCode + ':' + mutation
			mutations.append(extendedMutationCode)
			
		return mutations

	def __init__(self, code=None, name=None, reducedDefault=True):
	
		'''Initialises a new instance
		
		When initialised you can set type of name returned by defaultFilename().
		Parameters:
			reducedDefault: If true the defaultFilename() method uses the reduced naming scheme.
			Otherwise it uses the normal naming scheming. 
			code: An optional codeString defining a set of mutation
			
		Errors:
			If the supplied codeString is invalid an Exception.MutationCodeFormatError is raised.'''
	
		self.mutationCodeList = []
		self.reducedDefault = reducedDefault
		if name is None:
			name = 'None'
		
		self.name = name
		
		if code is not None:
			if code.find("+") != -1:
				components = code.split("+")
			elif code.find(MutationSet.mutationSeparator) != -1:
				#Could be taken from a filename	
				components = code.split(MutationSet.mutationSeparator)
			elif code.find('_') != -1:
				#Old filename separator
				components = code.split('_')	
			else:
				#Assume its a single mutation code
				components = [code]
				
			for component in components:
				self.addMutationCode(component)

	def __cmp__(self, object):
	
		'''Returns true if object describes the same set of mutations as the receiver'''
	
		if hasattr(object, 'mutationCodeList'):
			return cmp(self.mutationCodeList, object.mutationCodeList)
		else:
			return -1

	def __str__(self):
	
		'''Prints out a string description of the receiver'''
	
		if len(self.mutationCodeList) == 0:
			return 'WildType'
	
		string = 'Mutations:'
		for mutation in self.mutationCodeList:
			string = string + ' %s,' % mutation
		
		return string
		
	def __getattr__(self, name):
	
		#Update for pre reducedDefault instances
		if name == 'reducedDefault':
			self.reducedDefault = False
			return self.reducedDefault
		else:
			raise AttributeError, name

	def addMutation(self, chain, residueIndex, mutation):
	
		'''Adds the change of residue at position residueIndex in chain to mutation
		
		Errors:
			Raises an Exceptions.MutationCodeFormatError if there is a problem with any of the parameters'''
		
		if len(mutation) == 1:
			try:
				mutation = Utilities.InvertedCodeDictionary()[mutation]
			except KeyError:
				raise Exceptions.MutationCodeFormatError, "Unknown single letter residue code %s" % mutation	
	
		if Utilities.allResidues.count(mutation.upper()) is not 1:
			raise Exceptions.MutationCodeFormatError, "Unknown residue code %s" % mutation
	
		if len(chain) > 1:
			raise Exceptions.MutationCodeFormatError, "Invalid value for chain id (%s) - must be a single letter, number of blank" % chain
	
		residueCode = Utilities.CreateResidueCode(chain, int(residueIndex))
		self.mutationCodeList.append(residueCode + ':' + mutation)
		
	def addMutationCode(self, mutationCode):
	
		'''Adds the mutation represented by mutationCode to the receiver.
		
		The mutation can be in any format
		
		Errors:
			Raises an Exceptions.MutationCodeFormatError if there is a problem with the code'''

		extended = Utilities.IsExtendedMutationFormat(mutationCode)
	
		try:
			if Utilities.IsReducedCode(mutationCode):
				mutationCode = Utilities.ConvertMutationCodeFromReducedFormat(mutationCode, extended)
		
			data = Utilities.ParseMutationCode(mutationCode)
			self.addMutation(chain=data[0], residueIndex=data[1], mutation=data[-1])
		except Exception, data:
			raise Exceptions.MutationCodeFormatError, "Unable to parse mutation code - %s." % mutationCode
				
			
	def mutationCodes(self, pdb=None, reduced=True):
	
		'''Returns a list of the mutation codes in the receiver in standard format
		
		Parameters:
			pdb - An optional protool instance.
			If this is provided the mutations are returned in extended format.
			reduced - If True the codes are returned in reduced format.
			
		Exceptions:
			Raises an AttributeError if pdb does not respond to resname() and is not None.
			Raises an MutationCodeApplicationError if the mutation code cannot be 
			applied to the supplied pdb to generate an extended code.'''
	
		extended = False
		if pdb is not None:
			extended = True
			#This method will raise MutationCodeApplicationError if the extended code can't be created
			mutationCodes = self._extendedMutationCodes(pdb)
		else:
			mutationCodes = copy.deepcopy(self.mutationCodeList)
		
		if reduced is True:
			mutationCodes = [Utilities.ConvertMutationCodeToReducedFormat(code, extended=extended) for code in mutationCodes]	
					
		return mutationCodes
		
	def residueCodes(self, pdb=None, reduced=True):
	
		'''Returns the residue codes of all the residues mutated in the receiver in standard format
		
		Parameters:
			pdb - An optional protool instance.
			If this is provided the filename uses extended residue code format.
			reduced - If True the codes are returned in reduced format.

			
		Exceptions:
			Raises an AttributeError if pdb does not respond to resname() and is not None'''
		
		residueCodes = []
		extended = False
		for mutation in self.mutationCodeList:
			chain, residueIndex, mutation = Utilities.ParseMutationCode(mutation)
			residueCode = Utilities.CreateResidueCode(chain, residueIndex)
			if pdb is not None:
				extended=True
				residueCode = Utilities.CreateResidueCode(chain=chain,
								number=residueIndex,
								name=pdb.resname(residueCode))
			if reduced is True:
				residueCode = Utilities.ConvertResidueCodeToReducedFormat(residueCode, extended=extended)					

			residueCodes.append(residueCode)
			
		return residueCodes	

	def chainCodes(self, pdb=None, reduced=False):

		'''Returns individual code strings describing the mutations in each chain of the receiver'''

		chains = {}

		for mutationCode in self.mutationCodeList:
			chain, residueIndex, mutation = Utilities.ParseMutationCode(mutationCode)
			if chains.has_key(chain):
				chains[chain].addMutationCode(mutationCode)
			else:
				chains[chain] = MutationSet(code=mutationCode)
		
		return [set.codeString(pdb, reduced=reduced) for set in chains.values()]	
	
	def reducedMutationCodes(self, pdb=None):

		'''Returns a list of the mutation codes in the receiver in reduced format
		
		Deprecated: Use mutationCodes passing True for reduced
		
		Parameters:
			pdb - An optional protool instance.
			If this is provided the mutations are returned in extended format.
			
		Exceptions:
			Raises an AttributeError if pdb does not respond to resname() and is not None'''
		
		return self.mutationCodes(pdb=pdb, reduced=True)
		
	def codeString(self, pdb, reduced=True):
	
		'''Generates a string by concatentate all extended residue code with +.
		
		If reduced is true, the reduced format is used, otherwise the default format is used'''
		
		mutations = self.mutationCodes(pdb=pdb, reduced=reduced)
	
		return "+".join(mutations)
	
	def filename(self, reduced=True, pdb=None):	
	
		'''Returns a filename for a pdb containin the mutant represented by the receiver
		
		This is formed by concatenating all the mutation codes with _
		
		Parameters:
			pdb - An optional protool instance.
			If this is provided the filename uses extended residue code format
			
		Exceptions:
			Raises an AttributeError if pdb does not respond to resname() and is not None'''	
	
		mutations = self.mutationCodes(reduced=reduced, pdb=pdb)
		return MutationSet.mutationSeparator.join(mutations) + '.pdb'	
			
	def defaultFilename(self, pdb=None):
		
		return self.filename(reduced=self.reducedDefault, pdb=pdb)
		
		
	def applyToSequence(self, wildTypeSequence, id='A', pdb=None):

		'''Applies the mutations defined by the receiver to wildTypeSequence
		
		Note: wildTypeSequence must correspond to a single chain
		Also sequence are assumed to have a 1 offset.
		Thus the mutation A1A is applied to element 0 of wildTypeSequence
		
		Parameters:
			wildTypeSequence: A string of one letter amino-acid codes
			
			chain: The chain the sequence corresponds to.
				Defauts to A.
			
			pdb: Optional. A protool instance of a pdb containing the chain
			corresponding to wild-type sequence. If this is not None the wild-type
			residue from the pdb is compared to that in the supplied sequence
			and an exception is raised if the two are not the same.
			
		Returns: 
			The mutant sequence'''
			
		seq = list(wildTypeSequence)	
		for code in self.mutationCodes(pdb=pdb):
			components = Utilities.ParseMutationCode(code)
			if components[0] == id:
				if pdb is not None:
					chain, index, wt, mut = components
					#Check the wt res defined by the sequence is the
					#same as that defined by the pdb
					if seq[index - 1] == Utilities.aminoAcidCodes[wt]:
						seq[index - 1] =  Utilities.aminoAcidCodes[mut]
					else:
						message = 'Residue in pdb does not correspond to reside in sequence (%s)' % seq[index - 1]
						raise Exceptions.MutationCodeApplicationError(mutationCode=code, message=message)
				else:
					chain, index, mut = components
					seq[index - 1] = Utilities.aminoAcidDict[mut]
				
		return "".join(seq)		
			
			
		
		
				
