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
'''Contains classes and functions for processing a submission from the WebApp's main web page'''
import subprocess, os, urlparse, urllib2, string, time, StringIO, sys
import PEATSA.Core as Core
import UtilityFunctions, Exceptions, Data

def GetPDB(pdbCode, dict={}):

	'''Fetches a pdb from the pdb. Returns a stream to the pdb contents.
	
	If dict is provided it contains the key 'stream' on success whose value is the stream.
	On fail it contains two keys 'error' and 'description'''
	
	url = "http://www.rcsb.org//pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=" + pdbCode
	 
	try:
		stream = urllib2.urlopen(url, None, 10)		
		#Check for an error
		info = stream.info()
		status = info.status
		if status is not "":
			stream = None
			dict['error'] = 'Error status %s' % str(status)
		elif not info.has_key('content-disposition'):
			stream = None
			dict['error'] = 'Request for %s returned nothing' % pdbCode
		else:
			lines = stream.readlines()
			string = "".join(lines)
			stream = StringIO.StringIO(string)
			dict['stream'] = stream
	except urllib2.HTTPError, data:
		dict['error'] = 'Unable to retrive structure for pdb code %s from the Protein Data Bank. ' % pdbCode
		dict['description'] = 'Reason: %s' % data
		stream = None
	except urllib2.URLError, data:
		dict['error'] = 'Unable to retrive structure for pdb code %s from the Protein Data Bank. ' % pdbCode
		dict['description'] = 'Reason: %s' % data.reason
		stream = None
	except Exception, data:
		dict['error'] = 'Unable to retrive structure for pdb code %s from the Protein Data Bank. ' % pdbCode
		dict['description'] = 'Reason: Encountered unexpected exception %s' % data
		stream = None
	
	return stream
	

def CreateMutationString(mutationData, webStorageDirectory, jobInputDirectory, job):

	'''Creates the mutation part of the command line arguments based on the mutationData dictionary'''

	mutationType = mutationData['Type']
	mutationString = "--%s=%s"

	if mutationType == 'mutationList':
		#Write out the file to the web storage directory
		filename = os.path.join(webStorageDirectory, 'mutationList')
		stream = open(filename, 'w+')
		stream.write(mutationData['Data'].read())
		
		try:
			#Fixme - catch mutation list exception
			stream.seek(0)
			mutationList = Core.Data.mutationListFileFromStream(stream)
			job.setMutationListFile(mutationList)
			
			stream.close()
		except:	
			pass
		
		#Create the string
		#Note we must give the path to the mutationFile as it will 
		#be when the job is run
		filename = os.path.join(jobInputDirectory, 'mutationList')
		mutationString = mutationString % (mutationType, filename)
	elif mutationType == 'mutation':
		mutationString = mutationString % (mutationType, mutationData['Data'])
		job.setMutation(mutationData['Data'])

	return mutationString	
	
def CreateCalculationString(calculations, ligandFile):
	
	'''Creates the calculation parts of the command line arguments bassed on the calculations list'''
	
	calculationString = ""
	if calculations.count('scan') != 0:
		calculationString = calculationString + "--scan "
	
	if calculations.count('stability') != 0:
		calculationString = calculationString + " --stability"
		
	if ligandFile is not None and calculations.count('binding') != 0:
		calculationString = calculationString + " --ligand=%s" % ligandFile
		
	return calculationString																					

def ConstructErrorURL(domain="PDT.UnknownDomain", 
		description="Unexpected Error", 
		detailedDescription="Unknown error while processing job",
		recoverySuggestion="Please contact the developers with details of what you were doing"):

	'''Constructs an URL for the WebApps error page
	
	Automatically determines the server name and port from os.envrion
	
	Parameters
		domain - The domain of the error - Indicates what part of the program failed
		description - A short string which very briefly describes what happened
		detailedDescription - A longer string elaborating on the error
		recoverySuggestion - A string explaining what to do, if known
		
	Return 
		A string containing an URL'''
		
	#Add error data as GET data - Possibly a better way?
	#e.g. generate html here?
	dataDict = {"domain":domain,
			"description":description, 
			"detailedDescription":detailedDescription, 
			"recoverySuggestion":recoverySuggestion}
	data = ["%s=%s" % (key, dataDict[key]) for key in dataDict.keys()]
	query = "&".join(data)
	
	location = os.environ['SERVER_NAME'] + ":" + os.environ['SERVER_PORT']
	components = ("http", location, "PEATSA/Pages/Error.php", "", query, "")	
	return urlparse.urlunparse(components)	
	
def ProcessPDBFilename(filename):

	'''Processes the input string to be of the form "PDBID.pdb"
	
	This is done in the following way:
	- filename is split into base + extension
	- Trailing whitespace and underlines are removed
	- All punctuation (except for underlines) is removed
	- All spaces are replaces with underlines
	- The base is lower-cased and appended with .pdb
	
	If the file does not have the '.pdb' extension (in any mixture of cases)
	an empty string is returned'''
		
	filename = os.path.basename(filename)	
	extension = os.path.splitext(filename)[1]
	if extension.lower() != ".pdb":
		return ""
		
	pdbId = os.path.splitext(filename)[0]
	
	#Strip stuff
	pdbId = pdbId.strip()
	pdbId = pdbId.strip("_")
	
	#Replace underlines with spaces
	#This is so these aren't removed in the next step
	pdbId = pdbId.replace("_", " ")
	
	#Remove all punctuation characters
	for character in string.punctuation:
		pdbId = pdbId.replace(character, "")
	
	#Put the underlines back - this also replaces any
	#preexisting spaces with underlines
	pdbId = pdbId.replace(" ", "_")
	
	pdbId = pdbId.lower()
	return pdbId + '.pdb'

class JobConstructor:

	'''Setup a WebApp job using data submitted by the user.
	
	This includes creating output directories, placing files in the
	correct location, and forming the command line to be executed.
	
	Many of the options controlling the JobConstructor instances behaviour
	are given in the WebApps configuration file. 
	
	Principal attributes
		- options - The web app options, read from the configuration file
		- job - A WebApp.Data.Job instance representing the job
		- connection - A connection to the webapp db
		- runString - The command line for the job'''
	
	def __init__(self, formData, construct=False):
	
		'''Constructor
		
		Parameters:
			formData - A FormData instance
			construct - If True this method calls construct immediately,
				otherwise construct must be called at a later stage.
				This allows parameters to be modified.'''
	
		self.formData = formData
		self.runString = None
		self.errorData = None
		self.job = None
		self.jobManager = None
	
		try:
			self.options = UtilityFunctions.DefaultWebAppConfiguration()
			self.connection = UtilityFunctions.ConnectionFromConfiguration(self.options)
			self.webServerDirectory = self.options.get('WEB APPLICATION', 'webServerDirectory')
			self.uploadLimit = int(self.options.get('WEB APPLICATION', 'uploadLimit'))
		except Core.Exceptions.EnvironmentError, data:
			self._setErrorData(data)
			return	
	
		#Check if the uploaded files exceed the uploadLimit
		if self._checkUploadedFiles():
			#Connect to the db
			self.jobManager = Data.JobManager(self.connection)
			if construct is True:
				self.construct()
		
	def __del__(self):
	
		self.connection.close()	
		
	def __str__(self):
	
		if self.runString != None:
			return "JobConstructor - Job not created yet"
				
		else:
			return "JobConstructor - Run string:\n\t%s" % self.runString
		
	def _checkUploadedFiles(self):
	
		'''Checks all the uploaded files to see if they are within limits'''
	
		#If a PKA code wasn't provided a 
		#pdb file must be present - check the file size
		if not self.formData.isPKACodePresent():
			content = self.formData.pdbFileStream().read()
			if len(content) > self.uploadLimit:
				self.errorData = {"domain":"PDT.SubmissionDomain",
						"description":"Error with pdb.",
						"detailedDescription":"Filesize exceeds size limit (%.2lfMB)" % (self.uploadLimit/(1024.0*1024.0)),
						"recoverySuggestion":"Unfortunately we only can accept files under this limit."}			
			
			if len(content) == 0:
				self.errorData = {"domain":"PDT.SubmissionDomain",
						"description":"Error with pdb.",
						"detailedDescription":"No data in file",
						"recoverySuggestion":"Check that the correct file was provided."}
									
		if self.errorData is not None and self.formData.isLigandPresent():
			content = self.formData.ligandFileStream().read()
			if len(content) > self.uploadLimit:
				self.errorData = {"domain":"PDT.SubmissionDomain",
						"description":"Error with uploaded ligand.",
						"detailedDescription":"Filesize exceeds upload limit (%.2lfMB)" % (self.uploadLimit/(1024.0*1024.0)),
						"recoverySuggestion":"Unfortunately we only can accept files under this limit."}
						
		if self.errorData is not None and self.formData.isMutationList():
			content = self.formData. mutationListFileStream().read()
			if len(content) > self.uploadLimit:
				self.errorData = {"domain":"PDT.SubmissionDomain",
						"description":"Error with uploaded ligand.",
						"detailedDescription":"Filesize exceeds upload limit (%.2lfMB)" % (self.uploadLimit/(1024.0*1024.0)),
						"recoverySuggestion":"Unfortunately we only can accept files under this limit."}				
						
		if self.errorData is not None:
			return False
		else:
			return True
			
	def _setErrorData(self, data):
	
		'''Convenience method for creating errorData due to a configuration error.
		
		The exact reason for the error is provided by data.'''
	
		self.errorData = {"domain":"PDT.ConfigurationDomain",
				"description":"Error initialising job submission environment.",
				"detailedDescription":data,
				"recoverySuggestion":"This is a bug - please contact the developers."}
				
	def _writeStringToFile(self, string, filename):
	
		'''Convenience method for writing to a file'''
		
		stream = open(filename, "w+")
		stream.write(string)
		stream.close()
	
	def construct(self):
	
		'''Performs all necessary steps for setting up a WebApp job based on the data submitted by a user via the WebApp main page.
		
		This basically involves three steps.
		
			- Creation of a entry for the Job in the WebApp database
			- Outputing the job files uploaded by the user to correct locations
			- Assembing the command line string that will be exectuted when runBackend() is called.
			
		Check the result of errorData()	to see if there were any problems with construction'''
				
		#For each job run a instance of the Job class is created using the JobManager object.
		#This creates the necessary entries in the MySQL database.
		#The Job instance contains info on the job and allows the data stored for the job to be modified.
		self.job = self.jobManager.createJob(self.formData.pdbId(), self.formData.calculations())
		
		#Create the input/output directory names
		try: 
			jobOutputDirectory = self.options.get('WEB APPLICATION', 'jobOutputDirectory')
			jobInputDirectory = self.options.get('WEB APPLICATION', 'jobInputDirectory')
			pKaDataDirectory = self.options.get('WEB APPLICATION', 'pKaDataDirectory')
			backendExecutable = self.options.get('WEB APPLICATION', 'launchScript')
		except 	Core.Exceptions.EnvironmentError, data:
			self._setErrorData(data)
			return
			
		#Get the various IO directories that will be used by the job
		#If the job is to be run on the local host then jobInputDirectory must
		#be the same as webServerDirectory
		webStorageDir = os.path.join(self.webServerDirectory, self.job.identification)
		outputDir = os.path.join(jobOutputDirectory, self.job.identification + '_Out')
		
		if self.formData.isPKACodePresent():
			workingDir = os.path.join(pKaDataDirectory, self.formData.pKaCode())
		else:
			workingDir = os.path.join(jobOutputDirectory, self.job.identification + '_Work')
			
		inputDir = os.path.join(jobInputDirectory, self.job.identification + '_Out')
		os.mkdir(webStorageDir)	
	
		#If this is not a delta pKa calculation we have to write
		#Otherwise write the uploaded/download one to a file
		if not self.formData.isPKACodePresent():
			#To be deprecated
			filename = os.path.join(webStorageDir, '%s.pdb' % self.formData.pdbId())
			stream = self.formData.pdbFileStream()
			self._writeStringToFile(stream.read(), filename)
			pdbFilename = self.formData.pdbFilename()
			#Add structure to db
			stream.seek(0)
			self.job.setStructure(stream.read())
			
			#Check the provided pdb
			structure = self.job.protoolStructure()
			if structure.hasMissingMainChainAtoms():
				suggestion = "PEAT-SA requires that all main-chain heavy atoms are present in the structure.<br>"
				suggestion = suggestion + "You could try submitting a fragment of the structure that meets this requirement."
				self.errorData = {"domain":"PDT.SubmissionDomain",
					"description":"Error in submitted pdb structure.",
					"detailedDescription":"The supplied structure is missing main-chain heavy atoms",
					"recoverySuggestion": suggestion}
				return
			elif structure.hasChainBreak():
				suggestion = "PEAT-SA requires that all chains in submitted structures are complete.<br>"
				suggestion = suggestion + "You could try submitting a fragment of the structure that meets this requirement."
				self.errorData = {"domain":"PDT.SubmissionDomain",
					"description":"Error in submitted pdb structure.",
					"detailedDescription":"The supplied structure contains at least one chain-break.",
					"recoverySuggestion": suggestion}
				return
		else:
			pdbFilename = os.path.join(workingDir, self.formData.pdbFilename())	
	
	 
	
	
		if self.formData.isLigandPresent():
			filename = os.path.join(webStorageDir, '%s' % self.formData.ligandFilename())
			stream = self.formData.ligandFileStream()
			self._writeStringToFile(stream.read(), filename)
			stream.seek(0)
			self.job.setLigand(stream.read())
			
		#Add email address - This could be just 'Unknown' if none was provided
		self.job.setEmail(self.formData.email())	
	
		#Create the mutation string.
		#This also writes out the mutation file if neccessary
		mutationString = CreateMutationString(self.formData.mutationData, webStorageDir, inputDir, self.job)
		calculationString = CreateCalculationString(self.formData.calculations(), self.formData.ligandFilename())
		
		#Create the run string	
		self.runString = "%s -p %s -w %s -o %s -j %s -v %s %s" % (backendExecutable, os.path.join(inputDir, pdbFilename), 
									workingDir, outputDir, self.job.identification, 
									calculationString, mutationString) 
									
		if self.formData.isIonisableGroupsPresent():
			self.runString += " --ionisableGroups=%s" % self.formData.ionisableGroups()							
									
	def runBackend(self):
	
		'''Executes the command string for the job via Popen
		
		Returns:
			An URL which specifies a page giving information about the Job
			or None if construct() has not been called.'''
	
		if self.runString == None:
			return
	
		try:
			#Put selected calculations into state Queued regardless of whether this following works or not.
			#This avoid a possible a race condition if the backed is launched quickly between who
			#modified the jobs state first
			states = self.job.calculationStates()
			for calculation in states.keys():
				if states[calculation] != "NotSelected":
					self.job.setCalculationState(calculation, "Queued")
		
			#Start the job running
			#FIXME - Wont really be able to use files for a log
			standardOut = open(os.path.join(self.webServerDirectory, "PDTWebRuns.log"), "a")
			standardOut.write("\n----------------------------------------------------------\n")
			standardOut.write("\nRunning Job %s\nDate %s\nWebScript command line %s\n\n" % (self.job.identification, self.job.date, self.runString))
			standardOut.flush()
			standardError = open(os.path.join(self.webServerDirectory, "PDTWebErrors.log"), "a")
			process = subprocess.Popen(self.runString, shell=True, stdout=standardOut, stderr=standardError)
			
			#Wait until the job is running
			time.sleep(1.0)
			process.poll()
				
			if process.returncode != None and process.returncode != 0:
				string = "Unable to launch job - launch script exited with error %d" % process.returncode
				standardError.write(string)
				raise Exceptions.SubmissionException, string

			standardOut.close()
			standardError.close()	
				
			#Constrtuct the url for the processing job page
			#Pass the information on what is being calculated on aswell
			#The elements here are scheme, location, hierarchical path, parameters, query, fragment
			location = os.environ['SERVER_NAME'] + ":" + os.environ['SERVER_PORT']
			components = ("http", location, "PEATSA/Pages/Results.php", "", "jobId=%s" % self.job.identification, "")
			resultURL = urlparse.urlunparse(components)
	
		except BaseException, data:
			if hasattr(data, "child_traceback"):
				errorString = "Exception %s. \n Child traceback %s" % (data, data.child_traceback)
			else:
				errorString = "Exception - %s" % data
				
			#Delete job information from the db if it exists
			if self.job is not None:
				self.jobManager.deleteJob(self.job)

			self.errorData = {"domain":"PDT.SubmissionDomain",
					  "description":"Error when attempting to run job.",
					  "detailedDescription":errorString,
					   "recoverySuggestion":"This is a bug - please contact the developers."}

			resultURL = self.errorURL()
	
		return resultURL
	
			
	def errorURL(self):
	
		'''Returns an URL for the WebApp error page if there was an error with the form data.
		
		On loading this URL the user will be presented with information regarding what went wrong.
		
		If there was no problem this methods return None'''
	
		if self.errorData is None:
			return None
	
		return ConstructErrorURL(domain=self.errorData["domain"],
			description=self.errorData["description"],
			detailedDescription=self.errorData["detailedDescription"],
			recoverySuggestion=self.errorData["recoverySuggestion"])

	def error(self):
	
		'''See FormData.errorURL docs for information'''
	
		return self.errorData
		
		
class FormData:

	'''Class representing the form data submitted from the WebApp main page'''
	
	def __init__(self, formData):
	
		'''Initialises the Data class.
		
		formData must be an instance of the cgi.FieldStorage class'''
		
		self.errorData = None
		self.formData = formData
		self.pdbStream = None
		self._processMutationData()
		self._checkSubmittedData()
		if self.errorData is None:
			self._setPDBStream()
	
	def _setPDBStream(self):
	
		'''Assigns a stream to the pdb data to the pdbStream ivar.
		
		If the stream cannot be created it sets an error.
		Note, a stream is only created if a delta-pKa calculation is not requested.
		In this case the pdb file to be used is already available'''
		
		if self.isPDBFilePresent():
			self.pdbStream = self.formData["pdbFile"].file
			self.pdbStream.seek(0)
		elif self.isPDBCodePresent():
			data = {}
			self.pdbStream = GetPDB(self.pdbId(), dict=data)
			if data.has_key('error'):
				self.errorData = {"domain":"PDT.SubmissionDomain",
					  "description":data['error'],
					  "detailedDescription": data['description'],
					  "recoverySuggestion":"Check the supplied code is valid"}
			else:	
				self.pdbStream.seek(0)

	def _checkSubmittedData(self):
	
		'''Performs a series of checks on the submitted data'''
		if not self.isPDBFilePresent() and not (self.isPDBCodePresent() or self.isPKACodePresent()):
			#No pdb supplied - error
			self.errorData = {"domain":"PDT.SubmissionDomain",
					  "description":"Error in submitted form data.",
					  "detailedDescription":"No PDB file was uploaded and no PDB code was provided. Hard to do a calculation on nothing!",
					  "recoverySuggestion":"Head back to the main page and upload a PDB or provide a PDB code."}
		elif not self.isCalculationDataPresent():
			#No calculation specified
			self.errorData = {"domain":"PDT.SubmissionDomain",
					"description":"Error in submitted form data.",
					"detailedDescription":"At least one calculation type must be selected.",
					"recoverySuggestion":"Head back to the main page and choose some calculations."}
		elif not self.isMutationDataPresent():
			self.errorData = {"domain":"PDT.SubmissionDomain",
					"description":"Error in submitted form data.",
					"detailedDescription":"The mutations to perform must be specified.",
					"recoverySuggestion":"Head back to the main page and choose some mutations or upload a mutation file."}
		elif self.calculations().count('binding') == 1 and not self.isLigandPresent():
			self.errorData = {"domain":"PDT.SubmissionDomain",
					"description":"Error in submitted form data.",
					"detailedDescription":"Binding selected but no ligand provided.",
					"recoverySuggestion":"Head back to the main page and upload a ligand."}
		elif self.calculations().count('scan') == 1 and not self.isPKACodePresent():
			self.errorData = {"domain":"PDT.SubmissionDomain",
					"description":"Error in submitted form data.",
					"detailedDescription":"pKa Scan selected but no pKa calculation code provided.",
					"recoverySuggestion":"In order to perform a scan you must have previously completed a pKa calculation."}
		elif self.calculations().count('scan') == 0 and self.isPKACodePresent():
			self.errorData = {"domain":"PDT.SubmissionDomain",
					"description":"Error in submitted form data.",
					"detailedDescription":"pKa calculation code provided but pKa scan not selected.",
					"recoverySuggestion":"Please select delta pKa option if you want to perform a scan."}				

		#If theres been an error at this stage return now
		if self.errorData is not None:
			return

		#Check the submitted PDB filename
		# In order to standardize the names of the directories 
		# (so each can be identified with a specific id) the following operations are required
		# - The filename must be of the form PDBID.pdb
		# - The PDBID must be all lowercase
		# - The PDBID must not containing any punctuation marks except for underscores
		# - No spaces allowed
		
		if self.isPDBFilePresent():
			self.standardizedPDBFilename = ProcessPDBFilename(self.formData["pdbFile"].filename)
		elif self.isPDBCodePresent():
			self.standardizedPDBFilename = ProcessPDBFilename(self.formData.getvalue("pdbCode") + ".pdb")
		elif self.isPKACodePresent():
			self.standardizedPDBFilename = self.pKaCode() + ".pdb"	
		
		#Check that after the processing pdbFilename is not just an extension
		if self.standardizedPDBFilename == "" or self.standardizedPDBFilename[0] == ".":
			self.errorData = {"domain":"PDT.SubmissionDomain",
					"description":"Error in submitted form data.",
					"detailedDescription":"The filename of the uploaded pdb is invalid.",
					"recoverySuggestion":"Go back to the main page and check the naming guidelines for uploaded pdb files."}
					
		#Check the ligand file extension is mol2 (if it exists).
		if self.isLigandPresent():
			components = os.path.splitext(self.ligandFilename())
			if len(components) != 2:
				self.errorData = {"domain":"PDT.SubmissionDomain",
					"description":"Error in submitted form data.",
					"detailedDescription":"The filename of the uploaded ligand is invalid (missing extension).",
					"recoverySuggestion":"Go back to the main page and check the naming guidelines for uploaded ligand files."}
			elif components[1].lower() != ".mol2":
				self.errorData = {"domain":"PDT.SubmissionDomain",
					"description":"Error in submitted form data.",
					"detailedDescription":"The filename of the uploaded ligand is invalid - %s." % components[1],
					"recoverySuggestion":"The filename extension must be mol2"}
		
	def _processMutationData(self):
	
		self.mutationData = {"Type":"Unknown", "Data":"Unknown"}
		if self.isMutationDataPresent():
			#Process the mutations
			if self.isMutationList():
				self.mutationData['Type']='mutationList'
				#See if they uploaded a file or typed one in
				if self.formData["mutationFile"].filename:
					#The file element is a stream
					self.mutationData['Data']=self.formData["mutationFile"].file
				else:
					list = self.formData.getvalue('mutationFileBox')
					#Create a file-like stream object for the list
					self.mutationData['Data'] = StringIO.StringIO(list)
			else:
				#Must be Resiude Scan since otherwise we wouldn't be here
				self.mutationData['Type']='mutation'
				self.mutationData['Data']=self.formData["residue"].value

	def error(self):
	
		'''Returns a dictionary containing details on any errors with the form data
		
		The dictionary has the following keys
		
		- domain
		- description
		- detailedDescription
		- recoverySuggestion
		
		The method returns None if there is no problem with the form data'''
	
		return self.errorData

	def errorURL(self):
	
		'''Returns an URL for the WebApp error page if there was an error with the form data.
		
		On loading this URL the user will be presented with information regarding what went wrong.
		
		If there was no problem this methods return None'''
	
		if self.errorData is None:
			return None
	
		return ConstructErrorURL(domain=self.errorData["domain"],
			description=self.errorData["description"],
			detailedDescription=self.errorData["detailedDescription"],
			recoverySuggestion=self.errorData["recoverySuggestion"])

	def isPDBFilePresent(self):

		'''Returns True if the form data contains a PDB file'''

		pdbProvided = False
		if self.formData.has_key("pdbFile"):
			if self.formData["pdbFile"].filename != "":
				pdbProvided = True;
		
		return pdbProvided
	
	def isLigandPresent(self):

		'''Returns True if formData contains a PDB file'''

		provided = False
		if self.formData.has_key("ligandFile"):
			if self.formData["ligandFile"].filename != "":
				provided = True;
		
		return provided
				
	def isCodePresent(self):
	
		'''Returns True if the code field contains text'''
		
		codeProvided = False
		if self.formData.getvalue("pdbCode") != "":
			codeProvided = True		
		
		return codeProvided
		
	def isPDBCodePresent(self):
			
		'''Returns True if the form data contains a PDB code.
		
		The form data contains a PDB code is their is an entry in the 
		code field and "Scan" is not selected'''
			
		pdbProvided = False
		if self.isCodePresent() and not self.isDeltaPKACalculationRequested():
			pdbProvided = True
		
		return pdbProvided
		
	
	def isPKACodePresent(self):
	
		'''Returns True of formData contain a pKa code
		
		The form data contains a pKa code if "Scan" is selected
		and their is an entry in the code field'''
		
		pKaProvided = False
		if self.isCodePresent() and self.isDeltaPKACalculationRequested():
			pKaProvided = True
		
		return pKaProvided
		
	def isDeltaPKACalculationRequested(self):
	
		'''Returns True if a delta pKa calculation was requested'''
		
		provided = False
		if self.calculations().count('scan') == 1:
			provided = True
			
		return provided	
		
	def isMutationDataPresent(self):

		'''Returns True if formData contains mutation information'''

		mutationProvided = False
		mutationChoice = self.formData.getlist("mutation")
		if len(mutationChoice) == 1:
			mutationProvided = True
			
		return mutationProvided
		
	def isCalculationDataPresent(self):
	
		'''Returns True if data on what calculations to perform is present'''
	
		present = False
		if len(self.calculations()) != 0:
			present = True
			
		return present	
		
		
	def isIonisableGroupsPresent(self):
	
		'''Returns True if ionisable groups were specified AND a dpKa calculation was requested'''
	
		present = False
		if self.isDeltaPKACalculationRequested():
			string = self.ionisableGroups()
			if string != "" and string.lower() != "all":
				present = True
		
		return present
		
	def isMutationList(self):

		'''Returns True if the formData contains a mutationList file'''

		mutationList = False
		mutationChoice = self.formData.getlist("mutation")
		if mutationChoice.count('mutationFile') != 0:
			mutationList=True
		
		return mutationList	
			
	def calculations(self):
	
		'''Returns a list containing the names of the calculations requested'''
	
		return self.formData.getlist("calculation")	
		
	def pdbFilename(self):
	
		'''Returns the filename of the pdb - note this does not include a path'''
	
		return self.standardizedPDBFilename
		
	def pdbId(self):
	
		return os.path.splitext(self.pdbFilename())[0]
		
	def pdbFileStream(self):
	
		'''Returns an opened stream to the pdbFile'''
		
		self.pdbStream.seek(0)
		return self.pdbStream	
		
	def ligandFilename(self):
	
		filename = self.formData["ligandFile"].filename
		filename = os.path.basename(filename)
		return self.formData["ligandFile"].filename
		
	def ligandFileStream(self):
	
		'''Returns an opened stream to the ligand file'''
	
		stream = self.formData["ligandFile"].file
		stream.seek(0)
		return stream	
								
	def mutationListFileStream(self):
	
		'''Returns an opened stream to the mutationList file'''
	
		if self.mutationData['Type'] == 'mutationList':
			stream = self.mutationData["Data"]
			stream.seek(0)
			return stream	
			
	def pKaCode(self):
	
		'''Returns the pKa code.
		
		If none is present returns None.
		A pKa code is deemed present if the code field is filled
		and a delta-pKa calculation is requested'''
		
		retval = None
		if self.isPKACodePresent():
			retval = self.formData.getvalue("pdbCode")
		
		return retval
		
	def email(self):
	
		'''Returns the email address if their is one'''
		
		retval = self.formData.getvalue("email")
		if retval == "":
			retval = "Unknown"	
		
		return retval
		
	def ionisableGroups(self):
	
		return self.formData.getvalue("ionisableGroups")
				
				
		
