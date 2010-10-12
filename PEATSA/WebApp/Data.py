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

'''Contains classes representing the WebApp's data and for creating and managing it'''

import MySQLdb, os.path, StringIO, pickle, threading
import UtilityFunctions, Exceptions
import PEATSA.Core as Core

def SerializeDictionary(aDict):

	'''Serializes a dictionary as a string'''

	stream = StringIO.StringIO()
	pickler = pickle.Pickler(stream)
	pickler.dump(aDict)
	stream.flush()
	data = stream.getvalue()
	stream.close()

	return data
	
def UnserializeDictionary(string):

	'''Unserializes a dictionary stored as a string'''

	stream = StringIO.StringIO(string)
	dict = pickle.load(stream)
	stream.close()

	return dict

def DataIDsInTable(dataTable, connection):

	'''Returns all the values of the PDBID column from the given table in the database.'''

	cursor = connection.cursor()
	query = "SELECT PDBID FROM %s" % (dataTable)
	cursor.execute(query)
	data = cursor.fetchall()
	cursor.close()
	
	if len(data) is not 0:
		return [element[0] for element in data]
	else:
		return None
		
def JobIDsInTable(dataTable, connection):

	'''Returns all the values of the JobID column from the given table in the database.'''

	cursor = connection.cursor()
	query = "SELECT JobID FROM %s" % (dataTable)
	cursor.execute(query)
	data = cursor.fetchall()
	cursor.close()
	
	if len(data) is not 0:
		return [element[0] for element in data]
	else:
		return None	

def jobManagerWithDefaultConnection():

	return JobManager(UtilityFunctions.ConnectionFromDefaultConfiguration())
	
class RunManager:

	'''Factory class for creating new Run instances'''

	def __init__(self, connection, groupId):
		
		'''Returns a new RunManager instance.
		
		RunManagers create new Run objects along with the corresponding entries in an associated database.
		As such they need to be supplied with a connection to the associated database as a user who has INSERT privilages
		
		Parameters: 
			groupId - A string used to uniquely name the tables in the database this RunManager will use.
			
		(TODO:) If not already present the RunManager will create a set of tables with the folowing names
			- $(groupID)Runs
			- $(groupID)Jobs
			- $(groupID)RunJobs
			- $(groupID)JobData'''
	
		self.connection = connection
		self.cursor = connection.cursor();
		self.groupId = groupId
		self.runTable = groupId + 'Runs'
		self.jobTable = groupId + 'Jobs'
		self.jobAssocTable = groupId + 'RunJobs'
		self.jobDataTable = groupId + 'JobData'
		self.manager = None
		
	def __del__(self):
		
		'''Close database connection'''
		
		try: 
			#Inside a try incase conneciton is already closed
			self.cursor.close()
			self.connection.close()	
		except Exception:
			pass

	def __str__(self):

		return 'RunManager managing %d runs for group %s on %s' % (self.numberOfRuns(), self.groupId, self.connection.get_host_info())

	def __getitem__(self, index):
	
		'''runManager[i] returns the i'th run (ordered by id)'''
	
		id = self.allRunIds()[index]
		run = Run(id, self.connection, self.groupId)
		
		return run
		
	def numberOfRuns(self):
		
		'''Returns the number of runs managed by the receiver'''

		query = """SELECT COUNT(*) FROM %s""" % self.runTable
		self.cursor.execute(query)
		data = self.cursor.fetchall()

		return int(data[0][0])
	
	def createRun(self, name=None, calculations=[], metadata={}):

		'''Creates a new Run.
		
		Parameters 
			name - A string used to identify the run
			If not provided defaults to Run%d where %d is the number of the run in the run table	
		
		Returns:
			A instance of the Job class'''

		#Create a unique id for the job
		runId = UtilityFunctions.GetNewSID(self.groupId)
		
		if name is None:
			name = "Run%d" % (len(self.allRuns()))

		#Serialize metadata
		data = SerializeDictionary(metadata)
		data = MySQLdb.escape_string(data)	

		#Create the insert SQL command 		
		dateString = UtilityFunctions.CreateDateString()
		data = (self.runTable, runId, name, dateString, data)
		insertQuery = """INSERT INTO %s (RunID, Name, Date, Metadata) VALUES ('%s', '%s', '%s', '%s')""" % data
				
		res = self.cursor.execute(insertQuery)
		self.connection.commit()
		run = Run(runId, self.connection, self.groupId)
		
		return run
		
	def deleteRun(self, run):
	
		'''Deletes a run from the database
		
		Parameters
			 run - a instance of the Run class'''
		
		deleteQuery = """DELETE FROM %s WHERE RunId='%s'""" % (self.runTable, run.identification)
		self.cursor.execute(deleteQuery)
		self.connection.commit()
		
	def allRunIds(self):
	
		'''Returns all the run ids managed by the receiver'''	
		
		query = "SELECT RunID FROM %s" % (self.runTable)
		self.cursor.execute(query)
		data = self.cursor.fetchall()
		
		if len(data) is not 0:
			return [element[0] for element in data]
		else:
			return None	

	def jobManager(self):

		'''Returns a JobManager instance that manages the Job table'''

		if self.manager is None:
			self.manager = JobManager(self.connection, self.jobTable)

		return self.manager
	
class Run:

	'''Instances of this class represent a set of PEAT-SA calculations . 
	The actual Run data, for example the calculations being run, is stored in a SQL database. 
	The methods of this class allow this data to be retrived and set, hiding all access and query details.
	
	The main attributes of a Run object are a unique id and and a name.
	Run objects must be initialy created using the createRun() method of a RunManager class'''

	def __init__(self, runID, connection, groupId=""):
		
		'''Returns a Run object representing a Run in a SQL database.
		
		Note: The Run must have been created previously using the RunManager.
		
		Parameters:
			runID - The unique ID of the job.
			connection - A MySQLDB connection object representing a connection
			to the database containing the job.
			groupId - A string identifying the tables in the SQL database where this Runs data is stored
			'''
		
		self.identification = runID
		self.connection = connection
		self.cursor = self.connection.cursor();
		self.runTable = groupId + 'Runs'
		self.jobTable = groupId + 'Jobs'
		self.jobAssocTable = groupId + 'RunJobs'
		self.jobDataTable = groupId + 'JobData'
				
		#Check the job exists.
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "RunId %s does not exist in the database" % self.identification
			
		#Get date
		self.cursor.execute("SELECT Date FROM %s WHERE RunID = '%s'" % (self.runTable, self.identification))
		rows = self.cursor.fetchall()
		self.date = rows[0][0]
			
	def __del__(self):
				
		'''Close database connection'''

		self.cursor.close()

	def __str__(self):

		if self.exists():
			data = (self.name(), self.runTable, self.connection.get_host_info(), self.numberOfJobs()) 
			return 'Run %s in %s on %s. Contains %d jobs' % data
		else:
			return 'Run %s no longer exists in table %s of database %s' % (self.identification, self.runTable, self.connection.get_host_info())
		
					
	def __getitem__(self, index):
	
		'''run[i] returns the i'th job in the run (ordered by id)'''
	
		id = self.allJobIds()[index]
		job = Job(id, self.connection, self.jobTable, self.jobDataTable)
		
		return job							
														
	def exists(self):
		
		'''Returns True if data corresponding to the receiver exists in the database'''
		
		#Check the job entry exists in the database
		query = "SELECT * FROM %s WHERE RunID = '%s'" % (self.runTable, self.identification)
		self.cursor.execute(query)
		rows = self.cursor.fetchall()

		if len(rows) == 0:
			return False
		else:
			return True	
		
	def allJobIds(self):
	
		'''Returns the identifications of all the jobs in this run'''
		
		#Check the job entry exists in the database
		query = "SELECT JobID FROM %s WHERE RunID = '%s'" % (self.jobAssocTable, self.identification)
		self.cursor.execute(query)
		rows = self.cursor.fetchall()
		
		if len(rows) != 0:
			return [row[0] for row in rows]
		else:
			return []	
		
	def allJobs(self):
	
		'''Returns Job instances representing all the jobs in this run'''
		
		#Check the job entry exists in the database
		query = "SELECT * FROM %s WHERE RunID = '%s'" % (self.jobAssocTable, self.identification)
		self.cursor.execute(query)
		rows = self.cursor.fetchall()
		
		ids = self.allJobIds()
		jobs = []
		for id in ids:
			job = Job(id, self.connection, self.jobTable, self.jobDataTable)
			jobs.append(job)
		
		return jobs
		
	def addJob(self, job):
	
		'''Adds a job to the set of jobs in this run
		
		Parameters:
			job - A Job instance'''
	
		values = (self.jobAssocTable, self.identification, job.identification)
		statement = """INSERT IGNORE INTO %s (RunID, JobID) VALUES ('%s', '%s')""" % values
		self.cursor.execute(statement)
		self.connection.commit()
		
	def removeJob(self, job):
	
		'''Removes a job from the set of jobs in the run
		
		Parameters: 
			job - A job instance'''
			
		deleteQuery = """DELETE FROM %s WHERE RunId='%s'""" % (self.jobAssocTable, self.identification)
		self.cursor.execute(deleteQuery)
		self.connection.commit()	
			
	def numberOfJobs(self):
	
		'''Returns the number of jobs in this run'''
		
		#Check the job entry exists in the database
		query = "SELECT COUNT(*) FROM %s WHERE RunID = '%s'" % (self.jobAssocTable, self.identification)
		self.cursor.execute(query)
		data = self.cursor.fetchall()
		
		return int(data[0][0])
	
	def name(self):
	
		'''Returns the name of the receiver'''
		
		query = "SELECT Name FROM %s WHERE RunID = '%s'" % (self.runTable, self.identification)
		self.cursor.execute(query)
		rows = self.cursor.fetchall()
		
		return rows[0][0]
		
	def setName(self, name):
	
		'''Sets the name of the receiver'''
		
		statement = """UPDATE %s SET Name='%s' WHERE RunID='%s' """ % (self.runTable, name, self.identification)
		self.cursor.execute(statement)
		self.connection.commit()

class JobManager:

	'''Factory class for creating new Job instances'''

	def __init__(self, connection, jobTable='Jobs'):
		
		'''Returns a new JobManager instance.
		
		JobManagers create new Job objects along with the corresponding entries in the PDT database.
		As such they need to be supplied with a connection to a valid PDT database as a user who has INSERT privilages'''
	
		self.connection = connection
		self.cursor = connection.cursor();
		self.jobTable = jobTable
		
		#Logging ivars
		self.logging = False
		self.logFile = None
		self.logInterval = None
		self.logThread = None
		
	def __del__(self):
		
		'''Close database connection'''
		
		try: 
			#Inside a try incase conneciton is already closed
			self.cursor.close()
			self.connection.close()	
		except Exception:
			pass
	
	def __createInsertString__(self, jobID, pdbID, dateString, scan, stability, binding, metadata):

		'''Creates an SQL insert statement based on the provided data.
		
		Note: The table has five columns - JobID, Date, Stability, Scan and Binding.
		The last three are bools which are set to one if the corresponding entry is
		present in the calculations list.
		
		Parameters:
			jobID: The id for the job
			date: The date for the entry as a string in the format yyyy-mm-dd hh:mm:ss'''

		data = SerializeDictionary(metadata)
		data = MySQLdb.escape_string(data)
		
		data = [self.jobTable, jobID, pdbID, dateString, scan, stability, binding, data]
		
		insertData = """INSERT INTO %s (JobID, PDBID, Date, Scan, Stability, Binding, Metadata) VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s')""" % tuple(data)
	
		return insertData
	
	def createJob(self, pdbId, calculations=[], dataTable='Data', metadata={}):

		'''Creates a new job entry in the database.
		
		Parameters
			Calculations is a list containing one or more of 'scan', 'stability' or 'binding'.
			dataTable - The table the job will store its data sets to
		
		Returns:
			A instance of the Job class'''


		if len(calculations) is 0:
			raise Core.Exceptions.ArgumentError, 'At least one calculation must be specified'
		else:
			aDict = {'scan':'NotSelected', 'stability':'NotSelected', 'binding':'NotSelected'}
			for calculation in calculations:
				if aDict.has_key(calculation):
					aDict[calculation] = 'Selected'
				else:
					raise Core.Exceptions.ArgumentError, 'Invalid calculation types requested - %s' % calculations

		#Create a unique id for the job
		jobId = UtilityFunctions.GetNewSID("PDT")

		#Create the insert SQL command 		
		dateString = UtilityFunctions.CreateDateString()
		insertQuery = self.__createInsertString__(jobId, pdbId, dateString, 
					aDict['scan'], aDict['stability'], aDict['binding'], metadata)
		
		res = self.cursor.execute(insertQuery)
		self.connection.commit()
		job = Job(jobId, self.connection, self.jobTable, dataTable=dataTable)
		
		return job
		
	def deleteJob(self, job):
	
		'''Deletes a job from the database
		
		Parameters
			job - a instance of the Job class'''
		
		#Check job is not running
		if job.state() != 'Running' and job.state() != 'Launched':
			deleteQuery = """DELETE FROM %s WHERE JobId='%s'""" % (job.jobTable, job.identification)
			self.cursor.execute(deleteQuery)
			self.connection.commit()
		else:
			raise Core.Exceptions.ArgumentError, 'Cannot delete job - is running'	
		
	def allJobs(self):
	
		'''Returns all the job ids stored in the database'''	
		
		return JobIDsInTable(self.jobTable, self.connection)
		
	def jobStates(self):
	
		'''Returns a dict describing the state of each job
		
		The dicts keys are job ids. Each values is a dict containing the the jobs state and date'''
		
		self.connection.commit()
		selectQuery = """SELECT JobID, State, Date FROM %s""" % (self.jobTable)
		self.cursor.execute(selectQuery)		
		rows = self.cursor.fetchall()
		
		stateDict = {}
		for row in rows:
			stateDict[row[0]] = {'State':row[1], 'Date':row[2]}
				
		return stateDict
		
		
	def jobsWithCalculationsInState(self, state):
	
		'''Returns any job which has a calculation in \e state
		
		See Job.calculationState() for more info'''
	
		self.connection.commit()
		selectQuery = """SELECT JobID, Stability, Scan, Binding FROM %s""" % (self.jobTable)
		self.cursor.execute(selectQuery)		
		rows = self.cursor.fetchall()
		print rows
		
		ids = []
		for row in rows:
			if state in row:
				ids.append(row[0])
				
		return ids
		
	def jobsInState(self, state):
	
		'''Returns any job which is in state
		
		Valid states:
			UnderConstruction
			Ready
			Launched
			Running
			Finished'''
	
		#Refresh connection
		self.connection.commit()
		selectQuery = """SELECT JobID, State FROM %s""" % (self.jobTable)
		self.cursor.execute(selectQuery)		
		rows = self.cursor.fetchall()
		
		ids = []
		for row in rows:
			if state in row:
				ids.append(row[0])
				
		return ids
		
	def logJobStates(self, file):
	
		'''Writes the state of all jobs to file as a pickled dictionary'''
	
		stateDict = self.jobStates()
		stream = open(file, 'w+')
		pickler = pickle.Pickler(stream)
		pickler.dump(stateDict)
		stream.flush()
		stream.close()
		
	def setJobStateLogging(self, file, interval=60):	
	
		'''Sets the receiver to write the state of jobs to file every interval secs'''
	
		self.logging = True
		self.logFile = file
		self.logInterval = interval
		#Create initial log.
		self.logJobStates(file)
		self.logThread = threading.Timer(interval, JobManager.logJobStates, [self, file])
		self.logThread.start()
	
	def stopLogging(self):
	
		'''Stops automatic logging'''
	
		if self.isLogging():
			self.logThread.cancel()
			self.logThread = None
			self.logFile = None
			self.logInterval = None
			self.logging = False		
										
	def isLogging(self):
	
		'''Returns True if automatic logging is turned on'''

		return self.logging																			
																			
	
class Job:

	'''Instances of this class represent a PDT web-server job . 
	The actual Job data, for example the calculations being run, is stored in a SQL database. 
	The methods of this class allow this data to be retrived and set, hiding all access and query details.
	
	The main attributes of a Job object are a unique id and the pdb code of the protein the job is working on.
	Job objects must be initialy created using the createJob() method of a JobManager class
	
	Note: JobData must be set through this class'''

	def __init__(self, jobID, connection, jobTable="Jobs", dataTable="Data"):
		
		'''Returns a Job object representing a job in a SQL database.
		
		Note: The job must have been created previously using the JobManager.
		
		Parameters:
			jobID - The unique ID of the job.
			connection - A MySQLDB connection object representing a connection
			to the database containing the job.
			jobTable - The name of the table in the database where the job will be added.
			Note the table must have the correct fields etc.
			dataTable - The name of the table where the jobs data is/willbe stored
			'''
		
		self.identification = jobID
		self.connection = connection
		self.cursor = self.connection.cursor();
		self.jobTable = jobTable
		self.dataTable = dataTable
		self.validStates = ['UnderConstruction', 'Ready', 'Launched', 'Running', 'Finished']
		
		#Check the job exists.
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
			
		#Get date
		self.cursor.execute("SELECT Date FROM %s WHERE JobID = '%s'" % (self.jobTable, self.identification))
		rows = self.cursor.fetchall()
		self.date = rows[0][0]
		
		#Get PDBID and Create SQLDataSet instance
		self.cursor.execute("SELECT PDBID FROM %s WHERE JobID = '%s'" % (self.jobTable, self.identification))
		rows = self.cursor.fetchall()
		self.pdbID = rows[0][0]
		self.data = SQLDataSet(jobIdentification=self.identification, location=self.connection, table=self.dataTable)	
			
	def __del__(self):
				
		'''Close database connection'''

		self.cursor.close()
			
	def exists(self):
		
		'''Returns True if data corresponding to the receiver exists in the database'''
		
		#Check the job entry exists in the database
		self.connection.commit()
		query = "SELECT * FROM %s WHERE JobID = '%s'" % (self.jobTable, self.identification)
		self.cursor.execute(query)
		rows = self.cursor.fetchall()

		if len(rows) == 0:
			return False
		else:
			return True											
	
	def setState(self, aState):
	
		'''Sets the state of the receiver to aState.
		
		aState must be one of the following
		
		UnderConstruction - Not all data for the job has been set yet. It should not be launched.
		Ready - The job is ready to go whenever
		Launched - The job has been launched by some service. Note: It may not be actually running.
		Running	- The job is running
		Finished - The job has finished'''	
		
		if self.exists():
			if aState in self.validStates:
				statement = """UPDATE %s SET State='%s' WHERE JobID='%s'""" % (self.jobTable, aState, self.identification)
				self.cursor.execute(statement)
				self.connection.commit()
			else:
				raise Core.Exceptions.ArgumentError, 'Provided state %s not valid' % aState	
		else:
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification								
		
	def state(self):
	
		'''Returns the jobs state
		
		See setState() for more information.
		Use calculationState() to get more detailed information
		on  the status of the component calculations'''
		
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
		
		self.connection.commit()
		selectQuery = """SELECT State FROM %s WHERE JobID='%s'""" % (self.jobTable, self.identification)
		self.cursor.execute(selectQuery)		
		rows = self.cursor.fetchall()
		
		return rows[0][0]

	def setQueueStatusMessage(self, message):

		'''Sets the jobs queue status message.

		The status message should described the jobs current position in the queue etc.'''

		if self.exists():
			statement = """UPDATE %s SET QueueStatusMessage='%s' WHERE JobID='%s'""" % (self.jobTable, message, self.identification)
			self.cursor.execute(statement)
			self.connection.commit()
		else:
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification								

	def queueStatusMessage(self):
	
		'''Returns the jobs queue status message (a string)
		
		The message has no predefined format. It is used to provide some informative
		details about the current position of a job in a queue to users'''
		
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
		
		self.connection.commit()
		selectQuery = """SELECT QueueStatusMessage FROM %s WHERE JobID='%s'""" % (self.jobTable, self.identification)
		self.cursor.execute(selectQuery)		
		rows = self.cursor.fetchall()
		
		return rows[0][0]
	
	def addStabilityResults(self, matrix):
	
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
			
		#Add the data to the SQL database
		self.data.addMatrix(matrix, 'StabilityResults')
	
	def addScanResults(self, matrix):
	
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
			
		#Add the data to the SQL database
		self.data.addMatrix(matrix, 'ScanResults')
	
	def addResults(self, matrix, name):
	
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
			
		#Add the data to the SQL database
		self.data.addMatrix(matrix, name)	

	def addDataSet(self, dataSet):
	
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
			
		for name in dataSet.matrixNames():
			self.data.addMatrix(dataSet.__getattr__(name), name)	
		
	def addImageData(self, data, name):
	
		'''Adds the image data to the database calling it name'''
		
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification

		values =  (self.identification, name, MySQLdb.escape_string(data), len(data), len(data), MySQLdb.escape_string(data))
		statement = """INSERT INTO Images (JobID, Name, Content, Size) VALUES ('%s', '%s', '%s', '%d') ON DUPLICATE KEY UPDATE Size='%d', Content='%s'""" % values
		self.cursor.execute(statement)
		self.connection.commit()
	
	def setCalculationState(self, calculation, state):
	
		'''Sets the state of a calculation in the job.
		
		Parameters:
			calculation - One of Scan, Binding or Stability.
			state - Queued, Waiting, Running or Finished.
			Only applicable if initial state is selected'''
	
		if self.exists():
			if self.calculationStates()[calculation] != 'NotSelected':
				statement = """UPDATE %s SET %s='%s' WHERE JobID='%s'""" % (self.jobTable, calculation, state, self.identification)
				self.cursor.execute(statement)
				self.connection.commit()
			else:
				raise Core.Exceptions.ArgumentError, 'Specified calculation, %s, not selected. Cannot set state' % calculation	
		else:
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
			
	def calculationStates(self):
	
		'''Returns the state of the calculations.
		
		This is a dictionary of calculationName:JobState pairs.
		Valid States are 
			- NotSelected
			- Selected
			- Queued
			- Waiting
			- Running
			- Finished'''	
			
		#FIXME add test for valid job states									
		if self.exists():
			self.connection.commit()
			statement = """SELECT Stability, Scan, Binding FROM %s WHERE JobID='%s' """ % (self.jobTable, self.identification)
			self.cursor.execute(statement)
			rows = self.cursor.fetchall()
			return dict(zip(['Stability', 'Scan', 'Binding'], rows[0]))
		else:
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
	
			
	
	def setError(self, description, detailedDescription):

		if self.exists():
			description = MySQLdb.escape_string(description)
			detailedDescription = MySQLdb.escape_string(detailedDescription)
			values =  (self.jobTable, description, detailedDescription, self.identification)
			statement = """UPDATE %s SET Error='1', ErrorDescription='%s', DetailedDescription='%s' WHERE JobID='%s' """ % values
			self.cursor.execute(statement)
			self.connection.commit()
		else:
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
			
	def error(self):
	
		'''If there was an error during the job returns a dictionary detailing the error, otherwise returns None
		
		The dictionary has two keys:
			ErrorDescription - A short description of the type of error
			DetailedDescription - A longer description of what went wrong'''
	
		if self.exists():
			self.connection.commit()
			statement = """SELECT Error, ErrorDescription, DetailedDescription FROM %s WHERE JobID='%s' """ % (self.jobTable, self.identification)
			self.cursor.execute(statement)
			rows = self.cursor.fetchall()
			if bool(rows[0][0]) == False:
				return None
			else:
				return dict(zip(['ErrorDescription', 'DetailedDescription'], rows[0][1:]))
		else:
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification		
	
		
	def setMutation(self, mutationCode):
	
		'''Sets the job to do a scan by mutating each residue in structure to mutation.
		
		Parameters:
			mutationCode - A Three letter amino acid code.
		
		Exceptions:
			Raises an IOError if structureFile can't be opened'''
		
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification

		#Set mutation command
		statement = """UPDATE %s SET MutationCommand='%s' WHERE JobID='%s' """ % (self.jobTable, '--mutation', self.identification)
		self.cursor.execute(statement)
		self.connection.commit()	
		
		#Set the mutation type
		statement = """UPDATE %s SET MutationData='%s' WHERE JobID='%s' """ % (self.jobTable, mutationCode, self.identification)
		self.cursor.execute(statement)
		self.connection.commit()
		
	def mutation(self):
	
		'''Retrieves the scan mutation the job should perform.
		
		If no mutation has been set, or a mutation list has been set instead, this method return None'''
		
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
	
		mutationData = None
		if not self.isMutationList():
			self.connection.commit()
			statement = """SELECT MutationData FROM %s WHERE JobID='%s' """ % (self.jobTable, self.identification)
			self.cursor.execute(statement)
			rows = self.cursor.fetchall()
			mutationData  = rows[0][0]

		return mutationData

		
	def setMutationListFile(self, mutationList):
	
		'''Adds a mutation list data to the job.
		
		Note: This overwrites mutation data if it exists
		
		Parameters:
			mutationList - A Core.Data.MutationListFile instance
		
		Exceptions:
			Raises an IOError if structureFile can't be opened'''
		
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification

		stream = StringIO.StringIO()
		mutationList.writeToStream(stream)
		stream.flush()
		data = stream.getvalue()
		stream.close()
		
		#Set mutation command
		statement = """UPDATE %s SET MutationCommand='%s' WHERE JobID='%s' """ % (self.jobTable, '--mutationList', self.identification)
		self.cursor.execute(statement)
		self.connection.commit()	
		
		#Set the mutation type
		statement = """UPDATE %s SET MutationData='%s' WHERE JobID='%s' """ % (self.jobTable, MySQLdb.escape_string(data), self.identification)
		self.cursor.execute(statement)
		self.connection.commit()	
	
	def mutationListFile(self):
	
		'''Returns a Core.Data.MutationListFile instance containing the mutations the job is to run on
		
		If the no mutationListFile has been set for the Job this method returns None'''
	
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
	
		mutationList = None
		if self.isMutationList():
			self.connection.commit()
			statement = """SELECT MutationData FROM %s WHERE JobID='%s' """ % (self.jobTable, self.identification)
			self.cursor.execute(statement)
			rows = self.cursor.fetchall()
			mutationData = rows[0][0]
			mutationDataStream = StringIO.StringIO(mutationData)
			mutationList = Core.Data.mutationListFileFromStream(stream=mutationDataStream)

		return mutationList
		
	def isMutationList(self):
	
		'''Returns True if the jobs mutations are specified by a mutation list'''
	
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
			
		#FIXME: Catch index errors if nothing is returned - (if there are any)
		#It may be that always [[]] is returned if nothing is found	
		self.connection.commit()
		statement = """SELECT MutationCommand FROM %s WHERE JobID='%s' """ % (self.jobTable, self.identification)
		self.cursor.execute(statement)
		rows = self.cursor.fetchall()
		mutationCommand = rows[0][0]
		
		isMutationList = False
		if mutationCommand == '--mutationList':
			isMutationList = True
			
		return isMutationList					
				
	def setStructure(self, structure):
	
		'''Adds the structure of a protein to the job data.
		
		Parameters:
			structure - A string containing a pdb structure
		'''
		
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
	
		#Set mutation command
		statement = """UPDATE %s SET Structure='%s' WHERE JobID='%s' """ % (self.jobTable, MySQLdb.escape_string(structure), self.identification)
		self.cursor.execute(statement)
		self.connection.commit()	
		
	def setStructureFromFile(self, structureFile):
	
		'''Adds the structure of the protein from a PDB file to the job data
		
		Parameters: 
			structureFile - The path to a file containing a PDB structure
			
		Exceptions:
			Raises an IOError if structureFile can't be opened'''
			
		stream = open(structureFile, 'r')
		data = stream.read()
		stream.close()
		self.setStructure(data)			
		
	def structure(self):
	
		'''Returns a string containing the structure of the protein in PDB format'''

		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification

		self.connection.commit()	
		statement = """SELECT Structure FROM %s WHERE JobID='%s' """ % (self.jobTable, self.identification)
		self.cursor.execute(statement)
		rows = self.cursor.fetchall()
		return rows[0][0]
		
	def protoolStructure(self):
	
		'''Returns a protool structure instance initialised using the structure data stored in the database'''
		
		import tempfile, Protool
		
		structureData = self.structure()
		temp = tempfile.NamedTemporaryFile()
		temp.write(structureData)
		temp.flush()

		try:
			object = Protool.structureIO()
			object.readpdb(temp.name)
		except Exception, data:	
			raise Exceptions.FileFormatError, 'Format of stored PDB file for job %s not valid.\nUnderlying error - %s' % (self.identification, data)
		finally:
			temp.close()
			
		return object				 	
		
	def setLigand(self, ligand):
	
		'''Adds the structure of a protein to the job data.
		
		Parameters:
			ligand - A string containing a pdb structure
		'''
		
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
	
		#Set mutation command
		statement = """UPDATE %s SET Ligand='%s' WHERE JobID='%s' """ % (self.jobTable, MySQLdb.escape_string(ligand), self.identification)
		self.cursor.execute(statement)
		self.connection.commit()	
		
	def setLigandFromFile(self, ligandFile):
	
		'''Adds the structure of the ligand from a mol2 file to the job data
		
		Parameters: 
			ligandFile - The path to a file containing a mol2 structure of a ligand
			
		Exceptions:
			Raises an IOError if ligandFile can't be opened'''
			
		stream = open(ligandFile, 'r')
		data = stream.read()
		stream.close()
		self.setLigand(data)			
		
	def ligand(self):
	
		'''Returns a string containing the structure of the ligand in mol2 format.
		
		If no ligand is present the method returns None'''

		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification

		self.connection.commit()
		statement = """SELECT Ligand FROM %s WHERE JobID='%s' """ % (self.jobTable, self.identification)
		self.cursor.execute(statement)
		rows = self.cursor.fetchall()
		
		return rows[0][0]

	def setEmail(self, email):
	
		'''Adds an email to job'''
		
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification

		statement = """UPDATE %s SET Email='%s' WHERE JobID='%s' """ % (self.jobTable, email, self.identification)
		self.cursor.execute(statement)
		self.connection.commit()
		
	def email(self):
	
		'''Returns the email related to the job (if any)'''
		
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
		
		self.connection.commit()
		selectQuery = """SELECT Email FROM %s WHERE JobID='%s'""" % (self.jobTable, self.identification)
		self.cursor.execute(selectQuery)		
		rows = self.cursor.fetchall()
		
		return rows[0][0]
		
	def setEmailSent(self):
	
		'''Updates the stored Job data to indicate an email has already been sent about this job'''	
		
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification

		statement = """UPDATE %s SET SentMail='1' WHERE JobID='%s' """ % (self.jobTable, self.identification)
		self.cursor.execute(statement)
		self.connection.commit()	
		
	def isEmailSent(self):
	
		'''Returns True if an email has been sent in relation to this job already'''
		
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
		
		self.connection.commit()
		selectQuery = """SELECT SentMail FROM %s WHERE JobID='%s'""" % (self.jobTable, self.identification)
		self.cursor.execute(selectQuery)		
		rows = self.cursor.fetchall()
		
		return bool(rows[0][0])
		
	def metadata(self):
	
		'''Returns the metadata associated with the job.
		
		If no metadata was associated with the job this method returns None'''
		
		if not self.exists():
			raise Exceptions.DatabaseRetrievalError, "Job Id %s does not exist in the database" % self.identification
		
		self.connection.commit()
		selectQuery = """SELECT Metadata FROM %s WHERE JobID='%s'""" % (self.jobTable, self.identification)
		self.cursor.execute(selectQuery)		
		rows = self.cursor.fetchall()
		
		data = rows[0][0]
		dict = UnserializeDictionary(data)
		
		return dict
		

class SQLDataSet(Core.Data.DataSet):

	'''Class representing the output of a PEATSA run stored in a SQL database.
	
	The Database must have a table called Data with the correct format.
	
	Attributes:
		name - The name of the results
		connection - MySQLdb connection object'''

	def __addMatrix__(self, matrix, dataType):
	
		'''Adds the matrix for dataType to the SQL database.
		
		If the data for dataType is already present it is over-written
		Note: This is different to the filesystem case.'''
		
		data = matrix.csvRepresentation()
		data = MySQLdb.escape_string(data)
		values =  (self.table, self.jobId, self.name, dataType, data, len(data), len(data), data)
		statement = """INSERT INTO %s (JobID, DataSetName, MatrixName, Content, Size) VALUES ('%s', '%s', '%s', '%s', '%d') ON DUPLICATE KEY UPDATE Size='%d', Content='%s'""" % values
		self.cursor.execute(statement)
		self.connection.commit()
		
		self.data[dataType] = matrix

	def __loadMatrix__(self, name):

		'''Unarchives the matrix called name from the database and returns it.
		
		Returns None if no matrix called name exists or if the data could not be read.'''

		self.connection.commit()
		query = "SELECT Size, Content FROM %s WHERE JobID = '%s' AND DataSetName = '%s' and MatrixName = '%s'" % (self.table, self.jobId, self.name, name)
		self.cursor.execute(query)
		rows = self.cursor.fetchall()
		if len(rows) == 0:
			matrix = None
		else:	
			content = rows[0][1]
			size = rows[0][0]
			matrix = Core.Matrix.matrixFromCSVRepresentation(content)	
			
		return matrix	

	def __readData__(self):
	
		'''Reads the available data'''
		
		self.data = {}
		for name in Core.Data.peatMatrices:
			self.data[name] = self.__loadMatrix__(name)
		
	def __init__(self, jobIdentification, name="Output",  location=None, table="Data", overwrite=True):
	
		'''Initialise a new SQLDataSet instance. The instances is created if it doesn't exist.
		
		Parameters:
			jobIdentification - The identification of the job that generated this data
			name - The name of the PEATSA results.
			location - A MySQLdb connection object connecting to the database containing the results.
			If its none a connection to the default database, defined by the 
			default configuration file, is created
			table - The dataset table in the database to add the data to'''
			
		if location is not None:	
			self.connection = location
			self.ownsConnection = False
			self.overwrite=overwrite
		else:
			self.connection = UtilityFunctions.ConnectionFromDefaultConfiguration()
			self.ownsConnection = True
			if self.connection is None:
				raise Exceptions.DatabaseRetrievalError, 'Unable to connect to default database'
				
		self.cursor = self.connection.cursor()	
	
		if name is None:
			raise TypeError, 'Name must be a string'
			
		self.name = name
		self.table = table
		self.jobId = jobIdentification
		self.__readData__()	
		
	def __str__(self):
	
		'''Returns a description string'''
	
		description = 'PEAT-SA dataset %s in SQL database on %s. Contents :\n' %(self.name, self.connection.get_host_info())
		contents = ""
		for name in self.data.keys():
			if self.data[name] is not None:
				contents = contents + '\t%s\n' % name	
		
		if contents == "":
			contents = 'None'
		
		return description + contents
		
	def __del__(self):
				
		'''Close database connection'''

		self.cursor.close()
		if self.ownsConnection is True:
			self.connection.close()
			
	def addMatrix(self, matrix, name):
	
		'''Adds a matrix called name to the data dir.
		
		Note: Name is capitalized if its not already
		
		If data for name is already present, and overwriting is disabled,
		an Exceptions.DataAlreadyExistsError is raised.
		Otherwise the data is added to the self.data dictionary.'''

		if self.data.has_key(name) and self.data[name] is not None:
			if not self.overwrite:
				raise Core.Exceptions.DataAlreadyExistsError, 'A matrix called %s already present' % name
		
		name = name.capitalize()[0] + name[1:]

		self.__addMatrix__(matrix, name)
		self.data[name] = matrix
		
