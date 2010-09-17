#! /usr/bin/env python
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

'''WebScript - Wrapper around ProteinDesignTool.py containing web related extensions. 

	Written by Michael Johnston
	This programs copyright - Michael Johnston and Jens Nielsen 2008.
	For copyright of consituent programs check each

	NOTE: For detailed instructions see the programs README file.
	This tool should only be run from the JobSubmission.RunJob function or the RemoteLaunch script.
	It requires access to the PDT web-apps MySQL database

	Usage: WebScript -p [pdbname] -w [directory] -j [jobid] [optional arguments]
	
	-p - The name of the pdb the tool is to be run on
	
	-w - The directory where the program will be run. It MUST contain the pdb called [pdbname]
	     For certain calculations this directory must contain various other data files.
	     See README for more.
	     
	-j - The JobID for the web-app job     
	
	All other options are the same as for ProteinDesignTool.
	Except that the directory containing the output will be outputDir/jobid.
	The value for outputDir is determined by the entry in the webApplication.conf file
	and is not the directory from which the program is run.'''
 
import getopt, sys, os, traceback, math
import Data
import UtilityFunctions, Exceptions
import PEATSA.Core as Core

class CommandLineParser(Core.Utilities.CommandLineParser):

	'''Class for parsing the WebScripts command line and providing access to the data
	
	Adds handling of a couple of special features of the WebScript cl'''
	
	def parseCommandLine(self):

		'''Parses the command line
		
		Execeptions:
			Raises Core.Exceptions.ArgumentError if there was a problem with the arguments supplied'''
	
		try:
			self.opt = getopt.getopt(sys.argv[1:], 
					"vw:p:o:j:h", 
					["pKa", "scan", "path", "modes", "stability", 
					"mutationList=", "ionisableGroups=","mutation=", "mutants=", "ligand=", "configurationFile=", "create-configuration"])
			self.opt = dict(self.opt[0])
		except getopt.GetoptError, data:
			raise Core.Exceptions.ArgumentError, data
		
		#Create configuration instance 
		#This is needed for reading the default data and output directories
		if self.opt.has_key("--configurationFile"):
			configurationFile = self.opt["--configurationFile"]
		else:
			#Use default configuration for the web application
			configurationFile = os.path.join(UtilityFunctions.ResourceDirectory(), "webApplication.conf") 
			self.opt["--configurationFile"] = configurationFile

		configuration = Core.Environment.Configuration(filename=configurationFile)

		#
		#Check required args
		#

		if not self.opt.has_key("-w") and not self.createConfigurationRequested():
			raise Core.Exceptions.ArgumentError, "Working directory (-w) must be specified"
		
		if not self.opt.has_key("-p") and not self.createConfigurationRequested():	
			raise Core.Exceptions.ArgumentError, "PDB (-p) must be specified"
		
		if not self.opt.has_key("-j"):
			raise Core.Exceptions.ArgumentError, "A job id must be specified"
		
		if not self.opt.has_key("-o"):
			try:
				 self.opt["-o"] = os.path.join(os.getcwd(), self.opt["-j"])
				 self.environment.output('[PEAT-SA] No output directory specified - Defaulting to %s' % self.opt["-o"])
			except:
				Core.Exceptions.ArgumentError, "Output directory (-w) not specified and default value not available in configuration file"
		else:
			self.environment.output('[PEAT-SA] Output directory is %s' % self.opt["-o"])

		#
		#Check if any of the optional flags were specified	
		#
			
		counts = [self.opt.keys().count(option) for option in self.optionalArguments]
		value = reduce(lambda x,y: x+y, counts)
		if value == 0 and  not self.createConfigurationRequested():
			self.environment.output('[PEAT-SA] Performing pka scan, stability analysis, ENM analysis and transition path finding')
			#Add the options arguments
			for option in self.optionalArguments:
				if option != '--mutationList':
					self.opt[option] = ""	

		#Make any supplied paths absolute
		self.__fixPaths__()	
		
	def jobID(self):
	
		if not self.opt.has_key('-j'):
			return None
		
		return self.opt['-j']
		
	
class WebScript:

	def __init__(self):
	
		#Get environment
		self.environment = Core.Environment.Environment()
	
		#Parse command line
		self.parser = CommandLineParser()
		self.parser.parseCommandLine()
		
		#Read Configuration, connect to database and create Job object.
		self.configuration = Core.Environment.Configuration(filename=self.parser.configurationFile())
		self.connection = None
		self.job = None
		if self.environment.isRoot():
			self.connection = UtilityFunctions.ConnectionFromConfiguration(self.configuration)
			print '[WEBAPP] Connection %s' % self.connection
			print '[WEBAPP] Retrieving job data'
			self.job = Data.Job(jobID=self.parser.jobID(), connection=self.connection)
			if not self.job.exists():
				raise Exceptions.DatabaseRetrievalError, 'Data for job %d is not present in database' % self.job.identification
			else:
				print '[WEBAPP] Obtained job data'
				self.job.setState('Running')
				self.job.setQueueStatusMessage('Your job is running')
		
		self.environment.wait()
		
		#Set the $HOME environment variable for apache to the path to the .WHAT_IF file
		#specifying the path to the whatif executable (used by pKarun.WHATIF.py) 
		os.environ['HOME'] = self.configuration.get('PATHS', 'whatif')
		
		#Create the output directory
		if self.environment.isRoot() and self.parser.outputDirectory() is not None:
			if not os.path.exists(self.parser.outputDirectory()):
				os.mkdir(self.parser.outputDirectory())
		
		self.environment.wait()				

	def main(self):
		
		#Set all job states to waiting
		if self.environment.isRoot():
			calculationStates = self.job.calculationStates()
			for calculation in calculationStates.keys():
				if calculationStates[calculation] == 'Queued':
					self.job.setCalculationState(calculation, 'Waiting')
					
		#If a delta-pka calculation is requested check that the directory exists
		#If it doesn't - raise an error
		if self.parser.scan() is True and not os.path.exists(self.parser.workingDirectory()):
			raise Core.Exceptions.MissingPKADataError, "Delta pKa scan requested but specified working dir does not exist"
				
		#Create the ProteinDesignTool instance - this does all cleaning etc.
		tool = Core.ProteinDesignTool.ProteinDesignTool(self.parser.configurationFile(), 
				workingDirectory=self.parser.workingDirectory(),
				pdbFile=self.parser.pdbFile(), 
				outputDirectory=self.parser.outputDirectory())
		
		#Check if a previously created MutantCollection was specified
		mutantCollection = self.parser.mutants()
		if mutantCollection is None:	
			self.environment.output('[PEAT-SA] No mutants available from previous run.\n')
			#Note we use the ProteinDesignTool instances pdbFile attribute.
			#This points to the cleaned pdbFile the tool operates on.
			#Therefore we must use this to generate the mutants.
			mutantCollection = Core.ProteinDesignTool.CreateMutants(tool.pdbFile, tool.configuration, self.parser, self.environment)			
			
										
		#Check if any mutants were succesfully created
		if len(mutantCollection.mutantFiles()) == 0:
			raise Exceptions.MutantCollectionError, 'None of the specified mutants could be modelled due to excessive clashes'															
		
		#Do a scan if requested
		if self.parser.scan() is True:
			if self.environment.isRoot():
				self.job.setCalculationState('Scan', 'Running')
				#Add PDB file to job data since it wasn't uploaded by the user
				self.job.setStructureFromFile(self.parser.pdbFile())
				
			tool.runScan(mutantCollection=mutantCollection, verbose=self.parser.verbose())
			if self.environment.isRoot():
				self.job.addResults(matrix=tool.dataDirectory.scanResults, name='ScanResults')	
				self.job.setCalculationState('Scan', 'Finished')
			self.environment.wait()
	
		#Get the mutants for the stability and binding calculations if these were requested
		if self.parser.stability() is True:
			if self.environment.isRoot():
				self.job.setCalculationState('Stability', 'Running')
				
			tool.runStabilityCalculation(mutantFiles=mutantCollection.mutantFiles(), verbose=self.parser.verbose())			
			if self.environment.isRoot():
				self.job.addResults(matrix=tool.dataDirectory.stabilityResults, name='StabilityResults')
				#Create graph
				points = len(tool.dataDirectory.stabilityResults.total)
				data = Core.Utilities.CreateBarChartData(matrix=tool.dataDirectory.stabilityResults, column="Total", 
						xlabel="Mutant Index", ylabel="KJ/mol", title="Stability", 
						xticIncrement=int(math.ceil(points/20.0)))
				self.job.addImageData(data=data, name="StabilityImage")
				self.job.setCalculationState('Stability', 'Finished')
			self.environment.wait()
			
		if self.parser.binding() is True:
			if self.environment.isRoot():
				self.job.setCalculationState('Binding', 'Running')
				
			tool.runBindingCalculation(mutantFiles=mutantCollection.mutantFiles(), 
					ligandFile=self.parser.ligandFile(), 
					verbose=self.parser.verbose())
			if self.environment.isRoot():
				self.job.addResults(matrix=tool.dataDirectory.deltaBindingResults, name='BindingResults')
				#Create graph
				points = len(tool.dataDirectory.deltaBindingResults.total)
				data = Core.Utilities.CreateBarChartData(matrix=tool.dataDirectory.deltaBindingResults, column="Total", 
						xlabel="Mutant Index", ylabel="KJ/mol", title="Binding", 
						xticIncrement=int(math.ceil(points/20.0)))
				self.job.addImageData(data=data, name="BindingImage")
				self.job.setCalculationState('Binding', 'Finished')
				
			self.environment.wait()	
			
		#FIXME: Delete temp data??
		if self.environment.isRoot():
			self.job.setState('Finished')
			self.connection.close()
			
		#Necessary in case something goes wrong in above root only section
		#to force all processes to quit	
		self.environment.wait()	
			
		
if __name__ == "__main__":

	tool = None
	environment = None
	
	try: 
		environment = Core.Environment.Environment()
		tool = WebScript()
		tool.main()
	except Core.Exceptions.MissingPKADataError, data:
		print "\nError - %s" % data.errorString
		print "Description - %s\n" % data
		if tool is not None and environment.isRoot():
			code = os.path.split(tool.parser.workingDirectory())[1]
			tool.job.setError("No pKa data corresponding to supplied code %s is present" % code, 
				"You must have run a previous pKa calculation on this pdb using the pKD server")
			tool.job.setState("Finished")	
		#In serial this does nothing
		#In parallel if all processes are here they exit
		#In addition if only one process raised the exception and is here, 
		#the others will be waiting elsewhere in the code
		#In this case, passing True here causes all the processes to exit.
		environment.wait(error=True)
		#This is only called if we are in serial
		environment.exit(1)	
	except Core.Exceptions.ProteinDesignToolException, data:
		print "\nError - %s" % data.errorString
		print "Description - %s\n" % data
		if tool is not None and environment.isRoot():
			tool.job.setError(data.errorString,data.message)
			tool.job.setState("Finished")	

		if data.__class__ == Core.Exceptions.ArgumentError:
			print __doc__

		environment.wait(error=True)
		environment.exit(1)
	except SystemExit, data:
		#This is raised by sys.exit - check was the arg 0
		if data.message != 0:
			if tool is not None and environment.isRoot():
				tool.job.setError("Encountered an unexpected error - This is likely a bug", 
						"Program exited with exit code %d" % data.message)
				tool.job.setState("Finished")	
		
			print "\nEncountered an unexpected error - This is likely a bug"
			print "Please contact the developers with the following traceback information",
			print "along with details on the calculation being run (arguments, configuration file etc.)\n"
			print "Program exited with exit code %d" % data.message
	
		environment.wait(error=True)
		environment.exit(1)
	except BaseException, data:
		#Get the traceback as a string = convoluted ....
		trace = traceback.format_list(traceback.extract_tb(sys.exc_info()[2]))
		trace = reduce(lambda x, y: x + y, trace)
		#Need to remove single quotes from the traceback	
		trace = trace.replace("'", " ")
		if tool is not None and environment.isRoot():
			tool.job.setError("Encountered an unexpected error - This is likely a bug", 
					"Check error logs for more information.")
			tool.job.setState("Finished")	
				
		print "\nEncountered an unexpected error - This is likely a bug"
		print "Please contact the developers with the following traceback information",
		print "along with details on the calculation being run (arguments, configuration file etc.)\n"
		print "%s" % data
		environment.wait(error=True)
		environment.exit(1)
			
	
