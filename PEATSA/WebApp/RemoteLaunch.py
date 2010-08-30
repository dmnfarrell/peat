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
'''Wrapper for launching WebScript.py on a remote cluster via PBS

The arguments to this script must be exactly those that you want to pass to the WebScript.
This script relies on the default webapp configuration file containing a PARALLEL and REMOTE sections.
For the content of the PARALLEL section see Core.PEATSAParallel.
The REMOTE section must contain the following two options
	identityFile - An file containing the private key for logging in to the cluster
	user - The id of the user to log in as.'''

import sys, os, subprocess
import Data
import WebScript, UtilityFunctions
import PEAT_SA.Core.PEATSAParallel as PEATSAParallel

def Copy(source, destination, identityFile, user, host):

	'''Copies a file to a remote host using scp
	
	Parameters
		source: The file to copy
		destination: Where to put it on the remote host
		identityFile: Location of ssh identityFile
		user: User to connect to remote host as (related to identityFile)
		host: The remote host name
		
	Returns:
		0 if successful. Otherwise anyother integer.
		The integer is the return code of scp'''

	runString = "scp -i %s %s %s@%s:%s" % (identityFile, source, user, host, destination) 
	print 'Running %s' % runString
	process = subprocess.Popen(runString, shell=True)
	process.wait()
	
	return process.returncode
	
def Command(command, argumentString, identityFile, user, host):

	'''Runs a command on a remote host using ssh
	
	Parameters
		command: The command to run
		argumentString: The commands arguments
		identityFile: Location of ssh identityFile
		user: User to connect to remote host as (related to identityFile)
		host: The remote host name
		
	Returns:
		0 if successful. Otherwise anyother integer.
		The integer is the return code of ssh'''

	runString = "ssh -i %s -l %s %s %s %s" % (identityFile, user, host, command, argumentString)
	print 'Running %s' % runString
	process = subprocess.Popen(runString, shell=True)
	process.wait()
	
	return process.returncode

def main(parser, options, job):

	'''Launches the job on a remote cluster'''
	
	#Get information for running the job remotely
	identityFile = options.get('CLUSTER', 'identityFile')
	user = options.get('CLUSTER', 'user')
	host = options.get('CLUSTER', 'host')
	webServerDirectory =  options.get('WEB APPLICATION', 'webServerDirectory')
	outputDirectory = parser.outputDirectory()	
		
	#Create PBS file
	pbsScript = PEATSAParallel.PBSScriptFromConfiguration(options, arguments=sys.argv[1:], 
					 runDirectory=outputDirectory, logDirectory=outputDirectory)
	
	#Write the PBS file to the web server directory for the job	
	filename = os.path.join(webServerDirectory, "%s.pbs" % parser.jobID())
	stream = open(filename, "w+")
	stream.write(pbsScript)
	stream.close()
	
	#Create the output directory for the job on the remote machine
	retval = Command(command="mkdir", 
			argumentString=outputDirectory, 
			identityFile=identityFile, 
			user=user, 
			host=host)
	
	if(retval != 0):
		job.setError(description="Error during remote launch",
			detailedDescription="Unable to create job output directory")
		sys.exit(retval)	
		
	#Copy the PBS script to the run dir
	Copy(source=filename, 
		destination=outputDirectory, 
		identityFile=identityFile, 
		user=user, 
		host=host) 
		
	if(retval != 0):
		job.setError(description="Error during remote launch",
			detailedDescription="Unable to copy PBS script to output directory")
		sys.exit(retval)	
	
	#Copy the pdb file to the run dir
	#Split because pdbFile is full file name
	#Note: A pdb file is not present if scan was selected
	if not parser.scan():
		pdbCode = os.path.split(parser.pdbFile())[1]
		localPdbFile = os.path.join(webServerDirectory, parser.jobID())
		localPdbFile = os.path.join(localPdbFile, pdbCode)
		Copy(source=localPdbFile, 
			destination=outputDirectory, 
			identityFile=identityFile, 
			user=user, 
			host=host) 
		
	if(retval != 0):
		job.setError(description="Error during remote launch",
			detailedDescription="Unable to copy pdb to output directory")
		sys.exit(retval)	
	
	#If there is a ligand copy it to the run dir
	if parser.ligandFile() is not None:
		localFile = os.path.join(webServerDirectory, parser.jobID())
		localFile = os.path.join(localFile, os.path.split(parser.ligandFile())[1])
		Copy(source=localFile, 
			destination=outputDirectory, 
			identityFile=identityFile, 
			user=user, 
			host=host)
			
		if(retval != 0):
			job.setError(description="Error during remote launch",
				detailedDescription="Unable to copy ligand to output directory")
			sys.exit(retval)	
			
	#If there is a mutation list copy it to the run dir
	if parser.mutationListFile() is not None:
		localFile = os.path.join(webServerDirectory, parser.jobID())
		localFile = os.path.join(localFile, os.path.split(parser.mutationListFile())[1])
		Copy(source=localFile, 
			destination=outputDirectory, 
			identityFile=identityFile, 
			user=user, 
			host=host)
			
		if(retval != 0):
			job.setError(description="Error during remote launch",
				detailedDescription="Unable to copy mutation file to output directory")
			sys.exit(retval)	

	#Submit the job
	Command(command="qsub", 
		argumentString=os.path.join(outputDirectory, "%s.pbs" % parser.jobID()), 
		identityFile=identityFile, 
		user=user, 
		host=host)
		
	if(retval != 0):
		job.setError(description="Error during remote launch",
			detailedDescription="Failed to submit job to queue")
		sys.exit(retval)

if __name__ == "__main__":
		
	connection = None
	job = None

	try: 
		#Get WebServer configuration
		options = UtilityFunctions.DefaultWebAppConfiguration()
		
		#Parse the command line for more info
		parser = WebScript.CommandLineParser()
		parser.parseCommandLine()
		
		#Connect to the database in case there is an error
		connection = UtilityFunctions.ConnectionFromDefaultConfiguration()
		job = Data.Job(jobID=parser.jobID(), connection=connection)
		
		#call main
		print 'Setting up remote launch for job %s' % parser.jobID()
		main(parser, options, job)
		print 'Job sent'
		
	except SystemExit, data:
		connection.close()
		#This is raised by sys.exit 
		#Should have already been logged by the script so ignore it
		print "\nEncountered an error with remote launch"
		print "Check database entry for job for more info"
		raise
	except BaseException, data:
		#Catch any other exception here
		if job is not None:
			job.setError(description="Encountered an unexpected exception - This is likely a bug", 
					detailedDescription="%s" % data)
		raise
		
	if connection is not None:
		connection.close()
	sys.exit(0)
		
				
