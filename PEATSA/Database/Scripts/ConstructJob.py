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
'''Exports all the data corresponding to a job in PEAT-SA webserver database and creates all files neccessary to run that job

Note: Requires a correct webApplication.conf file to be present in PEATSA/WebApp/Resources/'''
import math,sys, optparse, os
import PEATSA.WebApp as WebApp
import PEATSA.Core as Core
import PEATSA.Core.PEATSAParallel as Parallel

def ConstructJob(jobId, outputDir, configuration, connection, jobConfigurationFile=None):

	job = WebApp.Data.Job(jobID=jobId, connection=connection)

	if outputDir is None:
		print 'No output directory specified - default to current directory'
		outputDir = os.getcwd()
	elif not os.path.exists(outputDir):
		os.mkdir(outputDir)

	outputDir = os.path.abspath(outputDir)

	#
	# In the following we
	# 1. retrieve the data from the db
	# 2. write the data to files (if necessary)
	# 3. construct the PEATsSA command line arguments
	#

	#Get the mutation data
	if job.isMutationList() is True:
		mutationList = job.mutationListFile()
		outputFile = os.path.join(outputDir, 'mutationList')
		mutationList.writeToFile(outputFile)
		mutationCommand = '--mutationList=%s' % outputFile
	else:
		mutationCommand = '--mutation=%s' % job.mutation()

	#Get the protein structure and write it out
	structureString = job.structure()
	filename = os.path.join(outputDir, 'protein.pdb')
	stream = open(filename, 'w+')
	stream.write(structureString)
	stream.close
	structureCommand = '-p %s' % filename

	#Get ligand structure and write it out (if present)
	bindingCommand = ""
	if job.ligand() is not None: 
		structureString = job.ligand()
		filename = os.path.join(outputDir, 'ligand.mol2')
		stream = open(filename, 'w+')
		stream.write(structureString)
		stream.close
		bindingCommand = '--ligand=%s' % filename

	#Get information on the requested calculations
	calculationStates = job.calculationStates()	
	calculationString = ""
	for calculation in calculationStates.keys():
		state = calculationStates[calculation]
		#Skip binding because it was handled when getting ligand data
		if state != 'NotSelected' and calculation != 'Binding':	
			print state, calculation
			calculationString = calculationString + '--%s ' % calculation.lower()


	#Create a list containing all the PEAT-SA command line args
	args = [structureCommand, '-j', jobId, '-w', os.path.join(outputDir, 'Work'), calculationString, bindingCommand, mutationCommand, '-v']
	if jobConfigurationFile != None:
		args.append('--configurationFile=%s' % jobConfigurationFile)

	#Create a pbs file and write it out
	pbsScript = Parallel.PBSScriptFromConfiguration(configuration, args, runDirectory=outputDir, logDirectory=outputDir)
	pbsFile = os.path.join(outputDir, 'PEATSA.pbs')
	stream = open(pbsFile, 'w+')
	stream.write(pbsScript)
	stream.close()

	return pbsFile

if __name__ == '__main__':

	usage = "usage: %prog [options]"

	parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
	parser.add_option("-j", "--job-id", dest="jobID",
			  help="A valid job id.", metavar="JOBID")
	parser.add_option("-d", "--directory", dest="dir",
			  help="Where to write to job data.", metavar="DIR")
			  
	(options, args) = parser.parse_args()

	if options.jobID is None:
		print 'Job id must be provided'
		sys.exit(1)

	outputDir = options.dir

	#Connect to the database as get the job
	connection = WebApp.UtilityFunctions.ConnectionFromDefaultConfiguration()
	configuration = WebApp.UtilityFunctions.DefaultWebAppConfiguration()
	ConstructJob(options.jobID, outputDir, configuration, connection)
	connection.close()
