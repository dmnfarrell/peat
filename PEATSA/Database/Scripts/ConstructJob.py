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
'''Exports all the data corresponding to a job in PEATSA webserver database and creates all files neccessary to run that job

Note: Requires a correct webApplication.conf file to be present in PEATSA/WebApp/Resources/'''
import math,sys, optparse, os
import PEATSA.WebApp as WebApp
import PEATSA.Core as Core
import PEATSA.Core.PEATSAParallel as Parallel

#FIXME - Construct DeltapKA calculations
def ConstructJob(jobId, outputDir, configuration, connection, jobConfigurationFile=None):

	job = WebApp.Data.Job(jobID=jobId, connection=connection)

	if outputDir is None:
		print 'No output directory specified - default to current directory'
		outputDir = os.getcwd()
	elif not os.path.exists(outputDir):
		os.mkdir(outputDir)

	outputDir = os.path.abspath(outputDir)
	
	args = ['-j', jobId, '-w', os.path.join(outputDir, 'Work'), '-v']
	
	#Add file args
	fileArgs = job.fileArguments()
	count = 0
	for option in fileArgs.keys():
		contentsAttribute = fileArgs[option]['contentsAttribute']
		if contentsAttribute is None:
			print 'Error - no contents attribute specified for file option %s' % option
			print 'Skipping'
			continue
	
		fileName = fileArgs[option]['fileName']
		if fileName is None:
			fileName = 'file%d.txt' % count
		
		outputFile = os.path.join(outputDir, fileName)		
					
		contents = getattr(job, contentsAttribute)()		
		f = open(outputFile, 'w+')
		f.write(contents)
		f.close()
		
		if option[:2] == '--':
			command = option+'='+outputFile
		elif option[:1] == '-':
			command = option+' '+outputFile			
							
		args.append(command)
																
		count += 1

	#Add option args
	optionArgs = job.optionArguments()
	for key in optionArgs.keys():
		if key[:2] == '--':
			if optionArgs[key] == "":
				arg = [key]
			else:
				arg = [key+'='+optionArgs[key]]
		elif key[:1] == '-':	
			arg = [key, optionArgs[key]]
		else:
			print 'Option %s not valid - must be long or short option' % key
			arg = []	
				
		args.extend(arg)		
	
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
	parser.add_option("-c", "--configurationFile", dest="configurationFile",
			  help="Configuration file containing DATABASE section", metavar="CONF")
			  
	(options, args) = parser.parse_args()

	if options.jobID is None:
		print 'Job id must be provided'
		sys.exit(1)

	if options.configurationFile is None:
		print 'PEAT-SA configuration file must be provided'
		sys.exit(1)
	else:
		configurationFile = os.path.abspath(options.configurationFile)
		print 'Using configuration file ', configurationFile

	outputDir = options.dir

	#Connect to the database as get the job
	configuration = Core.Environment.Configuration(filename=configurationFile)
	connection = WebApp.UtilityFunctions.ConnectionFromConfiguration(configuration)
	ConstructJob(options.jobID, outputDir, configuration, connection)
	connection.close()
