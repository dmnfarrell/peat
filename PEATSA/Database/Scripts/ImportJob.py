#!/usr/bin/env python
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
'''Imports the results of an existing peatsa calculation to the database
'''
import optparse, os, sys
import PEATSA.WebApp as WebApp
import PEATSA.Core as Core

usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-r", "--results", dest="resultsDir",
                  help="The calculations *.peatsa directory.", metavar="RESULTS")
parser.add_option("-m", "--mutantCollection",
                  dest="mutantDir",
                  help="The calculations *.mutants directory",
                  metavar='MUTANTS')		  
parser.add_option("-d", "--database-configuration",
                   dest="configuration",
                  help="Configuration file containing database information",
                  metavar='CONF')		  		  		  			  		  		  		  		  		  
		  
(options, args) = parser.parse_args()

if options.resultsDir is None:
        print 'The calcultions results output directory must be provided'
        sys.exit(1)
else:
	options.resultsDir = os.path.abspath(options.resultsDir)
	components = os.path.split(options.resultsDir)
	data = Core.Data.DataSet(name=components[1], location=components[0])
	name = os.path.splitext(data.name)[0]
			
if options.mutantDir is None:
	print 'The calculations mutants output directory must be provided'
	sys.exit(1)
else:
	options.resultsDir = os.path.abspath(options.mutantDir)
	components = os.path.split(options.mutantDir)
	mutantCollection = Core.Data.MutantCollection(name=components[1], location=components[0])	
	
if options.configuration is None:
        print 'A database configuration file must be supplied'
        sys.exit(1)	

configuration = Core.Environment.Configuration(filename=options.configuration)
connection = WebApp.UtilityFunctions.ConnectionFromConfiguration(configuration)

#Hack: Have to pass something for calculations ....
jobTable = configuration.get("DATABASE", "jobTable")
dataTable = configuration.get("DATABASE", "dataTable")
jobManager = WebApp.Data.JobManager(connection, jobTable=jobTable)
job = jobManager.createJob(name, calculations=['scan'], dataTable=dataTable)

try:
	job.setStructureFromFile(mutantCollection.pdbFile)
	mutationList = mutantCollection.mutationListFile('temp')
	job.setMutationListFile(mutationList)
	job.addDataSet(data)

	if len(mutantCollection.ligandFiles()) > 0:
		job.setLigandFromFile(mutantCollection.ligandFiles()[0])

	print job.identification	
	
except Exception, data:
	print 'Encountered a %s exception' % Exception
	print 'Reason: %s' % data
	jobManager.deleteJob(job)
	raise

