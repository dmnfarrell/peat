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
'''Creates a job entry in the database using the provided information
The database information must be given in a configuration-file
'''
import optparse, os, sys
import PEAT_SA.WebApp as WebApp
import PEAT_SA.Core as Core

usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-s", "--structure", dest="structureFile",
                  help="The pdb file the job was run on.", metavar="STRUCTURE")
parser.add_option("-d", "--database-configuration",
                   dest="configuration",
                  help="Configuration file containing database information",
                  metavar='CONF')		  		  		  
parser.add_option("-l", "--ligand",
                  dest="ligandFile",
                  help="A mol2 file of a ligand",
                  metavar='LIGAND')
parser.add_option("" "--mutation",
                  dest="mutation",
                  help="The mutation run [optional]",
                  metavar='MUTATION')
parser.add_option("" "--mutationList",
                  dest="mutationList",
                  help="The mutation list used [optional]",
                  metavar='MUTATION')			  		  		  		  		  		  
		  
(options, args) = parser.parse_args()

if options.structureFile is None:
        print 'Job id must be provided'
        sys.exit(1)
else:
	structureName = os.path.split(options.structureFile)[1]
	structureName = os.path.splitext(structureName)[0]	
	
if options.configuration is None:
        print 'A database configuration file must be supplied'
        sys.exit(1)	

if options.mutation is None and options.mutationList is None:
	print 'One of mutation and mutationList must be supplied'
	sys.exit(1)

configuration = Core.Environment.Configuration(filename=options.configuration)
connection = WebApp.UtilityFunctions.ConnectionFromConfiguration(configuration)

jobTable = configuration.get("DATABASE", "jobTable")
dataTable = configuration.get("DATABASE", "dataTable")
jobManager = WebApp.Data.JobManager(connection, jobTable=jobTable)
#Hack: Have to pass something for calculations ....
job = jobManager.createJob(structureName, calculations=['scan'], dataTable=dataTable)

try:
	job.setStructureFromFile(options.structureFile)

	if options.ligandFile is not None:
		job.setLigandFromFile(options.ligandFile)

	if options.mutationList is not None:
		mutationList = Core.Data.MutationListFile(filename=options.mutationList)
		job.setMutationListFile(mutationList)
	else:
		job.setMutation(options.mutation)
		
	print job.identification	
except Exception, data:
	print 'Encountered a %s exception' % Exception
	print 'Reason: %s' % data
	jobManager.deleteJob(job)
	raise

