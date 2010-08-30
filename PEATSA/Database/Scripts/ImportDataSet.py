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
'''Script for adding a data-set to the database.
Note: The data-set must be associated with an existing job'''
import math,sys, optparse
import PEATSA.WebApp as WebApp
import PEATSA.Core as Core

usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-j", "--job-id", dest="jobID",
                  help="A valid job id.", metavar="JOBID")
parser.add_option("-d", "--data",
                   dest="data",
                  help="The peatsa output data directory",
                  metavar='DATA')		  		  		  
parser.add_option("-c", "--configuration",
                   dest="configurationFile",
                  help="Configuration file containing database information",
                  metavar='CONF')		  		  		  
		  
(options, args) = parser.parse_args()

if options.jobID is None:
        print 'Job id must be provided'
        sys.exit(1)
	
if options.data is None:
        print 'Valid data directory must be provided'
        sys.exit(1)	

if options.configurationFile is None:
        print 'A configuration file must  be provided'
        sys.exit(1)	
		  
data = Core.Data.DataSet(options.data)
configuration = Core.Environment.Configuration(options.configurationFile)
connection = WebApp.UtilityFunctions.ConnectionFromConfiguration(configuration)

jobTable = configuration.get("DATABASE", "jobTable")
dataTable = configuration.get("DATABASE", "dataTable")
job = WebApp.Data.Job(jobID=options.jobID, connection=connection, jobTable=jobTable, dataTable=dataTable)
job.addDataSet(data)
connection.close()
