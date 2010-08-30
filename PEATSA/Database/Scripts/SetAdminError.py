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
'''Sets an error for a PEAT-SA webserver job from the command line.

Requires that a valid webApplication.conf file is present.'''
import PEAT_SA.WebApp.Data
import PEAT_SA.WebApp.UtilityFunctions
import optparse, sys

usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-j", "--job-id", dest="jobID",
                  help="A valid job id.", metavar="JOBID")
parser.add_option("-e", "--error",
                   dest="error",
                  help="A short error description",
                  metavar='ERROR')		  		  		  
parser.add_option("-d", "--description",
                  dest="description",
                  help="Description of the error",
                  metavar='MIN')		  		  		  		  		  		  
		  
(options, args) = parser.parse_args()

if options.jobID is None:
        print 'Job id must be provided'
        sys.exit(1)
		  
if options.error is None:
        print 'Error must be provided'
        sys.exit(1)
	
if options.description is None:
	print 'Description  must be provided'
        sys.exit(1)


connection = PEAT_SA.WebApp.UtilityFunctions.ConnectionFromDefaultConfiguration()
job = PEAT_SA.WebApp.Data.Job(options.jobID, connection)
job.setError(description=options.error, detailedDescription=options.description) 


