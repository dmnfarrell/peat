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
'''Script for sending PEAT-SA webserver data to the database from the command line.
The data is a csv file containing the results of a calculation'''
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
parser.add_option("-c", "--calculation",
                  dest="calculation",
                  help="The calculation data to send",
                  metavar='CALC')		  		  		  		  		  		  
		  
(options, args) = parser.parse_args()

if options.jobID is None:
        print 'Job id must be provided'
        sys.exit(1)
	
if options.data is None:
        print 'Valid data directory must be provided'
        sys.exit(1)	
		  
calculations = ['Binding', 'Scan', 'Stability']		  
try:
	calculations.index(options.calculation)
except ValueError, data:	
	print '%s is an invalid calculation option. Valid options %s' % (options.calculation, calculations)
	sys.exit(1)

connection = WebApp.UtilityFunctions.ConnectionFromDefaultConfiguration()
job = WebApp.Data.Job(jobID=options.jobID, connection=connection)

data = Core.Data.DataSet(options.data)

if options.calculation == 'Stability':
	matrix = data.stabilityResults
elif options.calculation == 'Binding':
	matrix = data.bindingResults
elif options.calculation == 'Scan':
	matrix = data.scanResults
	
if matrix == None:
	print 'No data for calculation %s present in data directory %s' % (options.calculation, options.data)					

job.addDataSet(matrix=matrix, name='%sResults' % options.calculation)

if options.calculation != 'Scan':
	points = len(matrix.total)
	chart = Core.Utilities.CreateBarChartData(matrix=matrix, column="Total",
			xlabel="Mutant Index", ylabel="KJ/mol", title=options.calculation,
			xticIncrement=int(math.ceil(points/20.0)))
	job.addImageData(data=chart, name='%s.Image' % options.calculation)

job.setCalculationState(options.calculation, 'Finished')
connection.close()
