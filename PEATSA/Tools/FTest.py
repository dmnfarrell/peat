'''Performs an FTest for two fits performed using different numbers of parameters to the 
same data based on chisquared fit values (which are supplied to this script).
Note: The fit with least parameters must be passed as the first matrix.
The chisquared data for each fit must be in separate matrices - 
this enables calculations of multiple ftests with one call'''
#!/bin/env python
#
#  FTest.py
#  ProteinDesignTool
#
#  Created by Michael Johnston on 21/06/2010.
#  Copyright (c) 2010 UCD. All rights reserved.
#
from PEATSA import Core
import sys, optparse
from scipy.stats import distributions


usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-a", "--matrixOne",
		dest="matrixOne",
		help="A csv file with a column containing chisquared data for a least-squares fit",
		metavar='MATRIX-ONE')
parser.add_option("-b", "--matrixTwo",
		  dest="matrixTwo",
		  help="A csv file with a column containing chisquared data for a least-squares fit",
		  metavar='MATRIX-TWO')
parser.add_option("-n", "--numberSamples",
		  type="int",
		  dest="numberSamples",
		  help="The number of samples used for the fits (must be the same in both cases)",
		  metavar='NUMBER-SAMPLES')
parser.add_option("-p1", "--parameters1",
		  type="int",
		  dest="parametersOne",
		  help="The number of parameters used for the first fit",
		  metavar='PARAMETERS-ONE')
parser.add_option("-p1", "--parameters2",
		  type="int",
		  dest="parametersTwo",
		  help="The number of parameters used in the second fit",
		  metavar='PARAMETERS-TWO')		
parser.add_option("-c", "--columnName",
		  dest="columnName",
		  default="Chisquared",
		  help="The name of the column containing the chisquared values. Defauts to %default",
		  metavar='CHI-COLUMN')		  	  		  
parser.add_option("", "--row", type="int",
		  dest="row",
		  help="The row in each file to use for the test")	


(options, args) = parser.parse_args()

if options.matrixOne is None or options.matrixTwo is None:
	print 'Two matrices must be specified'
	sys.exit(1)
	
if options.numberSamples is None:
	print 'Number of samples used must be specified'
	sys.exit(1)
	
if options.parametersOne is None or options.parameters2 is None:
	print 'Number of parameters for both fits must be specified'
	sys.exit(1)
	
(options, args) = parser.parse_args()	
	
m1 = Core.Matrix.matrixFromCSVFile(options.matrixOne)
m2 = Core.Matrix.matrixFromCSVFile(options.matrixTwo)

data1 = m1.columnWithHeader(options.columnName)
data2 = m2.columnWithHeader(options.columnName)

pDiff = options.parametersTwo - options.parametersOne
dof = options.numberSamples - options.parametersTwo

p95 = 1 - distributions.f.ppf(0.95, pDiff, dof)
p99 = 1 - distributions.f.ppf(0.99, pDiff, dof)

print '95% critical value is %lf' % p95
print '99% critical value is %lf' % p99

rows = []
if options.row is None:
	pairs = zip(data1, data2)
	for pair in pairs:
		fValue = (pair[0] - pair[1])/pDiff
		fValue /= pair[1]/dof
		pValue = 1 - distributions.f.cdf(fValue, pDiff, dof)
		rows.append([fValue, pValue])
		
results = Core.Matrix.Matrix(rows=rows, headers=['F-Value', 'P-Value'])	
		
		
	





				







