#! /bin/env python
'''Reweights the terms of an uffbaps calculation and recalcualtes the totals.

The components of the UFFBAPS calculation are assumed to be in contiguous columns in the provided file.'''
from PEATSA import Core
import sys, operator, optparse

usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-f", "--file", dest="file",
		  help="A csv file containing unweighted UFFBAPS calculation results", metavar="FILE")
parser.add_option("-w", "--weightsFile", dest="weightsFile",
		  help="A csv file each row of which contains weights", metavar="WEIGHTS")
parser.add_option("-s", "--start",
		  type="int",
		  default=0,
		  dest="start",
		  help="The column in the file containing the first component to apply weights to",
		  metavar='START')		  		  		  		  		  		  
parser.add_option("-e", "--end",
		  type="int",
		  dest="end",
		  help="The column in the file containing the last component to apply weights to",
		  metavar='START')		  		  		  		  		  		  
parser.add_option("-t", "--totals",
		  type="int",
		  dest="totalIndex",
		  help="The column in the file containing the unweighted total",
		  metavar='TOTALS')		  		  		  		  		  		  
parser.add_option("-r", "--weightsRow",
		  type="int",
		  default=0,
		  dest="rowIndex",
		  help="The row in the weights file containing the weights to use",
		  metavar='ROW')		  		  		  		  		  		  

(options, args) = parser.parse_args()

if not sys.stdin.isatty():
	csv = sys.stdin.read()
	m = Core.Matrix.matrixFromCSVRepresentation(csv)
elif options.file is None:
	print 'CSV file must be provided'
	sys.exit(1)
else:
	m = Core.Matrix.matrixFromCSVFile(options.file)

if options.weightsFile is None:
	print 'Weights file not specified'
	sys.exit(1)

if options.totalIndex is None:
	options.totalIndex = options.start + 6
	print >>sys.stderr, 'Index of total column not specified - defaulting to', options.totalIndex

if options.end is None:
	options.end = options.start + 5
	print >>sys.stderr, 'Index of column containing last component not specified - defaulting to', options.end
else:
	#Offset
	options.end = options.end + 1

weightsMatrix = Core.Matrix.matrixFromCSVFile(options.weightsFile)
weights = weightsMatrix.row(options.rowIndex)

#Check if the number of weights corresponds to the number of components
if len(weights) != options.end - options.start:
	string  = 'Number of weights, %d, does not equal number of components, %d' % (len(weights), options.end - options.start)
	print >>sys.stderr, string
	sys.exit(1)

rows = []
for row in m:
	row = list(row)
	components = row[options.start:options.end]
	weighted = map(lambda x,y: x*y, components, weights)
	row[options.start:options.end] = weighted
	total = reduce(operator.add, weighted)
	row[options.totalIndex] = total
	rows.append(row)

results = Core.Matrix.Matrix(rows=rows, headers=m.columnHeaders())
print results.csvRepresentation(),
