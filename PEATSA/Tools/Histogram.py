#! /usr/bin/env python
#
#   Project: ProteinDesignTool
#
#   Copyright (C) 2008 Michael Johnston & Jen Nielsen
#
#   Author: Michael Johnston
#
'''Creates a histogram from a column of a csv file'''
import PEATSA.Core.Matrix as Matrix
import PEATSA.Core.Utilities as Utilities
import PEATSA.Core.Environment as Environment
import optparse, sys

usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-m", "--matrix", dest="matrix",
                  help="Csv file to create histogram from.", metavar="MATRIX")
parser.add_option("-c", "--column",
                   dest="column",
                  help="The header of the column to create the histogram from",
                  metavar='COLUMN')
parser.add_option("-b", "--binSize",
                  dest="binSize", default="0.1",
                  help="The size of the histogram bins.\nDefault value %default",
                  metavar='BIN_SIZE')
parser.add_option("-o", "--outputFile",
                  dest="outputFile", default="output.csv",
                  help="The file to write the histogram to.\nDefault value %default",
                  metavar='OUTPUT')
parser.add_option("", "--chart",
                  dest="chart",
                  help="If supplied a bar chart of the histogram will be written to the specified file.",
                  metavar='CHART')		  		  		  
parser.add_option("", "--min",
                  dest="min",
                  help="Sets the minimum value on the x-axis of the histogram. If not supplied this is obtained from the data",
                  metavar='MIN')		  		  		  
parser.add_option("", "--max",
                  dest="max",
                  help="Sets the maximum value on the x-axis of the histogram. If not supplied this is obtained from the data",
                  metavar='MAX')		  		  		  
		  
(options, args) = parser.parse_args()

if options.matrix is None:
        print 'Matrix must be provided'
        sys.exit(1)
		  
if options.column is None:
        print 'Column must be provided'
        sys.exit(1)
	
if options.outputFile is None:
	print 'Output file must be provided'
        sys.exit(1)

if options.max == None:
	maxLimit = None
else:
	maxLimit = float(options.max)
	print 'Setting max limit to %d' % maxLimit

if options.min == None:
	minLimit = None
else:
	minLimit = float(options.min)
					
#Create histogram
matrix = Matrix.matrixFromCSVFile(options.matrix)
matrix = Matrix.PEATSAMatrix(rows=matrix.matrix, headers=matrix.columnHeaders())

print 'Creating histogram from column %s of matrix %s - %s data points' % (options.column, options.matrix, matrix.numberOfRows())
print minLimit, maxLimit
histogram = matrix.createHistogram(column=options.column, binSize=float(options.binSize), minLimit=minLimit, maxLimit=maxLimit, verbose=True)

stream = open(options.outputFile , 'w+')
stream.write(histogram.csvRepresentation())
stream.close()

print 'Histogram has %d rows' % histogram.numberOfRows()
if options.chart != None:
	conf = Environment.Configuration()
	Utilities.CreateBarChart(matrix=histogram, column='count', outputFile=options.chart, 
				xaxis=options.column, title=options.matrix, xlabel=options.column,
				ylabel='count', xticIncrement=0.01)

