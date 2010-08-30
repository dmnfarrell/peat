#! /bin/env python
'''Creates a sub-matrix from selected columns of a csv file.

The cvs file can be passed to stdin'''
from PEAT_SA import Core
import sys, optparse

usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-f", "--file", dest="file",
                  help="A csv file", metavar="FILE")
parser.add_option("-c", "--columns",
                  dest="columns",
                  help="A comma separated list of column indices",
                  metavar='COLUMN')		  		  		  		  		  		  

(options, args) = parser.parse_args()

if not sys.stdin.isatty():
	csv = sys.stdin.read()
	m = Core.Matrix.matrixFromCSVRepresentation(csv)
elif options.file is None:
	print 'CSV file must be provided'
	sys.exit(1)
else:
	m = Core.Matrix.matrixFromCSVFile(options.file)

if options.columns is None:
	print 'Column indexes must be provided'
	sys.exit(1)
else:
	columns = [int(value) for value in options.columns.split(',')]

subm = m[:, columns[0]:columns[0] + 1]
for index in columns[1:]:
	subm.addColumn(m.column(index))
	subm.setColumnHeader(subm.numberOfColumns() -1, m.headerForColumn(index))

print subm.csvRepresentation(),
