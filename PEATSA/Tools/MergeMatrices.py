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
'''Merges rows corresponding to the same key e.g. PDB name, in multiple csv files.

A key is created for each row in each file using the data in one or more of the matrix columns.
Rows with the same key in all supplied csv files are merged into one row.
The resulting file contains only entries whose key is present in all csv files.
Note each rows key should be unique within that matrix.'''
import sys, optparse
from PEATSA import Core

usage = "usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="% 0.1", description=__doc__)
parser.add_option("-f", "--files", dest="file",
                  help="A comma separated list of csv files", metavar="FILE")
parser.add_option("-k", "--key-columns", dest="keyColumns", default='Mutations',
                  help="A comma separated list of column headers common to all csv files.\n"
		  "Each rows key is created by concatentating the data in these columns.\nDefaults to %default", metavar="COLUMNS")
parser.add_option("-o", "--outputFile", dest="outputFile", default='Merged.csv',
                  help="Where to write the results. Defaults to %default", metavar="OUTPUT")
parser.add_option("-v", "--verbose", action='store_true', default=False, dest="verbose",
                  help="Output more info.\n", metavar="VERBOSE")

(options, args) = parser.parse_args()

if options.file is None:
	print 'CSV file must be provided'
	sys.exit(1)
else:
	files = [file.strip() for file in options.file.split(',')]
	print 'Reading matrices ...'
	matrices = [Core.Matrix.matrixFromCSVFile(file) for file in files]

keyColumns = [header.strip() for header in options.keyColumns.split(',')]

print 'Checking key-columns are present in all matrices'
#Check all matrices contain the keyColumns
for matrix in matrices:
	for header in keyColumns:
		try:
			matrix.indexOfColumnWithHeader(header)
		except ValueError:
			print 'Column %s not present in all matrices' % header	
			sys.exit(1)

#Create the id columns for each matrix

print 'Creating row ids for each matrix'
id = []
ids = []
counter = 0
for matrix in matrices:
	print '\tMatrix %s' % files[counter]
	columns = [matrix.columnWithHeader(header) for header in keyColumns]
	for j in range(matrix.numberOfRows()):
		for i in range(len(keyColumns)):
			id.append(str(columns[i][j]))
		
		ids.append("+".join(id))
		id = []

	if options.verbose is True:
		print '\tKey-Code Example: %s' % ids[0]

	print '\tHashing key-column'
	matrix.addColumn(ids)
	matrix.setColumnHeader(-1, 'Key')
	matrix.setKeyColumn('Key')
	print '\tDone'
	ids = []
	counter += 1

keys = matrices[0].keys()

#mutants = [Core.Data.MutationSet(code) for code in matrices[0].mutations]
print 'Combining rows ...'
combinedData = []
i=0
row = []
skipped = 0
matched = 0
for key in keys:
	i = 0
	for matrix in matrices:
		try:
			data = matrix.rowsForKey(key)
			if len(data) > 1:
				print 'WARNING: Multiple rows in matrix %s match key %s - only using first' % (files[i], key)
		except KeyError: 
			if options.verbose is True:
				print 'No rows matching %s in matrix %s - Skipping all entries with this key' % (key, files[i])
			row = []
			break
		else:
			row.extend(data[0])

		i += 1
	
	if len(row) != 0:
		matched += 1
		combinedData.append(row)
		row = []
	else:
		skipped += 1

print 'Matched %d rows. Skipped %d rows' % (matched, skipped)

#Create combined headers array
headers = []
for matrix in matrices:
	headers.extend(matrix.columnHeaders())

print 'Creating results'
results = Core.Matrix.PEATSAMatrix(rows=combinedData, headers=headers)
	
print 'Removing redundant columns'
#Remove redundant key columns
indexes = []
keyColumns.append('Key')
for column in keyColumns:
	data = [m.indexOfColumnWithHeader(column) for m in matrices]
	for i in range(1, len(matrices)):
		offset = matrices[i-1].numberOfColumns()
		data[i] += offset

	indexes.extend(data)

indexes.sort()
keyColumns.remove('Key')

indexes = indexes[len(keyColumns):]
indexes.reverse()
for index in indexes:
	results.removeColumn(index)

print 'Writing results to %s' % options.outputFile
f = open(options.outputFile , 'w+')
f.write(results.csvRepresentation())
f.close()
print 'Merge complete'

