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
'''Creates a new matrix from a PEAT-SA results file where the information in the mutation code column is expanded across four new columns

Usage: 
	ExpandResults.py [CSVFile]

Description:

Each mutation code in the original matrix must refer to a single point mutation.

For example an entry A54AG in the old matrix would be expanded to  A, 54, Ala, Gly in the new matrix.
Each of these four components would be in separate columns called 'Chain', 'Residue Index', 'Residue Name' and 'Mutation respectively'''

import sys, os
import PEATSA.Core as Core

fileName = sys.argv[1]

#Read in the csv file
matrix = Core.Matrix.matrixFromCSVFile(fileName)
#Expand all the SPM codes - assumes there in a column called Mutations.
expandedCodes = [Core.Utilities.ParseMutationCode(code) for code in matrix.mutations]
#Create a new matrix from the expanded data
newMatrix = Core.Matrix.Matrix(expandedCodes, ['Chain', 'Residue Index', 'Residue Name', 'Mutation'])

#Append all the columns in the original matrix to the new one
for i in range(matrix.numberOfColumns()):
	newMatrix.addColumn(matrix.column(i))
	newMatrix.setColumnHeader(newMatrix.numberOfColumns() - 1, matrix.headerForColumn(i))
	
#Write the new matrix out	
stream = open('%sExtended.csv' % os.path.splitext(fileName)[0], 'w+')
stream.write(newMatrix.csvRepresentation())
stream.close()	
