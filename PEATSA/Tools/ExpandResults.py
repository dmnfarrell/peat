#!/usr/bin/env python
#
#  ExpandResults.py
#  ProteinDesignTool
#
#  Created by Michael Johnston on 25/02/2009.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
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
