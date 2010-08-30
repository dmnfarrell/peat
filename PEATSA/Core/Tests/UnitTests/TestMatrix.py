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
'''Contains unittests for the ProteinDesignTool.Matrix class'''

import sys, os
sys.path.append(".")
sys.path.insert(0, "../../")

import Matrix
import unittest 

#Simple data for testing.
lists = [[1,2], ["michael", 34.66]]
headers = ("Header1", "Header2")
testCSVString = "Header1	Header2\n1	2\nmichael	34.66\n"
testCSVString2 = "1	2\nmichael	34.66\n"

class TestMatrix(unittest.TestCase):

	'''Test class for Matrix'''

	def setUp(self):

		self.matrix = Matrix.Matrix(rows=lists, headers=headers)
	
	def testDefautInitialisation(self):
	
		'''Tests that instantiating the matrix with no arguments works'''
		
		matrix = Matrix.Matrix()

	def testSetupInitInitialisation(self):
	
		'''Tests that the rows and columns are correctly set'''
		
		self.assertEqual(2, self.matrix.numberOfRows())
		self.assertEqual(2, self.matrix.numberOfColumns())

	def testInitialisationWithoutHeaders(self):
	
		'''Test that the matrix is initialised correctly if no headers are supplied'''
	
		matrix = Matrix.Matrix(rows=lists)
		self.assertEqual(2, matrix.numberOfRows())
		self.assertEqual(2, matrix.numberOfColumns())
		self.assertEqual(("None", "None"), matrix.columnHeaders())

	def testIncorrectInitialisation(self):
	
		'''Test that initialisation fails if the length of headers is incompatible with number of columns provided'''
		
		try:
			Matrix.Matrix(headers=headers)
			self.fail("No exception raised")
		except IndexError:
			self.assert_(True)
		except:	
			self.fail("Unexpected exception raised")

	def testSetName(self):
		
		'''Tests that setting the matrices name works'''	
	
		self.matrix.setName("someName")
		self.assertEqual(self.matrix.name(), "someName")
		
	def testSetInvalidName(self):
	
		'''Test that the matrix name can only be a string'''
	
		self.assertRaises(TypeError, self.matrix.setName, [])

	def testGetElement(self):

		'''Tests that element method works'''
		
		element = self.matrix.element(1,1)
		self.assertEqual(element,34.66)
		
		self.assertRaises(IndexError, self.matrix.element, 2,2)

	def testRowAddition(self):
	
		'''Tests adding a new row
		
		Tests default addition and addition to a specified index.
		Tests that the correct error is raised if index is out of bounds.'''
	
		#Default add
		newRow = ["Adding", "Row"]
		self.matrix.addRow(newRow)
		row = self.matrix.row(2)
		self.assertEqual(row, tuple(newRow))
		self.assertEqual(3, self.matrix.numberOfRows())
				
		#Specified index add
		self.matrix.addRow(newRow, index=1)
		row = self.matrix.row(1)
		self.assertEqual(row, tuple(newRow))
		
		#Test raises an error if index is incorrect
		self.assertRaises(IndexError, self.matrix.addRow, row=newRow, index=10)
		
		#Test raises an error if number of elements in row is incorrect
		self.assertRaises(IndexError, self.matrix.addRow, row=['Wrong'])
		
	def testColumnAddition(self):
	
		'''Tests adding a new column
		
		Tests default addition and addition to a specified index.
		Tests that the correct error is raised if index is out of bounds.
		Tests that the correct headers are set for the new columns'''
	
		#Default add
		newColumn = ["Adding", "Column"]
		self.matrix.addColumn(newColumn)
		column = self.matrix.column(2)
		self.assertEqual(column, tuple(newColumn))
		self.assertEqual(3, self.matrix.numberOfColumns())
		
		#Check column header is correct
		header = self.matrix.headerForColumn(2)
		self.assertEqual("None", header)
		
		#Specified index add
		self.matrix.addColumn(newColumn, index=1)
		column = self.matrix.column(1)
		self.assertEqual(column, tuple(newColumn))
		
		#Check column header is correct
		header = self.matrix.headerForColumn(1)
		self.assertEqual("None", header)
		
		#Test range
		self.assertRaises(IndexError, self.matrix.addColumn, newColumn, 10)
	
	def testGetRow(self):
	
		'''Tests the correct row is returned.
		Tests that the correct error is raised if the requested index is out of bounds'''
	
		#Test getting row
		row = self.matrix.row(1)
		self.assertEqual(("michael", 34.66), row)
			
		#Check its a tuple
		self.assertEqual(type(row), type(()))	
		
		#Test fail
		self.assertRaises(IndexError, self.matrix.row, 4)
		
	def testGetColumn(self):
	
		'''Tests the correct column is returned.
		Tests that the correct error is raised if the requested index is out of bounds'''
	
		#Test getting column
		column = self.matrix.column(1)
		self.assertEqual( (2, 34.66), column)	
		
		#Test fail
		self.assertRaises(IndexError, self.matrix.column, 4)
		
		#Check its a tuple
		self.assertEqual(type(column), type(()))
	
	def testIteration(self):
	
		'''Test that iteration works'''
	
		#Check that the following accumulates all matrix
		#rows in a new list - which should have identical contents
		#to lists global variable except that the entries are tuples.
		rows = [row for row in self.matrix]
		
		for i in range(len(rows)):
			self.assertEqual(rows[i], tuple(lists[i]))
			
	def testCompare(self):
				
		'''Test that matrices with identical data compare as equal'''
		
		newMatrix =  Matrix.Matrix(rows=lists, headers=headers)
		self.assertEqual(self.matrix, newMatrix)
		
		#Not equal test
		
		self.assertNotEqual(self.matrix, None)																
	
	def testGetHeaders(self):				
		
		'''Tests returning of column headers'''
		
		retval = self.matrix.columnHeaders()
		self.assertEqual(retval, headers)	
		
		#Check its a tuple
		self.assertEqual(type(retval), type(()))
			
	def testGetColumnForHeader(self):
	
		'''Tests correct column is returned for specified header
		Tests reponse when a non-existant header is supplied'''
	
		column = self.matrix.columnWithHeader("Header2")
		self.assertEqual( (2, 34.66), column)	
		
		#Test if the header doesn't exist
		column = self.matrix.columnWithHeader("Non-existant header")
		self.assertEqual(None, column)				
	
	def testSetColumnHeaders(self):
	
		'''Tests setting new headers for all columns.'''
	
		newHeaders = ("Col 1", "Col 2")
		self.matrix.setColumnHeaders(newHeaders)
		retval = self.matrix.columnHeaders()
		self.assertEqual(retval, newHeaders)
		
		#Check if fails if the incorrect number of headers are supplied
		newHeaders = ("Col 1", "Col 2", "Extra")
		self.assertRaises(IndexError, self.matrix.setColumnHeaders, newHeaders)		
	
	def testSetColumnHeader(self):
	
		'''Tests setting a column header to a new value'''
		
		newHeader = "New Head"
		self.matrix.setColumnHeader(header=newHeader, index=1)
		header = self.matrix.headerForColumn(1)
		self.assertEqual(newHeader, header)
		
		#Check fail
		self.assertRaises(IndexError, self.matrix.setColumnHeader, header=newHeader, index=10)
	
	def testIndexForColumn(self):
	
		'''Tests getting the index of the column associated with a given header'''
	
		index = self.matrix.indexOfColumnWithHeader("Header2")
		self.assertEqual(index, 1)
		
		self.assertRaises(ValueError, self.matrix.indexOfColumnWithHeader, "Non existant header")
	
	def testReadCSV(self):
	
		'''Tests if a matrix can be properly initialised from a csv file'''
	
		#Test can read csv files produce by the class itself
		string = self.matrix.csvRepresentation()
		file = open('temp', 'w')
		file.write(string)
		file.close()
		
		matrix = Matrix.matrixFromCSVFile('temp')
		
		self.assertEqual(matrix, self.matrix)
		os.remove('temp')
		
		#Test initialisation for a csv file not in the same format
		#as the one created by the Matrix class
		file = open('temp', 'w')
		file.write(testCSVString)
		file.close()
		
		matrix = Matrix.matrixFromCSVFile('temp')
		self.assertEqual(matrix, self.matrix)
		
		#Test it works when told not to read header
		file = open('temp', 'w')
		file.write(testCSVString2)
		file.close()
		
		matrix = Matrix.matrixFromCSVFile(filename='temp', readHeaders=False)
		matrix.setColumnHeaders(["Header1", "Header2"])
		self.assertEqual(matrix, self.matrix)
		
		#Clean up
		os.remove('temp')
		
#	def testReadingCSVWithKnownFormat(self):
#	
#	'''Test reading a csv file which has a defined format'''
#
#	self.fail("Add test for reading known format csv file e.g. excel")

	def testCSVRepresentation(self):
	
		'''Tests the csv representation returned by the object is whats expected'''
	
		string = "Header1, Header2\n1, 2\nmichael, 34.66\n"
		self.assertEqual(string, self.matrix.csvRepresentation())
	
	def testArchiving(self):
		
		'''Tests that archiving a matrix and then reading it back produces the same matrix'''
		
		filename = 'temp'
		self.matrix.writeToFile('temp')
		newMatrix = Matrix.matrixFromFile(filename)
		os.remove(filename)
		
		self.assertEqual(newMatrix, self.matrix)
		
	def testDescription(self):
	
		'''Tests the matrix description string is whats expected'''
	
		string = "Matrix with 2 rows and 2 columns - 4 elements in total"
		self.assertEqual(string, self.matrix.description())
		
	def testRowMutability(self):
			
		'''Checks that matrix copies row additions so they can't be changed externally'''
		
		myRow = [1,12342]
		self.matrix.addRow(myRow)
		myRow.pop()
		self.assertEqual(self.matrix.row(2), tuple([1,12342]))
		
	def testInitialisationMutability(self):
			
		'''Checks that matrix copies the data passed on initialisation so it can't be changed externally'''
		
		#Remove the last row from the list
		lists.pop()
		
		#Check the last row still exists in the matrix and check it has correct data
		self.assertEqual(self.matrix.row(1), tuple(["michael", 34.66]))
	
		#Add it again for other tests
		#NOTE: If the above test fails then so will others since the readdition wont occur.
		lists.append(["michael", 34.66])
		
       #def testReadAdMatrix(self):
       #
       #	'''Tests initialisation of a Matrix instance from an archived AdMatrix instance'''
       #	
       #	matrix = Matrix.matrixFromArchivedAdunMatrix('testAdMatrix.ad')
       #	self.assertEqual(matrix, self.matrix)
       #	
       #def testWriteAdMatrix(self):
       #
       #	'''Test can the matrix write itself as an AdMatrix'''
       #	
       #	self.fail("Test not implemented")
		
	def testColumnAttributeAccess(self):
	
		'''Test we can access a column as an attribute using its header'''
		
		column = self.matrix.Header1
		self.assertEqual(self.matrix.columnWithHeader('Header1'), column)
	
if __name__ == "__main__":
	unittest.main()
	
		


