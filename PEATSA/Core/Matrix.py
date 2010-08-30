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

'''Module containing a Matrix class used by the various Program classes to store data

Also contains some functions which act like class initialisers.
These allow Matrix instances to be created from csv and archived AdMatrix instances for example'''

import pickle, csv, StringIO, operator, Utilities, Exceptions

charged = ['ARG', 'LYS', 'HIS', 'ASP', 'GLU']
polar =  ['SER', 'THR', 'ASN', 'GLN']
nonpolar = ['ALA', 'LEU', 'ILE', 'VAL', 'MET', 'TYR', 'TRP', 'PHE']
special = ['PRO', 'CYS', 'GLY']

groupToResidue = {'charged':charged, 'polar':polar, 'nonpolar':nonpolar, 'special':special}

def residueIsInGroup(code, group):

	'''Returns YES if the residue defined by code is in the given group.
	
	group can be - charged, polar, nonpolar, special.'''
	
	retval = True
	members = groupToResidue[group]
	if Utilities.IsExtendedResidueFormat(code):
		chain, number, residueName = Utilities.ParseResidueCode(code)
		try: 
			members.index(residueName)
		except ValueError:
			retval = False
	else:
		raise ValueError, 'Residue codes passed to this method must be extended - (passed %s)' % code
		
	return retval

def isPolarResidue(code):

	return residueIsInGroup(code, 'polar')
	
def isChargedResidue(code):

	return residueIsInGroup(code, 'charged')
		
def isSpecialResidue(code):

	return residueIsInGroup(code, 'special')

def isNonPolarResidue(code):

	return residueIsInGroup(code, 'nonpolar')

def residueToGroupDict():

	newDict = {}
	for group in groupToResidue.keys():
		residues = groupToResidue[group]
		for residue in residues:
			newDict[residue] = group

	return newDict

def matrixFromFile(filename):

	'''Creates a matrix instance from one stored in a file.
	
	Parameters
		filename: The name of a file with a archived Matrix instance in it
		
	Return
		Returns a Matrix instance
		
	Exceptions:
		Raises an IOError if no file called filename exists
		Raises a TypeError if filename does not contain an archived matrix'''
		
	file = open(filename)					
	unpickler = pickle.Unpickler(file)
	matrix = unpickler.load()
	file.close()
	
	if not isinstance(matrix, Matrix):
		raise TypeError, "Unarchived object is not an instance of the Matrix class"
	
	return matrix

def matrixFromArchivedAdunMatrix(filename):

	'''Creates a matrix instance from one stored in file in Adun format.
	
	Parameters
		filename: The name of a file with a archived AdDataMatrix instance in it
		
	Return
		Returns a Matrix instance
		
	Exceptions:
		Raises an IOError if no file called filename exists
		Raises a TypeError if filename does not contain an archived AdDataMatrix matrix'''
		
	#Check we can run the conversion script
	
	return None	
	
def matrixFromCSVRepresentation(string, readHeaders=True, dialect=None):

	'''Creates a matrix instance from a csv representation of a matrix (a string).
	
	Note: Any lines beginning with a hash (#) are ignored.
	
	Parameters
		string: The csv reprensentation
		readHeaders: True if the first (non-comment) line of the string should be interpreted
			as giving the column headers. Default is True.
		dialect: An instance of the csv.Dialect class containing information
		on the dialect of the csv representation. Defaults to None.
		In this case the function tries to deduce the dialect itself.
		
	Return
		Returns a Matrix instance
		If the first column header is Mutations returns a PEATSAData matrix
		
	Exceptions:
		Raises a TypeError if a matrix could not be read from the file'''

	#Find the first line that isn't a # into a string.
	#Then use StringIO to treat this string like a file
	#for use with the csv reader object.
	file = StringIO.StringIO(string)
	
	#Seek to the first line that does not begin with a hash
	position = 0
	flag = 0
	while flag == 0:
		line = file.readline()
		if line[0] != '#':
			flag = 1
		else:
			position = file.tell()
	
	file.seek(position)
	fileContent = file.read()
	file.close()
	
	#Create the file-like string object
	csvFile = StringIO.StringIO(fileContent)
	
	#Check the file dialect before processing
	if dialect == None:
		line = csvFile.readline()
		sniffer = csv.Sniffer()
		dialect = sniffer.sniff(line)
		csvFile.seek(0)	
	
	#The reader reads the csv file and converts each
	#line of the file into a list containing the csv elements.
	reader = csv.reader(csvFile, dialect=dialect)
	
	#Read all the rows
	rows = [row for row in reader]
	
	#Get the column headers if specified
	if readHeaders:
		headers = rows.pop(0)
	else:
		headers = None
		
	csvFile.close()	

	#Convert any floats or ints stored as strings
	convertedRows = []
	for row in rows:
		newRow = []
		for element in row:
			try: 
				#Check if its an int
				#Do this first since all ints can be converted to floats
				element = int(element)
				newRow.append(element)
				continue
			except ValueError:
				pass
				
			try:
				#Check if its a float
				element = float(element)
				newRow.append(element)
				continue
			except ValueError:
				#Its not an int or float
				newRow.append(element)
		
		convertedRows.append(newRow)	
	
	if headers is not None:
		headers = [header.strip() for header in headers]
	
	#If the matrix contains this first columnsreturn a PEATSA matrix
	if headers is not None and 'Mutations' in headers:
		matrix = PEATSAMatrix(convertedRows, headers)
	else:	
		matrix = Matrix(convertedRows, headers)	
					
	return matrix
	
def matrixFromCSVFile(filename, readHeaders=True, dialect=None):

	'''Creates a matrix instance from a csv file.
	
	Note: Any lines beginning with a hash  (#) are ignored
	
	Parameters
		filename: The name of a file with a matrix in csv format
		readHeaders: True if the first (non-comment) line of the csv file should be interpreted
			as giving the column headers. Default is True.
		dialect: An instance of the csv.Dialect class containing information
		on the dialect of the csv file. Defaults to None.
		In this case the function tries to deduce the dialect itself.
		
	Return
		Returns a Matrix instance
		
	Exceptions:
		Raises an IOError if no file called filename exists
		Raises a TypeError if a matrix could not be read from the file'''
	
	#Open the file and read the contents
	#Then pass it to matrixFromCSVRepresentation
	fileObject = open(filename, 'r')
	string = fileObject.read()
	fileObject.close()
	
	return matrixFromCSVRepresentation(string, readHeaders, dialect)
	
class Matrix:

	'''Basic matrix class
	
	Can store heterogenous data and contains functionality for handling named columns.
	Also supports comparison and iteration and can be converted to other formats'''
	
	isInitialised = False
	
	def initialiseClass(thisClass):
	
		'''Class method that sets up various things.
		
		Basically required to see if we can locate and run the script to convert AdDataMatrices
		to csv format'''
	
		if not Matrix.isInitialised:
			Matrix.isInitialised = True
	
	initialiseClass = classmethod(initialiseClass)

	def __init__(self, rows=[], headers=None, name='No Name'):
	
		'''Initialise a new matrix
		
		Parameters -
			rows: A list of lists with the contents of the matrix.
				All lists must be the same length otherwise a TypeError exception is raised
			headers: A list containing the headers for each column
			name: A name for the matrix
			
		Exceptions -	
			The number of headers must be equal to the number of elements in each list in rows.
			If not a TypeError is raised
			A TypeError is also raised is the supplied matrix name is not a string'''		
				
		self.numberRows = len(rows)
		self.setName(name)
		self.columnHash = None
		self.currentKeyColumn = None
		
		if len(rows) != 0:
			self.numberColumns = len(rows[0])
		else:
			self.numberColumns = 0
		
		#Check length of each row is the same
		for row in rows:
			if len(row) != self.numberColumns:
				raise IndexError, 'All supplied rows not of same length (%d, %d, %s)' % (self.numberColumns, len(row), str(row))
		
		#Check header length
		if headers != None:
			if len(headers) != self.numberColumns:
				raise IndexError, 'Incorrect number of column headers supplied (%d, %d)' % (self.numberColumns, len(headers))
			else:
				self.headers = list(headers)
		else:
			self.headers = []
			for i in range(self.numberOfColumns()):
				self.headers.append('None')
		
		#Copy the rows
		self.matrix = []
		for row in rows:
			newRow = []
			for element in row:
				newRow.append(element)
			self.matrix.append(newRow)	
		
		#Call the class initialisation method	
		self.__class__.initialiseClass()
	
	def __cmp__(self, object):
	
		'''Compares two matrix instances.
		
		Return
			Returns 0 if the data in the matrices ns identical and they both have identical headers.
			Otherwise returns -1'''
			
		if isinstance(object, Matrix):
			#Cmp returns 0 if two objects are equal
			value = cmp(self.matrix, object.matrix) + cmp(self.columnHeaders(), object.columnHeaders())
			value = value + cmp(self.name(), object.name())
			return value
		else:
			return -1
		
	def __getRowSlice(self, obj):
	
		'''Processes the row slice'''
		
		start, end, step = obj.indices(self.numberOfRows())
		rows = []
		for i in range(start, end):
			rows.append(self.row(i))
			i = i + step
		return self.__class__(rows=rows, headers=self.columnHeaders())
			
	def __getColumnSlice(self, obj):
	
		'''Processes the slice'''
			
		start, end, step = obj.indices(self.numberOfColumns())
		matrix = self.__class__(rows=[])
		for i in range(start, end):
			matrix.addColumn(self.column(i))
			i = i + step
			
		matrix.setColumnHeaders(self.columnHeaders()[start:end:step])

		return matrix				
		
	def __getitem__(self, obj):
	
		'''For matrix[i] returns row i
		For matrix[i:j] return another Core.Matrix instance with the specified rows'''

		if type(obj) is slice:
			return self.__getRowSlice(obj)
		elif type(obj) is tuple and len(obj) == 2:
			#Assume first slice is row slice and second is column slice
			#First get the row slice
			matrix = self.__getRowSlice(obj[0])
			#Now get a column slice of the row slice
			return matrix.__getColumnSlice(obj[1])	
		else:	
			return self.row(obj)
		
	def __str__(self):
	
		'''Returns an information string'''
		
		string = ("Matrix with %d rows and %d columns - %d elements in total"  %
				(self.numberOfRows(), self.numberOfColumns(),
				self.numberOfRows()*self.numberOfColumns()))
		return string
		
	def __getattr__(self, name):
	
		'''Allows column headers to be used as attributes
		
		Check if name is a column header - if it is return it.
		Otherwise raise AttributeError.
		
		If the column name has space you can access it by converting it to camelCase.
		In this case the name is assumed to have all capitals'''
				
		#Don't search for headers attribute if it hasn't been defined!
		#Otherwise start an infinite recursion
		if name == 'headers':
			raise AttributeError, name

		if name in self.headers:
			return self.columnWithHeader(name)
		elif name[:2] != "__":
			#try converting the name from camel case
			indexes = [0]
			for i in range(len(name)):
				if name[i].isupper() and i != 0:
					indexes.append(i)
					
			indexes.append(len(name))
			substrings = []
			for i in range(len(indexes) - 1):
				substring = name[indexes[i]:indexes[i+1]]
				substring = substring.capitalize()
				substrings.append(substring)
									
			convertedNameSpaces = " ".join(substrings)	
			convertedNameNoSpaces = "".join(substrings)
			if convertedNameSpaces in self.headers:
				return self.columnWithHeader(convertedNameSpaces)
			elif convertedNameNoSpaces in self.headers:
				return self.columnWithHeader(convertedNameNoSpaces)
			else:					
				raise AttributeError, name
		else:
			raise AttributeError, name	
		
	def name(self):
		
		'''Returns the matrices names'''
		
		return self.matrixName
		
	def setName(self, name):
	
		'''Sets the matrices name
		
		Parameters: 
			name - a string
			
		Exceptions:
			Raises a TypeError if name is not a string'''
			
		if not isinstance(name, str):
			raise TypeError, "Matrix name must be a string object"
			
		self.matrixName = name	
	
	def numberOfRows(self):
	
		'''Returns the number of rows in the matrix'''
		
		return self.numberRows
		
	def numberOfColumns(self):
	
		'''Returns the number of columns in the receiver'''
	
		return self.numberColumns
	
	def element(self, rowIndex, columnIndex):
	
		'''Returns a matrix element
		
		Parameters
			rowIndex: The elements rowIndex
			columnIndex: The elements columnIndex
		
		Exceptions:
			Raises an IndexError if either rowIndex or columnIndex is out of the range of the reciever'''
		
		return self.matrix[rowIndex][columnIndex]		
	
	def row(self, index):
	
		'''Returns a row of the matrix
		
		Parameters 
			index: an int indicating the row to return
		
		Return
			A tuple containing the row elements.
		
		Exceptions:
			Raises an IndexError of index is out of range of the receiver'''
		
		return tuple(self.matrix[index])
		
	def column(self, index):
	
		'''Returns a row of the matrix
		
		Parameters 
			index: an int indicating the row to return
			
		Return
			A tuple containing the column elements
		
		Exceptions:
			Raises an IndexError of index is out of range of the receiver'''
	
		if index >= self.numberOfColumns():
			raise IndexError, 'Column index %d is out of range of the receiver (%d)' % (index, self.numberOfColumns())
			
		column = [row[index] for row in self.matrix]
		
		return tuple(column)
	
	def addRow(self, row, index=-1):
	
		'''Adds a row to the receiver
		
		Parameters:
			aList: A list containing the elements to add
				Must containing the correct number of elements.
			index: The index at which to add the new row.
				-1 indicates it should be placed at the end.
			
		Exceptions
			Raises an IndexError if index is beyond the range of the receiver and is not -1
			Raises an IndexError if aList does not containing the correct number of elements.'''
		
		if len(row) != self.numberOfColumns():
			raise IndexError, 'Incorrect number of elements (%d) in supplied row %d' % (len(row), self.numberOfRows())
		
		if index != -1:
			if index > self.numberOfRows():
				raise IndexError, 'Supplied index %d is beyond the range of the receiver %d' % (index, self.numberOfRows())
			else:
				self.matrix.insert(index, list(row))
		else:
			self.matrix.append(list(row))
		
		self.numberRows = self.numberRows + 1
		
		#Update column hash if present
		if self.keyColumn() is not None:
			self._createColumnHash(self.keyColumn())
		
	def addColumn(self, column, index=-1):
	
		'''Adds a column to the receiver
		
		The new column has header 'None'
		
		Parameters:
			aList: A list containing the elements to add
				Must containing the correct number of elements.
			index: The index at which to add the new column.
				-1 indicates it should be placed at the end.
			
		Exceptions
			Raises an IndexError if index is beyond the range of the receiver.
			Raises an IndexError if aList does not containing the correct number of elements.'''
	
		if len(column) != self.numberOfRows() and self.numberOfRows() != 0:
			raise IndexError, 'Incorrect number of elements (%d) in supplied column (%d)' % (len(column), self.numberOfColumns())
		
		#If index is -1 change to the number of columns in the receiver
		if index == -1:
			index = self.numberOfColumns()
	
		if index > self.numberOfColumns():
			raise IndexError, 'Supplied index %d is beyond the range of the receiver %d' % (index, self.numberOfColumns())
		elif self.numberOfRows() == 0:
			for i in range(len(column)):
				self.matrix.append([])
			self.numberRows = len(column)
		
		for i in range(len(column)):
			self.matrix[i].insert(index, column[i])
	
		#Add a new header 
		self.headers.insert(index, 'None')
		self.numberColumns = self.numberColumns + 1		
	
	def removeRows(self, index):
	
		'''Removes the specified row
		
		Parameters:
			index - The index of the row to remove.
			A negative index means remove from the end'''
		
		if index < 0:
			index = self.numberOfRows() - index
			
		if abs(index) > self.numberOfRows():
			raise IndexError, 'Supplied index %d is beyond the range of the receiver %d' % (index, self.numberOfColumns())
		
		
		self.matrix.pop(index)
		self.numberRows = self.numberRows - 1
		
		#Update column hash if present
		if self.keyColumn() is not None:
			self._createColumnHash(self.keyColumn())
			
	def removeColumn(self, index):
	
		'''Removes the specified column
		
		Parameters:
			index - The index of the column to remove.
			A negative index means remove from the end'''
		
		if index < 0:
			index = self.numberOfColumns() - index
			
		if abs(index) > self.numberOfColumns():
			raise IndexError, 'Supplied index %d is beyond the range of the receiver %d' % (index, self.numberOfColumns())
		
		self.headers.pop(index)
		for i in range(self.numberOfRows()):
			self.matrix[i].pop(index)		
			
		self.numberColumns = self.numberColumns - 1					
									
	def columnWithHeader(self, aString):
		
		'''Returns the column whose header is given by aString.
		
		Parameters
			header: A aString 
		
		Return
			A tuple containing the column elements.
			If more than one column has the same header the first is returned.
			None if no column has aString as a header'''
			
		try:
			index = self.indexOfColumnWithHeader(aString)
			column = self.column(index)
		except ValueError:
			#No header called aString exists
			column = None	
			
		return column
				
	def headerForColumn(self, index):	
	
		'''Returns the header for a specified column.
		
		Parameters
			index: The index of a column.
		
		Return
			The column header - a aString.
			
		Exceptions
			Raises an IndexError if index is out of range of the receiver'''
	
		return self.headers[index]
			
	def columnHeaders(self):
	
		'''Returns the headers of all the columns
		
		Return
			A tuple containing the column headers'''
			
		return tuple(self.headers)			
	
	def indexOfColumnWithHeader(self, aString):
	
		'''Returns the index of the column whose header is equal to aString
		
		Parameters
			aString: A aString correpsonding to a column header
			
		Exceptions: 
			Raises an ValueError if aString is not the header of any column'''
		
		return self.headers.index(aString)
	
	def setColumnHeader(self, index, header):
		
		'''Sets a new header for the specified column
		
		Parameters
			index: The index of the column whose header is to be changed.
			header: A string - The new header
			
		Exceptions	
			Raises an IndexError if index is out of range of the receiver'''
			
		self.headers[index] = header
		
	def setColumnHeaders(self, headers):
	
		'''Sets new headers for all the columns.
		
		Parameters
			headers - A list or tuple containing the new values.
				The first value in this variable if set as the header for the first column and so on.
				
		Exceptions
			Raises an IndexError if headers does not have the correct number of elements.
			It must have the same number of elements as rows'''
	
		if len(headers) != self.numberOfColumns():
			raise IndexError, "Incorrect number of headers supplied %d. Requires %d" % (len(headers), self.numberOfColumns())
			
		self.headers = list(headers)
	
	def description(self):
	
		'''Returns a description of the receiver (a string)'''
		
		return self.__str__()
	
	def writeToFile(self, filename):
	
		'''Writes the matrix to a file.
		
		The matrix is written in a binary representation.
		Use matrixFromFile() to read the matrix.
		
		Parameters
			filename: Name of the file to write to.
		
		Exceptions
			Raises an IOError if the file cannot be created'''
	
		file = open(filename, 'w')
		pickler = pickle.Pickler(file)
		pickler.dump(self)
		file.close()
		
	def writeAsAdMatrixArchive(self, filename):
		
		'''Writes the matrix to a file as an archived AdMatrix
		
		Parameters: 
			filename: Name of the file to write to.
			
		Exceptions:
			Raises an IOError if the file cannot be created
			Rasies a EnvironmentError if the scripts necessary for
			creating the archive cannot be exectuted'''
			
		pass	
			
	def csvRepresentation(self, comments=None, includeHeaders=True):
	
		'''Returns a csv representation of the receiver.
		
		This is a aString with row elements seperated by commas
		and rows separated by newlines
		
		Parameters:
			comments: A list of strings which will be inserted as comments in the returned string.
			The strings should not contain new line characters.
			includeHeaders: If false the headers are not present in the CSV representation.
			Default: True'''
		
		csvRep = ""
		if includeHeaders is True:
			csvRep = reduce(lambda x, y: str(x) +", " + str(y), self.columnHeaders())
			csvRep = csvRep + "\n"
		
		rowStrings = []
		for row in self:
			rowStrings.append(reduce(lambda x,y: str(x) + ", " + str(y), row))

		csvRep = csvRep + '\n'.join(rowStrings) + '\n'	
			
		if comments is not None and len(comments) is not 0:
			commentString = ""
			for comment in comments:
				#Get rid of any newlines as these could affect the formatting
				comment = comment.replace("\n", " ")
				commentString = commentString + "# " + comment + '\n'	
			csvRep = commentString + csvRep	
				
		return csvRep	
		
	def sort(self, columnHeader='Total', descending=True, function=cmp):
	
		'''Sorts the entire matrix in place by the data the column with header'columnHeader'.
		
		Params:
			columnHeader: The name of the column to sort
			
		Returns: Nothing'''
		
		index = self.indexOfColumnWithHeader(columnHeader)
		self.matrix.sort(function, lambda x: x[index], descending)	

	def filter(self, columnHeader='Total', filterFunction=None, sort=False):
	
		'''Returns a matrix filtered using filterFunction on the given column
		
		Parameters:
			columnHeader - One of the matrices in the receiver. Defaults to 'Total'
			filterFunction - The function to filter using. Must return True or False for each value in the choosen column
			
		Returns:
			A Matrix instance
			
		Exceptions:
			Raises a ValueError if there is no column with the given header'''	

		filteredRows = []
		index = self.indexOfColumnWithHeader(columnHeader)
		
		for row in self:
			value = row[index]
			if filterFunction(value) is True:
				filteredRows.append(row)
		if sort is True:
			filteredRows.sort(cmp, lambda x: x[index], True)
		
		return Matrix(rows=filteredRows, headers=self.columnHeaders())

	def _createColumnHash(self, columnHeader):

		'''Creates a general hash from a given matrix column'''

		hashColumnIndex = self.indexOfColumnWithHeader(columnHeader)
		self.columnHash = {}
		index = 0
		for row in self.matrix:
			key = row[hashColumnIndex]
			if not self.columnHash.has_key(key):
				self.columnHash[key] = []
			
			self.columnHash[key].append(index)
			index += 1
	
	def keys(self):
	
		'''Returns the current hash keys'''
		
		return self.columnHash.keys()										
							
	def setKeyColumn(self, columnHeader):
	
		'''Sets the hash column to \e columnHeader
		
		Raises a ValueError if no column has the supplied header'''
		
		self.currentKeyColumn = columnHeader
		self._createColumnHash(columnHeader)
		
	def keyColumn(self):
	
		'''Returns the header of the current key column or None if none has been set'''
		
		return self.currentKeyColumn
		
	def rowIndexesForKey(self, key):
	
		'''Returns the indexes of the rows containing key in the current hash column.
		
		Raises KeyError if no column contains key
		Returns None if no key column has been set'''
		
		if self.columnHash is not None:
			return tuple(self.columnHash[key])
		else:
			return None	
		
	def rowsForKey(self, key):
	
		'''Returns the rows corresponding to the given key
		
		Raises KeyError if no column contains key
		Returns None if no key column has been set'''
		
		if self.columnHash is not None:
			indexes = self.columnHash[key]
			return [self.row(index) for index in indexes]
		else:
			return None	

class PEATSAMatrix(Matrix):

	'''Adds some extra methods handy for working with PEATSA results matrices
	
	PEATSAMatrices have a column called Mutations - which associates the data in each 
	row with a particular mutant'''

	def __init__(self, rows=[], headers=None, name='No Name'):

		Matrix.__init__(self, rows, headers, name)
		self.mutationHash = None

	def addRow(self, row, index=-1):

		Matrix.addRow(self,row, index)
		self._createMutationHash()

	def removeRows(self, index):
		
		Matrix.removeRows(self, index)
		self._createMutationHash()
	
	def duplicateEntries(self):
	
		'''Identifies duplicate rows in the reciever
		
		An entry is considered a duplicate of a previous one if they both refer to the same mutation
		
		Returns:
			A list of [index, mutationCode] pairs.
			index is the row index where the duplicate occurs (0 offset)'''
		
		import Data
		
		mutationSets = [Data.MutationSet(code) for code in self.mutations]
		seen = {} 
		duplicates = [] 
		for i in range(len(mutationSets)): 
			set = mutationSets[i]
			#Check if set is in the seen dictionary				
			if set.filename(reduced=True) in seen:
				duplicates.append([i, "+".join(set.mutationCodes(reduced=True))])
				 
			#First time seeing this set - add it to seen	
			seen[set.filename(reduced=True)] = 1 
			
		return duplicates	
			
	def mutatedResidues(self, pdb=None):
	
		'''Returns the residue codes corresponding to positions where there are mutations
		
		Returns: a list or reduced residue codes'''
		
		import Data
		
		mutationSets = [Data.MutationSet(code) for code in self.mutations]
		residueCodes = []
		codes = []
		
		for mutationSet in mutationSets:
			for code in mutationSet.residueCodes(reduced=True, pdb=pdb):
				codes.append(code)
			residueCodes.extend(codes)
			codes[:] = []

		residueCodes = list(set(residueCodes))

		return residueCodes
		
	def typeSubstitutions(self, pdb=None):
	
		'''Returns a list of the type subsititutions observed in the data
		
		There are 16 possible type subs.
		
		Returns: a list of tuples'''
		
		mutationsColumn = self.indexOfColumnWithHeader('Mutations')			
		groupMap = residueToGroupDict()
		exchangeDict = {}
		
		for row in self:
			codes = row[mutationsColumn].split('+')
			subs = [Utilities.ParseMutationCode(code)[3] for code in codes]
			wt = [Utilities.ParseMutationCode(code)[2] for code in codes]
			
			pairs = zip(wt, subs)
			for pair in pairs:
				exchange = [groupMap[pair[0]], groupMap[pair[1]]]
				name = '-'.join(exchange)
				if not exchangeDict.has_key(name):
					exchangeDict[name] = 1
				
		return exchangeDict.keys()		
		
	def mutationCodes(self, pdb=None):
	
		'''Returns the mutations codes corresponding the mutations present in the data-set
		
		Returns: a list of reduced mutation codes'''
		
		import Data
		
		mutationSets = [Data.MutationSet(code) for code in self.mutations]
		mutationCodes = []
		codes = []
		
		for mutationSet in mutationSets:
			codes = mutationSet.reducedMutationCodes(pdb=pdb)
			mutationCodes.extend(codes)
			codes[:] = []

		mutationCodes = list(set(mutationCodes))

		return mutationCodes	

	def  pairedMutatedResidues(self, pdb=None):
	
		'''Returns pairs of residue code which occur together
		
		Returns a list of 2 element lists residue codes'''
		
		import Data
		
		mutationSets = [Data.MutationSet(code) for code in self.mutations]
		residueCodes = []
		codes = []
		for mutationSet in mutationSets:
			for code in mutationSet.residueCodes(reduced=True, pdb=pdb):
				codes.append(code)
			
			for i in range(len(codes)):
				for j in range(i+1, len(codes)):
					residueCodes.append('%s+%s' % (codes[i], codes[j]))
			
			codes[:] = []

		residueCodes = list(set(residueCodes))

		return [code.split('+') for code in residueCodes]
	
	def substitutions(self):
	
		'''Returns a list of the amino acid substitutions present'''
		
		import Data
		
		mutationSets = [Data.MutationSet(code) for code in self.mutations]
		substitutions = []
		for mutationSet in mutationSets:
			codes = mutationSet.mutationCodes()
			subs = [Utilities.ParseMutationCode(code)[2] for code in codes]
			substitutions.extend(subs)		

		substitutions = list(set(substitutions))

		return substitutions
	
	def entriesWithMutatedResidues(self, positions, matchMode="All"):
	
		'''Returns a filtered matrix based on the specified positions and matchMode
	
		Params: 
			positions - a list of residue codes
			matchMode - Specified how matches are determined
				All - Returns mutants that contain mutations at all the specified positions
				Any - Returns mutants that contain a mutation at any of the specified posisions
				Exclusive - Returns mutants that only contain mutations at the specified positions.
				
		Example:
			Positions [A15,A127]
			Example mutant - A15A+A145A+A127L
			"All" - Mutants that contain the double mutation A15+A127 - Example matched
			"Any" - Mutants that contain a mutation at either A15 or A127 - Example matched
			"Exclusive" - Mutants that only contain mutations at A15 and/or A127.
				Example not matched since it also contains A145'''
		
		import Data	
					
		#Convert positions to be in full standard format
		filteredPositions = []
		for code in positions:
			components = Utilities.ParseResidueCode(code)
			code = Utilities.CreateResidueCode(components[0], components[1])
			filteredPositions.append(code)
			
		rows = []
		mutationsColumn = self.indexOfColumnWithHeader('Mutations')													
		for row in self:
			mutationSet = Data.MutationSet(row[mutationsColumn])
			codes = mutationSet.residueCodes(reduced=False)
			
			if matchMode == "All":
				matched = True
				for position in filteredPositions:
					if position not in codes:
						matched = False
						break
			elif matchMode == "Any":
				matched = False
				for position in filteredPositions:
					if position in codes:
						matched = True
						break
			elif matchMode == "Exclusive":
				matched = True
				for code in codes:
					if code not in filteredPositions:
						matched = False
						break			
			
			if matched:
				rows.append(row)				
			
		if len(rows) == 0:
			return None
		else:	
			return PEATSAMatrix(rows=rows, headers=self.columnHeaders())		
	
	def entriesWithMutations(self, mutations, matchMode="All"):
	
		'''Returns a filtered matrix based on the specified mutations and matchMode
	
		Params: 
			mutations - a list of mutations codes
			matchMode - Specified how matches are determined
				All - Returns mutants that contain all the specified mutations
				Any - Returns mutants that contain any of the specified mutations
				Exclusive - Returns mutants that only contain the specified mutations.
				
		'''
		
		import Data	
					
		#Convert positions to be in full standard format
		filteredMutations = []
		for code in mutations:
			code = Utilities.ConvertMutationCodeFromReducedFormat(code)
			filteredMutations.append(code)
			
		rows = []
		mutationsColumn = self.indexOfColumnWithHeader('Mutations')													
		for row in self:
			mutationSet = Data.MutationSet(row[mutationsColumn])
			codes = mutationSet.mutationCodes(reduced=False)
			
			if matchMode == "All":
				matched = True
				for mutation in filteredMutations:
					if mutation not in codes:
						matched = False
						break
			elif matchMode == "Any":
				matched = False
				for mutation in filteredMutations:
					if mutation in codes:
						matched = True
						break
			elif matchMode == "Exclusive":
				matched = True
				for mutation in codes:
					if mutation not in filteredMutations:
						matched = False
						break			
			
			if matched:
				rows.append(row)				
			
		if len(rows) == 0:
			return None
		else:	
			return PEATSAMatrix(rows=rows, headers=self.columnHeaders())	
	
	
	def entriesWithSubstitutions(self, substitutions, matchMode="All"):
	
		'''Returns the rows corresponding to mutants containing the specified substitutions.
		
		Params: 
			substitutions: A list of single or three-letter amino-acid codes
			matchMode - see entriesWithMutatedResidues for explanation'''
		
		import Data	
		
		#Convert any single letter codes to triple letter codes
		codeDict = Utilities.InvertedCodeDictionary()	
		temp = substitutions
		substitutions = []
		for substitution in temp:
			if len(substitution) == 1:
				substitution = codeDict[substitution]
			substitutions.append(substitution)	
							
		rows = []	
		mutationsColumn = self.indexOfColumnWithHeader('Mutations')													
		for row in self:
			mutationSet = Data.MutationSet(row[mutationsColumn])
			codes = mutationSet.mutationCodes()
			subs = [Utilities.ParseMutationCode(code)[2] for code in codes]
			
			if matchMode == "All":
				matched = True
				for substitution in substitutions:
					if substitution not in subs:
						matched = False
						break
			elif matchMode == "Any":
				matched = False
				for substitution in substitutions:
					if substitution in subs:
						matched = True
						break
			elif matchMode == "Exclusive":
				matched = True
				for sub in subs:
					if sub not in substitutions:
						matched = False
						break 				
										
			if matched:
				rows.append(row)		
					
		if len(rows) == 0:
			return None
		else:	
			return PEATSAMatrix(rows=rows, headers=self.columnHeaders())	
	
	def entriesWithTypeSubstitutions(self, typeSubs, matchMode="All"):
	
		'''Returns the rows containing mutations involving the specified
		amino-acid type substitutions e.g. charged->polar etc
		
		Note the mutation codes in each row must contain the wild-type residue
		
		Params: 
			typeSubs: A list of tuples. Each tuple contains two strings, one describing
			the wild-type residues type and the other the mutant residue type.
			The types are 'charged', 'polar', 'special', 'nonpolar'
			matchMode - see entriesWithMutatedResidues for explanation'''
				
		rows = []	
		mutationsColumn = self.indexOfColumnWithHeader('Mutations')													
		groupMap = residueToGroupDict()
		for row in self:
			codes = row[mutationsColumn].split('+')
			subs = [Utilities.ParseMutationCode(code)[3] for code in codes]
			wt = [Utilities.ParseMutationCode(code)[2] for code in codes]
			
			pairs = zip(wt, subs)
			exchanges = []
			for pair in pairs:
				exchange = [groupMap[pair[0]], groupMap[pair[1]]]
				exchanges.append(exchange)
				
			if matchMode == "All":
				matched = True
				for substitution in typeSubs:
					if substitution not in exchanges:
						matched = False
						break
						
			elif matchMode == "Any":
				matched = False
				for substitution in typeSubs:
					if substitution in exchanges:
						matched = True
						break
						
			elif matchMode == "Exclusive":
				matched = True
				for sub in typeSubs:
					if sub not in exchanges:
						matched = False
						break 				
										
			if matched:
				rows.append(row)		
					
		if len(rows) == 0:
			return None
		else:	
			return PEATSAMatrix(rows=rows, headers=self.columnHeaders())

	def _createMutationHash(self):

		import Data

		mutationColumn = self.indexOfColumnWithHeader('Mutations')
		self.mutationHash = {}
		index = 0
		for row in self.matrix:
			try:
				code = row[mutationColumn]
				set = Data.MutationSet(code=code)
				self.mutationHash['+'.join(set.mutationCodes())] = index
			
			except Exceptions.MutationCodeFormatError, data:
				print 'Unable to parse mutation code', code
				print 'Reason -', data

			index += 1
	
	def dataForMutationSet(self, mutationSet, skipUnknownCodes=False):
	
		'''Returns the row corresponding to mutant defined by mutationSet

		If skipUnknownCodes is True, rows whose mutation codes cannot be
		interpreted are skipped.
		
		Raises an IndexError if now such row exists.'''

		if self.mutationHash is None:
			self._createMutationHash()

		try:
			index = self.mutationHash["+".join(mutationSet.mutationCodes())]
		except KeyError:
			raise IndexError, 'No data for mutation %s' % mutationSet

		return self.row(index)
				
	def dataForMutation(self, mutation):
		
		'''Returns the row corresponding to mutation
		
		Raises an IndexError is mutation doesn't exist.
		Raises an Exceptions.MutationCodeFormatError if mutation can't be parsed'''
		
		import Data
		
		if self.mutationHash is None:
			self._createMutationHash()

		#Have to ensure the key used to get index is in the correct format
		try:
			set = Data.MutationSet(code=mutation)
			index = self.mutationHash["+".join(set.mutationCodes())]
		except KeyError:
			raise IndexError, 'No data for mutation %s' % mutation

		return self.row(index)
		
	def filterColumnByValue(self, columnHeader='Total', value='5', comparisonFunction=None, isAbsolute=False, sort=False):
	
		'''Filters a column of a the receiver by comparing the values in it to a set value. 
		
		Parameters:
			columnHeader - One of the matrices in the receiver. Defaults to 'Total'
			value - The value to filter using
			comparisonFunction - The function to filter using. Each value in the list will be compared to value
			using this function - defaults to 'greater than or equal to'.
			The function should take two arguments - the first value is compared to the second.
			The function must result True or False depending on the comparison.
			isAbsolute - If true the absolute value of each element of data is used
			
		Returns:
			A list of  (mutationCode, value) pairs containing only the value that pass the filter'''	

		value = float(value)
		data = self.columnWithHeader(columnHeader)
		data = zip(self.mutations, data)
		filterData = []
		
		if comparisonFunction == None:
			comparisonFunction = operator.ge
		
		for pair in data:
			result = pair[1]
			if isAbsolute is True:
				result = abs(result)
			if comparisonFunction(result, value) is True:
				filterData.append(pair)
		
		if sort is True:
			filterData.sort(cmp, lambda x: x[1], True)
		
		return filterData
		
	def filterColumn(self, columnHeader='Total', filterFunction=None, sort=False):
	
		'''Filters a column of the receiver using filterFunction
		
		Parameters:
			columnHeader - One of the matrices in the receiver. Defaults to 'Total'
			filterFunction - The function to filter using. Must return True or False for each value in the choosen column
			
		Returns:
			A list of  (mutationCode, value) pairs containing only the value that pass the filter'''	

		data = self.columnWithHeader(columnHeader)
		data = zip(self.mutations, data)
		filterData = []
		
		for pair in data:
			value = pair[1]
			if filterFunction(value) is True:
				filterData.append(pair)

		if sort is True:
			filterData.sort(cmp, lambda x: x[1], True)
		
		return filterData
		
	def mapColumn(self, columnHeader='Total', mapFunction=None, sort=False):	
	
		'''Transforms a column of the receiver using mapFunction
		
		Parameters:
			columnHeader - One of the matrices in the receiver. Defaults to 'Total'
			mapFunction - The function to map using. Must return a new value for each value in the choosen column
			
		Returns:
			A list of  (mutationCode, value) pairs containing only the mapped values for the column'''
	
		data = self.columnWithHeader(columnHeader)
		mappedData = map(mapFunction, data)
		data = zip(self.mutations, mappedData)

		if sort is True:
			data.sort(cmp, lambda x: x[1], True)
		
		return data

	def sortColumn(self, columnHeader='Total', idColumnHeader='mutations', descending=True):
	
		'''Sorts the column given by 'columnHeader'.
		
		 If idColumnHeader is supplied this method returns a list of (idValue, value) pairs
		 Otherwise it returns a list of values.
		 The use of idValue is to easily identify which records the sorted values correspond to
		
		Params:
			columnHeader: The name of the column to sort
			idColumnHeader: The name of a column containing a string identifiying each row
				This defaults to 'mutations'. In this case the sorted list will contain values like
				[(A10V, -0.1), (A128L, -0.2) ...]'''
		
		data = self.columnWithHeader(columnHeader)
		if idColumnHeader is not None:
			data = zip(self.__getattr__(idColumnHeader), data)	
			data.sort(cmp, lambda x: x[1], descending)
		else:
			data = list(data)	
			data.sort(reverse=descending)	
		
		return data
			
			
	def createHistogram(self, column, binSize, minLimit=None, maxLimit=None, verbose=False):
	
		'''Creates a histogram from column
		
		Returns:
			A Matrix.PEATSAMatrix instance containing the histogram'''
	
		import math
		
		data = self.sortColumn(column, idColumnHeader=None, descending=False)

		min = float(data[0])
		max = float(data[-1])

		if minLimit is not None and min < minLimit:
			print "Over-riding %lf with %lf" % (min, minLimit)
			min = minLimit
		
		if maxLimit is not None and max > maxLimit:
			print "Over-riding %lf with %lf" % (max, maxLimit)
			max = maxLimit

		if verbose:
			print "Minimum %lf. Maximum %lf" % (min, max)

		histRange = max - min
		bins = int(math.ceil(histRange/binSize))

		if verbose:
			print "Range %lf. Bins %lf" % (histRange, bins)

		matrix = []
		minOutliers = 0
		maxOutliers = 0
		
		for i in range(bins):
			row = [min + (i + 0.5)*binSize, 0, 0]
			matrix.append(row)

		for element in data:
			index = int(math.floor((element - min)/binSize))	
			if index < 0:
				print "Min outlier ", element
				minOutliers = minOutliers + 1
				continue
			elif index > bins:
				print "Max outlier ", element
				maxOutliers = maxOutliers + 1
				continue
				
			matrix[index][1] = matrix[index][1] + 1

		for index in range(len(matrix)):
			if index > 0:
				matrix[index][2] = matrix[index-1][2] + matrix[index][1]
			else:
				matrix[index][2] = matrix[index][1]
				
		print "%d points. %d min outliers %d max outliers" % (len(data), minOutliers, maxOutliers)

		return PEATSAMatrix(rows=matrix, headers=[column, 'count', 'cumulative'], name='Histogram')	

