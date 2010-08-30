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
'''Contains unittests for the ProteinDesignTool.Data.MutationList class'''

import sys, os
sys.path.append(".")
sys.path.insert(0, "../../")

import Data, Exceptions
import unittest 

#Simple data for testing.
simpleList = "A:0001, ALA, GLY\nA:0053, HIS, PRO"
standardList = "A:0001:ALA+A:0053:GLY+A:0123:PRO\nA:0004:HIS+B:0012:ALA+A:0032:GLY"
#Has incorrect mutant residue name
invalidSimpleListOne = "A:0001, ALA, GLY\nA:0053, HIST, PRO"
#Has incorrect residue id
invalidSimpleListTwo = "A:0001, ALA, GLY\n0053, HIS, PRO"
#Has incorrect line format (misssing commas)
invalidSimpleListThree = "A:0001, ALA, GLY\nA:0053, HIS PRO"
#Has incorrect mutant residue name
invalidStandardListOne = "A:0001:ALA+A:0053:GLY+A:0123:PRO\nA:0004:HIS+B:0012:ALAS+A:0032:G"
#Has incorrect residue id
invalidStandardListTwo = "A:0001:ALA+A:0053:GLY+A:0123:PRO\nA:0004:HIS+A:0012:ALA+A::GLY"
#Has incorrect line format (incorrect separator - space instead of +)
invalidStandardListThree = "A:0001:ALA+A:0053:GLY+A:0123:PRO\nA:0004:HIS B:0012:ALA+A:0032:GLY"
#Has incorrect mutant residue name


class TestSinglePointMutationFormat(unittest.TestCase):

	'''Test class for Matrix'''

	def setUp(self):

		'''Creates a SPM format file and a Data.MutationList instance from it'''

		stream = open("SimpleMutationList", "w+")
		stream.write(simpleList)
		stream.close()
		self.list = Data.MutationListFile(filename="SimpleMutationList")
		
		#Create a real list of the mutations
		self.mutants = []
		
		set = Data.MutationSet()
		set.addMutation(chain="A", residueIndex=1, mutation='ALA')
		self.mutants.append(set)
		
		set = Data.MutationSet()
		set.addMutation(chain="A", residueIndex=1, mutation='GLY')
		self.mutants.append(set)
		
		set = Data.MutationSet()
		set.addMutation(chain="A", residueIndex=53, mutation='HIS')
		self.mutants.append(set)
		
		set = Data.MutationSet()
		set.addMutation(chain="A", residueIndex=53, mutation='PRO')
		self.mutants.append(set)
	
	def tearDown(self):
	
		os.remove("SimpleMutationList")
		
	def testSinglePointMutationFormatParse(self):
	
		'''Test that the results of parsing an SPM file are correct
		
		Also tests compare'''
		
		self.assertEqual(self.list.mutantList(), self.mutants)
		
	def testFormatDetection(self):
	
		'''Test that a single point mutation format is detected'''		
		
		self.assertEqual(self.list.isSinglePointMutationFormat(), True)	
		
	def testSinglePointMutationFormatWrite(self):
	 
		'''Test that we can write a single point mutation format'''
		
		self.list.writeToFile("testFile")
		newList = Data.MutationListFile("testFile")
		self.assertEqual(newList, self.list)
			
		os.remove("testFile")	
		
	def testSinglePointMutationFormatMutationError(self):
	
		'''Test that the class responds correctly to an incorrect mutation format'''
		
		stream = open("SimpleMutationList", "w+")
		stream.write(invalidSimpleListOne)
		stream.close()
		self.assertRaises(Exceptions.MutationListFileFormatError, Data.MutationListFile, "SimpleMutationList")
	
	def testSinglePointMutationFormatResidueError(self):
	
		'''Test that the class responds correctly to an incorrect residue format'''
		
		stream = open("SimpleMutationList", "w+")
		stream.write(invalidSimpleListTwo)
		stream.close()
		self.assertRaises(Exceptions.MutationListFileFormatError, Data.MutationListFile, "SimpleMutationList")
		
	def testSinglePointMutationLineError(self):
	
		'''Test that the class responds correctly to an incorrectly formatted line'''
		
		stream = open("SimpleMutationList", "w+")
		stream.write(invalidSimpleListThree)
		stream.close()
		self.assertRaises(Exceptions.MutationListFileFormatError, Data.MutationListFile, "SimpleMutationList")		
	
	def testInitialisationWithNonExistantFile(self):
	
		'''Test that the class raises the correct error when initialised with a non-existant file'''
		
		self.assertRaises(Exceptions.MutationListError, Data.MutationListFile, "NonexistantFile")
				
				
class TestStandardFormat(unittest.TestCase):

	def setUp(self):

		'''Creates a mutation list from a standard formatted file'''

		stream = open("StandardMutationList", "w+")
		stream.write(standardList)
		stream.close()
		self.list = Data.MutationListFile(filename="StandardMutationList")

		self.mutants = []
		
		set = Data.MutationSet(code="A:0001:ALA+A:0053:GLY+A:0123:PRO")
		self.mutants.append(set)
	
		set = Data.MutationSet(code="A:0004:HIS+B:0012:ALA+A:0032:GLY")
		self.mutants.append(set)

	def tearDown(self):
	
		os.remove("StandardMutationList")
		
	def testStandardParse(self):
	
		'''Test that the results of parsing a standard file are correct'''
		
		self.assertEqual(self.list.mutantList(), self.mutants)		

	def testFormatDetection(self):
	
		'''Test that standard format is detected'''		
		
		self.assertEqual(self.list.isStandardMutationFormat(), True)	
		
	def testStandardFormatWrite(self):
	 
		'''Test that we can write a standard format file'''
		
		self.list.writeToFile("testFile")
		newList = Data.MutationListFile("testFile")
		self.assertEqual(newList, self.list)
		
	def testStandardFormatMutationError(self):
	
		'''Test that the class responds correctly to an incorrect mutation format'''
		
		stream = open("StandardMutationList", "w+")
		stream.write(invalidStandardListOne)
		stream.close()
		self.assertRaises(Exceptions.MutationListFileFormatError, Data.MutationListFile, "StandardMutationList")
	
	def testStandardFormatResidueError(self):
	
		'''Test that the class responds correctly to an incorrect residue format'''
		
		stream = open("StandardMutationList", "w+")
		stream.write(invalidStandardListTwo)
		stream.close()
		self.assertRaises(Exceptions.MutationListFileFormatError, Data.MutationListFile, "StandardMutationList")
		
	def testStandardFormatLineError(self):
	
		'''Test that the class responds correctly to an incorrectly formatted line'''
		
		stream = open("StandardMutationList", "w+")
		stream.write(invalidStandardListThree)
		stream.close()
		self.assertRaises(Exceptions.MutationListFileFormatError, Data.MutationListFile, "StandardMutationList")		

class TestListModifications(unittest.TestCase):

	'''Tests that lists can be modified and written'''

	def setUp(self):

		'''Creates a simple format file and a Data.MutationList instance from it'''

		stream = open("SimpleMutationList", "w+")
		stream.write(simpleList)
		stream.close()
		self.list = Data.MutationListFile(filename="SimpleMutationList")
		
		#Create a real list of the mutations
		self.mutants = []
		
		set = Data.MutationSet()
		set.addMutation(chain="A", residueIndex=1, mutation='ALA')
		self.mutants.append(set)
		
		set = Data.MutationSet()
		set.addMutation(chain="A", residueIndex=1, mutation='GLY')
		self.mutants.append(set)
		
		set = Data.MutationSet()
		set.addMutation(chain="A", residueIndex=53, mutation='HIS')
		self.mutants.append(set)
		
		set = Data.MutationSet()
		set.addMutation(chain="A", residueIndex=53, mutation='PRO')
		self.mutants.append(set)

	def testAddSinglePointMutant(self):
	
		'''Test that we can add a mutant and write the resulting file.
		
		We add a SPM to an list of SPM's. The next test tests more complex possibilities.'''
				
		set = Data.MutationSet()
		set.addMutation(chain="B", residueIndex=120, mutation='CYS')
		self.list.addMutant(set)
		
		#Check the mutant is added
		self.assertEqual(self.list.mutants.count(set), 1)
		
		#Create a new list from the simple format file
		newList = Data.MutationListFile("SimpleMutationList")
		
		self.assertEqual(newList, self.list)
		
	def testAddComplexMutant(self):	
	
		'''Test can a change of format (Single to Standard) can be detected and dealt with'''
		
		set = Data.MutationSet()
		set.addMutation(chain="B", residueIndex=120, mutation='CYS')
		set.addMutation(chain="A", residueIndex=42, mutation='ALA')
		self.list.addMutant(set)
		
		#Test the format change is detected
		self.assertEqual(self.list.isStandardMutationFormat(), True)
		
		newList = Data.MutationListFile("SimpleMutationList")
		#Test that the file has been updated to the new format
		self.assertEqual(newList.isStandardMutationFormat(), True)
		#Test that the two list are the same
		self.assertEqual(newList, self.list)
		
	def testRemoveMutant(self):
	
		'''Test we can remove an entry from the list'''
		
		mutationSet = self.list.mutantList()[0]
		self.list.removeMutant(mutationSet)
		
		#See if its removed
		self.assertEqual(self.list.mutantList().count(mutationSet), 0)
		
		#Check the file is updated properly
		newList = Data.MutationListFile("SimpleMutationList")
		self.assertEqual(newList, self.list)

	
if __name__ == "__main__":
	unittest.main()
	
		


