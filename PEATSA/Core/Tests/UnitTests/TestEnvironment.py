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

'''Contains unittests for the ProteinDesignTool.Environment class'''

import sys, os
sys.path.append("../../")

import Environment, Exceptions
import unittest

#Locations of directories used in the tests
testDirectory = os.path.join(Environment.moduleLocation, 'Tests/UnitTests') 
validDirectory = os.path.join(testDirectory, 'TestDirs/Valid')
invalidDirectory = os.path.join(testDirectory, 'TestDirs/Invalid')
uffbapsRunDirectory = os.path.join(testDirectory, 'TestDirs/Valid/UFFBAPSRun')
testPDB = '1crn.pdb'

class TestConfguration(unittest.TestCase):	

	'''Tests the configuration returned on default initialisation'''

	def setUp(self):
	
		'''Create a configuration that by default reads the test configuration file'''
		
		#This will read the configuration in the Test directory
		self.configuration = Environment.Configuration()
	
	def tearDown(self):
		
		'''Removes the path added by the configuration object instantiation'''
	
		sys.path.remove(self.configuration.pkaModuleLocation())
	
	def testPKAModulePath(self):
	
		'''Checks the path given for PKA Tools is whats expected
		
		This is set to /usr/lib for test purposes'''
		
		path = '/usr/lib'
		
		#If it it should be in sys.path
		self.assertNotEqual(sys.path.count(path), 0)
		#Also test return value from the class itself
		self.assertEqual(self.configuration.pkaModuleLocation(), path)
		
	def testInvalidPKAModulePath(self):
	
		'''Checks that the object detects invalid paths
		
		Here we read an invalid configuration file'''
		
		self.assertRaises(Exceptions.ConfigurationError,
			Environment.Configuration, filename='invalid.conf')
		
	def testCompare(self):
	
		'''Checks comparing configurations'''
		
		newConfiguration = Environment.Configuration()
		self.assertEqual(newConfiguration, self.configuration)
		
		#Compare to some other object
		
		self.assertNotEqual(None, self.configuration)
	
	def testWriteToFile(self):
	
		'''Tests writing and reading'''
		
		self.configuration.writeToFile('temp')
		newConfiguration = Environment.Configuration(filename='temp')
		self.assertEqual(newConfiguration, self.configuration)
		
	def testDesignPKAOptions(self):
	
		'''Test the design pka options
		
		Requires that the design pka tools are installed and in the path'''
		
		import pKa
		
		returnedOptions = self.configuration.designPKAOptions()
		realOptions = pKa.Design_pKa.get_defaults()
		
		#Check they have the same keys
		self.assertEqual(realOptions.keys(), returnedOptions.keys())
		
		#The value for each pka design default key is an list with the form [optionValue, type] 
		#Check the types are all corrected
		
		returnedTypes = [returnedOptions[key][1] for key in returnedOptions.keys()]
		realTypes = [realOptions[key][1] for key in realOptions.keys()]

		self.assertEqual(returnedTypes, realTypes)
		
		
	def testOverides(self):
	
		'''Explicit test that information in a configuration file overides default ones'''
		
		#Create a default instance that doesn't read the configuration file
		newConfiguration = Environment.Configuration(searchDefaultLocations=False, writeFile=False)
		self.assertNotEqual(newConfiguration, self.configuration)
		
	def testFowarding(self):
	
		'''Test that any methods the Configuration doesn't respond to are fowarded to the underlying configuration object'''
		
		sections = ['PKA SCAN METHOD', 'PKA SCAN PARAMETERS', 'PKA SCAN OUTPUT', 'PATHS']
		self.assertEqual(sections, self.configuration.sections())
		
	def testUFFBAPSLocation(self):
	
		'''Test the path to UFFBAPS is correct'''
		
		#This is the path set in the configuration file
		path =  '/usr/bin'
		self.assertEqual(path, self.configuration.uffbapsLocation())
	
	
class TestDefaultConfguration(unittest.TestCase):	

	'''Tests the Configuration returned when no searching is done
	
	Basically the exact same tests but with no file reading'''

	def setUp(self):
	
		'''Create a configuration that by default reads the test configuration file'''
		
		#This will create a configuratio WITHOUT reading the default file
		self.configuration = Environment.Configuration(searchDefaultLocations=False, writeFile=False)
	
	def tearDown(self):
	
		'''Removes the path added by the configuration object instantiation'''
		
		sys.path.remove(self.configuration.pkaModuleLocation())
	
	def testPKAModulePath(self):
	
		'''Checks the path given for PKA Tools is whats expected
		
		This is set to $HOME/PEAT for test purposes'''
		
		path = os.path.join(os.path.expanduser('~'), 'PEAT')

		#If it worked the last entry in sys.path should be this path
		self.assertNotEqual(sys.path.count(path), 0)
		#Also test return value from the class itself
		self.assertEqual(self.configuration.pkaModuleLocation(), path)
		
	def testCompare(self):
	
		'''Checks comparing configurations'''
		
		newConfiguration = Environment.Configuration(searchDefaultLocations=False, writeFile=False)
		self.assertEqual(newConfiguration, self.configuration)
	
	def testWriteToFile(self):
	
		'''Tests writing and reading the configuration file'''
		
		self.configuration.writeToFile('temp')
		newConfiguration = Environment.Configuration(filename='temp')
		self.assertEqual(newConfiguration, self.configuration)
	
	def testUFFBAPSLocation(self):
	
		'''Test the default location'''

		#This is the default path
		path = os.path.join(os.path.expanduser('~'), 'PEAT/UFFBAPS')
		self.assertEqual(path, self.configuration.uffbapsLocation())

class TestWorkingDirectory(unittest.TestCase):
		
	def tearDown(self):
	
		'''Resets the current working directory to the initial directory'''
		
		os.chdir(self.initialWorkingDirectory)
		#Clean any created working directory instance
		try:
			self.workingDirectory.clean()
		except AttributeError:
			pass
		
	def setUp(self):
	
		self.initialWorkingDirectory = os.getcwd()	
	
	def testDefaultInitialisation(self):
	
		'''Tests initialisation on valid working directory
		
		See if no errors are raised'''
		
		try:
			self.workingDirectory = Environment.WorkingDirectory(location=validDirectory)
		except Exceptions.WorkingDirectoryError:
			self.fail()
		
		#If we are here we have passed the test
		self.assertEqual(1,1)
				
	def testMissingDesignFiles(self):
	
		'''Tests initialisation on a working directory not containing the design files
		
		If these files are missing the object should copy them there'''
		
		#Remove the design files from the test dir
		for file in Environment.requiredPKADesignFiles:
			path = os.path.join(validDirectory, file)
			os.remove(path)
		
		#The object should identify the files are missing and copy them back
		#Put stability=false to prevent copying of stability files
		try:
			object = Environment.WorkingDirectory(location=validDirectory)
			object.setup(pdbFile=testPDB, scan=True, stability=False)
		except Exception, data:
			print data
			self.fail()
			
		#Check it copied the files	
		files = os.listdir(validDirectory)
		error = ''
		for file in Environment.requiredPKADesignFiles:	
			try:
				files.index(file)
			except ValueError:
				error = error + 'File %s not copied\n' % file
				
		if len(error) != 0:
			self.fail(error)
		else:
			self.assertEqual(1,1)
	
	def testMissingProteinFiles(self):
	
		'''Tests initialisation on a working directory not containing the required protein files
		
		Should raise a WorkingDirectoryError.
		Also checks that it cleans up after itself before raising the error'''
		
		preContents = os.listdir(invalidDirectory)
		self.workingDirectory = Environment.WorkingDirectory(location=invalidDirectory) 
		self.assertRaises(Exceptions.WorkingDirectoryError, 
			self.workingDirectory.setup,
			pdbFile='Unknown.pdb',
			scan = True,
			stability = False)
		
		#Check it cleaned itself up properly - in this case it should remove anything it put there
		postContents = os.listdir(invalidDirectory)
		self.assertEqual(preContents, postContents)
	
	def testCurrentWorkingDir(self):
	
		'''Tests initialisation changes the current working directory'''
		
		self.workingDirectory = Environment.WorkingDirectory(location=validDirectory)
		self.assertEqual(self.workingDirectory.path(), os.getcwd())
	
	def testSkipScan(self):
	
		'''Tests that not specifying scan skips checking for scan files'''
		
		#Remove the design files from the test dir
		#They could be there if another test failed
		for file in Environment.requiredPKADesignFiles:
			path = os.path.join(invalidDirectory, file)
			try:
				os.remove(path)	
			except OSError:
				#The file isn't there
				pass
								
		try:
			self.workingDirectory = Environment.WorkingDirectory(location=invalidDirectory)
			self.workingDirectory.setup(pdbFile=testPDB, scan=False, stability=True)
		except BaseException, data: 
			print data
			self.fail("Raised exception when it should not have")
		
		#Check that none of the required scan files were copied
		contents = os.listdir(invalidDirectory)
		for file in Environment.requiredPKADesignFiles:
			try: 
				contents.index(file)
				self.fail("Copied scan files (%s) when asked not to" % file)
			except ValueError:
				pass
				
		self.assertEqual(1,1)		
			
	def testUFFBAPSPath(self):
	
		'''Tests that a valid path to the UFFBAPS run dir is returned and that everything is present'''
	
		self.workingDirectory = Environment.WorkingDirectory(location=validDirectory)
		self.assertEqual(self.workingDirectory.uffbapsRunDirectory(), uffbapsRunDirectory)
		self.workingDirectory.setupUFFBAPS()
	
		#Check everthing was copied
		contents = os.listdir(self.workingDirectory.uffbapsRunDirectory())
		expectedContents = os.listdir(Environment.UffbapsResourceDirectory())
		self.assertEqual(contents, expectedContents)
	
	def testClean(self):
	
		'''Tests that clean works
		
		Does this by checking the directory is in the same state afterwards as before
		Note this check only works does not setup a scan since this may copy the design files
		and these will not be removed by clean (as they could have been present already) 
		They are only removed if there was an error during the setup stage.'''
		
		contents = os.listdir(validDirectory)
		
		self.workingDirectory = Environment.WorkingDirectory(location=validDirectory)
		self.workingDirectory.setup(pdbFile=testPDB, 
					scan=False,
					stability=True)
		self.workingDirectory.clean()
		
		#Check the directory contents are the same as when started				
		newContents = os.listdir(validDirectory)
		
		self.assertEqual(contents, newContents)
	

if __name__ == "__main__":
	unittest.main()
	
	
	
	
