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
import subprocess, os, shutil

location = '../../'

def TestPEAT_SA(pdbFile, workingDir, outputDir, calculation, ligand=None):

	'''Runs a test for a calculation on a given pdbFile'''

	runString = "%sProteinDesignTool.py -p %s -w %s -o %s %s --mutationList=testList -v" % (location, pdbFile, workingDir, outputDir, calculation) 
	name = os.path.split(outputDir)[1]
	if os.path.exists(outputDir):
		shutil.rmtree(outputDir)

	os.mkdir(outputDir)
	retval = 0

	try:	
		print 'Testing %s' % name
		print 'Executing - %s' % runString
		standardOut = open("%s/%sTest.log" % (name, name), "w")
		standardError = open("%s/%sError.log" % (name,name), "w")
		process = subprocess.Popen(runString, shell=True, stdout=standardOut, stderr=standardError)
		standardOut.close()
		standardError.close()
		process.wait()
		if process.returncode is not 0:
			print '%s test failed - exit code %d' % (name, process.returncode)
			print 'Tail of test log follows\n\n'
			os.system('tail %s/%sTest.log' % (name, name))
			print 'Tail of error log follows\n\n'
			os.system('tail %s/%sError.log' % (name, name))
			print '\nCheck %s/%sErrors.log for details' % (name, name)
			retval = 1
		else:
			print '\nTest suceeded'
	except BaseException, data:
		print '%s test failed before scan process started' % name
		print 'Exception %s' % data
		if hasattr(data, "child_traceback"):
			print "Test traceback %s" % data.child_traceback
		retval = 1

	return retval

