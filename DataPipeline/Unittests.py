#!/usr/bin/env python
#
# DataPipeline - A data import and fitting tool
# Copyright (C) 2011 Damien Farrell
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
# Email: damien.farrell_at_ucd.ie
# Normal mail:
# Damien Farrell
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
#

"""Unit Tests for data pipeline. Uses the code in Testing.py to
   generate most of the tests."""

import unittest
import os, random
import DataPipeline
from Base import Pipeline
import Testing

basictests = Testing.basictests

class BaseTestCase(unittest.TestCase):
    """Basic importer testcase"""
    def setUp(self):
        self.p = Pipeline()
        modulepath = os.path.dirname(DataPipeline.__file__)
        self.confpath = os.path.join(self.p.defaultpath,'temp.conf')
        self.filepath = os.path.join(modulepath, 'testfiles')

class ImporterTestCase(BaseTestCase):
    def dotest(self, filename, conf):
        p = self.p
        p.createConfig(self.confpath,**conf)
        lines = p.openRaw(os.path.join(self.filepath,filename))
        data = p.doImport(lines)
        self.assertTrue(len(data)>0)
        if p.model1 != '':
            p.run()

'''class MultiFileTestCase(BaseTestCase):
    def runTest(self):
        Testing.multiFileTest()'''

class MultiFolderTestCase(BaseTestCase):
    def runTest(self):
        Testing.multiFolderTest()

class GroupedFilesTestCase(BaseTestCase):
    def runTest(self):
        Testing.groupedFilesTest()

class ReplicatesTestCase(BaseTestCase):
    def runTest(self):
        Testing.replicatesTest()

class FitPropagationTestCase(BaseTestCase):
    def runTest(self):
        Testing.fitPropagationTest()

class GroupbyFieldsTestCase(BaseTestCase):
    def runTest(self):
        Testing.groupbyFieldsTest()

class KineticsDataTestCase(BaseTestCase):
    def runTest(self):
        Testing.kineticsTest()

class PreProcessingTestCase(BaseTestCase):
    def runTest(self):
        Testing.preProcessingTest()

class PeakDetectionTestCase(BaseTestCase):
    def runTest(self):
        Testing.peakDetectionTest()

#this adds the basic tests dynamically
def _add_test(name, filename, conf):
    def testmethod(self):
        self.dotest(filename, conf)
    setattr(ImporterTestCase, 'test_'+name, testmethod)
    testmethod.__name__ = 'test_'+name

#create format tests from dictionary
for t in sorted(basictests.keys()):
    info = basictests[t]
    conf = info[0]
    filename = info[1]
    _add_test(str(t), filename, conf)

if __name__ == '__main__':
    unittest.main()

