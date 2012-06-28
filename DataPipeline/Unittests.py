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

"""Unit Tests for data pipeline."""

import unittest
import os, random
import DataPipeline
from Base import Pipeline

basictests = {1:({'format':'databyrow','model1':'Linear'},'databyrow1.txt'),
                2:({'format':'databycolumn'},'databycol1.txt'),
                #rows, multiple groups
                3:({'format':'databyrow','colrepeat':6},'databyrow2.txt'),
                #cols, multiple groups
                4:({'format':'databycolumn','rowrepeat':6}, 'databycol2.txt'),
                #paired x-y data in rows
                5:({'format':'paireddatabyrow'},'databyrow_paired.txt'),
                #paired x-y data in cols with colheader offset
                6:({'format':'paireddatabycolumn','colstart':1,'colheaderstart':2},
                         'databycol_paired.txt'),
                #paired x-y data in double rows
                7:({'format':'paireddatabydoublerow','rowstart':0,'rowheader':0},
                        'databyrow_paired_double.txt'),
                #paired x-y data in double cols
                8:({'format':'paireddatabydoublecolumn','rowstart':2,'colheader':0},
                  'databycol_paired_double.txt'),
                #various non-default formatting
                9:({'format':'databyrow','delimeter':'tab','decimalsymbol':',',
                        'colrepeat':6}, 'databyrow_errors.txt'),
                #fitting models included
                10:({'format':'groupeddatabyrow','rowrepeat':4,'rowheader':0,'rowstart':1,
                         'model1':'Linear'}, 'databyrow_grouped.txt'),
                11:({'format':'groupeddatabycolumn','colrepeat':4,'colheader':0,'colstart':1,
                        'model1':'Linear','model2':'sigmoid','variable1':'a','variable2':'tm',
                        'xerror':0.2,'yerror':0.3},
                        'databycol_grouped.txt'),
                12:({'format':'databyrow','rowheader':"aaa,bbb,ccc,ddd"},
                        'databyrow_noheader.txt')
                }

class ImporterTestCase(unittest.TestCase):
    """Basic importer testcase"""
    def setUp(self):
        self.p = Pipeline()
        modulepath = os.path.dirname(DataPipeline.__file__)
        self.confpath = os.path.join(self.p.defaultpath,'temp.conf')
        self.filepath = os.path.join(modulepath, 'testfiles')

    def dotest(self, filename, conf):
        p = self.p
        p.createConfig(self.confpath,**conf)
        lines = p.openRaw(os.path.join(self.filepath,filename))
        data = p.doImport(lines)
        self.assertTrue(len(data)>0)
        if p.model1 != '':
            p.run()

#this adds named tests to the unittests object
def _add_test(name, filename, conf):
    def testmethod(self):
        self.dotest(filename, conf)
    setattr(ImporterTestCase, 'test_'+name, testmethod)
    testmethod.__name__ = 'test_'+name

#dynamically create tests from dictionary
for t in sorted(basictests.keys()):
    info = basictests[t]
    conf = info[0]
    filename = info[1]
    _add_test(str(t), filename, conf)


if __name__ == '__main__':
    unittest.main()

