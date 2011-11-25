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

"""Tests for data pipeline."""
from Base import Pipeline
import os

basictests = {#'test1':({'format':'databyrow'},'databyrow1.txt'),
              #'test2':({'format':'databycolumn'},'databycol1.txt'),
              #rows, multiple groups
              #'test3':({'format':'databyrow','colrepeat':6},'databyrow2.txt'),
              #cols, multiple groups
              #'test4':({'format':'databycolumn','rowrepeat':6},'databycol2.txt'),
              #paired x-y data in rows         
              #'test5':({'format':'paireddatabyrow'},'databyrow_paired.txt'),
              #paired x-y data in cols    
              #'test6':({'format':'paireddatabycolumn'},'databycol_paired.txt'),                     
              #various non-default formatting
              #'test7':({'format':'databyrow','delimeter':'tab','decimalsymbol':','},
              #         'databyrow3.txt'),
              'test8':({'format':'groupeddatabyrow','rowrepeat':4,'rowheader':0,'rowstart':1,
                         'model1':'Linear'},
                         'databyrow_grouped.txt'),
              'test9':({'format':'groupeddatabycolumn','colrepeat':4,'colheader':0,'colstart':1,
                        'model1':'Linear'},
                        'databycol_grouped.txt'),              
              #'test10':({'format':'databyrow','rowheader':"aaa,bbb,ccc,ddd"},
              #          'databyrow_noheader.txt'),                        
              }
            
exceltests = {'test1':{'conf':{'format':'databyrow'},
                       'filename':'databyrow1.xls'}}
    
path = 'testfiles'

def formatTests(testinfo):
    """Test basic standard format handling"""
    
    for t in sorted(testinfo.keys()):
        print t
        p = Pipeline()
        conf=testinfo[t][0]
        filename=testinfo[t][1]
        p.createConfig('temp.conf',**conf)
        p.openRaw(os.path.join(path,filename))
        data = p.doImport()
        if data != None:
            print '%s import ok' %t
            print data.keys()
        print '-------------------'
    
def queueTests():
    return

def fittingTests():
    return
    
formatTests(basictests)

