#!/usr/bin/env python
#
# Protein Engineering Analysis Tool DataBase (PEATDB)
# Copyright (C) 2010 Damien Farrell & Jens Erik Nielsen
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

"""Tests for data pipeline."""
from Base import Pipeline
import os

basictests = {'test1':{'conf':{'format':'databyrow','delimeter':','},
                       'filename':'databyrow1.txt'},
              'test2':{'conf':{'format':'databycolumn','delimeter':','},
                       'filename':'databycol1.txt'},
              #rows, multiple groups
              'test3':{'conf':{'format':'databyrow','delimeter':',','colrepeat':6},
                     'filename':'databyrow2.txt'},
              #cols, multiple groups
              'test4':{'conf':{'format':'databycolumn','delimeter':',','rowrepeat':6},
                       'filename':'databycol2.txt'},
              #paired x-y data in rows         
              'test5':{'conf':{'format':'paireddatabyrow','delimeter':','},
                       'filename':'databyrow_paired.txt'},
              #paired x-y data in cols    
              'test6':{'conf':{'format':'paireddatabycolumn','delimeter':','},
                       'filename':'databycol_paired.txt'},                     
              #various non-default formatting
              #'test7':{'conf':{'format':'databyrow','delimeter':'tab','decimalsymbol':','},
              #         'filename':'databyrow3.txt'},
            }
            
exceltests = {'test1':{'conf':{'format':'databyrow'},
                       'filename':'databyrow1.xls'}}
    
path = 'testfiles'

def formatTests(testinfo):
    """Test basic standard format handling"""
    
    for t in testinfo:
        p = Pipeline()
        conf=testinfo[t]['conf']
        filename=testinfo[t]['filename']
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

