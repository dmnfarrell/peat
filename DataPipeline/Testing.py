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
import os, random
import numpy as np
import csv
from PEATDB.Ekin.Base import EkinProject

basictests = {'test1':({'format':'databyrow','model1':'Linear'},'databyrow1.txt'),
              'test2':({'format':'databycolumn'},'databycol1.txt'),
              #rows, multiple groups
              'test3':({'format':'databyrow','colrepeat':6},'databyrow2.txt'),
              #cols, multiple groups
              #'test4':({'format':'databycolumn','rowrepeat':6},'databycol2.txt'),
              #paired x-y data in rows
              #'test5':({'format':'paireddatabyrow'},'databyrow_paired.txt'),
              #paired x-y data in cols
              #'test6':({'format':'paireddatabycolumn'},'databycol_paired.txt'),
              #various non-default formatting
              #'test7':({'format':'databyrow','delimeter':'tab','decimalsymbol':',',
              #          'colrepeat':6}, 'databyrow_errors.txt'),
              #'test8':({'format':'groupeddatabyrow','rowrepeat':4,'rowheader':0,'rowstart':1,
              #           'model1':'Linear'}, 'databyrow_grouped.txt'),
              'test9':({'format':'groupeddatabycolumn','colrepeat':4,'colheader':0,'colstart':1,
                        'model1':'Linear','model2':'sigmoid','xerror':0.2,'yerror':0.3},
                        'databycol_grouped.txt'),
              #'test10':({'format':'databyrow','rowheader':"aaa,bbb,ccc,ddd"},
              #          'databyrow_noheader.txt')
              }

exceltests = {'test1':{'conf':{'format':'databyrow'},
                       'filename':'databyrow1.xls'}}

def doTest(info, name='test', path='testfiles'):
    p = Pipeline()
    conf = info[0]
    filename = info[1]
    p.createConfig('temp.conf',**conf)
    p.openRaw(os.path.join(path,filename))
    #data = p.doImport()
    p.run()
    print '%s completed' %name
    print '-------------------'

def formatTests(testinfo):
    """Test basic standard format handling"""
    for t in sorted(testinfo.keys()):
        doTest(testinfo[t], t)

def multiFileTests():
    """Tests the processing and grouping of multiple files"""
    pth = 'testfiles/group1'
    if not os.path.exists(pth):
        os.mkdir(pth)
    createFakeFiles(pth)
    conf = {'format':'databycolumn', 'model1':'Linear', 'groupbyfile':1}
    p = Pipeline()
    p.createConfig('temp.conf',**conf)
    p.addFolder(pth, ext='txt')
    p.run()
    return

def fitPropagationTest():
    """Tests the propagation of fit data direct from a dict - no importing"""
    import time
    start=time.time()    
    p = Pipeline()
    conf = {'model1':'linear','model2':'linear','model3':'sigmoid',
            'variable1':'a','variable2':'a','variable3':'tm','xerror':0.1,'yerror':0.2}    
    p.createConfig('temp.conf',**conf)   
    data = createNestedData()
    Em = EkinProject()    
    E,fits = p.processFits(data, Em=Em)
    print 'final fits', fits
    Em.saveProject('results.ekinprj')
    print 'completed fit propagation test'
    print 'took %s seconds' %round((time.time()-start),2)
    print '-------------------'
    
    return

def customTests():
    """Tests kinetics data for paper"""

    info = ({'format':'kineticsdata','colrepeat':4,'colheader':0,'colstart':1,
              'model1':'Linear'},'setG_110309_1_pH7,5.txt')
    doTest(info, 'kinetics test', 'novo_setG/rep1')

def createGroupedData1(path='testfiles'):
    """Create sets of fake data to test queuing and file grouping"""

    names = ['aaa','bbb','ccc','ddd']
    for i in np.arange(2,9,1.0):
        fname = os.path.join(path,'ph_'+str(i)+'.txt')
        cw = csv.writer(open(fname,'w'))
        cw.writerow(['temp']+names)
        for x in range(250,360,5):
            vals = [round(i*x/random.normalvariate(10,0.3),2) for j in range(len(names))]
            vals.insert(0,x)
            cw.writerow(vals)
    return
    
def createNestedData():
    """Create fake nested data similar to our kinetics data"""
    data={}
    names = ['aaa','bbb','ccc']
    labels = range(1,8) #e.g. a ph range
    sublabels1 = np.arange(0.2,0.6,0.1)
    x = range(20)    
    for n in names:
        data[n]={}
        for l in labels:
            data[n][l]={}
            for s in sublabels1:
                data[n][l][s]={}
                y=[round(i*s/l+random.normalvariate(1,0.1),3) for i in x]
                data[n][l][s] = (x,y) 
    return data  

#formatTests(basictests)
#customTests()
#multiFileTests()
fitPropagationTest()

