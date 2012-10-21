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
import os, random, time, string, datetime
import numpy as np
from math import *
import csv
from PEATDB.Ekin.Base import EkinProject
import Utilities

"""This module does any manual pipeline tests, the remainder
    are all in the Unittests module"""

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
              12:({'format':'databyrow','rowheaderlabels':"aaa,bbb,ccc,ddd"},
                        'databyrow_noheader.txt')
              }

exceltests = {1:({'format':'databyrow','colheaderstart':0,'sheet':1},
                  'databyrow1.xls')}

def doTest(info, name='test', path='testfiles'):
    print 'running test %s' %name
    p = Pipeline()
    conf = info[0]
    filename = info[1]
    confpath = os.path.join(p.defaultpath,'temp.conf')
    p.createConfig(confpath,**conf)
    lines = p.openRaw(os.path.join(path,filename))
    data = p.doImport(lines)
    if p.model1 != '':
        p.run()

    return

def formatTests( testinfo, names=None):
    """Test basic standard format handling and fitting"""

    total=0; success=0
    if names == None:
        names = testinfo.keys()
    for t in sorted(names):
        doTest(testinfo[t], t)

def multiFileTest():
    """Test handling of single datasets per file with grouping per filename"""
    path = 'testfiles/singlefiles'
    Utilities.createSingleFileData(path)
    conf = {'format':'databycolumn', 'groupbyname':1,  'parsenamesindex':0,
            'parsemethod':'numeric', 'filenameseparator':'_',
            'model1':'linear','variable1':'a','model2':'sigmoid','variable2':'tm'}
    p = Pipeline()
    p.createConfig('temp.conf',**conf)
    p.addFolder(path)
    p.run()
    return

def groupedFilesTest():
    """Tests the processing and grouping of multiple files with the same
         sets of datasets in all files"""
    path = 'testfiles/grouped'
    Utilities.createGroupedData(path)
    conf = {'format':'databycolumn','groupbyname':1, 'parsenamesindex':0, 'parsemethod':'numeric',
            'model1':'linear','variable1':'a','model2':'sigmoid','variable2':'tm'}
    p = Pipeline()
    p.createConfig('temp.conf',**conf)
    p.addFolder(path)
    p.run()
    return

def multiFolderTest():
    """Handling of multiple folders in a hierarchy with replicates"""
    p = Pipeline()
    conf = {'format':'databycolumn','groupbyname':1,'parsenamesindex':0, 'parsemethod':'numeric',
            'replicates':1, #'saveplots':1,
            'model1':'linear','variable1':'a','model2':'sigmoid','variable2':'tm'}
    p.createConfig('temp.conf',**conf)
    path = 'testfiles/multifolders'
    Utilities.createDirectory(path)
    phs = range(2,10)
    reps = range(1,4)   #replicates
    names = Utilities.createRandomStrings(3,6)
    today = str(datetime.date.today())
    for i in phs:
        #sigmoid dependence of the slopes on 'ph'
        #so we know we are getting the right results
        val = 1/(1+exp((i-4)/1.04))
        folder = os.path.join(path,'ph'+str(i))
        Utilities.createDirectory(folder)
        for r in reps:
            fname = os.path.join(folder,'r'+str(r)+'_'+today+'.txt')
            Utilities.createTempData(fname, names, val)
    p.addFolder(path)
    p.run()
    return

def replicatesTest():
    """Tests handling of replicates"""

    p = Pipeline()
    conf = {'format':'databycolumn','groupbyname':1, 'parsenamesindex':1, 'parsemethod':'numeric',
            'replicates':1,
            'model1':'linear','variable1':'a','model2':'sigmoid','variable2':'tm'}
    p.createConfig('temp.conf',**conf)
    reps = ['rep1','rep2','rep3']
    path = 'testfiles/replicates'
    Utilities.createDirectory(path)
    names = Utilities.createRandomStrings(3,6)
    for r in reps:
        rpath = os.path.join(path, r)
        Utilities.createGroupedData(rpath,names=names)
    p.addFolder(path)
    p.run()
    return

def fitPropagationTest():
    """Tests the propagation of fit data direct from a dict - no importing"""

    start=time.time()
    p = Pipeline()
    conf = {'model1':'linear','model2':'Michaelis-Menten','model3':'sigmoid',
            'variable1':'a','variable2':'Km','variable3':'tm',#'xerror':.1,'yerror':0.05,
            }
    p.createConfig('temp.conf',**conf)
    data = Utilities.createNestedData()
    Em = EkinProject()
    E,fits = p.processFits(data, Em=Em)
    print 'final fits', fits
    fname = os.path.join(p.workingdir,'results')
    Em.saveProject(fname)
    p.saveEkinPlotstoImages(Em, fname)
    print 'completed fit propagation test'
    print 'took %s seconds' %round((time.time()-start),2)
    print '-------------------'
    return

def groupbyFieldsTest():
    """Tests grouping by fields function using NMRdata"""

    conf = {'format':'databycolumn','colheaderlabels':'15N,1H',
            'parsenamesindex':0,'parsemethod':'numeric',
            'delimeter':' ', 'groupbyfields':1,'extension':'.inp'}
    path = 'testfiles/NMRdata'
    p = Pipeline()
    p.createConfig('temp.conf',**conf)
    p.addFolder(path)
    p.run()
    return

def kineticsTest():
    """Tests kinetics data for paper"""

    p = Pipeline()
    colheaderlabels = 'wt 5,wt 3,wt 2,68 5,68 3,68 2,138 5,138 3,138 2,248 5,248 3,248 2'
    rowheaderlabels = '3.2,1.6,0.8,0.4,0.2,0.1,0.05,0.025'
    conf = {'format':'kineticsdata', 'delimeter':'tab','rowstart':3,'colend':12,
            'rowrepeat':9,
            'colheaderlabels':colheaderlabels, 'rowheaderlabels':rowheaderlabels,
            'decimalsymbol':',','xformat':'%M:%S',
            'groupbyname':1, 'parsenamesindex':2,  'parsemethod':'numeric',
            'model1':'linear','model2':'Michaelis-Menten','model3':'1 pKa 2 Chemical shifts',
            'variable1':'a','variable2':'Km','variable3':'pKa',#'xerror':.1,'yerror':0.05,
            }
    p = Pipeline()
    p.createConfig('temp.conf',**conf)
    path = 'testfiles/kineticsdata/jan/setF/MM/rep1'
    p.addFolder(path)
    p.run()
    return

def preProcessingTest():
    """Test processing steps like differentation of the data"""
    path = "testfiles"
    names = Utilities.createRandomStrings(8,6)
    fname = os.path.join(path,'preprocessingtest.txt')
    Utilities.createCDData(fname, names, 300, .5)
    conf = {'format':'databycolumn','model1':'gaussian',
            'function1':'differentiate','function2':'gaussiansmooth',
            'iterations':100,
            'variable1':'a','saveplots':1}
    p = Pipeline()
    p.createConfig('temp.conf',**conf)
    p.openRaw(fname)
    p.run()
    return

def peakDetectionTest(path=None, noise=0.08):
    """Use pre-processing funcs to detect peaks"""
    if path == None:
        path = "testfiles"
    names = Utilities.createRandomStrings(8,6)
    fname = os.path.join(path,'spectraldata.txt')
    peaks = Utilities.createSimulatedSpectralData(fname, names, noise=noise)
    conf = {'format':'databycolumn', 'saveplots':1, 'marker':'-',
            'markers':'-,x','alpha':0.7,'normalise':1,
            'function1':'smooth', 'function2':'baselinecorrection',
            'function3':'detectpeaks'}
    p = Pipeline()
    p.createConfig('temp.conf',**conf)
    p.openRaw(fname)
    results = p.run()

    #compare predicted peaks
    successrates=[]
    res  = results[results.keys()[0]]
    for name in peaks:
        #print name, sorted(peaks[name]), results[name][0]
        orig = set(peaks[name])
        pred = set(res[name][0])
        s = float(len(orig.intersection(pred)))/(len(orig))
        successrates.append(s)
    return np.mean(successrates)

def renameFilesTest():
    import glob
    path = "testfiles/renametest"
    Utilities.createDirectory(path)
    for fname in Utilities.createRandomStrings(5,5):
        f=(open(os.path.join(path,fname)+'.txt','w'))
        f.close()
    print glob.glob(os.path.join(path,'*.txt'))
    import Rename
    B=Rename.doFindReplace(wildcard=os.path.join(path,'*.txt'),
                      find=".", replace='_',rename=True)
    print glob.glob(os.path.join(path,'*'))
    return

def main():

    #formatTests(basictests)
    #formatTests(exceltests)
    multiFileTest()
    #groupedFilesTest()
    #multiFolderTest()
    #replicatesTest()
    #fitPropagationTest()
    #groupbyFieldsTest()
    #kineticsTest()
    #preProcessingTest()
    #peakDetectionTest()
    #renameFilesTest()

if __name__ == '__main__':
    main()
