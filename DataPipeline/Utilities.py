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

'''Module containing utility classes and functions'''

import os, random, string, types
import re, glob
import csv
import ConfigParser
from math import *
import numpy as np
from PEATDB.Ekin.Base import EkinProject, EkinDataset

def setAttributesfromConfigParser(obj, cp):
    """A helper method that makes the options in a ConfigParser object
       attributes of obj"""

    for s in cp.sections():
        obj.__dict__[s] = cp.items(s)
        for f in cp.items(s):
            try: val=int(f[1])
            except: val=f[1]
            obj.__dict__[f[0]] = val

def createConfigParserfromDict(data, sections, **kwargs):
    """Helper method to create a ConfigParser from a dict and/or keywords"""

    cp = ConfigParser.ConfigParser()
    for s in sections:
        cp.add_section(s)
        if not data.has_key(s):
            continue
        for i in data[s]:
            name,val = i
            cp.set(s, name, val)

    #use kwargs to create specific settings in the appropriate section
    for s in cp.sections():
        opts = cp.options(s)
        for k in kwargs:
            if k in opts:
                cp.set(s, k, kwargs[k])
    #handle model and variable sections which can have zero or multiple
    #options
    for k in sorted(kwargs):
        if k.startswith('model'):
            cp.set('models', k, kwargs[k])
        elif k.startswith('variable'):
            cp.set('variables', k, kwargs[k])
        elif k.startswith('function'):
            cp.set('functions', k, kwargs[k])
    return cp

def getListFromConfigItems(items):
    """Get a list from a set of ConfigParser key-value pairs"""
    lst = [i[1] for i in sorted(items) if i[1] != '']
    return lst

def clearDirectory(path):
    """Remove all files in folder"""
    for f in os.listdir(path):
        filepath = os.path.join(path, f)
        try:
            if os.path.isfile(filepath):
                os.unlink(filepath)
        except Exception, e:
            print e
    return

def createDirectory(path):
    """Create or clear a directory"""
    if not os.path.exists(path):
        os.mkdir(path)
    else:
        clearDirectory(path)

def getdirectoryStructure(rootdir):
    """Creates a nested dictionary that represents the folder structure
       of rootdir
       taken from http://code.activestate.com/recipes/577879/"""

    D = {}
    rootdir = rootdir.rstrip(os.sep)
    start = rootdir.rfind(os.sep) + 1
    for path, dirs, files in os.walk(rootdir):
        folders = path[start:].split(os.sep)
        filenames = [os.path.join(path,f) for f in files]
        subdir = dict.fromkeys(filenames)
        parent = reduce(dict.get, folders[:-1], D)
        parent[folders[-1]] = subdir
    return D

def createRandomStrings(l,n):
    """create list of l random strings, each of length n"""

    names = []
    for i in range(l):
        val = ''.join(random.choice(string.ascii_uppercase) for x in range(n))
        names.append(val)
    return names

def createTempData(fname, names, slope, noise=0.2):
    """Create some simulated linear data with noise as a function of temp"""

    cw = csv.writer(open(fname,'w'))
    cw.writerow(['temp']+names)
    for x in range(250,360,2):
        vals = [round(slope*x/random.normalvariate(10,noise),2) for j in range(len(names))]
        vals.insert(0,x)
        cw.writerow(vals)
    return

def createCDData(fname, names, tm, noise=0.2):
    """Create some simulated cd data with noise as a function of temp"""

    cw = csv.writer(open(fname,'w'))
    cw.writerow(['temp']+names)
    for x in range(250,350,1):
        val = 10/(1+exp((tm-x)/2))
        vals = [round(val+random.normalvariate(0,noise),2) for j in range(len(names))]
        vals.insert(0,x)
        cw.writerow(vals)
    return

def createSimulatedSpectralData(fname, names, peaks=10, noise=0.08):
    """Create spectral type data with noise and peaks that you might find
       commonly in experimental data"""

    cw = csv.writer(open(fname,'w'))
    cw.writerow(['temp']+names)
    offset=100
    n=800
    pheight = offset*2
    percnoise = offset*noise
    baselinefunc = lambda x: offset + pheight*pow(.98/x,0.2)

    def gaussianpeaks(x, peaks):
        #add random gaussian-shaped peaks as signals
        res=0
        for i in peaks:
           height = random.normalvariate(pheight,pheight/5)
           res+= height*exp(-(i-x)**2/1.3)
        return res

    data={}
    freqs = [i+100 for i in range(0,n,1)]
    peakdata = {}
    for name in names:
        peakvals = [int(random.normalvariate(n/2,n/4)) for p in range(peaks)]
        print name, sorted(peakvals)
        vals=[]
        for i in freqs:
            val = gaussianpeaks(i, peakvals)
            val += baselinefunc(i)+random.normalvariate(0,percnoise)
            vals.append(val)
        data[name] = vals
        peakdata[name] = peakvals
    for x in range(len(freqs)):
        row=[]
        for name in names:
            row.append(data[name][x])
        row.insert(0,freqs[x])
        cw.writerow(row)
    return peakdata

def createSingleFileData(path='testfiles', clear=False):
    """Create sets of individual data files all in one folder,
       one per xy datasew with multiple copies for each label representing ph values"""

    createDirectory(path)
    names = createRandomStrings(3,6)
    reps = range(1,2)
    for name in names:
        for i in np.arange(2,10,1.0):
            for rep in reps:
                val = 1/(1+exp((i-5)/1.2))
                fname = os.path.join(path, name+'_ph'+str(i)+'_rep'+str(rep)+'.txt')
                createTempData(fname, ['ph'], val)
    return

def createGroupedData(path='testfiles', clear=False, names=None):
    """Create sets of grouped data to test queuing and file grouping"""

    createDirectory(path)
    if names == None:
        names = createRandomStrings(3,6)
    for i in np.arange(2,10,1.0):
        val = 1/(1+exp((i-5)/1.2))
        fname = os.path.join(path,'ph_'+str(i)+'__xx_.txt')
        createTempData(fname, names, val, noise=0.5)
    return

def createNestedData():
    """Create simulated nested data similar to our kinetics data"""
    data={}
    names = ['aaa','bbb','ccc']
    tms = [7.5,8.0,8.5]
    phs = range(2,14) #e.g. a ph range
    concs = [0.01,0.1,0.2,0.3,0.5,1.0,2.0,5.0]
    x = range(0,100,10)

    for n in names:
        i=names.index(n)
        tm=tms[i]#+random.normalvariate(0,0.1)
        data[n]={}
        for p in phs:
            km = 2/(1+exp((p-tm)/1.2))
            #print p,km
            data[n][p]={}
            for s in concs:
                vel = (10 * s)/(s + km)
                y = [round(i*vel,3) for i in x]
                y = [i * random.normalvariate(1,0.03) for i in y]
                data[n][p][s] = (x,y)
    return data

def getEkinProject(data, xerror=None, yerror=None, sep='__'):
    """Get an ekin project from a dict of the form
         {label:([x],[y]),..} or
         {label:([x],[y],[xerr],[yerr]),..}"""

    E = EkinProject(mode='General')
    for d in data.keys():
        if type(data[d]) is types.DictType:
            for lbl in data[d]:
                name = str(d)+sep+str(lbl)
                xy = data[d][lbl]
                ek=EkinDataset(xy=xy)
                E.insertDataset(ek, name)
        else:
            #print data[d]
            if len(data[d]) == 4:
                x,y,xerrs,yerrs = data[d]
            else:
                x,y = data[d]
                xerrs = []; yerrs=[]
                if xerror!=None:
                    xerrs=[xerror for i in x]
                if yerror!=None:
                    yerrs=[yerror for i in y]
            ek = EkinDataset(xy=[x,y], xerrs=xerrs, yerrs=yerrs)
            E.insertDataset(ek, d)
            #print ek.errors
    return E

def getFits(E, model, varname='a', filename=None,
            iterations=50, xerror=0, yerror=0):
    """Fit an Ekin project
       model: model to fit
       varname: variable to extract"""

    fits = []
    xerrors = []
    yerrors = []
    E.fitDatasets('ALL', models=[model], noiter=iterations,
                  conv=1e-9, grad=1e-8, silent=True)
    labels = E.datasets
    if xerror != 0 or yerror != 0:
        print 'getting exp uncert'
        expUncert(E, xerr=float(xerror))
        for d in labels:
            yerrors.append(E.getMeta(d,'exp_errors')[varname][1])

    for d in labels:
        fits.append(E.getMetaData(d)[varname])
    '''if self.saveplots == 1 and filename != None and filename != '':
        print 'plotting %s' %filename
        self.saveEkinPlotstoImages(E, filename)'''
    return E,(labels,fits,xerrors,yerrors)

def expUncert(E, xerr=0, yerr=0):
    """Get fitted variable errors using the iterative method"""

    for d in E.datasets:
        ferrs = E.estimateExpUncertainty(d, runs=10, xuncert=xerr, yuncert=yerr)
        E.addMeta(d, 'exp_errors', ferrs)
    E.saveProject()
    return

def parseFileNames(filenames, ind=0, sep='', match='numeric'):
    """Parse file names"""
    filenames = [os.path.basename(f) for f in filenames]
    if ind == '':
        labels = None
    else:
        labels = parseNames(filenames, ind, sep, match)
    return labels

def parseNames(filenames, ind=0, sep='', match='numeric'):
    """Parse a list of strings to extract a numerical value
       match: whether to match numeric, string or whole word values
       ind: extract the ith instance of a number in the filename
       returns: a dict of the original keys with labels as items
       """

    labels = {}
    if match == 'numeric':
        expr = r'([0-9.]*[0-9]+)'
    elif match == 'text':
        expr = r'[a-z]*[A-Z]+'
    elif match == 'both':
        expr = r'([a-z]*[A-Z]+?)|([0-9.]*[0-9]+?)'
    elif match == 'words':
        expr = r'^\w+'

    #print expr
    r = re.compile(expr,re.I)
    for f in filenames:
        vals = r.findall(f)
        if match == 'both':
            result = []
            for v in vals:
                for i in v:
                    if i != '':
                        result.append(i)
        else:
            result = vals
        labels[f] = result[ind]
        #print result
        #print ind, f, labels[f]
    return labels

def findNumeric(line):
    """Parse string to extract all numerical values"""
    return re.findall("([0-9.]*[0-9]+)", line)

def findWords(line):
    """Parse string to extract all non-numeric strings"""
    return re.findall(r'[a-z]*[A-Z]+', line, re.I)

def findAlphanumeric(line):
    """Parse string to extract all non-numeric strings"""
    return re.findall(r'^\w+', line)

def addReplicates(data, labels):
    """If the configuration specifies replicates than we average the
       corresponding data points in the raw dict """

    print 'processing replicates..'
    import operator
    from itertools import groupby
    newdata = {}

    sorteditems = sorted(labels.iteritems(), key=operator.itemgetter(1))
    for key, group in groupby(sorteditems, lambda x: x[1]):
        c = 0
        subdata = []
        for g in group:
            f = g[0]    #key in corresponding data dict
            c+=1
            subdata.append(data[f])
            #print f
        newdata[key] = averageDicts(subdata)
        print '%s replicates for label %s' %(c, key)

    #print newdata
    return newdata

def averageDicts(dictslist):
    """Average dicts of the form
       {label1: [[x1],[y1]],..],label2:[[x2],[y2]]...}"""

    newdata = {}
    names = dictslist[0].keys()
    for n in names:
        arrs = []
        for D in dictslist:
            arrs.append(np.array(D[n]))
        newdata[n] = list(sum(arrs)/len(arrs))
    return newdata

def arrangeDictbySecondaryKey(data, labels, sep='__'):
    """Re-arrange a dict of dicts by each records secondary keys"""

    newdata = {}
    if not type(data) is types.DictType:
        return data
    key1 = data.keys()[0]
    if type(data[key1]) is types.DictType:
        fields = data[key1].keys()
    for f in fields:
        newdata[f] = {}
    for f in fields:
        for d in data.keys():
            #print d
            lbl = labels[d]
            if data[d].has_key(f):
                xy = data[d][f]
                newdata[f][lbl] = xy
    return newdata

def extractSecondaryKeysFromDict(data):
    """Re-arrange a dict of the form {key1:([names],[values]),..
       into the form {name1:([keys],[values]),..}"""

    newdata = {}
    kys = data.keys()
    datasets = data[kys[0]][0]
    for d in datasets:
        newdata[d] = []
    for l in data:
        names, vals, xerr, yerr = data[l]
        #print zip(names, vals)
        for d,v in zip(names, vals):
            newdata[d].append((l,v))
    for d in newdata:
        newdata[d] = zip(*newdata[d])
    return newdata

def buildNestedStructure(data, labels1, labels2):
    """Rebuild a nested dictionary from a flat sructure according to labels"""
    newdata = {}
    for key in data:
        pri = labels1[key]
        sec = labels2[key]
        print key, pri, sec
        if not newdata.has_key(pri):
            newdata[pri] = {}
            if len(data[key])==1:
                item = data[key][data[key].keys()[0]]
            else:
                item = data[key]
            newdata[pri][sec] = item
    return newdata

def differentiate(self, x,y):
    dy = numpy.diff(y,1)
    dx = x[:len(dy)]
    return dx,dy

