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

import os, random, string
import re, glob
import csv
import ConfigParser
from math import *
import numpy as np

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
        #print x, val
        vals = [round(val+random.normalvariate(0,noise),2) for j in range(len(names))]
        vals.insert(0,x)
        cw.writerow(vals)
        
    return

def createSingleFileData(path='testfiles', clear=False):
    """Create sets of individual data files all in one folder,
       one per xy datasew with multiple copies for each label representing ph values"""

    createDirectory(path)
    names = Utilities.createRandomStrings(8,6)
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
        names = Utilities.createRandomStrings(3,6)
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

def differentiate(self, x,y):
    dy = numpy.diff(y,1)
    dx = x[:len(dy)]
    return dx,dy

