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
import csv

def setAttributesfromConfigParser(obj, cp):
    """A helper method that makes the options in a ConfigParser object
       attributes of obj"""
   
    for s in cp.sections():
        obj.__dict__[s] = cp.items(s)       
        for f in cp.items(s):        
            try: val=int(f[1])
            except: val=f[1]
            obj.__dict__[f[0]] = val

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
    
def createTempData(fname, names, slope):
    """Create some simulated linear data with noise as a function of temp"""
    
    cw = csv.writer(open(fname,'w'))
    cw.writerow(['temp']+names)
    for x in range(250,360,2):
        vals = [round(slope*x/random.normalvariate(10,0.5),2) for j in range(len(names))]
        vals.insert(0,x)
        cw.writerow(vals)
    return
