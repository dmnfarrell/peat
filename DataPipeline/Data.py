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

import os, sys, copy, types
import Utilities, random

class Dataset(object):
    """This class stores a multiply nested dictionary and it's properties such the
       levels of nesting. It also allows the dictionary to be re-arranged by secondary
       keys. The methods assume that the lowest level children are lists or tuples"""

    def __init__(self, data=None):
        if data != None:
            self.data = data
        else:
            self.data = {}
        return

    def add(self, rec, name):
        return

    def levels(self):
        """Level of nesting"""
        return self.getDictNesting(self.data)

    def getDictNesting(self, data, level=1):
        """Get level of nesting for a dictionary"""
        keys = data.keys()
        if len(keys) == 0:
            return 0
        if type(data[keys[0]]) is types.DictType:
            return self.getDictNesting(data[keys[0]], level+1)
        else:
            return level

    def __setitem__(self, key, data):
        self.data[key] = data

    def __repr__(self):
        if len(self.data)==0:
            s = 'dataset with no data'
        else:
            s = 'dataset with %s primary keys and ' %len(self.data.keys())
            s += '%s levels of nesting' %self.levels()
        return s

    def show(self):
        self.printStructure(self.data)

    def printStructure(self, data):
        if type(data) is types.DictType:
            keys = data.keys()
            print keys
            self.printStructure(data[keys[0]])
        else:
            print type(data)
        return

    def fit(self, model='linear', param='a'):
        """Fits the bottom level of data and returns the dataset"""
        import PEATDB.Ekin.Fitting as Fitting
        data = self.data
        keys = data.keys()
        #for k in keys:

        #E = Utilities.getEkinProject(rawdata, xerror=xerror, yerror=yerror)
        #E,fit = self.getFits(E, currmodel, currvariable, str(parentkey))
        #Em.addProject(E, label=parentkey)'''

def getKeys(keytypes=0):
    if keytypes == 0:
        keytypes = random.randint(1,2)
    if keytypes == 1:
        keys = Utilities.createRandomStrings(3,6)
    elif keytypes == 2:
        keys = [i for i in range(random.randint(1,4),random.randint(6,9))]
    return keys

def tempData():
    slope=0.2
    noise=0.01
    x = range(250,360,10)
    vals = [round(slope*i/random.normalvariate(10,noise),2) for i in x]
    return x,vals

def createNestedData(level=1, current=None, keytypes=0, keys={}):
    """create test nested data"""
    data = {}
    if level==0:
        data = tempData()
        print type(data),level
    else:
        if not keys.has_key(level):
            keys[level] = getKeys()
        print keys[level],level
        for k in keys[level]:
            data[k] = createNestedData(level=level-1, keytypes=keytypes, keys=keys)
    return data

def test():
    data = createNestedData(2)
    #print data
    D = Dataset(data)
    print D
    D.show()
    #print D.levels()

if __name__ == '__main__':
    test()
