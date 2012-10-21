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
from math import *
import numpy as np


class NestedData(object):
    """This class manages a multiply nested dictionary and derives properties such the
       levels of nesting. It also allows the dictionary to be re-arranged by secondary
       keys. The methods assume that the lowest level children are lists or tuples"""

    def __init__(self, data=None):
        if data != None:
            self.data = data
        else:
            self.data = {}
        self.initialdata = copy.deepcopy(data)
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
        """Preview the dict structure"""
        if type(data) is types.DictType:
            keys = data.keys()
            print keys
            self.printStructure(data[keys[0]])
        else:
            print type(data)
        return

    def restoreData(self):
        """Reset data to the original form"""
        self.data = copy.deepcopy(self.initialdata)
        return

    def arrangeDictbySecondaryKey(self):
        self.getLabels()
        self.data = Utilities.arrangeDictbySecondaryKey(self.data, self.namelabels)
        print self.data
        return

    def buildNestedStructure(self, parsenamesindex=[]):
        """Rebuild a nested dictionary from a flat sructure according to labels
           extracted from the key names, which must be separated by some symbol"""

        data = self.data
        keys = data.keys()
        labels = []
        #get set of labels in the correct order
        for i in parsenamesindex:
            labels.append(Utilities.parseNames(keys, ind=i, sep='', match='numeric'))

        self.data = self.recursiveBuild(labels, data, 0)
        return

    def recursiveBuild(self, labels, data, i):
        """Recursively build a new dictionary based on the labels extracted
           from a flat dict keys, this is order dependant"""

        names = labels[i]
        newdata = {}
        if i==len(labels)-1:
            for key in data:
                n = names[key]
                #print key, n
                #tuple expected
                if type(data[key]) is types.DictType:
                    item = data[key][data[key].keys()[0]]
                else:
                    item = data[key]
                newdata[n] = item
        else:
            for j in names:
                n = names[j]
                newdata[n] = self.recursiveBuild(labels, data, i+1)
        return newdata

    def flatten(self):
        """Flatten a nested dictionary by creating one key from each
           childs sub-keys"""
        return

    def averageReplicates(self, level=0):
        """Average replicates"""
        #self.data = Utilities.addReplicates(self.data)
        data = self.data
        newdata = {}
        #labels = self.parseNames(keys, ind=0, match='words')
        print 'processing replicates..'

        self.data = newdata
        return

def getKeys(keytypes=0):
    if keytypes == 0:
        keytypes = random.randint(1,2)
    if keytypes == 1:
        keys = Utilities.createRandomStrings(3,6)
    elif keytypes == 2:
        keys = [i for i in range(random.randint(1,4),random.randint(6,9))]
    return keys

def tempData(slope=0.2,noise=0.01):
    x = range(250,360,10)
    vals = [round(slope*i/random.normalvariate(10,noise),2) for i in x]
    return x,vals

def createNestedData(level=1, current=None, keytypes=0, keys={}):
    """create test nested data"""
    data = {}
    if level==0:
        data = tempData()
        #print type(data),level
    else:
        if not keys.has_key(level):
            keys[level] = getKeys()
        #print keys[level],level
        for k in keys[level]:
            data[k] = createNestedData(level=level-1, keytypes=keytypes, keys=keys)
    return data

def multiFileData():
    """Simulate a dict representing imported multiple file data with replicates"""
    data={}
    names = Utilities.createRandomStrings(3,6)
    reps = range(1,3)
    for name in names:
        for rep in reps:
            for i in np.arange(2,10,1.0):
                key = name+'_ph'+str(i)+'_rep'+str(rep)+'.txt'
                val = 1/(1+exp((i-5)/1.2))
                data[key] = {}
                data[key]['ph'] = tempData(val)
    #print data
    return data

def test():
    #data = createNestedData(2)
    data = multiFileData()
    #print data
    D = NestedData(data)
    D.buildNestedStructure([0],[0,1])
    print D
    D.show()
    D.averageReplicates()

if __name__ == '__main__':
    test()
