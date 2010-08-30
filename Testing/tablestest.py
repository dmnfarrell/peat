#!/bin/env python
#
# This file is part Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#
# Author: Damien Farrell May 2010

"""Table and labbook testing"""

import os, random, string
from PEATDB.Labbook import LabbookTableModel

def simulate_sequence(length):
    """simulate a dna seq"""
    dna = ['A', 'C', 'G', 'T']
    sequence = ''
    for i in range(length):
        sequence += random.choice(dna)
    return sequence

def createTable():
    L = LabbookTableModel()
    for r in range(5000):
        c=''
        for k in range(8):
            c += random.choice(string.letters)        
        L.addRecord(name=c,
                    seq=simulate_sequence(20),
                    stab=str(round(random.normalvariate(1,2),3)),
                    choice='a')
            
    print L
    #save as labbook
    saveasLabbook(L,'large.labbook')
    return

def saveasLabbook(model,filename):
    data={}
    data['test'] = model.getData()  
    fd=open(filename,'w')
    import pickle
    pickle.dump(data,fd)
    fd.close()
    return

def listTable():
    
    return

if __name__ == '__main__':
    path=os.getcwd()
    createTable()
