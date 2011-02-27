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
# Author: Damien Farrell 2011

import pickle, sys, os, copy, time, types
import numpy
from PEATDB.Base import PDatabase 
from PEATDB import Utils
from PEATDB.plugins.PEATSAplugin import PEATSAPlugin
from PEATDB.plugins.Correlation import CorrelationAnalyser
from PEATDB.PEATTables import PEATTableModel
import matplotlib.pyplot as plt

settings={'server':'peat.ucd.ie','username':'guest',
           'password':'123'}
path = '/home/farrell/Desktop/SADBPaperData/projects'
dbnames = ['1a2p.fs','1bf4.fs']

def createStabilityJobs():
    """do PEATSA runs for all projects"""
    return

def summarise(projects):
    summDB = PDatabase(local='summary.fs')
    C = CorrelationAnalyser()
    fig = plt.figure()
    for p in projects:
        DB = Utils.loadDB(os.path.join(path,p), remote=False)
        S = PEATTableModel(DB)
        print DB.meta.info
        print DB.meta.userfields
        exp,pre = S.getColumns(['Exp','predictions'],allowempty=False)
        print exp,pre
        ax,fr,mh = C.plotCorrelation(pre,exp,title=p)
        fig.add_subplot(111)
        fig.savefig('test.png')    
    #print S.projects
    #summDB.commit()
    return

if __name__ == '__main__':
    
    summarise(dbnames)
  
