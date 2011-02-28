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
from PEATDB.Actions import DBActions
from PEATDB.plugins.PEATSAplugin import PEATSAPlugin
from PEATDB.plugins.Correlation import CorrelationAnalyser
from PEATDB.PEATTables import PEATTableModel
import matplotlib.pyplot as plt

settings={'server':'peat.ucd.ie','username':'guest',
           'password':'123'}
path = '/home/people/farrell/Desktop/SADBPaperData/'
savepath = os.path.join(path,'projects')
dbnames = ['1a2p.fs','1bf4.fs']

def submitPEATSAJobs(prjs):
    """do PEATSA runs for all projects"""
    for name in prjs:
        DB = PDatabase(local=os.path.join(savepath,name))        
        PS = PEATSAPlugin()
        PS.main(DB=DB)
        pdb = DB['wt'].Structure
        mutlist = []
        for p in DB.getRecs():
            mutlist.append(DB.get(p).Mutations)
        #print mutlist
        pdbfile = PS.writetempPDB()    
        PS.submitJob(name='mycalc', pdbname=DB.meta.refprotein, pdbfile=pdbfile, 
                     mutations=mutlist, calcs=['stability'], meta={'protein':name})
        PS.jobManager.stopLogging()
        
    return

def createProjects():
    """Create multiple projects at once from csv files"""

    csvfiles = os.listdir(path)[:2]    
    for filename in csvfiles:
        print filename
        name = os.path.splitext(filename)[0]
        #create/open db
        DB = PDatabase(local=os.path.join(savepath,name))
        DB.add('wt')
        #add wt pdb
        stream = DBActions.fetchPDB(name)
        DBActions.addPDBFile(DB, 'wt', pdbdata=stream, pdbname=name, gui=False)        
        DB.meta.refprotein = 'wt'
        #import data from csv
        DB.importCSV(os.path.join(path,filename), namefield='Mutations')
        #DB.deleteField('PDB')
        DB.commit()
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
    #createProjects()
    #submitPEATSAJobs(['1a2p','2chf'])
    summarise(dbnames)
  
