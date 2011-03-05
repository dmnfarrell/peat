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

"""This script will create multiple projects from csv files and
add pdbs based on the csv names. It can also create peatsa jobs
and merge them back into the database"""

import pickle, sys, os, copy, time, types
import numpy
from PEATDB.Base import PDatabase 
from PEATDB import Utils
from PEATDB.Actions import DBActions
from PEATDB.plugins.PEATSAplugin import PEATSAPlugin
from PEATDB.plugins.Correlation import CorrelationAnalyser
from PEATDB.PEATTables import PEATTableModel
import matplotlib.pyplot as plt

#plt.rc('text',usetex=True)
plt.rc('font',size=8)
plt.rc('savefig',dpi=300)

settings={'server':'peat.ucd.ie','username':'guest',
           'password':'123'}
path = '/home/people/farrell/Desktop/SADBPaperData'
savepath = os.path.join(path,'projects')
cpath = os.path.join(path,'data')
csvfiles = os.listdir(cpath)#[:4]
dbnames = [os.path.splitext(i)[0] for i in csvfiles]
print dbnames

    
def PEATSAJobs(prjs):
    """Submit PEATSA runs for all projects or merge results if done"""
    for name in prjs:
        print name
        DB = PDatabase(local=os.path.join(savepath,name))
        pdb = DB['wt'].Structure
        PS = PEATSAPlugin()
        PS.main(DB=DB)        
        if hasattr(DB.meta,'peatsa_jobs'):
            if 'mycalc' in DB.meta.peatsa_jobs:
                print 'job already present'
                #try to merge results
                S = PEATTableModel(DB)
                job,n = PS.getJob('mycalc')
                PS.mergeResults(job, 'prediction', S)
                DB.commit()
                print 'merged results'
        else:
            mutlist = []
            for p in DB.getRecs():
                mutlist.append(DB.get(p).Mutations)
            #print mutlist
            pdbfile = PS.writetempPDB()
            PS.submitJob(name='mycalc', pdbname=DB.meta.refprotein, pdbfile=pdbfile, 
                         mutations=mutlist, calcs=['stability'], meta={'protein':name})
            #required to end process
        PS.jobManager.stopLogging()
        DB.close()
    return

def createProjects(files):
    """Create multiple projects at once from csv files"""

    for filename in files:
        print filename
        name = os.path.splitext(filename)[0]
        #create/open db
        DB = PDatabase(local=os.path.join(savepath,name))
        DB.add('wt')
        #add wt pdb
        stream = DBActions.fetchPDB(name)
        DBActions.addPDBFile(DB, 'wt', pdbdata=stream, pdbname=name, gui=False)
        DB.meta.refprotein = 'wt'
        DB.meta.info['protein'] = name
        #import data from csv
        DB.importCSV(os.path.join(cpath,filename), namefield='Mutations')
        print 'imported ok'
        DB.deleteField('PDB')
        DB.commit()
        DB.close()
        print 'done'
    return

def summarise(projects):
    summDB = PDatabase(local='summary.fs')
    C = CorrelationAnalyser()
    fig = plt.figure()
    from mpl_toolkits.axes_grid1 import AxesGrid
    grid = AxesGrid(fig, 111, 
                    nrows_ncols = (5, 5),
                    share_all=True,
                    #label_mode = "1",      
                    axes_pad = 0.3)
    i=0
    data=[]
    for p in projects:
        print p
        DB = PDatabase(local=os.path.join(savepath,p))
        S = PEATTableModel(DB)
        #print DB.meta.info
        try:
            exp,pre = S.getColumns(['Exp','prediction'],allowempty=False)
        except:
            print 'no results'
            continue
        #print exp,pre
        ax,fr,mh = C.plotCorrelation(pre,exp,title=p,ms=2,ax=grid[i])        
        cc,rmse = C.getStats(pre,exp)
        print 'rmse',rmse
        data.append({'name':p,'rmse':rmse,'cc':cc})
        #,'Structure':DB['wt'].Structure})
        i+=1        
    
    summDB.importDict(data)    
    summDB.commit()
    fig.savefig('test.png')
    plt.show()
    return

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option("-c", "--create", dest="create", action='store_true',
                       help="create/import", default=False)
    parser.add_option("-j", "--jobs", dest="jobs", action='store_true',
                       help="do/merge jobs", default=False)
    parser.add_option("-s", "--summary", dest="summary", action='store_true',
                       help="do summary", default=False)
    opts, remainder = parser.parse_args()
    if opts.create == True:
        createProjects(csvfiles)
    if opts.jobs == True:    
        PEATSAJobs(dbnames)
    if opts.summary == True:
        summarise(dbnames)
  
