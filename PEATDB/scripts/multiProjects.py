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

import pickle, sys, os, copy, time, types, math
import numpy
from PEATDB.Base import PDatabase 
from PEATDB import Utils
from PEATDB.Actions import DBActions
from PEATDB.plugins.PEATSAplugin import PEATSAPlugin
from PEATDB.plugins.Correlation import CorrelationAnalyser
from PEATDB.PEATTables import PEATTableModel
import PEATDB.Utils
from PEATDB.Parsers import PDBParser
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import stats

#plt.rc('text',usetex=True)
plt.rc('font',size=7)
plt.rc('legend',fontsize=6)
plt.rc('savefig',dpi=300)
plt.rc('axes',linewidth=.5)

settings={'server':'enzyme.ucd.ie','username':'guest',
           'password':'123'}
#path = '/home/people/farrell/Desktop/SADBPaperData'
path = os.getcwd()
savepath = os.path.join(path,'projects')
cpath = os.path.join(path,'data')
csvfiles = os.listdir(cpath)#[:4]
dbnames = [os.path.splitext(i)[0] for i in csvfiles]
print dbnames


def PEATSAJobs(prjs, resubmit=False):
    """Submit PEATSA runs for all projects or merge results if done"""
    for name in prjs:
        print name
        DB = PDatabase(local=os.path.join(savepath,name))
        pdb = DB['wt'].Structure
        PS = PEATSAPlugin()
        PS.main(DB=DB)
        if hasattr(DB.meta,'peatsa_jobs') and resubmit == False:
            if 'mycalc' in DB.meta.peatsa_jobs:
                print 'job is present'
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
            #we add source project data so exp data can be read from summary
            prjdata = {'server':'enzyme.ucd.ie','username':'guest',
                              'project':name,'password':'123','port':'8080'}
            PS.submitJob(name='mycalc', pdbname=DB.meta.refprotein, pdbfile=pdbfile, 
                         mutations=mutlist, calcs=['stability'],
                         meta={'protein':name,'expcol':'Exp','project':prjdata})
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
    figs = []
    for f in range(4):
        figs.append(plt.figure())
    
    gs = gridspec.GridSpec(5, 5, wspace=0.3, hspace=0.5)    
    i=0
    data=[]    

    for p in projects:
        print 'structure:',p
        DB = PDatabase(local=os.path.join(savepath,p))
        S = PEATTableModel(DB)           
        
        try:
            exp,pre = S.getColumns(['Exp','prediction'],allowempty=False)
            errs = [j[0]-j[1] for j in zip(exp,pre)]
        except:
            print 'no results'
            continue
            
        #DB.close()
        #add link to proj
        summDB.add(p)
        summDB.addField('project',fieldtype='Project')
        summDB[p]['project'] = {'server':'enzyme.ucd.ie','username':'guest',
                              'project':p,'password':'123','port':'8080'}              
        #stats
        cc,rmse,meanerr = C.getStats(pre,exp)
        #ttest for mean errs 0        
        ttp = round(stats.ttest_1samp(errs, 0)[1],2)
        #normality of errs
        w,swp = C.ShapiroWilk(errs)
        '''ax = figs[0].add_subplot(gs[0, i])
        C.plotCorrelation(pre,exp,title=p,ms=2,axeslabels=False,ax=ax)
        ax = figs[1].add_subplot(gs[0, i])
        C.showHistogram([pre,exp],title=p,labels=['pre','exp'],ax=ax)                
        ax = figs[2].add_subplot(gs[0, i])
        C.plotNorm(errs,title=p,lw=1,ax=ax)
        #qqplot
        ax = figs[3].add_subplot(gs[0, i])
        C.QQplot(errs,title=p,ax=ax)'''
        x={'name':p,'mutants':len(pre),'rmse':rmse,'corrcoef':cc,'meanerr':meanerr,
           'ttest':ttp,'shapirowilk':swp}
        #print x
        parser = PDBParser()
        descr = parser.getDescription(p)
        x.update(descr)
        data.append(x)
        i+=1
        print summDB.isChanged()
        print summDB[p]['project']
        
    summDB.importDict(data)    
    summDB.commit()

    #add all peatsa jobs to summary proj also
    '''print 'adding peatsa job info'
    PS = PEATSAPlugin()
    PS.main(DB=summDB)
    #summDB.meta.peatsa_jobs = None
    #from ZODB.PersistentMapping import PersistentMapping
    #summDB.meta.peatsa_jobs = PersistentMapping()    
    PS.checkJobsDict()
    PS.jobManager.stopLogging()
    for p in projects:
        #print summDB.meta
        DB = PDatabase(local=os.path.join(savepath,p))
        job = DB.meta.peatsa_jobs['mycalc']
        summDB.meta.peatsa_jobs[p] = job
        print job
        #DB.close()
    print summDB.isChanged()
    print summDB.meta.peatsa_jobs
    summDB.commit()'''

    #for i in range(len(figs)):
    #    figs[i].savefig('fig%s.png' %i)
    #plt.show()
        
    return

def findOutliers(data):
    """Outliers in all corr data"""
    C = CorrelationAnalyser()   
    
    return ax

def send2Server(projects):
    """Send all projects to remote versions"""
    settings={'server':'enzyme.ucd.ie','username':'guest',
               'password':'123','port':8080}
    adminsettings={'host':'localhost','user':'peatadmin',
               'passwd':'123','port':8080}    
    '''for p in projects:
        print p
        DB = PDatabase(local=os.path.join(savepath,p))        
        #Utils.createDBonServer(prj=p,settings=adminsettings,
        #                       access='guest')
        Utils.copyDBtoServer(DB,p,settings)'''
        
    DB = PDatabase(local='summary.fs')
    Utils.copyDBtoServer(DB,'PotapovDataset',settings)
    return
    
if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()    
    parser.add_option("-i", "--importcsv", dest="importcsv", action='store_true',
                       help="create/import", default=False)
    parser.add_option("-j", "--jobs", dest="jobs", action='store_true',
                       help="submit/merge jobs", default=False)
    parser.add_option("-s", "--summary", dest="summary", action='store_true',
                       help="do summary/stats", default=False)  
    parser.add_option("-p", "--path", dest="path",
                        help="Path with csv files")
    parser.add_option("-c", "--copy", dest="copy",action='store_true',
                        help="copy to server", default=False)    
    opts, remainder = parser.parse_args()

    if opts.path != None:
        print path
    if opts.importcsv == True:
        createProjects(csvfiles)
    if opts.jobs == True:    
        #PEATSAJobs(['2lzm'], resubmit=True)
        PEATSAJobs(dbnames, resubmit=True)
    if opts.summary == True:
        summarise(dbnames)
        #summarise(['2lzm'])
    if opts.copy == True:
        send2Server(dbnames)
  
