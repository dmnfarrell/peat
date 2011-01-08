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
# Author: Damien Farrell 2009

import pickle, sys, os, copy, time, types
import numpy
from PEATDB.Base import PDatabase 
from PEATDB.Ekin.Titration import TitrationAnalyser
from PEATDB.Ekin.Base import EkinProject, EkinDataset

path=os.environ['HOME']

t = TitrationAnalyser()
H = '1H NMR'
N = '15N NMR'
C = '13C NMR'

'''yuncerts = {H: 0.03, fields[1]:0.2, fields[2]: 0.1} #from lawrence
minspans = {fields[0]: 0.06, fields[1]:0.2, fields[2]: 0.2} '''

complete= ['HEWL', 'Bovine Beta-Lactoglobulin', 'Plastocyanin (Anabaena variabilis)',
           'Plastocyanin (Phormidium)', 'CexCD (Apo)', 'Protein G B1',
           'Glutaredoxin','Staphylococcal Nuclease D+PHS']

def loadDB():
    from PEATDB.Base import PDatabase
    DB = PDatabase(server='peat.ucd.ie', username='guest',
                             password='123', project='titration_db',port=8080)
    return DB

def refitAll():
    #refit whole db and and/or re-get all exp errors, then save the db
    for i in ekindicts:
        e=ekindicts[i]
        #t.fitAll(e)
        t.getExpErrors(e, xuncert=0.1, yuncert=yuncerts[i])
        t.returnData()
   
def doSummary(DB):
    ekindicts = t.getEkinDicts(DB)
    t.dotitDBStats(ekindicts)
    for e in ekindicts:
        #p, img1, img2 = t.analysepKas(ekindicts[e], silent=True)
        p = t.extractpKas(ekindicts[e], titratable=False)
        t.makepKasTable(p, outfile='pkastab_'+e+'.html')#, primary=True)

def miscanalysis(DB):
    
    #p=t.extractpKas(DB,'1H NMR',titratable=True)

    #t.compareProteinpKas(p, prot1='Plastocyanin (Phormidium)',  prot2='Plastocyanin (Anabaena variabilis)')
    t.publicationSetting()
    t.compareNuclei(DB, N, H, titratable=False, names=complete)

    #combekindata = t.makeCombined(h, n, names=complete)
    #pickle.dump(combekindata, open('1h15ncomb.obj', 'w'))
    #combekindata = pickle.load(open('1h15ncomb.obj', 'r'))
    #t.compareNuclei(h, combekindata, titratable=True)

    #p = t.analysepKas(h, silent=True)#, satoms=['H','HA']) #exclude=complete)
    #t.makepKasTable(p, outfile='pkas.html')

    #t.correlatewithStructure(h, ptDB)
    return

def transfer2meta():
    """send residue/atom info to meta_data (for web display)"""
    for i in ekindicts:
        e=ekindicts[i]
        t.showMeta(e)
        #t.transfer2meta(e)
        #t.returnData()

    return

def setpmidInfo():
    """Add extra pmid info to pmid fields using bio.Entrez module"""
    DB = PT.DB
    for p in PT.proteins:
        f = DB[p]['PMID_link']
        #print f
        try:
            auth, tit = t.fetchPMIDSummary(f['text'])
            #print 'got info', tit
        except:
            print 'no pmid'
        try:
            f['authors'] = auth
            f['title'] = tit
            print auth, tit
            #print DB[p]['PMID_link']
        except:
            print 'no dict'
   
    return

def refitandAnalyse(refit=True, usepickle=False, savedb=False):
    """Do everything in one go!"""

    models = ['Modified Hill']
    '''models = ['1 pKa 2 Chemical shifts', 'Modified Hill',
                    '2 pKas, 3 Chemical shifts',
                    '3 pKas, 4 Chemical shifts']'''

    for e in ekindicts:
        if usepickle == True:
            filepi = open('ekindict_'+e, 'r')
            ekindicts[e] = pickle.load(filepi)
            filepi.close()
        elif refit == True:
            t.fitAll(ekindicts[e], models, strictchecking=False)
            filepi = open('ekindict_'+e, 'w')
            pickle.dump(ekindicts[e], filepi)
            filepi.close()

        #p = t.extractpKas(ekindicts[e])
        saveout = sys.stdout
        fsock = open('pkastab_'+e+'.html', 'w')
        sys.stdout = fsock

        #p=t.extractpKas(ekindicts[e], silent=True)
        p, img1, img2 = t.analysepKas(ekindicts[e], silent=True, prefix=e)#,  satoms=['H','HB*'])
        t.makepKasTable(p, primary=True)
        #t.getExpErrors(e, xuncert=0.1, yuncert=yuncerts[i])
        #t.returnData()
        sys.stdout = saveout
    #analyseHill(ekindicts)

    #saveout = sys.stdout
    #fsock = open('fit_stats.html', 'w')
    #sys.stdout = fsock
    #t.dotitDBStats(ekindicts)
    #t.compareNuclei(ekindicts['15N NMR'], ekindicts['1H NMR'])
    #sys.stdout = saveout

    return

def analyseHill(ekindicts):
    """make hist of n coefficents for hill fits"""

    import pylab
    pylab.rc('text', usetex=True)
    f=pylab.figure()
    f.suptitle('n distributions- No linear (case 3)')
    i=1
    for e in ekindicts:
        ekindata = ekindicts[e]
        proteins = ekindata.keys()
        nvals = []
        for prot in proteins:
            edata = ekindata[prot]
            E = EkinProject(data=edata)
            for d in E.datasets:
                fdata = E.getMetaData(d)
                if fdata != None and fdata.has_key('model'):
                    if fdata['model'] == 'Modified Hill':
                        n=fdata['n']
                        if n<5 and n>-5:
                            nvals.append(n)
                            print 'n=', n

        ax = f.add_subplot(2,2,i)
        n, b, patches = pylab.hist(nvals, 30, histtype='bar', alpha=0.8)
        std = round(numpy.std(nvals), 2)
        ave = round(numpy.mean(nvals), 2)
        ax.set_title(e +' mean= '+str(ave)+r' $\sigma$= '+str(std))
        i+=1
    f.subplots_adjust(hspace=0.4)
    f.savefig('n_hist.png')
    return

def dofigure5():
    """make fig for paper, cheapskate Jens"""
    required = {'Plastocyanin (Phormidium)':[61,92],
                'Plastocyanin (Anabaena variabilis)':[62,92],
                'Bovine Beta-Lactoglobulin':[64,114]}
    
    for p in DB.getRecs():
        prot = DB[p]['name']
        if prot in required:
            Eh = data=DB[p]['1H NMR']
            En = data=DB[p]['15N NMR']
            for d in Eh.datasets:
                hdata = Eh.getDataset(d)
                num = int(t.getResidueFields(hdata)['res_num'])
                #print d, type(num), required[prot]
                if num in required[prot]:
                    print num
                    ndata = En.getDataset(d)
                    Eh.plotDatasets(d, filename=d+'.png')
    return

if __name__ == '__main__':
    
    DB = loadDB()
    #doSummary(DB)
    miscanalysis(DB)
    
    #transfer2meta()   
    #refitAll()
    #t.getExpErrors(c, xuncert=0.1, yuncert=0.2)
                    #names=['Xylanase (Bacillus agaradhaerens)'])
    #t.fitAll(h)
    #t.returnData()
    #t.showMeta(c)
    #setpmidInfo()
    #refitandAnalyse(refit=False, usepickle=False)
    #dofigure5()

