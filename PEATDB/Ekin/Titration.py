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

"""A class for doing some analyses on Ekin NMR titration data"""


import sys,os,math
import pickle
import Pmw
import numpy
import string
try:
    import matplotlib
    #matplotlib.use('TkAgg')
    from matplotlib.font_manager import FontProperties
    import pylab
except:
    pass

from NMR import NMR_data
import Utils
from PEATDB.Ekin.Base import EkinProject, EkinDataset
from PEATDB.Ekin.Fitting import Fitting
from PEATDB.Ekin.Pylab import Options
try:
    from Bio import Entrez
except:
    pass

class TitrationAnalyser():
    """Will be used for doing stats on fits etc """

    residue_list = ['ALA', 'ARG', 'ASN', 'ASP',
                            'CYS', 'GLU', 'GLN', 'GLY',
                            'HIS','ILE', 'LEU', 'LYS',
                            'MET', 'PHE', 'PRO', 'SER',
                            'THR', 'TRP', 'TYR', 'VAL']

    residue_letters = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN',
                        'D': 'ASP', 'C': 'CYS', 'E': 'GLU',
                        'Q': 'GLN', 'G': 'GLY', 'H': 'HIS',
                        'I': 'ILE', 'L': 'LEU', 'K':'LYS',
                        'M':'MET', 'F':'PHE', 'P':'PRO',
                        'S': 'SER', 'T': 'THR', 'W': 'TRP',
                        'Y':'TYR', 'V': 'VAL' }

    residue_names = {   'ALA': 'Alanine',
                        'ARG': 'Arginine',
                        'ASN': 'Asparagine',
                        'ASP': 'Aspartate',
                        'CYS': 'Cysteine',
                        'GLU': 'Glutamate',
                        'GLN': 'Glutamine',
                        'GLY': 'Glycine',
                        'HIS': 'Histidine',
                        'ILE': 'Isoleucine',
                        'LEU': 'Leucine',
                        'LYS': 'Lysine',
                        'MET': 'Methionine',
                        'PHE' : 'Phenylalanine',
                        'PRO' : 'Proline',
                        'SER' : 'Serine',
                        'THR' : 'Threonine',
                        'TRP' : 'Tryptophan',
                        'TYR' : 'Tyrosine',
                        'VAL' : 'Valine' }

    atom_types = ['C','CA','CB','CG','CD','CD1','CD2','CE','CE1','CE2','CE3',
                   'H','HA','HB','HB1','HB2','HB3','HD','HD1','HD2','HD21','HD22',
                   'HE','HE1','HE2','HG','HG2','HG3','HH11','HH21','HH12','HH22',
                   'N','NE','NE1','NE2','ND1','ND2','NZ']

    atom_greek = { 'H':'H$^N$', 'HA':r'H$\alpha$',
                  'HB1':r'H$\beta^1$', 'HB2':r'H$\beta^2$', 'HB3':r'H$\beta^3$',
                  'HB*':r'H$\beta^*$',
                  'HE1':'H$\epsilon^1$', 'HE2':'H$\epsilon^2$',
                  'HD1':'H$\delta^1$', 'HD2':'H$\delta^2$',
                  'HG1':'H$\gamma^1$', 'HG2':'H$\gamma^2$', 'HG3':'H$\gamma^3$',
                  'HG*':'H$\gamma^*$',
                  'CA':r'C$\alpha$', 'CB':r'C$\beta$',
                  'CG':'C$\gamma$', 'CD':'C$\delta$', 'CE':'C$\epsilon$',
                  'CD1':'C$\delta^1$', 'CD2':'C$\delta^2$',
                  'CE1':'C$\epsilon^1$', 'CE2':'C$\epsilon^2$',
                  'N':'N$^H$', 'ND1':'N$\delta^1$',
                  'NE':'N$\epsilon$', 'NE2':'N$\epsilon^2$' }

    atom_clrs = {'H':'#0000A0','HB*':'#FF0000','HG*':'#437C17'}

    titratable = ['ASP','GLU','HIS','LYS']#,'CYS','TYR']

    modelpkas = {'N-term': 8.0, 'C-term': 3.7,
                    'ASP': 3.9, 'GLU': 4.3,
                    'CYS': 8.6, 'TYR': 9.8,
                    'SER': 14.2, 'THR': 15.0,
                    'ARG': 13, 'LYS': 10.4,
                    'HIS': 6.5}

    models = ['Linear', '1 pKa 2 Chemical shifts',
                    '2 pKas, 3 Chemical shifts',
                    '3 pKas, 4 Chemical shifts']

    shapes = Options.shapes[:11]
    pylabcolors = Options.colors
    grays = ['1.0','0.8','0.5','0.1','0.6']

    def __init__(self, data=None, parent=None):
        if parent!=None:
            self.parent = parent
        self.currfile='/home/people/farrell/PEAT_projects/PROT_TITRA.PEAT'
        self.logfile='/tmp/nmr_titr_out.log'

        self.allmodels=['Linear', '1 pKa 2 Chemical shifts',
                        '3 pKas, 2 Chemical shifts',
                        '2 pKas, 3 Chemical shifts',
                        '3 pKas, 4 Chemical shifts',
                        '4 pKas, 5 Chemical shifts']

        #takes data in the form of dict for each protein, as stored in peat
        #so we have __datatabs_fits__ fields etc. use excluded to ignore them
        self.excluded=['__datatabs_fits__','__Exp_Meta_Dat__','__fit_matches__',
                  '__FITTER_fit__','__FITTER_model__','__distance_matches__',
                  '__datatab_structmapping__','__meta_data__']
        return

    def publicationSetting(self):
        pylab.rc("font", family="serif")
        pylab.rc("font", size=22)   
        pylab.rc('text', usetex=True)
        return
   
    def getEkinDicts(self, DB):
        '''reform peat db into sets of ekin data stored per nucleus type'''

        protnames={}
        ekindicts = {}
        for protein in DB.getRecs():
            name = DB[protein].name
            protnames[protein] = name
            for col in DB[protein].keys():
                if DB['userfields'].has_key(col):
                    if 'NMR' in DB['userfields'][col]['field_type']:
                        if not ekindicts.has_key(col):
                            ekindicts[col] = {}
                        ekindicts[col][protein] = DB[protein][col]
        self.protnames=protnames
        self.ekindicts = ekindicts
        return ekindicts
 
    def do_summary(self, DB, cols=None):
        """Print a summary of the current data loaded - no longer used"""
        summary=[]

        self.model_count={}
        self.residue_count={}
        self.residue_models={}
        self.atom_count={}
        for m in self.allmodels:
            self.model_count[m]=0
        #store per residue model data
        for res in self.residue_list:
            self.residue_count[res]=0
            self.residue_models[res]={}
            for m in self.allmodels:
                self.residue_models[res][m]=0
        for m in self.allmodels:
            self.atom_count[m]={}
            for a in self.atom_types:
                self.atom_count[m][a]=0
        for prot in ekindata:
            if ekindata[prot].has_key('__datatabs_fits__'):
                fits = ekindata[prot]['__datatabs_fits__']
            for dataset in ekindata[prot].keys():
                if dataset in self.excluded or ekindata[prot][dataset]==None:
                    continue
                if ekindata[prot][dataset].has_key('residue'):
                    resname = ekindata[prot][dataset]['residue']
                    self.residue_count[resname] += 1
                else:
                    resname = None
                if ekindata[prot][dataset].has_key('atom_type'):
                    atomtype = ekindata[prot][dataset]['atom_type']
                else:
                    atomtype = None
                if fits.has_key(dataset):
                    fdata = fits[dataset]

                    #print dataset, fdata
                    if fdata != None and len(fdata) != 0:
                        model = fdata['model']
                        self.model_count[model] += 1
                        if resname != None:
                            self.residue_models[resname][model] += 1
                        if atomtype != None:
                            self.atom_count[model][atomtype] += 1

        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        x=pp.pformat(self.model_count)
        summary.append('Models:')
        summary.append(x+'\n')
        x=pp.pformat(self.residue_count)
        summary.append('Residues:')
        summary.append(x+'\n')
        x=pp.pformat(self.residue_models)
        summary.append('Residue model breakdown:')
        summary.append(x+'\n')
        x=pp.pformat(self.atom_count)
        summary.append('Atom types breakdown:')
        summary.append(x+'\n')
        return summary

    def getProtNames(self, DB):
        self.protnames={}       
        for protein in DB.getRecs():
            name = DB[protein]['name']
            self.protnames[protein] = name
        return self.protnames              

    @classmethod
    def getFitStats(cls, ekindata):
        """Simple fit stats"""
        stats={}
        #for m in cls.models: stats[m]=0
        proteins = ekindata.keys()
        for prot in proteins:
            E = ekindata[prot]
            #E = EkinProject(data=edata)
            for d in E.datasets:
                fdata = E.getFitData(d)
                if fdata!=None and len(fdata)>1:
                    m = fdata['model']
                    if not m in stats:
                        stats[m]=0
                    stats[m]+=1

        return stats

    def dotitDBStats(self, ekindata, proteins=None, redirect=False):
        """Stats for titration datasets - takes a dict of ekindata, with the nmr
           column name as the key"""

        #stats on fits first
        stats={}

        for i in ekindata:
            s = self.getFitStats(ekindata[i])
            stats[i] = s
            modls=s.keys()

        print '<div>'
        print '<table id="mytable" width=60%>'
        print '<th></th>'
        for m in modls:
            print '<th>%s</th>' %m
        for x in stats:
            print '<tr><td>%s</td>' %x
            for m in modls:
                print '<td>%s</td>' %stats[x][m]
            #print '<td>%s</td>' %total
            print '</tr>'
        print '</table>'
        #plot bar chart of stats

        print '</div>'


        totaldatasets=0
        totalphpoints=[]
        minphpoints=[];maxphpoints=[]
        phranges=[]
        stdevphrange=0
        residuestats={}
        proteinstats={}
        lessXphstep=0

        totaltitres=0
        for t in self.residue_list: residuestats[t]=0
        names = {}
        for i in ekindata:
            names[i] = {}
            proteinstats[i] = {}

        for i in ekindata:
            ekinprjs = ekindata[i]
            avgchemshft=0
            names[i]['total'] = 0
            for prot in ekinprjs:
                proteinstats[prot]={}
                proteinstats[prot][i]={}
                E = ekinprjs[prot]                
                totaldatasets += len(E.datasets)
                names[i]['total'] += len(E.datasets)
                proteinstats[prot][i]['total'] = len(E.datasets)
                for ed in E.datasets:
                    try:
		                ek = EkinDataset(E.data[ed])
		                totalphpoints.append(ek.length())
		                maxphpoints.append(ek.maxX())
		                minphpoints.append(ek.minX())
		                phranges.append(ek.maxX()-ek.minX())
		                phstep = (ek.maxX()-ek.minX()) / ek.length()
		                if phstep <= 0.3:
		                    lessXphstep+=1
		                ek.avgY()
		                fit = E.getFitData(ed)
		                res = E.getField(ed, 'residue')
		                if res != None:
		                    residuestats[res] += 1
                    except:
                        pass


        print '<div id="left">'
        print '<h3>Global statistics for current data:</h3>'
        print '<a>'
        print 'total no of datasets: %s <br>' %totaldatasets
        print 'average/stdev. no. pH points: %.2f, %.2f <br>' %(numpy.mean(totalphpoints), numpy.std(totalphpoints))
        print 'average lowest/highest pH point: %.2f, %.2f, std. devs.(%s, %s) <br>' \
                %(numpy.mean(minphpoints), numpy.mean(maxphpoints), numpy.std(minphpoints), numpy.std(maxphpoints))
        print 'average/stdev. pH ranges: %.2f, %.2f <br>' %(numpy.mean(phranges), numpy.std(phranges))
        print 'total no. of pH points: %s <br>' %numpy.sum(totalphpoints)
        for t in self.titratable:
            totaltitres += residuestats[t]
        print 'total no. of titratable residues: %s (%.2f) <br>' %(totaltitres, float(totaltitres)/totaldatasets*100)
        print 'Per nucleus stats: <br>'
        for i in names:
            perc = round(float(names[i]['total'])/totaldatasets*100,2)
            print '%s: %s (%s) <br>' %(i,names[i]['total'], perc)

        perlessXcurves=0
        total=0
        for p in proteinstats:
            for i in names:
                if proteinstats[p].has_key(i):
                    n=proteinstats[p][i]['total']
                    total+=1.0
                    if n<=5:
                        perlessXcurves+=1

        print 'percent with less than 5 curves: %.2f <br>' %(perlessXcurves/total)
        print 'percent with lower than 0.2 pH step: %.2f <br>' %(float(lessXphstep)/totaldatasets)
        print '</a>'
        print '</div>'

        print '<div id="right">'
        print '<table id="summary" width=60% align=center cellspacing=0>'
        print '<h3>Breakdown by Residue:</h3>'

        print '<th>%s</th> <th>%s</th> <th>(%s)</th>' %('res','no.','% total')
        for r in residuestats:
            print '<tr>'
            if Utils.getDictSum(residuestats) == 0:
                perc=0.0
            else:
                perc = round(float(residuestats[r])/Utils.getDictSum(residuestats)*100,2)
            print '<td>%s</td> <td>%s</td> <td>(%s)</td>' %(r,residuestats[r], perc)
            print '</tr>'

        print '</td></tr>'
        print '</table>'
        print '</div>'
        print '</div>'
        return


    def fitAll(self, DB, col, proteins=None, models=None, strictchecking=False):
        """Fit selected recs/col in DB using find best model"""
        stats={}
        if strictchecking == True:
            checkfunc=self.dostrictchecking
        else:
            checkfunc=None

        for m in models: stats[m]=0
        if models == None:
            models = self.models
        if proteins == None:
            proteins = DB.getRecs()
        self.getProtNames(DB)
        
        c=0; f=0
        for prot in proteins:
            name = self.protnames[prot]
            print 'fitting', name
            E = DB[prot][col] 
            for d in E.datasets:
                print name, d
                if len(models)==1:
                    #we don't do f-test if only one model
                    fdata = E.fitDataset(d, model=models[0])
                else:
                    fdata, p = E.findBestModel(d, models=models, checkfunc=checkfunc)
                if fdata==None:
                    f+=1
                c+=1
                if fdata!=None:
                    stats[fdata['model']]+=1           
            DB[prot][col] = E
        print 'done. fitted %s datasets with %s failed' %(c, f)
        return stats

    def findBest(self, E, models, geterrs=False, xuncert=0.1, yuncert=0.03):
        """Find best model for all datasets in an ekin project"""
        for d in E.datasets:
            fdata, p = E.findBestModel(d, models=models)#, checkfunc=self.dostrictchecking)
            if geterrs==True:
                ferrs = E.estimateExpUncertainty(d, runs=20, xuncert=xuncert, yuncert=yuncert)
                if ferrs == None:
                    continue                
                E.addMeta(d, 'exp_errors', ferrs)
                    
        return E
        
    def getExpErrs(self, E, xuncert=0.1, yuncert=0.03, runs=20):
        """Get exp uncertainties on current fits"""  
        for d in E.datasets:
            print d
            ferrs = E.estimateExpUncertainty(d, runs=runs, 
                                    xuncert=xuncert, yuncert=yuncert)
            if ferrs == None:
                continue                
            E.addMeta(d, 'exp_errors', ferrs)        
        return E
        
    def showMeta(cls, ekindata):
        """Print out meta data for all recs - debug"""
        proteins = ekindata.keys()
        for prot in proteins:
            print prot
            edata = ekindata[prot]
            E = EkinProject(data=edata)
            E.printMeta()

        return

    @classmethod
    def dostrictchecking(cls, datapoints, fitdata):
        """Do some checking of pKa models as per McIntosh suggestions
           1. check pKa values are not too close to the end of the data points"""
        
        ph = datapoints
        phmin = min(ph)
        phmax = max(ph)        
        model = fitdata['model']
        vrs = Fitting.getFitVars(fitdata)
        X = Fitting.getFitter(model, vrs=vrs)
        fit = X.getFitDict()

        print fit
        pkas=[]
        for p in fit:
            if 'pK' in p:
                pka=float(fit[p])
                pkas.append(pka)
                diff1=phmax-pka
                diff2=pka-phmin
                print 'pka', pka, 'is', diff1, 'away from endpoint'
                print 'and', diff2, 'away from startpoint'
                if pka < phmin or pka > phmax:
                    print 'this pKa is off the pH scale!'
                    return 0
                if diff1 < 0.1 or diff2 < 0.1:
                    print 'this pKa is too close to an end point'
                    return 0
        return 1

    def getExpErrors(cls, ekindata, xuncert=0.1, yuncert=0.01, names=None, residues=None,
                        guess=True):
        """Get exp error on all pKa fits and store in meta_data
            - may take some time... """

        print 'getting errs from exp uncertainty..'
        print 'using x uncert: %s and yuncert: %s' %(xuncert, yuncert)
        proteins = ekindata.keys()
        for prot in proteins:
            name = cls.protnames[prot]
            print name
            print '----------------------------------'
            if names != None and name not in names:
                continue
            edata = ekindata[prot]
            E = EkinProject(data=edata)
            for d in E.datasets:
                if residues != None:
                    ed = E.getDataset(d)
                    #resdata = cls.getResidueFields(ed)
                    resdata = E.getMetadata(d)
                    resnum = int(resdata['res_num'])
                    if resnum not in residues:
                        continue
                print d
                fit = E.getFitData(d)
                if fit != None and len(fit)>0 and fit['model'] != 'Linear':
                    ferrs = E.estimateExpUncertainty(d, runs=20, xuncert=xuncert, yuncert=yuncert, guess=guess)
                    if ferrs == None:
                        continue
                    #print d, ferrs
                    E.addMeta(d, 'exp_errors', ferrs)
            #we have added meta_data, so overwrite dict with ekinproject
            ekindata[prot] = E.prepare_data()

        return E

    def findpKas(cls, E, titratable=True, reliable=True, minspan=0.06):
        """Get pkas for an ekin project"""
        pkasdict = {}
        for d in E.datasets:   
            resdata = E.getMetaData(d)            
            if not resdata.has_key('residue'): continue
            res = resdata['residue']
            try:
                atom = resdata['atom']
            except:
                atom = resdata['atom_type']
            resnum = resdata['res_num']          
            fitdata = E.getFitData(d)
            
            if titratable == True and not res in cls.titratable:
                continue
            if fitdata == None or len(fitdata)<2:
                continue
            if reliable == True:
                p,pval,sp = cls.getMainpKa(fitdata, res=res, minspan=minspan, maxspan=5)
                if pval == None:
                    continue
                pkas = [(p,pval,sp)]
            else:
                pnames, pvals, spans = cls.getallpKas(fitdata, minspan=minspan)
                if pnames == None:
                    continue                    
                pkas = zip(pnames, pvals, spans)

            if len(pkas)>0 and pkas!=None:
                pkasdict[d]={}
                pkasdict[d]['res']=res
                pkasdict[d]['resnum']=resnum
                pkasdict[d]['atom']=atom               
                pkasdict[d]['model']=fitdata['model']
                for plst in pkas:
                    (p, pval, sp) = plst
                    s = 'span' + p.strip('pKa')
                    if E.getMeta(d, 'exp_errors')==None:
                        #normally we should have stored this already
                        ferrs = E.estimateExpUncertainty(d, runs=10, xuncert=0.2, yuncert=0.1)
                        E.addMeta(d, 'exp_errors', ferrs)
                    try:
                        pkaerr = round(E.getMeta(d, 'exp_errors')[p][1],4)
                        spanerr = round(E.getMeta(d, 'exp_errors')[s][1],4)
                    except:
                        #print 'failed to get exp errs'
                        pkaerr = 0.0
                        spanerr = 0.0
                    #print  d, plst#, pkaerr, spanerr
                    if reliable==True and (pkaerr > 1.5 or spanerr > sp):                        
                        #print d, 'exp error too big, omitting: %s %s' %(pkaerr, spanerr)
                        continue
                    pkasdict[d][p]={}
                    pkasdict[d][p][p]=pval
                    pkasdict[d][p]['span']=sp                       
                    pkasdict[d][p]['pkaerror']=pkaerr
                    pkasdict[d][p]['spanerror']=spanerr

        return pkasdict        
                        
    def extractpKas(cls, DB, col, names=None, minspan=0.06, reliable=True,
                        titratable=True, silent=False):
        """Extract all or just reliable pKas from a set of ekin projects in the
           DB.
           Returns a dictionary with record/dataset keys storing pka info"""

        cls.getProtNames(DB)
        if silent == False:
            print 'extracting pka values..'
        pkainfo={}
        total=0
        for prot in DB.getRecs():           
            name = cls.protnames[prot]
            if names != None and name not in names:
                continue
            if not DB[prot].has_key(col):
                continue
            if silent == False:
                print 'processing %s' %name
                print '------------------------------------'            
            E = DB[prot][col]
            pkainfo[name] = cls.findpKas(E, titratable=titratable,
                                         reliable=reliable, minspan=minspan)
            
            if silent == False: print 'found %s' %len(pkainfo[name])
            total+=len(pkainfo[name])
        print 'found %s total pKa values' %total
        return pkainfo

    def analysepKas(cls, pkainfo, satoms='all', path=None, exclude=[],
                        silent=False, prefix=''):
        """Read in the pkainfo dict from extractpkas method and do 
          additional analyses, makes plots for titdb"""
    
        if path==None:
            path=os.getcwd()
        aminoacids = cls.residue_names.keys()
        titratable = cls.titratable        
       
        titrspans = {}      #histogram of reliable spans, for titr groups
        primaryspans = {}   #primary spans, for titr groups
        allspans = {}       #histogram of all spans for all groups
        titravgshfts = {}   #histogram of average shifts
        atomshfts = {}      #store per atom
        spansbyres = {}     #sort out spans by residue

        for t in aminoacids:
            titravgshfts[t]=[]
        for t in titratable:
            titrspans[t]=[]
            primaryspans[t]=[]
            spansbyres[t]={}
            atomshfts[t]={}
        for t in aminoacids:
            allspans[t] = []
            
        if silent == False:
            print
            print 'Analysing pKa info..'

        count=0
        titrcount = 0       
        for name in pkainfo.keys():            
            for d in pkainfo[name].keys():
                PI = pkainfo[name][d]
                res = PI['res'] 
                atom = PI['atom']
                model = PI['model']  
                for p in PI:
                    #print p, d, '<p>'                    
                    if 'pK' in p:
                        count += 1
                        pval = PI[p][p]
                        sp = PI[p]['span']                    	                      
                   
                        #print  d, p, pval,  sp
                        #change atom if handling HB2/3 in ASP, HIS and GLU, so
                        #we group them together as HB* for the stats
                        
                        #move to extract??
                        if res in ['ASP','GLU','HIS'] and atom in ['HB2', 'HB3']:
                            atom = 'HB*'
                        elif res == 'GLU' and atom in ['HG2','HG3']:
                            atom = 'HG*'
                        
                        if model == '1 pKa 2 Chemical shifts' and res in titratable:
                            primaryspans[res].append(sp)
                        elif res in titratable:
                            titrspans[res].append(sp)
                            titrcount+=1
                        else:
                            allspans[res].append(sp)
                       
                        #assign spans by atom
                        if not atomshfts[res].has_key(atom):
                            atomshfts[res][atom] = []
                        if satoms != 'all' and not atom in satoms:
                            if silent == False:
                                print 'excluding atom', atom
                                continue
                        else:
                            atomshfts[res][atom].append(sp)

            '''if allpkas != None:
                for s in allsp:
                    if s >= maxspan:
                        continue
                    allspans[res].append(s)
                    if res in titratable:
                        spansbyres[res][s] = prot+'__'+d'''

        #lenallspans=Utils.getlenDict(allspans)

        lentitrspans=Utils.getlenDict(titrspans)
        lenprimaryspans=Utils.getlenDict(primaryspans)
        #lentitrallspans=Utils.getlenDict(allspans, keys=titratable)

        #remove primary from allspans to get leftover ones
        otherspans = {}
        for k in titrspans:
            otherspans[k] = list( set(allspans[k]) - set(titrspans[k]) )
            #otherspans[k] = list( set(allspans[k]) - (set(titrspans[k]) | set(primaryspans[k])))

        print '<div>'
        print 'Total residues: %s <br>' %count
        print 'of which titratable grps: %s <br>' %titrcount
        #print 'All pkas: %s <br>' %lenallspans
        print 'Titratable "reliable" pkas: %s <br>' %(lentitrspans + lenprimaryspans)
        print 'Titratable "primary" pkas: %s <br>' %lenprimaryspans
        #print 'Titratable all pkas: %s <br>' %lentitrallspans
        print '</div>'
        print '<br>'

        for r in atomshfts:
            print r
            for a in atomshfts[r]:
                print a, len(atomshfts[r][a])

        pylab.rc("font", family="serif")
        pylab.rc("font", size=10)
        pylab.rc('text', usetex=True)

        fig1 = pylab.figure(1,figsize=(10,6))
        fig1.text(0.5, 0.95, 'Distribution: $\Delta\delta$ of extracted pKas: '+prefix, horizontalalignment='center',
                            fontproperties=FontProperties(size=18))

        #plot multiple histograms together in multiple subplots
        hdata = {}
        for k in titrspans:
            hdata[k] = {'primary':primaryspans[k], 'reliable':titrspans[k], 'other':otherspans[k]}

        cls.doMultipleHistograms(fig1, hdata, bins=30, xlabel='$\Delta\delta$ (ppm)',
                                ylabel='no. curves', xticks=False, yticks=False)
        img1 = prefix+'deltashifts_dist.png'
        fig1.savefig(os.path.join(path,img1), dpi=100)
        pylab.clf()
        fig6 = pylab.figure(6,figsize=(10,6))
        fig6.text(0.5, 0.95, 'Distribution: $\Delta\delta$ by atom: '+prefix, horizontalalignment='center',
                            fontproperties=FontProperties(size=18))

        cls.doMultipleHistograms(fig6, atomshfts, bins=25, xlabel='$\Delta\delta$ (ppm)',
                                xticks=True, yticks=True)
        img2 = prefix+'deltashifts_byatom.png'
        fig6.savefig(os.path.join(path,img2) ,dpi=100)

        '''fig7 = pylab.figure(7,figsize=(10,6))
        cls.do3DHistograms(fig7, atomshfts, bins=25, xlabel='$\Delta\delta$ (ppm)',
                                xticks=True, yticks=True)
        img3 = prefix+'deltashifts_3d.png'
        fig7.savefig(os.path.join(path,img3) ,dpi=100)'''

        #pylab.show()
        pylab.clf()

        return img1, img2

    def makepKasTable(cls, pkainfo=None, outfile=None, primary=False):
        """make table of pka fit/analysis results"""
        if outfile != None:
            saveout = sys.stdout
            fsock = open(outfile, 'w')
            sys.stdout = fsock

        kys = sorted(pkainfo.keys())
        resinfo = {}
        #for r in cls.titratable:
        for r in cls.residue_names.keys():    
            resinfo[r] = {}

        print '<table id="summary" width=80% align=center cellspacing="0">'
        row=1;c=1
        print '<th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th> <th>%s</th>'\
                    % ('protein', 'dataset', 'res type', 'resnum', 'atom', '&Delta&delta', '&Delta&delta err', 'pKa', 'pKa err', 'model')
        print '<tr>'

        for name in kys:
            for d in pkainfo[name]:
                PI = pkainfo[name][d]
                resnum = PI['resnum']
                res = PI['res']
                atom = PI['atom']
                model = PI['model']
                for p in PI:
                    if not 'pK' in p: continue
                    pka = PI[p][p]
                    span = PI[p]['span']
      
                    if primary == True and not '1 pKa' in model:
                        continue
                    if pkainfo[name][d][p].has_key('pkaerror'):
                        pkaerr = pkainfo[name][d][p]['pkaerror']
                        spanerr = pkainfo[name][d][p]['spanerror']
                    else:
                        pkaerr = 0
                        spanerr=0
                                              
                    if not resinfo[res].has_key(atom):
                        resinfo[res][atom]={}
                        resinfo[res][atom]['pkas']=[]
                        resinfo[res][atom]['pkaerrs']=[]
                        resinfo[res][atom]['spans']=[]
                        resinfo[res][atom]['spanerrs']=[]
                    resinfo[res][atom]['pkas'].append(pka)
                    resinfo[res][atom]['pkaerrs'].append(pkaerr)
                    resinfo[res][atom]['spans'].append(span)
                    resinfo[res][atom]['spanerrs'].append(spanerr)
                    #print name, res, resnum, atom, span, spanerr, pka, pkaerr, model
                    if pka != None:
                        print '<td>%s</td> <td>%s</td> <td>%s</td> <td>%s</td> <td>%s</td> <td>%.2f</td> <td>%.2f</td> <td>%.2f</td> <td>%.2f</td> <td>%s</td>'  \
                                    % (name, d, res, resnum, atom, span, spanerr, pka, pkaerr, model)
                        print '</tr>'
        print '</table>'
        print '<br>'
        #summary table for average pKas per residue
        print '<h2>Average pKa and &Delta&delta per atom/residue</h2><br\Delta'
        print '<table id="summary" width=60% align=center cellspacing="0">'
        print '<tr><td>Residue</td> <td>atom</td> <td>average &Delta&delta</td> <td>average &Delta&delta err</td> <td>average pKa</td> <td>av pKa error</td> <td>no. curves</td></tr>'
        for r in resinfo:
            allpkas=[]
            allspans=[]
            allpkaerrs=[]
            for a in resinfo[r]:
                numcurves = len(resinfo[r][a]['pkas'])
                averagepka = numpy.mean(resinfo[r][a]['pkas'])
                averagepkaerr = numpy.mean(resinfo[r][a]['pkaerrs'])
                averagespan = numpy.mean(resinfo[r][a]['spans'])
                averagespanerr = numpy.mean(resinfo[r][a]['spanerrs'])
                minpka = min(resinfo[r][a]['pkas'])
                maxpka = min(resinfo[r][a]['pkas'])
                print '<tr><td>%s</td> <td>%s</td> <td>%.2f</td> <td>%.2f</td> <td>%.1f</td> <td>%.2f</td> <td>%s</td></tr>' \
                    %(r, a, averagespan, averagespanerr, averagepka, averagepkaerr, numcurves)
                allpkas += resinfo[r][a]['pkas']
                allpkaerrs += resinfo[r][a]['pkaerrs']
                allspans += resinfo[r][a]['spans']

            finalpkaerr = numpy.sqrt(sum(numpy.square(allpkaerrs)))
            print '<tr><td>%s</td> <td>%s</td> <td>%.2f</td> <td>%.2f</td> <td>%.1f</td> <td>%.1f</td> <td>%s</td></tr>' \
                %(r, 'All', numpy.average(allspans), numpy.std(allspans), numpy.average(allpkas), numpy.std(allpkas), len(allpkas))

        print '</table>'
        print '<br>'
        if outfile != None:
            sys.stdout.flush()
            sys.stdout = saveout
            fsock.close()
        return

    def combineCurves(cls, d1, d2, factor=5):
        """Combine 2 sets of chem shifts to get one pH titr"""     
        ph1,ch1 = d1.getxy()
        ph2,ch2 = d2.getxy()
        
        def getnearest(x, v):
            #get index of nearest val
            n=100
            c=0
            for val in x:
                i = abs(val-v)
                if v < n:
                    n = v
                    i = c
                c=c+1
                return i

        print ph1
        print ph2
        #get ref vals for both , ie. min val
        ref1 = min(ch1)
        ref2 = min(ch2)
        comb = []
        if ph1==ph2:
            for v in range(len(ph1)):
                comb.append(sqrt(pow(ch1[v]-ref1,2)+pow((ch2[v]-ref2)/factor,2)))
        else:
            print 'not equal'
            #combine nearest ph vals instead
            for v in range(len(ph1)):
                #get nearest in ph2 to v1
                v2 = getnearest(ph2, v)
                comb.append(sqrt(pow(ch1[v]-ref1,2)+pow((ch2[v2]-ref2)/factor,2)))
        
        dc = EkinDataset(xy=[ph1, comb])
        return dc

    def makeCombined(cls, ekindata1, ekindata2, names=None):
        """Make combined chem shft and fit,then compare to originals"""

        proteins = ekindata1.keys()
        ekindatac = {}
        E1 = EkinProject(ekindata1)
        E2 = EkinProject(ekindata2)
        for prot in proteins:
            print
            if cls.protnames != None:
                name = cls.protnames[prot]
                print 'processing protein', name
            print '-------------------------'
            if names != None and name not in names:
                continue
            if not ekindata2.has_key(prot):
                print 'this protein is not in the 2nd ekindata'
                continue
            edata1 = ekindata1[prot]; E1 = EkinProject(data=edata1)
            edata2 = ekindata2[prot]; E2 = EkinProject(data=edata2)
            Ec = EkinProject()
            for d in E1.datasets:
                if d in cls.excluded or not d in E2.datasets:
                    continue
                d1 = E1.getDataset(d)
                d2 = E2.getDataset(d)
                print name
                dc = cls.combineCurves(d1, d2)
                Ec.insertDataset(dc, d)
            Ec.saveProject(name[:5]+'_combined')
            for d in Ec.datasets:
                f, p = Ec.findBestModel(d, models=cls.models, strictchecking=True, alpha=0.05)
            Ec.saveProject(name[:5]+'_combined')
            ekindatac[prot] = Ec.prepare_data()
        return 

    def getClosest(cls, val, pkas):
        #print val ,pkas
        c=100
        ind=0
        for j in pkas:
            d=abs(val - j)
            if d < c:
                c=d
                ind=pkas.index(j)
        return pkas[ind]

    def comparepKas(cls, E1, E2, d, titratable=False, errcutoff=1e2):
        """Compare pKas from fits in two ekin datasets
           returns:
               None if fails to find meaningful pKas from fits
               otherwise returns tuple of results"""
        reliable=True
        try:    
            res1 = E1.getMetaData(d)['residue']
        except:
            res1 = cls.getResidue(E1.getDataset(d))             
        if titratable == True and not res1 in cls.titratable:
            return None
        fitdata1 = E1.getFitData(d)
        fitdata2 = E2.getFitData(d)
        if len(fitdata1)<1:
            return None
        model1 = fitdata1['model']
        model2 = fitdata2['model']
        if model1 != model2:
            print 'different models'
            #return None
        
        p1,p1val,sp1 = cls.getMainpKa(fitdata1, res=res1)#,strict=False)
        p2,p2val,sp2 = cls.getMainpKa(fitdata2, res=res1)#,strict=False)
        allpkas1,allpkvals1,allsp1 = cls.getallpKas(fitdata1,minspan=0.03)
        allpkas2,allpkvals2,allsp2 = cls.getallpKas(fitdata2,minspan=0.03)

        if p1val != None and p2val != None and len(allpkas1) > 0 and len(allpkas2) > 0:
            pchk = cls.getClosest(p1val, allpkvals2)
            try:
                perr1 = round(E1.getMeta(d, 'exp_errors')[p1][1],4)
                perr2 = round(E2.getMeta(d, 'exp_errors')[p2][1],4)
            except:
                return None
            print d, p1val, p2val, perr1, perr2
            #print errcutoff
            if perr1>errcutoff or perr2>errcutoff: return None
            if p2val != pchk:
                p2val = pchk
        else:            
            if p1val != None and p2val == None:
                if len(allpkas2) > 0:
                    p2val = cls.getClosest(p1val, allpkvals2)                   
                else:                   
                    return None
            elif p2val != None and p1val == None:
                if len(allpkas1) > 0:
                    p1val = cls.getClosest(p2val, allpkvals1)                   
                else:
                    return None
            elif p1val==None and p2val==None:
                return None
            perr1=perr2=None
            reliable = False
            
        return p1val, p2val, sp1, sp2, perr1, perr2, reliable

    def compareAllpKas(cls, E1, E2, titratable=False, exclude=None, errcutoff=1e2):
        """Compare all pKas from 2 ekin projects"""
        relpkas1=[]
        relpkas2=[]
        relspans1=[]
        relspans2=[]        
        otherpkas1=[]
        otherpkas2=[]
        errs1=[]
        errs2=[]
        names=[]
        
        for d in E1.datasets:
            if exclude!=None and d in exclude: continue
            if not d in E2.datasets: continue
            X = cls.comparepKas(E1, E2, d, titratable, errcutoff)            
            if X == None:
                continue
            p1val, p2val, sp1, sp2, err1, err2, rel = X
            if rel == True:
                relpkas1.append(p1val)
                relpkas2.append(p2val)     
                relspans1.append(sp1)
                relspans2.append(sp2)
                errs1.append(err1)
                errs2.append(err2)
                names.append(d)
            else:                
                otherpkas1.append(p1val)
                otherpkas2.append(p2val)        
            
        print 'reliable pkas matched:', len(relpkas1)
        print 'others:', len(otherpkas1)

        return relpkas1, relpkas2, otherpkas1, otherpkas2, relspans1, relspans2, errs1, errs2, names

    def compareNuclei(cls, DB, col1, col2, names=None, titratable=True):
        """Compare corresponding datasets for proteins that have data
          for 2 different NMR nuclei e.g. 1H vs 15N over entire DB"""

        relpkas1=[]
        relpkas2=[]
        relspans1=[]
        relspans2=[]        
        otherpkas1=[]
        otherpkas2=[]
        #names=[]
        pkasbyres1={}
        pkasbyres2={}
        for t in cls.residue_list:
            pkasbyres1[t]=[]
            pkasbyres2[t]=[]
     
        cls.getProtNames(DB)

        for prot in DB.getRecs():            
            name = cls.protnames[prot]            
            if names != None and name not in names:
                continue
            print 'processing protein', name
            print '-------------------------'
            if not DB[prot].has_key(col1) or not DB[prot].has_key(col2):
                continue
            E1 = DB[prot][col1]
            E2 = DB[prot][col2]            
            
            X = cls.compareAllpKas(E1, E2, exclude=cls.excluded)
            print len(X)
            if X == None:
                continue            
            rp1, rp2, op1, op2, rsp1, rsp2, errs1, errs2, n = X
            relpkas1.extend(rp1)
            relpkas2.extend(rp2)
            otherpkas1.extend(op1)
            otherpkas2.extend(op2)
            #names.extend(n)
            
        print 'reliable pkas matched:', len(relpkas1)
        print 'others:', len(otherpkas1)
   
        f=pylab.figure(figsize=(10,20))
        ax1=f.add_subplot(211)
        ax2=f.add_subplot(212)
        cc = cls.doXYPlot(ax1, relpkas1, relpkas2, #names=names,
                        title='15N vs 1H : reliable pKas', xlabel='15N', ylabel='1H')
        print 'reliable pKas, correl coeff:', cc
        cc = cls.doXYPlot(ax2, otherpkas1, otherpkas2, color='r',
                        title='15N vs 1H : other pKas', xlabel='15N', ylabel='1H')
        print 'other pKas, correl coeff:', cc
        f.savefig('comparenuclei.png', dpi=300)
        
        '''f=pylab.figure(figsize=(10,10))
        ax=f.add_subplot(111)
        i=0
        leglines = [];series=[]
        for r in pkasbyres1.keys():
            if len(pkasbyres1[r])==0:
                continue            
            if i >= len(cls.shapes) or i>=len(cls.pylabcolors):
                i=0
            cls.doXYPlot(ax, pkasbyres1[r],pkasbyres2[r], symbol=cls.shapes[i], color=cls.pylabcolors[i],
                            markersize=40,
                            title='1H vs 15N: reliable pKa values', xlabel='15N', ylabel='1H')
            l = pylab.Line2D(range(10), range(10), color=cls.pylabcolors[i], alpha=0.7,
                                marker=cls.shapes[i])
            leglines.append(l)
            series.append(r)
            i+=1
        leg = ax.legend(leglines, series, markerscale=2,
                             numpoints=1,loc='lower right',
                             prop=FontProperties(size="small"))
        leg.draw_frame(False)
        f.savefig('comparenuclei_relbyres.png', dpi=300)'''
        #pylab.show()
        return

    def compareExtractedpKas(cls, DB, col, prot1, prot2):
        """Compare extracted pKas across similar proteins and plot correlations per res"""

        pkas = cls.extractpKas(DB, col)
        p1={};p2={};p1errs={};p2errs={}
        labels={}
        tp=['conserved','non-conserved']
        for t in tp:
            p1[t]=[];p2[t]=[]
            p1errs[t]=[];p2errs[t]=[]
            labels[t]=[]
        s1=[];s2=[]
        pkas1 = pkas[prot1]
        pkas2 = pkas[prot2]

        for n in pkas1.keys():
            if n in pkas1 and n in pkas2:
                res1 = pkas1[n]['res']; res2 = pkas2[n]['res']
                #this part requires alignment info!!
                if res1 == res2:# and res1 in cls.titratable:
                    t='conserved'
                else:
                    t='non-conserved'
                p1[t].append(pkas1[n]['pka']); p2[t].append(pkas2[n]['pka'])
                p1errs[t].append(pkas1[n]['pkaerror']); p2errs[t].append(pkas2[n]['pkaerror'])
                #s1.append(pkas1[n]['span']); s2.append(pkas2[n]['span']);
                labels[t].append(n)

        fig1=pylab.figure()
        i=0; leg=[]
        for t in tp:
            print p1[t], labels[t]
            if len(p1[t])==0 or len(p2[t])==0:
                continue
            cls.doXYPlot(fig1, p1[t], p2[t], xerrs=p1errs[t], yerrs=p2errs[t],
                            symbol=cls.shapes[i], color=cls.pylabcolors[i],
                            title='compared pKas', annotate=labels[t],
                            xlabel=prot1, ylabel=prot2)
            l = pylab.Line2D(range(10), range(10), color=cls.pylabcolors[i], alpha=0.8,
                                marker=cls.shapes[i])
            leg.append(l)
            i+=1
        pylab.legend(leg,tp,numpoints=1)
        a=fig1.get_axes()[0]
        a.set_xlim(4,5.5)
        a.set_ylim(4,5.5)
        fig1.savefig('compareprotpkas.png', dpi=200)

        '''fig2=pylab.figure()
        cc = cls.doXYPlot(fig2, s1, s2, xerrs=p1errs, yerrs=p2errs,
                            title='compared $\Delta\delta$', annotate=labels[t],
                            xlabel=prot1, ylabel=prot2)

        fig2.savefig('compareprotspans.png', dpi=200)'''
        return

    def mappKas(self, DB, col, pkainfo,names=None,spansfile=None,
                        nucleus='H',calculatespans=True):
        """Create mapping of pKas to nearby titratable groups, we use the
           dict of extracted pKas, but also requires the DB for structural info"""
        
        self.getProtNames(DB)
        path = os.getcwd()
        for prot in DB.getRecs():           
            name = self.protnames[prot]
            fname = name.replace(' ','').replace('(','').replace(')','')
            if names != None and name not in names:
                continue
            if not DB[prot].has_key(col):
                continue
            #structure    
            struct = DB[prot].Structure
            pdbname = fname+'.pdb'
            from PEATDB.Actions import DBActions
            DBActions.writePDB(struct, pdbname)
            #ekin project with exp data
            E = DB[prot][col]
            
            #actual titr grp pKa values we get from labbook
            try:
                titrpkas = DB.getLabbookSheet(name+'.pKas').data
                print 'got pka values for %s' %name
            except:
                print 'failed to get pKa values for %s' %name
            titrpkas = self.convertpkadict(titrpkas)
            print titrpkas
            
            #call CSP analysis
            if calculatespans==True:
                import pKaTool.Ghost.maps.Chemical_shift as ChemShift
                C = ChemShift.CSP(pdbname, method='Coulomb', pKaDict=titrpkas, nucleus=nucleus)
                calcspans = C.getSpanDict()
                self.savePickle(calcspans, os.path.join(path,fname+'_'+nucleus+'_spans.pickle'))
            else: 
                calcspans = self.loadPickle(fname+'_'+nucleus+'_spans.pickle')
                
            #assign ghosts - we use the spans dict to match our exp fits
            G = self.assignGhosts(pkainfo[name], calcspans, titrpkas, pdbname)            
            self.savePickle(G, fname+'_'+nucleus+'_ghosts.pickle')
            self.ghostStats(G)            
            '''#save csp curves to an ekinproject
            Ecsp = self.dict2Ekin(csps)
            Ecsp.saveProject('csps')
            #compare CSP curves to our actual curves
            rmsds = self.compareCSPCurves(E, Ecsp)'''
            
        return 
  
    def assignGhosts(self, pkainfo, calculatedspans, titrpkas, pdbfile):
        """Rough attempt to assign titrations from using calculated
           spans for all titr grps and exp fits"""
        ghosts = {}
        
        from pKaTool.Ghost.maps.CSP_explicit import CSP_coulomb
        import pKaTool.Ghost.maps.utilities_CSP as CSPutils
        CC = CSP_coulomb(None, None)
        coordmol, chargemol = CC.parseCoordCharges(pdbfile, pdb=True)
        
        #get titr residue coords first
        titrcoords = {}
        for res in titrpkas:            
            ch,resnum,rescode = res.split(':')
            r = int(resnum)            
            try:
                atom = CSPutils.distanceAtoms[rescode]
                titrcoords[res] = coordmol[ch][r][atom]
            except:
                pass
        
        #iterate over spans and get matching exp fits
        for res in calculatedspans:
            #calculated titr spans for this residue
            calcspans = calculatedspans[res]
            #get largest span
            mx = max(calcspans.values())
            tmax = calcspans.keys()[calcspans.values().index(mx)]
            #print res, tmax, mx
            ch,resnum = res.split(':')
            resnum = int(resnum)
            atom='N'
            coord = coordmol[ch][resnum][atom]

            #find exp fits for this residue            
            for d in pkainfo:
                F = pkainfo[d]
                r = int(F['resnum'])                
                rescode = F['res']
                if r != resnum:
                    continue
                pks = [i for i in F.keys() if 'pK' in i]
                ghosts[d]={}
                ghosts[d]['resnum']=r
                ghosts[d]['res']=F['res']
                for p in pks:
                    pka=F[p][p]
                    expspan=F[p]['span']
                    assignedspan=0
                    for t in calcspans:
                        cspan = calcspans[t]             
                        if abs(cspan) < assignedspan:
                            continue
                        if not titrcoords.has_key(t): continue                            
                        titrpka = titrpkas[t]                        
                        xx,distance = CSPutils.getDistance(coord, titrcoords[t])
                        if self.filterGhost(cspan, expspan, titrpka, pka, distance) == True:
                            print distance
                            print 'assigned %s to %s' %(res,t)
                            ghosts[d][p]=t
                            s = 'span' + p.strip('pKa')
                            ghosts[d][s] = expspan
                            assignedspan = abs(cspan)                            
        
        return ghosts
        
    def filterGhost(self, calcspan, expspan, titrpka, exppka, distance):
        """Apply criteria for assignment of ghost"""
        pkacutoff=1
        spancutoff=0.01
        #print calcspan, expspan, titrpka, exppka
        if abs(calcspan) < spancutoff:
            return False
        if calcspan>0 and expspan<0:
            return False
        if calcspan<0 and expspan>0:
            return False 
        if abs(exppka-titrpka) > pkacutoff:
            return False
        #if distance > 20:
        #    return False
        return True
        
    def ghostStats(self,G):
        """Some stats on ghost titrations"""
        i=0
        total = len(G)
        assigned=[]
        for r in G:
            for p in G[r].keys():
                if 'pK' in p:  
                    s = G[r][p]
                    if not r in assigned:
                        assigned.append(r)
                   
        print '%s of %s residues assigned' %(len(assigned),total)
 
        return
            
    def compareGhosts(self, g1, g2):
        """Compare two ghosts dicts"""
        
        return
        
    def compareCSPCurves(self, E1, E2):
        """Get normalised rmsd between 2 ekin proj"""
        rmsdvals={}
        for d1 in E1.datasets:
            ed1 = E1.getDataset(d1)                    
            r = ed1['res_num']
            if not r in E2.datasets:
                continue
            ed2 = E2.getDataset(r)
            x1,y1 = ed1.getxy()
            x2,y2 = ed2.getxy()
            #offset y1 to zero
            
            #get rmsd, we assume x pts. match            
            errs = sum([math.pow(i[0]-i[1],2) for i in zip(y1,y2)])           
            nrmsd = math.sqrt(errs/len(y1)) / (max(y1)-min(y1))            
            rmsdvals[r] = nrmsd
        return rmsdvals
        
    def dict2Ekin(self, data):
        E=EkinProject(mode='NMR titration')
        for d in data:
            ch, resnum = d.split(':')
            resnum = str(int(resnum))            
            x,y = (data[d].keys(), data[d].values())                  
            E.insertDataset(newname=resnum,xydata=(x,y))
            #edata = E.getDataset(resnum)
            #edata['res_num'] = resnum 
        return E
        
    def convertpkadict(self, t):
        """Convert pka info from table into format for CSP func"""        
        P={}
        s=':'
        for r in t:
            resnum = string.zfill(t[r]['resnum'], 4)   
            name = 'A'+s+resnum+s+t[r]['resname']
            if t[r]['pka'] != '':
                P[name] = float(t[r]['pka'])
           
        return P
    
    def getkey(self, d, value):
        """find the key(s) as a list given a value"""
        return [item[0] for item in d.items() if item[1] == value]

    def savePickle(self, data, filename):       
        f=open(filename,'w')       
        pickle.dump(data, f)        
        f.close()
        return
        
    def loadPickle(self, filename):
        f=open(filename,'r')
        data = pickle.load(f)
        f.close()
        return data        
    
    def correlatewithStructure(cls, ekindata, peatdb, pdbs=None):
        """Correlate a structural or other single property with each residue pKa/span"""
        #example here is for %surface accessibility

        import random

        count=0
        pkas=[]
        spans={}
        other={}
        for t in cls.titratable:
            spans[t]=[]
            other[t]=[]
        failedpdb=[]

        proteins = ekindata.keys()
        for prot in proteins:
            print prot
            print
            print 'doing protein', prot
            print '---------------------'
            #for structure stuff, we need pdb file here
            if peatdb != None:
                #get pdb file
                pdbfile = cls.save_PEAT_Structure(peatdb, prot)
                print 'trying to get pdb from PEAT db'
            elif pdbs.has_key(prot):
                pdbfile = pdbs[prot]
                #pdblines = cls.getPDBfromFile(pdbfile)
            if pdbfile == None:
                print 'no pdb file, skipping'
                failedpdb.append(protname)
                continue
            else:
                print 'got pdb'
            #get another attribute here eg. % rel. access.
            relacc = cls.getRelAccessibility(pdbfile)
            print relacc.keys()
            print

            E = EkinProject(ekindata[prot])
            for dataset in E.datasets:
                #resdata = cls.getResidueFields(ekindata[prot][dataset])
                resdata = E.getMetaData(dataset)
                res = resdata['residue']
                atom = resdata['atom']

                rk = cls.getPDBResidueName(resdata)
                #print 'trying', rk
                if res in cls.titratable:
                    if dataset in self.excluded:
                        continue
                    if relacc.has_key(rk):
                        accval = relacc[rk]['sum']['rel']
                        print 'found relative acc:', rk, accval
                        count+=1
                    else:
                        continue

                    if ekindata[prot]['__datatabs_fits__'].has_key(dataset):
                        fitdata = ekindata[prot]['__datatabs_fits__'][dataset]
                        #print 'fitdata', fitdata
                    pk,pkval,sp1 = cls.getMainpKa(fitdata, res=res)
                    print dataset, pk, pkval
                    if pkval != None:
                        pkas.append(pkval)
                        spans[res].append(sp1)
                        #other.append(random.uniform(pkval-0.4,pkval+0.4)*1.3)   #test
                        other[res].append(accval)
                    else:
                        print 'no primary pKa'
        print spans, other
        print 'failed to get pdbs for:', failedpdb
        print 'acc values found for %s residues'  % count
        pylab.figure(4)
        #cls.doXYPlot(pkas, other, title='compare x', xlabel='pkas', ylabel='y')
        #cls.doXYPlot(spans, other, title='Rel Access. vs $\delta$ shift', xlabel='$\delta$ chem. shift', ylabel='y')
        pylab.gcf().text(0.5, 0.95, 'Rel Access. vs $\delta$ shift', horizontalalignment='center',
                        fontproperties=FontProperties(size=18))

        i=1
        for r in spans.keys():
            pylab.subplot(2,2,i)
            cls.doXYPlot(spans[r], other[r], title=r) # xlabel='$\delta$ chem. shift', ylabel='% Surface Access.')
            i+=1

        pylab.show()
        return

    #
    # Utility classes for ekin titration data
    #

    @classmethod
    def getResidues(cls, E, titratable=False):
        """Get titratable residues in an ekin project
           Returns a tuple of dataset name, residue code, residue no."""
        t=[]
        for d in E.datasets:
            edata = E.getDataset(d)            
            if edata.has_key('residue'):
                if not edata.has_key('residue') or not edata.has_key('res_num'):
                    continue
                res = edata['residue']
                res_num = edata['res_num']
                if titratable == True and not res in cls.titratable:
                    continue
                else:                    
                    t.append((d,res,res_num))
        return t
    
    @classmethod
    def getResidue(cls, recdata):
        """get residue info if present"""
        if recdata.has_key('residue'):
            res =  recdata['residue']
        else:
            res = None
        return res

    @classmethod
    def getResidueFields(cls, recdata):
        """Get the residue info from the ekin dataset, if present"""
        residue = cls.getResidue(recdata)
        if recdata.has_key('res_num'):
            res_num = recdata['res_num']
        else:
            res_num = None
        if recdata.has_key('atom_type'):
            atom = recdata['atom_type']
        else:
            atom = None
        if recdata.has_key('chain_id'):
            chainid = recdata['chain_id']
        else:
            chainid = None
        return {'chain_id': chainid,'residue':residue, 'res_num':res_num, 'atom':atom}

    @classmethod
    def getPDBResidueName(cls, resdata):
        """Convert ekin residue info fields to Protool/PDB format"""
        import Utils
        s =':'
        resnum = Utils.leadingZeros(resdata['res_num'],4)
        if resdata['chain_id'] != None:
            resname = resdata['chain_id'] + s + resnum + s + resdata['residue']
        else:
            resname = s + resnum + s + resdata['residue']
        return resname

    @classmethod
    def setResidueNames(cls, E):
        """Try to set residue names based on dataset labels"""
        import string 
        for d in E.datasets:
            #edata = E.getDataset(d)
            name = string.upper(d[0:3])            
            found = 0
            if name in cls.residue_list:
                #edata['residue'] = name
                E.addMeta(d, 'residue', name)
                found = 1
            #if not found from 1st three letters try one-letter code
            if found == 0:
                if name[0] in cls.residue_letters.keys():
                    res = cls.residue_letters[name[0]]
                    #edata['residue'] = res
                    E.addMeta(d, 'residue', res)
        return E

    @classmethod
    def setResidueNumbers(cls, E):
        """Try to set residue numbers from dataset names"""
        import re
        r=re.compile('\d+')  
        for d in E.datasets:
            #edata = E.getDataset(d)
            nums = r.findall(d)
            #print nums
            if len(nums)>0:
                #edata['res_num'] = nums[0]
                E.addMeta(d, 'res_num', nums[0])
        return E
        
    @classmethod
    def setAtomTypes(cls, E, atom=None):
        d_titatoms = {
            'GLU':'CD',
            'ASP':'CG',
            'GLN':'CD',
            }     
        for d in E.datasets:    
            if atom!=None:              
                E.addMeta(d, 'atom_type', atom)                       
        return E

    @classmethod
    def setMetaInfo(self, E, atom=None):
        """Auto add Atom, residue info to ekin NMR data"""
        t = TitrationAnalyser()
        t.setResidueNames(E)
        t.setResidueNumbers(E)
        if atom!=None:
            t.setAtomTypes(E, atom=atom)        
        return E
        
    @classmethod
    def setMetaInfoFromDB(self, DB, col, proteins=None, atom=None):
        """Auto add Atom, residue info to ekin NMR data"""
        if proteins == None:
            proteins = DB.getRecs()
        for prot in proteins:
            E = DB[prot][col]
            DB[prot][col] = self.setMetaInfo(E, atom=atom)
        return
       
    @classmethod
    def checkStructMapProp(cls, structmap, prop):
        """Check if a particular property is in the struct map"""
        found=0
        for i in structmap.keys():
            for props in structmap[i]['Data property']:
                if prop in props:
                    found=1
                    break
        return found

    @classmethod
    def getMainpKa(cls, fitdata, datapoints=None, res=None, strict=True,
                     minspan=0.06, maxspan=5):
        """Check fitdata for the most reliable pKa using criteria
           and return its key name ie. 'pka1', pka value and span
           if strict is true(default), we use our criteria to filter them
        """
        accept = True

        try:
            model = fitdata['model']
        except:
            return None, None, None
        if len(fitdata) < 2:
            return None, None, None
        params = cls.get_parameter_names(model)
        pkas={}
        spans={}
        for p in params:
            if 'pK' in p:
                i = cls.get_param_index(p,model)
                pkas[p] = fitdata[i]
                #get span val for that pka in the fitdata
                num=p.strip('pKa')
                j = cls.get_param_index('span'+num,model)
                spans[p] = abs(float(fitdata[j]))

        #print 'pkas',pkas
        #print 'spans',spans
        if cls.modelpkas.has_key(res):
            modelpka = cls.modelpkas[res]
        else:
            modelpka = None

        if len(spans) == 0:
            return None, None, None
        elif len(spans) == 1:
            pkmain = 'pKa'
        else:
            pkmain = Utils.getDictMax(spans)

        errors = ['span is less than cutoff','too close to another pKa','other spans comparable scale',
                    'too far from the model pKa','too close to edge of data']

        span = spans[pkmain]
        pkval = pkas[pkmain]

        if span <=minspan or span >= maxspan:
            #print errors[0]
            accept = False
        if strict == True and accept == True:
            #apply our criteria here
            for p in spans.keys():
                if p != pkmain:
                    #too close to another pKa?
                    if abs(pkas[p]-pkval) < 0.2:
                        accept=False
                        #print errors[1]
                    #other spans comparable scale?
                    if span/spans[p] < 1.2:
                        accept=False
                        #print errors[2]
                    #too far from the model pKa?
                    if modelpka != None:
                        if abs(modelpka-pkval) > 3:
                            accept=False
                            #print errors[3]
                    #too close to edge of ph range - titration not finished
                    if datapoints != None:
                        if (pkval - phmin) < 0.3 or (phmax - pkval) < 0.3:
                            acccept=False
                            #print errors[4]

        if accept == True:
            return pkmain, pkval, span
        else:
            #return 'NA', pkval, span
            return None, None, None

    @classmethod
    def getallpKas(cls, fitdata, minspan=0):
        """Get all the pkas, spans in the fit as 2 lists"""
        try:
            model = fitdata['model']
        except:
            return None, None, None
        params = cls.get_parameter_names(model)
        pnames=[]
        pkas=[]
        spans=[]
        if fitdata==None or len(fitdata)<2:
            return None, None, None
        for p in params:
            if 'pK' in p:
                i = cls.get_param_index(p,model)
                #get span val for that pka in the fitdata
                num=p.strip('pKa')
                j = cls.get_param_index('span'+num,model)
                sp = abs(float(fitdata[j]))
                if sp < minspan:
                    continue
                pkas.append(fitdata[i])
                pnames.append(p)                
                spans.append(sp)

        #print pkas, spans
        return pnames, pkas, spans

    @classmethod
    def getPDBfromPEAT(cls, peatDB, protein):
        """If within PEAT we can try to get the structure from the DB"""
        if peatDB == None:
            print 'No peat data'
            return
        import types
        print peatDB[protein].keys()
        if type(peatDB[protein]['Structure']) is types.ListType:
            # Real X-ray structure
            import Protool
            X = Protool.structureIO_fast()
            pdblines = peatDB[protein]['Structure']
            fd = open('junkjens.pdb','w')
            for line in pdblines:
                fd.write(line)
            fd.close()
            try:
                X.parsepdb(pdblines)
                return pdblines, X
            except:
                print 'Protool has fafiled to parse PDB file'
                print
                return None, None
        else:
            return None, None

    @classmethod
    def save_PEAT_Structure(cls, peatdata, protein):
        """Save a PDB file from a Structure cell"""
        pdblines,X = cls.getPDBfromPEAT(peatdata, protein)
        if pdblines == None:
            return None
        import os
        filename=os.path.join(os.getcwd(),'temp')
        if filename:
            if filename[-4:]!='.pdb':
                filename=filename+'.pdb'
            fd=open(filename,'w')
            for line in pdblines:
                fd.write(line)
            fd.close()
        return filename

    @classmethod
    def getPDBfromFile(cls, filename):
        """Get pdblines from a PDB file"""
        import Protool
        X=Protool.structureIO_fast()
        X.readpdb(filename)
        pdblines = X.writepdb('dummy',nowrite=1)
        return pdblines, X

    @classmethod
    def fetchPMIDSummary(cls, pmid):
        """Get summary of paper from pmid"""
        h=Entrez.esummary(db="pubmed", id=pmid)
        r=Entrez.read(h)
        author = r[0]['AuthorList']
        title = r[0]['Title']
        return author, title

    @classmethod
    def getRelAccessibility(cls, pdbfile=None, pdblines=None):
        """Get relative accessibilities for structure"""
        import pKarun.WI_tools as wi
        print 'Getting rel acc. using', pdbfile
        if pdbfile != None:
            r = wi.relative_accessibility(pdbfile)
            #print 'rel_acc',r
        return r

    @classmethod
    def plotHistograms(cls, ax, data, bins=25, title='', xlabel=None, ylabel=None, color=None,
                          colors=False, xticks=True, yticks=True, xlim=None):
        """Do a pylab histogram overlaying 1 or more dicts on the same axis
            and using the same bins, so that they are clearly visible"""
        import Utils
        pyclrs = cls.pylabcolors
        if color == None:
            clr = 'blue'
        else:
            clr = color

        patches={}
        for h in data.keys():
            if len(data[h])==0:
                del data[h]
                continue
        if len(data)==0:
            return

        handles=[]
        if None in data.keys():
            del data[None]
        names = data.keys()

        legnames=[]
        for key in names:
            if key in cls.atom_greek:
                legnames.append(cls.atom_greek[key])
            else:
                legnames.append(key)


        #if len(data) > 1:
        mybins = Utils.get_bins(data, bins)
        w= max(mybins)/len(mybins)
        if mybins == None:
            return
        #else:
        #    mybins = bins
        #    w=1.0/bins

        h=0
        if len(data) > 0:
            for key in names:
                if len(data[key]) == 0:
                    h+=1
                    continue
                n, b, patches[h] = pylab.hist(data[key], bins=mybins, facecolor=pyclrs[h])#, alpha=0.9)
                hd, x = numpy.histogram(data[key], bins=mybins)
                #patches[h] = ax.bar(mybins[:len(mybins)-1], hd, width=w, color=pyclrs[h],
                #                    alpha=0.8)
                handles.append(patches[h][0])
                h+=1

            leg = ax.legend(handles, legnames, shadow=True, numpoints=1, loc='best',
                            markerscale = 0.8, prop = FontProperties(size='smaller'),
                            handletextsep = 0.04)

            leg.draw_frame(False)
        i=0
        '''for pa in patches:
            for rect in patches[pa]:
                #print patches[pa]
                #rect.set_x(rect.get_x() + i)
                #print p.get_zorder()
                #p.set_linewidth(1.1)
                rect.set_edgecolor('#737CA1')
                rect.set_zorder(0)
                #p.set_antialiased(False)
            i+=w'''
        #reorder bars, so that smallest are at front
        Utils.adjust_bars(patches)

        ax.set_title(title)

        if xlabel != None:
            ax.set_xlabel(xlabel)
        if ylabel != None:
            ax.set_ylabel(ylabel)
        if xlim != None:
            ax.set_xlim(xlim)
        if xticks != True:
            ax.set_xticks(pylab.xlim())
        if yticks != True:
            ax.set_yticks(pylab.ylim())

        return

    @classmethod
    def doMultipleHistograms(cls, fig, recs, bins=20, title='', xlabel=None, ylabel=None, color=None,
                         subplots=True, colors=False, xticks=True, yticks=True, xlim=None):
        """Do a pylab histogram of a dict of 1 or more dicts """		
        subplots=[]
        dim=int(math.ceil(len(recs)/2.0))
        i=1

        for r in recs:
            if len(recs[r])==0:
                i=i+1
                continue
            ax = fig.add_subplot(dim,2,i)
            cls.plotHistograms(ax, recs[r], title=r, bins=bins, xlabel=xlabel,ylabel=ylabel,
                                xticks=xticks, yticks=yticks, xlim=xlim)
            i=i+1

        fig.subplots_adjust(hspace=0.4)
        return fig

    @classmethod
    def do3DHistograms(cls, fig, recs, bins=20, title='', xlabel=None, ylabel=None,
                        xticks=True, yticks=True, xlim=None):
        """3d histograms"""
        from mpl_toolkits.mplot3d import Axes3D
        i=0
        print recs
        for r in recs:
            ax = Axes3D(fig)
            if len(recs[r])==0:
                i=i+1
                continue
            z=0
            for k in recs[r]:
                data = recs[r][k]
                if len(data) == 0:
                    continue
                #mybins = Utils.get_bins(data, bins)
                hd, x = numpy.histogram(data, bins=10)#, bins=mybins)
                ax.bar(range(len(hd)), hd, zs=z, zdir='y', alpha=0.8)
                z+=10

            i=i+1

        return
    
    @classmethod
    def rmse(cls, ar1, ar2):
        """Mean squared error"""
        ar1 = numpy.asarray(ar1)
        ar2 = numpy.asarray(ar2)
        dif = ar1 - ar2
        dif *= dif
        return numpy.sqrt(dif.sum()/len(ar1))

    @classmethod
    def doXYPlot(cls, ax, x, y, names=None, err=0.5,
                    title=None, xerrs=None, yerrs=None,
                    xlabel=None, ylabel=None,
                    xaxislabels=None, color=None, symbol=None, markersize=30):
        """Do xy plot of 2 lists and show correlation
           annotate is a list of tuples with x,y coord and text to print at that pt"""
           
        if len(x) == 0:
            return None
        if color == None:
            clr = 'b'
        else:
            clr = color
        if symbol == None:
            symb = 'o'
        else:
            symb = symbol

        if min(x)<min(y): a=min(x)-1
        else: a=min(y)-1
        if max(x)>max(y): b=max(x)+1
        else: b=max(y)+1    
        
        ax.plot((a,b),(a,b),color='black')
        ax.scatter(x, y, facecolor=clr, marker=symb, s=markersize, 
                            picker=4, alpha=0.6)

        if names!=None:
            c=0
            for i in zip(x,y):                
                if abs(i[0]-i[1])>err:  
                    #z = pylab.Circle((i[0], i[1]), 0.2,fill=False,alpha=0.7)
                    #ax.add_patch(z)
                    ax.annotate(names[c], (i[0]+0.1, i[1]),
                                xytext=None, textcoords='data',
                                fontsize=10)
                c+=1
                    
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if xaxislabels != None:         
            ax.set_xticklabels(xaxislabels, rotation='vertical', fontproperties=ft)
        if xerrs!=None or yerrs!=None:
            errline = ax.errorbar(x, y, xerr=xerrs, yerr=yerrs, fmt=None,
                                        elinewidth=.5, ecolor=clr, alpha=0.7)
        ax.set_xlim(a,b); ax.set_ylim(a,b)    
        cc = round(numpy.corrcoef(numpy.array([x,y]))[0][1],2)
        rmse = round(cls.rmse(x,y),2)
        print 'corr. coeff:', cc
        print 'RMSE:', rmse
        return cc, rmse

    @classmethod
    def get_param_index(self, paramname, model):
        """Returns an integer, the index of the input variable in dict for a given model"""
        index=None
        index_list = self.get_parameter_names(model)
        i=0
        for p in index_list:
            if p == paramname:
                index = i
            i=i+1
        return index
    
    @classmethod
    def get_parameter_names(self, model):
        X= Fitting.getFitter(model)
        names = X.getVarNames()
        return names
