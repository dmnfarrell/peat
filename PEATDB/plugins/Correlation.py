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

"""Feb 2010"""

try:
    from Plugins import Plugin
except:
    from PEATDB.Plugins import Plugin
from Tkinter import *
import Pmw
import os, sys, math, string, re
import math, numpy
from scipy.stats import stats
try:
    import matplotlib
    #matplotlib.use('Agg')
    from matplotlib.font_manager import FontProperties
    import matplotlib.pyplot as plt
    #plt.rc('text', usetex=True)
except:
    pass
from PEATDB.TableModels import TableModel
from PEATDB.Ekin.Base import EkinDataset
from PEATDB.Ekin.Fitting import Fitting

class CorrelationAnalyser(Plugin):
    """A class for processing specific Kinetics exp data"""
    """Author: Damien Farrell"""

    capabilities = ['gui','uses_sidepane']
    requires = ['pylab','scipy']
    menuentry = 'Correlation Analysis'
    #gui_methods = ['plotdeltatmpH', 'analyseTempResults',
    #              'ploterrvsMutations', 'analyseActivityResults']
    gui_methods = { 'doCorrelation':'show Correlation',
                    'ploterrvsMutations':'plot err vs mutations (10R)',
                    'plotdeltatmpH':'plot delta Tms (10R)'}
    about = 'This plugin is designed for some correlation analysis and for 10R project'

    def __init__(self):
        self.phlist = ['6','6.5','7','7.5','8','8.5','9','10']
        return

    def main(self, parent):
        if parent==None:
            return
        self.parent = parent
        self.DB = parent.DB
        if self.DB == None:
            return
        self._doFrame()
        return

    def _doFrame(self):
        if 'uses_sidepane' in self.capabilities:
            self.mainwin = self.parent.createChildFrame()
        else:
            self.mainwin=Toplevel()
            self.mainwin.title('Corerlation Analysis')
            self.mainwin.geometry('800x600+200+100')

        methods = self._getmethods()        
        methods = [m for m in methods if m[0] in self.gui_methods.keys()]
        self._createButtons(methods)        
        self.smenu = Pmw.OptionMenu(self.mainwin,
                labelpos = 'w',
                label_text = 'Sheet:',
                items = self.DB.meta.labbook.keys(),
                command=self.updateColsMenu,                
                menubutton_width = 8)
        self.smenu.pack(side=TOP,fill=BOTH)     
        self.updateColsMenu()
        return

    def _getmethods(self):
        """Get a list of all available public methods"""
        import inspect
        mems = inspect.getmembers(self, inspect.ismethod)
        methods = [m for m in mems if not m[0].startswith('_')]
        return methods

    def _createButtons(self, methods):
        """Dynamically create buttons for supplied methods, which is a tuple
            of (method name, label)"""     
        for m in methods:           
            b=Button(self.mainwin,text=self.gui_methods[m[0]],command=m[1])
            b.pack(side=TOP,fill=BOTH)
        return

    def updateColsMenu(self, evt=None):
        """create drop menus of available cols in chosen sheet"""
        sheet = self.smenu.getcurselection()
        S=self.DB.meta.labbook[sheet]
        model = TableModel(S)
        names = model.columnNames
        if hasattr(self, 'xcolsmenu'):
            self.xcolsmenu.destroy()
            self.ycolsmenu.destroy()
            self.filterbymenu.destroy()
            self.filtervalentry.destroy()
            self.bf.destroy()
        self.xcolsmenu = Pmw.OptionMenu(self.mainwin,
                labelpos = 'w',
                label_text = 'X-Col:',
                items = names,                
                menubutton_width = 10)
        self.xcolsmenu.pack(side=TOP,fill=BOTH) 
        self.ycolsmenu = Pmw.OptionMenu(self.mainwin,
                labelpos = 'w',
                label_text = 'Y-Col:',
                items = names,                
                menubutton_width = 10)
        self.ycolsmenu.pack(side=TOP,fill=BOTH)
        self.bf=Frame(self.mainwin); self.bf.pack(side=TOP,fill=BOTH)
        names.extend('-')
        self.filterbymenu = Pmw.OptionMenu(self.bf,
                labelpos = 'w',
                label_text = 'Filter by:',
                items = names,
                initialitem='-',
                menubutton_width = 5)
        self.filterbymenu.pack(side=LEFT,fill=BOTH)
        self.filtervalentry = Pmw.EntryField(self.bf,
                labelpos = 'w',
                label_text = 'Value:')
        self.filtervalentry.pack(side=LEFT,fill=BOTH)
        return
        
    def doCorrelation(self):
        """From gui"""
        if self.filterbymenu.getcurselection() != '-':
            filt = (self.filterbymenu.getcurselection(), self.filtervalentry.getvalue())
        else:
            filt = None
        self.plotCorrelation(sheet=self.smenu.getcurselection(), 
                             xcol=self.xcolsmenu.getcurselection(),
                             ycol=self.ycolsmenu.getcurselection(),
                             filterby=filt)        
        return        
     
    @classmethod 
    def simpleCorrelation(self, ax, x, y):
        if len(x) == 0 :
            return      
        line = ax.scatter(x, y, marker='o',alpha=0.8)
        cl = numpy.arange(-10,25)
        ax.plot(cl, cl, 'g', alpha=0.5)
        ax.set_xlabel('Predicted')
        ax.set_ylabel('Exp')
        ax.set_xlim(-10,20); ax.set_ylim(-10,20)
        ax.axhline(y=0, color='grey'); ax.axvline(x=0,color='grey')
        ax.set_title(g)
        n+=1
        return
     
    def showHistogram(self, x):
        f=plt.figure()
        ax=f.add_subplot(111)
        ax.hist(x)
        ax.set_title('%s points' %len(x))
        return
        
    def plotCorrelation(self, sheet=None, filterby=None, 
                          xcol='exp', ycol='Total', labelcol=None,
                          labeloutliers=False,
                          plottitle='Predicted vs Experimental',
                          xlabel='Predicted',
                          ylabel='Experimental',
                          limit=None,
                          filename='correlation.png'):
        """Show exp vs pred. from a Labook in DB - General case,
            filterby is a tuple of colname, value"""

        #plt.rc('text', usetex=True)
        model = self.DB.getLabbookSheet(sheet)
        colors = ['b','g','r','y','m','c']
        
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)        
        legs=[]; legnames=[]
        bad=[]; good=[]       
       
        x=self.tofloats(model.getColumnData(columnName=xcol, filterby=filterby))
        y=self.tofloats(model.getColumnData(columnName=ycol, filterby=filterby))
        if len(x) ==0:
            return
        errs = [a - b for a, b in zip(x, y)]
        line = ax.scatter(x, y, color='red', alpha=0.7)
                                  
        if labelcol != None:
            labels = model.getColumnData(columnName=labelcol,filterby=filterby)      
            #user peat_sa stats tools to get outliers here
            #outliers = self.findOutliers(x,y,errs)
            '''for l in labels:
                i=labels.index(l)
                #label and list outliers
                if l in outliers:
                    #c = plt.Circle((x[i], y[i]), 0.4,fill=False,alpha=0.7)
                    ax.annotate(l, (x[i]+0.4, y[i]), xytext=None, textcoords='data',
                                        fontsize=8)
                    #ax.add_patch(c)'''

            for e in errs:
                i=errs.index(e)
                #label bad ones
                if abs(e)>5:# and y[i]>8:
                    ax.annotate(labels[i], (x[i]+0.4, y[i]), xytext=None, textcoords='data',
                                            fontsize=8)

        #legs.append(line)
        #legnames.append(g)

        #draw expected correlation line with slope x
        slope=1
        #set range of axes
        if max(y)>max(x): lims = min(y)-5,max(y)+5
        else: lims = min(x)-5,max(x)+5
        cx = numpy.arange(lims[0],lims[1])
        cy = [slope*i for i in cx]
        ax.plot(cx, cy, 'g', alpha=0.7)
        ax.plot(cx-4, cy, 'g', alpha=0.5,linestyle='--')
        ax.plot(cx+4, cy, 'g', alpha=0.5,linestyle='--');
        ax.axhline(y=0, color='grey'); ax.axvline(x=0,color='grey')
        #ax.set_xlim(min(x)-2,max(x)+2)
        if limit==None:
            limit=lims[1]  
        ax.set_xlim(lims[0],limit)
        ax.set_ylim(lims[0],limit)           
        ax.set_title('%s points' %len(x))
        #cc = str(round(pow(stats.pearsonr(x,y)[0],2),2))
        #ax.text(1,16, r'$r^2= %s$' %cc, fontsize=16)
        '''finalx=[x[i] for i in range(len(labels)) if labels[i] not in outliers]
        finaly=[y[i] for i in range(len(labels)) if labels[i] not in outliers]
        cc1 = str(round(pow(stats.pearsonr(finalx,finaly)[0],2),2))
        #ax.text(1,14, r'$r^2_{final}= %s$' %cc1, fontsize=16)'''

        #fig.legend(legs, legnames)
        fig.suptitle(plottitle)
        fig.savefig(filename)
        
        #self.showHistogram(x); self.showHistogram(y)
        plt.show()
        return ax
                  
    def analyseCorrelation(self, sheet, filterby=None,
                            xcol='exp', ycol='Total', labelcol=None):
        """Analyse correlation"""
        model = self.DB.getLabbookSheet(sheet)
        x=self.tofloats(model.getColumnData(columnName=xcol, filterby=filterby))
        y=self.tofloats(model.getColumnData(columnName=ycol, filterby=filterby))
        if len(x) ==0:
            return
        labels = model.getColumnData(columnName=labelcol,filterby=filterby)
        names = model.getColumnData(columnName='name',filterby=filterby)
        errs = [a - b for a, b in zip(x, y)]
        print 'points with largest errors:\n'
        for e in errs:
            i=errs.index(e)            
            if abs(e)>5:
                print names[i],labels[i], e   
        print '\n'
        #get outliers
        
        print 'no. of points:', len(x)
        print 'correlation co-eff:', numpy.corrcoef(x,y)[0][1]
        print 'stdev of errors:', numpy.std(errs)

        return   

    def plotdeltatmpH(self):
        """Plot deltatms vs pH for variants with ok fits"""
        S = self.DB.meta.labbook['ok']
        model = TableModel(S)

        grps = ['ok','rt']
        names = []
        labels=model.getColumnData(columnName='name')
        names={}
        for l in labels: names[l]=[]

        for ph in self.phlist:
            dtms= self.tofloats(model.getColumnData(columnName='avdeltatm_'+ph))
            fit= model.getColumnData(columnName='fit_'+ph)
            for i in range(len(labels)):
                name = labels[i]
                if name in names:
                    if fit[i] not in grps:
                        del names[name]
                    else:
                        names[name].append(dtms[i])

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('pH')
        ax.set_ylabel(r'$\delta$Tm')
        legs=[]
        phs = self.tofloats(phlist)
        for n in names:
            line = ax.plot(phs, names[n], marker='x',linewidth=3,alpha=0.8)
            legs.append(line)
        fig.legend(legs, names)
        fig.suptitle(r'$\delta$Tm vs pH',fontsize=20)
        plt.show()
        fig.savefig('dtmvspH.png')

        #analyse stats of ph vs dtm matrix
        x=[names[n] for n in names]
        a = numpy.array(x).transpose()
        for i in a:
            print numpy.std(i)
        return names
    
    def ploterrvsMutations(self):
        """Plot errors vs. no. mutations"""

        S = self.DB.meta.labbook['ok']
        model = TableModel(S)
        nummutations = {}
        for ph in self.phlist:
            x=self.tofloats(model.getColumnData(columnName='Total',filterby=('fit_'+ph,'ok')))
            y=self.tofloats(model.getColumnData(columnName='avdeltatm_'+ph,filterby=('fit_'+ph,'ok')))
            nmuts=model.getColumnData(columnName='numberMutations',filterby=('fit_'+ph,'ok'))

            if len(x) == 0:
                continue
            errs = [pow(a - b,2) for a, b in zip(x, y)]
            for e in errs:
                i=errs.index(e)
                n=nmuts[i]
                if not nummutations.has_key(n):
                    nummutations[n]=[]
                nummutations[n].append(e)

        averrs=[]; stds=[]
        for n in nummutations:
            averrs.append(numpy.mean(nummutations[n]))
            stds.append(numpy.std(nummutations[n]))
            print n, len(nummutations[n])
        print averrs,stds, nummutations.keys()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('number of mutations')
        ax.set_ylabel('error')
        ind = numpy.array([int(i) for i in nummutations.keys()])
        ax.bar(ind, averrs, yerr=stds, linewidth=2,alpha=0.8)
        fig.suptitle('Number of Mutations vs Average Error',fontsize=20)
        fig.savefig('ploterrvsMutations.png')
        plt.show()
        return

    def plotSequenceInfo(self):
        """Plot some attribute as a func of redisue"""
        pdb = open('/enc/Novozymes/Parse_mutations/10R.pka.pdb')
        S = self.DB.meta.labbook['tempvspredicted']
        model = TableModel(S)
        mutations = {}
        x=model.getColumnData(columnName='mut')
        import re
        for m in x:
            vals = m.split('+')
            res = [int(re.sub("[\D]", "", i)) for i in vals]
            #print res
            for r in res:
                if not r in mutations:
                    mutations[r]=1
                else:
                    mutations[r]+=1
        ind=[];y=[]

        for m in sorted(mutations.keys()):
            ind.append(m)
            y.append(mutations[m])

        fig = plt.figure(figsize=(10,5))
        ax = fig.add_subplot(111)
        ax.bar(ind, y, linewidth=0,facecolor='black',alpha=0.8)
        ax.set_xlabel('Residue')
        ax.set_ylabel('number of mutations')
        #ax.set_xticks(numpy.array(ind)+0.5)
        #ax.set_xticklabels(ind,size=10)
        fig.autofmt_xdate()
        fig.suptitle('Number of Mutations vs Residue',fontsize=20)
        fig.savefig('seqvsNumberMutations.png')

        bfacs = []
        for i in range(1,189):
            if mutations.has_key(i): 
                bfacs.append(mutations[i])
            else:
                bfacs.append(0)

        return bfacs
               
    def QQplot(self, data=None, labels=None, n=3):
        """Do Q-Q plot to determine if sample is normally distributed"""

        from scipy.stats import norm
        import random
        outliers = []
        #test data
        if data==None:
            data = numpy.random.normal(loc=0, scale=1.0, size=1000)
            numpy.append(data,[45,173,201]) # add some outliers
            labels = ([random.choice(string.letters) for x in range(len(data))])
        if labels!=None:
            lmap = sorted(zip(data, labels))
        m=numpy.mean(data)
        sd=numpy.std(data)
        y = sorted([(i-m)/sd for i in data])
        #expected normal scores as ranked (z scores)
        x = [(i-0.5)/len(data) for i in range(1,len(data)+1)]
        #inverse cumul. distribution
        x = norm.ppf(x)

        def markoutliers():
            r=range(n)+range(len(x)-n, len(x))
            for i in r:
                ax.text(x[i]+0.2, y[i], lmap[i][1], size=6)
                outliers.append(lmap[i][1])
            return

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('normal distr')
        ax.set_ylabel('data')
        ax.plot(x, y, 'o', alpha=0.7)
        ax.plot(x, x, '-', alpha=0.7)
        if labels!=None:
            markoutliers()
        fig.suptitle('Q-Q plot')
        fig.savefig('qq.png')
        return outliers

    def tofloats(self, l):
        x=[]
        for i in l:
            x.append(float(i))
        return x

    def loadDB(self, local=None):
        from PEATDB.Base import PDatabase
        if local != None:
            self.DB = PDatabase(local=local)
        return    

    def mapStructure(self, positions=None):
        """Map some of the results on the structure"""
        
        s=self.DB['10R WT'].Structure
        #bf = self.plotSequenceInfo()
        from PEATDB.Actions import DBActions
        DBActions.DB = self.DB
        path = '/local/bin/yasara/'
        DBActions.displayStructure('10R WT', 'Structure', 
                                    molapp='yasara', path=path, color='green')       
        Y = DBActions.yasara       
        Y.SaveBMP('struct.bmp')
        if positions!=None:
            for p in positions:
                Y.ColorRes(p,'red')
        while Y != None:
            wait()
        return
        
    def extract10RStability(self, DB=None, ph='7', jobname='stab', save=False):
        """10R specific method.
           This gets the temp stability results from the 10R DB and creates
           sheets at this ph, combining them with PEATSA predictions.
           We can then run the correlations in those sheets"""
           
        if DB == None:
            DB=self.DB
        symb = unichr(177).encode("UTF-8") 
        S = DB.createLabbookSheet('stability'+ph)
        #fit quality
        F = self.DB.getLabbookSheet('fits')
        fits = dict(zip(F.getColumnData(columnName='name'),F.getColumnData(columnName='fit_'+ph)))
        #factors to get deltaG
        P = self.DB.getLabbookSheet('slopes') 
        #confirmed sequences
        C = self.DB.getLabbookSheet('confirmedseq')
        confirmed = dict(zip(C.getColumnData(columnName='name'), C.getColumnData(columnName='Tapo.check')))

        #get wt tm
        wtlist = ['wt 1a','wt 2a','wt 4a','wt 2','wt 3','wt 5']
        wtvals = []
        phc=ph.replace('.',',')
        for name in wtlist:
            E = DB[name].temp
            try:                
                wtvals.append(E.getMetaData('1_'+phc)['t50'])
                wtvals.append(E.getMetaData('2_'+phc)['t50'])
            except:
                continue
        wtTm = numpy.mean(wtvals)
        wtTmerr = round(numpy.std(wtvals),2)
        print ph, 'wt Tm:', wtTm, symb, wtTmerr       
        
        #we get delta tms and insert into sheet and also create deltaG vals
        #by dividing deltaS/deltaTm

        for name in DB.getRecs():            
            if not hasattr(DB[name], 'temp') or not hasattr(DB[name], 'Mutations'):
                continue
            if not name in confirmed or confirmed[name] == '':
                continue
                
            Et = DB[name].temp
            #average of 2 tms
            try:
                tm = round(numpy.mean(Et.getMetaData('1_'+phc)['t50'],Et.getMetaData('2_'+phc)['t50']),2)
                tmerr = round(numpy.mean(Et.getMeta('1_'+phc, 'exp_errors')['t50'][1], E.getMeta('2_'+phc, 'exp_errors')['t50'][1]),2) 
                deltatm = round(wtTm-tm,3)
                ddg = deltatm * 1.3
            except:
                continue
            #print tm, deltatm, tmerr            
            S.addRecord(name, Mutations=DB[name].Mutations,
                        tm=tm, tmerr=tmerr, deltatm=deltatm, ddg=ddg, fit=fits[name])
        
        #merge predictions
        from PEATSA import PEATSAPlugin
        PSA = PEATSAPlugin()
        PSA.main(DB=DB)  
        job, name = PSA.getJob(jobname)
        if job.error() != None or job.state() != 'Finished':
            return                                    
        stabmatrix = job.data.stabilityResults      
        PSA.mergeMatrix(stabmatrix, S)
        #add error column
        S.addColumn('error')
        S.addColumn('sqerror')
        for r in S.reclist:
            S.data[r]['error'] = float(S.data[r]['Total']) - float(S.data[r]['ddg'])
            S.data[r]['sqerror'] = math.pow(S.data[r]['error'],2)
        self.DB.saveLabbook('stability'+ph, S)
        self.DB.saveLabbooktoFile('10R_new.labbook')
        if save==True:
            self.DB.commit()
        return 

    def testResActFit(self, DB=None):
        """Test how fitting the Tm to res activity data varies with 
           actual Tm using simulated data with noise added"""
           
        from mpl_toolkits.axes_grid import AxesGrid
        if DB == None:
            DB=self.DB

        f=plt.figure(figsize=(10,6))
        grid = AxesGrid(f,111,nrows_ncols=(4,5),axes_pad=0.3,
                    share_all=True,aspect=False)                
           
        def getdata(tm=60,noise=False):
            #create fake data, noise is roughly based on errors in velocities from 10R assays
            x=[25.0,60.0,65.0,70.0,75.0]
            y=[];m=20
            n=[]
            for T in x:                
                val = 1.0/(1+math.pow((T/tm),m))
                y.append(val)               
                noise=numpy.random.normal(0, 0.03)
                val=val+noise
                n.append(val)
            return x,y,n
        
        tms=numpy.arange(40,78,2.0)
        c=0
        errors=[]; errors1=[]
        for i in tms:
            x,y,n=getdata(i)
            A,X=Fitting.doFit(expdata=zip(x,y),model='Residual Activity',silent=True)
            A,Xn=Fitting.doFit(expdata=zip(x,n),model='Residual Activity',silent=True)
            fitx=numpy.arange(0,100)
            fity = X.getFitLine(fitx)
            fitn = Xn.getFitLine(fitx)
            tm=X.getFitDict()['t50']      
            tm1=Xn.getFitDict()['t50']
            ax=grid[c]
            p=ax.plot(x,y,'x',mew=3)            
            p=ax.plot(fitx,fitn,'r',lw=2)
            ax.set_title('Tm=%s fit=%s' %(i,round(tm,1)))
            c+=1 
            errors.append(abs(i-tm))
            errors1.append(abs(i-tm1))

        legs=[]; legnames=['ideal','3% noise added']
        f1=plt.figure()
        ax=f1.add_subplot(111)
        p=ax.plot(tms,errors,lw=2); legs.append(p)
        ax.plot([55,55],[0,max(errors)], 'b', alpha=0.5,linestyle='--',lw=2)
        ax2 = plt.twinx()
        p=ax2.plot(tms,errors1,lw=2,color='green'); legs.append(p)
        ax.set_xlabel('True Tm')
        ax.set_ylabel('error (true-fit)')
        ax2.set_ylabel('error (noise added)')
        ax.set_title('normally distr noise added to test data')
        ax.legend(legs, legnames)
        
        f1.suptitle('Variation of error with true Tm (min pt at 60C)')
        f.savefig('tmfits.png')
        f1.savefig('tmrange.png')
        plt.show()
        return

    def testSN(self):
        """test signal-noise for absorbance->velocity measurements"""

        from mpl_toolkits.axes_grid import AxesGrid       
            
        f=plt.figure(figsize=(8,8))  
        DB=self.DB
        Et = DB['6 c9'].temp       
        names=[]
        ph='9'
        for i in ['60','65','70','75']:
            n='1_'+i+'_'+ph; names.append(n)
            fs = Et.estimateExpUncertainty(n, runs=50, guess=True)
            print n, fs['a'], fs['a'][1]/fs['a'][0]
        names.append('1_'+ph)        
        Et.plotDatasets(names,figure=f,legend=True,size=(6,6),showerrorbars=True,
                        plotoption=2)
        f.suptitle('')
        plt.show()
        #tm = round(numpy.mean(Et.getMetaData('1_'+ph)['t50'],Et.getMetaData('2_'+ph)['t50']),2)
        #tmerr = round(numpy.mean(Et.getMeta('1_'+ph, 'exp_errors')['t50'][1], E.getMeta('2_'+phc, 'exp_errors')['t50'][1]),2) 

        return        

    def fitarrhenius(self):                             
        x=range(270,400,5)
        y=[]
        A=1e10; Ea=1e-3;R=8.3144e-3
        for T in x:            
            val = A * math.exp(-Ea/R*T)            
            y.append(val)
            print T,val
        A,X=Fitting.doFit(expdata=zip(x,y),model='Arrhenius',silent=True)     
        fitx=numpy.arange(0,100)
        fity = X.getFitLine(fitx)
        f=plt.figure()
        ax=f.add_subplot(111)
        p=ax.plot(x,y,lw=2)
        plt.show()
        return
    
    def analyseMutations(self, sheet, filterby=None,
                            xcol='exp', ycol='Total'):
        """Analyse mutations"""
        from PEATDB.Ekin.Titration import TitrationAnalyser
        
        model = self.DB.getLabbookSheet(sheet)
        x=self.tofloats(model.getColumnData(columnName=xcol, filterby=filterby))
        mutations = model.getColumnData(columnName='Mutations',filterby=filterby)       
        print 'using %s points' %len(x)
        mutslist = [i[0].split('+') for i in zip(mutations,x)]
        #frequency of residue numbers
        
        freq = {}                
        import re        
        for muts in mutslist:            
            for m in muts:
                res = int(re.sub('\D','',m))
                if not freq.has_key(res):
                    freq[res]=0
                freq[res]+=1    
  
        #stats      
        mutlens = [len(i) for i in mutslist]
        singlemuts=0
        for i in mutlens:
            if i==1: singlemuts+=1
        rescodefreq={}
        for muts in mutslist:
            for m in muts:
                code = m[-1:]               
                if not rescodefreq.has_key(code): rescodefreq[code]=0
                rescodefreq[code]+=1                
  
        fr=freq.values()     
        mostfreqres=freq.keys()[fr.index(max(fr))]
        fr=rescodefreq.values()        
        mostmutated=rescodefreq.keys()[fr.index(max(fr))]
        #print stats
        print 'positions mutated: %s' %len(freq)
        print 'average no. of mutations per clone: %s' %numpy.mean(mutlens)
        print 'no. of single mutant clones: %s' %singlemuts
        print 'most frequently mutated residues: %s' %mostfreqres
        print 'most frequently inserted res type: %s' %TitrationAnalyser.residue_letters[mostmutated]        
        print 'no. of charge altering clones'
        
        #look at residues associated with high/low stab
       
        samples={}
        for i in zip(mutations,x):
            muts = i[0].split('+')
            stab=i[1]
            for m in muts:                
                res = int(re.sub('\D','',m))
                if not samples.has_key(res): samples[res]=[]
                samples[res].append(stab)
                
        meanstabs=[]      
        for s in samples:
            if len(samples[s]) > 2:
                l=samples[s]
                meanstabs.append((s,numpy.mean(l),numpy.std(l),len(l)))
        print meanstabs        
        f=plt.figure(figsize=(12,6))
        ax=f.add_subplot(111)
        ind=range(len(meanstabs))
        vals=[i[1] for i in meanstabs]
        names=[i[0] for i in meanstabs]
        stds=[i[2] for i in meanstabs]
        ax.bar(numpy.array(ind)+0.5,vals,yerr=stds)
        ax.set_title('Mean stab per residue')
        ax.set_xlabel('residue')
        ax.set_ylabel('mean deltaTm')
        ax.set_xticklabels(names)
        from matplotlib.ticker import MaxNLocator  
        ax.xaxis.set_major_locator(MaxNLocator(len(names)))         
        plt.show()
        
        self.plotfreq(freq)
        return freq

    def plotfreq(self, freq):
        #plot freq
        x=sorted(freq.keys())
        y=[freq[f] for f in x]
        f=plt.figure(figsize=(12,6))
        ax=f.add_subplot(111)
        ind=range(len(x))
        ax.bar(x,y)
        ax.set_xlabel('residue')
        ax.set_ylabel('freq')
        ax.set_title('Freq of each residue in all mutations - 10R')
        #ax.set_xticklabels(x)
        for l in ax.get_xticklabels():
            l.set_rotation(90)          
        from matplotlib.ticker import MaxNLocator  
        ax.xaxis.set_major_locator(MaxNLocator(12)) 
        yl = ax.get_ylim()
        legs=[]; legnames=[]
        def addrect(xy,w,h,clr,name):
            p = plt.Rectangle(xy,w,h,fc=clr,ec='none',alpha=0.5)        
            ax.add_patch(p)
            legs.append(p) 
            legnames.append(name)
        addrect((9,0),6,yl[1],'green','loop1')     
        addrect( (120,0),8,yl[1],'red','loop2')        
        addrect((160,0),10,yl[1],'yellow','loop3')       
        ax.legend(legs, legnames)
        f.savefig('mutsfreq.png')
        
    def classifyMutations(self, sheet, filterby=None):
        """Fit weights for positions using classifier"""
        model = self.DB.getLabbookSheet(sheet)
        x=self.tofloats(model.getColumnData(columnName='deltatm', filterby=filterby))
        mutations = model.getColumnData(columnName='Mutations',filterby=filterby)       
     
        from PEAT_SA.Core.Classify import StructureClassifier
        import PEAT_SA.Core.Matrix as Matrix
        SC = StructureClassifier()
        matrix = Matrix.matrixFromCSVFile('stabs.csv')
        w,s,o = SC.doRun(matrix, 'deltatm')
        #reconstruct stabilities using weights
        weights = {}
        for i in w:
            res,stab,num=i
            res=int(res[1:])
            weights[res] = stab
        y = []; lens=[]
        for i in zip(mutations,x):
            muts = i[0].split('+')
            val=0
            for m in muts:
                res = int(re.sub('\D','',m))               
                if res in weights:
                    val+=weights[res]                    
            y.append(val)
            lens.append(len(muts))
            print i[0], i[1], val
            
        f=plt.figure(figsize=(8,8))
        ax=f.add_subplot(111)
        ax.plot(y,x,'x',mew=2)
        ax.axhline(y=0, color='grey'); ax.axvline(x=0,color='grey')
        ax.set_xlim(-5,20); ax.set_ylim(-5,20)
        ax.set_ylabel('exp  stability')
        ax.set_xlabel('fitted')
        ax.set_title('Stabilities: exp vs linearly added from LARS fit')
        for i in range(len(lens)):            
            ax.annotate(lens[i], (y[i], x[i]+0.4), xytext=None, textcoords='data',
                             fontsize=8)        
        plt.show()
        return     

def main():
    """Run some analysis"""

    from optparse import OptionParser
    parser = OptionParser()
    app = CorrelationAnalyser()
    parser.add_option("-f", "--file", dest="file",
                        help="Open a local db")
    parser.add_option("-m", "--map structure", dest="mapstruct", action='store_true',
                       help="Map some properties on structure", default=False)
    parser.add_option("-c", "--classify", dest="classify", action='store_true',
                       help="Classify features using correl errors", default=False)
    parser.add_option("-t", "--tests", dest="tests", action='store_true',
                       help="Misc tests", default=False)  
    parser.add_option("-r", "--createstability10R", dest="createstability10R", action='store_true',
                        help=" -create temp data - 10R option", default=False)
    parser.add_option("-n", "--novo", dest="novo", action='store_true',
                       help="10R plots", default=False)      
    
    opts, remainder = parser.parse_args()

    if opts.file != None and os.path.exists(opts.file):
        app.loadDB(opts.file)
        print app.DB
    else:
        from PEATDB.Base import PDatabase
        app.DB = PDatabase(server='peat', username='farrell',
                             password='tafa', project='novo', port=8080)
        #print app.DB
    ph='7'
    if opts.createstability10R == True:
        app.extract10RStability(ph=ph)
        
    if opts.mapstruct == True:
        app.mapStructure(positions=[123,11])
      
    if opts.tests == True:
        #app.testSN()
        #app.testResActFit()
        app.fitarrhenius()
        
    #10R plot/analysis
    if opts.novo == True:
        '''ax = app.plotCorrelation(sheet='stability'+ph, filterby=('fit',['ok']), 
                              xcol='Total', ycol='ddg', labelcol='name',
                              labeloutliers=True,
                              plottitle='10R : PEAT-SA Predicted vs Experimental, ph '+ph,
                              xlabel='ddG Predicted (kJ/mol)',
                              ylabel='ddG Experimental',
                              #limit=40,
                              filename='10Rstability'+ph+'.png')
        
        app.analyseCorrelation(sheet='stability'+ph, filterby=('fit',['ok']),
                                xcol='Total', ycol='ddg', labelcol='Mutations')'''
        
        app.analyseMutations(sheet='stability'+ph, filterby=('fit',['ok','?']),
                                xcol='Total')
        app.classifyMutations(sheet='stability'+ph, filterby=('fit',['ok','?']))         
            
            
    #limit of valid exp vals 
    #ax.plot([-10,20],[16,16], 'b', alpha=0.5,linestyle='--')    
    #r = plt.Rectangle((-4, -4), 8, 8,fill=False,alpha=0.5,
    #                      edgecolor='gray',linestyle='dotted')        
    #ax.add_patch(r)    
 

if __name__ == '__main__':
    main()

