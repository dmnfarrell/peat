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

def main():
    """Run some analysis"""

    from optparse import OptionParser
    parser = OptionParser()
    app = CorrelationAnalyser()
    parser.add_option("-f", "--file", dest="file",
                        help="Open a local db")
    parser.add_option("-t", "--tests", dest="tests", action='store_true',
                       help="Misc tests", default=False)  
   
    
    opts, remainder = parser.parse_args()

    if opts.file != None and os.path.exists(opts.file):
        app.loadDB(opts.file)   

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
     
 

if __name__ == '__main__':
    main()

