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

try:
    from Plugins import Plugin
except:
    from PEATDB.Plugins import Plugin
from Tkinter import *
import Pmw
import os, sys, math, string
import csv
import math, numpy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from PEATDB.Ekin.Base import EkinProject,EkinDataset
from PEATDB.Ekin.Fitting import Fitting

class VantHoff(Plugin):
    """A plugin to do Van't Hoff Analysis of temperature melting curves"""
    """Author: Damien Farrell"""

    capabilities = ['gui','uses_sidepane']
    requires = ['pylab','numpy']
    menuentry = "Van't Hoff Analysis"

    gui_methods = {'getCSV': 'Import CSV',
                    'loadEkin':'Load Ekin Proj',
                    'doAnalysis':"Do Analysis",
                    'benchmark': 'Do Benchmark',                  
                    'close':'Close' }
    about = "A plugin to do Van't Hoff Analysis of temperature melting curves"
    R = 8.3144
    
    def __init__(self):
        
        return

    def main(self, parent):
        if parent==None:
            return
        self.parent = parent
        self.DB = parent.DB
        if self.DB == None:
            return
        self.xydata = None
        self.E = None
        self._doFrame()
        return

    def _doFrame(self):
        if 'uses_sidepane' in self.capabilities:
            self.mainwin = self.parent.createChildFrame()
        else:
            self.mainwin=Toplevel()
            self.mainwin.title(self.menuentry)
            self.mainwin.geometry('800x600+200+100')
         
        methods = self._getmethods()        
        methods = [m for m in methods if m[0] in self.gui_methods.keys()]
        self._createButtons(methods)
        self.showDatasetSelector()
        self.doall = Pmw.RadioSelect(self.mainwin,
                buttontype = 'checkbutton',
                orient = 'horizontal',
                labelpos = 'w')
        self.doall.add('Process All')
        self.doall.pack()
        self.methods = Pmw.RadioSelect(self.mainwin,
                buttontype = 'checkbutton',
                orient = 'horizontal',
                labelpos = 'w',               
                label_text = 'Methods:')        
        for m in ['method 1','method 2','method 3']:
            self.methods.add(m)
        self.methods.invoke('method 1')    
        self.methods.pack() 
        self.sm = Pmw.EntryField(self.mainwin,
                labelpos = 'w',
                value = 5,
                label_text = 'Smoothing:')        
        self.sm.pack()
        self.tw = Pmw.EntryField(self.mainwin,
                labelpos = 'w',
                value = 60,
                label_text = 'Width of transition:')
        self.tw.pack()       
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
    
    def close(self):
        self.mainwin.destroy()
        self.plotframe = None
        return
        
    def showDatasetSelector(self):
        if self.E==None:
            return
        if hasattr(self, 'dmenu'):
            self.dmenu.destroy()
        self.dmenu = Pmw.OptionMenu(self.mainwin,
                labelpos = 'w',
                label_text = 'Dataset:',
                items = sorted(self.E.datasets),
                command=self.showPreview,
                menubutton_width = 8)
        self.dmenu.pack(side=TOP,fill=BOTH)        
        return
        
    def showPreview(self,event=None):
        if self.E == None:
            return
        if not hasattr(self, 'plotframe') or self.plotframe == None:               
            from Ekin.Ekin_main import PlotPanel
            self.plotframe = PlotPanel(parent=self.mainwin, side=BOTTOM)         
        self.plotframe.setProject(self.E)
        d = self.dmenu.getcurselection()
        self.plotframe.plotCurrent(d)
        plt.close(1)
        return
        
    def getCSV(self):
        """Import a csv file"""        
        self.E = EkinProject()
        from PEATDB.Ekin.IO import Importer
        importer = Importer(self,parent_win=self.mainwin)
        newdata = importer.import_multiple()
        if newdata == None: return
        for n in newdata.keys():
            self.E.insertDataset(newdata[n], n, update=None)
        print 'imported %s datasets' %len(self.E.datasets)            
        self.showDatasetSelector()
        self.showPreview()
        return

    def loadEkin(self):
        """Load the ekin prj"""
        import tkFileDialog
        filename=tkFileDialog.askopenfilename(defaultextension='.ekinprj',
                                                  initialdir=os.getcwd(),
                                                  filetypes=[("ekinprj","*.ekinprj"),
                                                             ("All files","*.*")],
                                                  parent=self.mainwin)
        if not os.path.isfile(filename):
            return         
        self.E = EkinProject()
        self.E.openProject(filename)
        self.showDatasetSelector()
        self.showPreview()
        return
        
    def doAnalysis(self):
        """Execute from GUI"""
        if self.E == None:            
            return        
        methods = self.methods.getcurselection()
        if  'Process All' in self.doall.getcurselection():        
            self.doAll(methods=methods)
        else:            
            if 'method 1' in methods:
                self.fitVantHoff(E=self.E,d=self.dmenu.getcurselection(),
                        transwidth=int(self.tw.getvalue()))
            if 'method 2' in methods:              
                self.fitDifferentialCurve(E=self.E,d=self.dmenu.getcurselection(),
                                            smooth=int(self.sm.getvalue()))
            if 'method 3' in methods:              
                self.fitElwellSchellman(E=self.E,d=self.dmenu.getcurselection(),
                                            transwidth=int(self.tw.getvalue()))
                
        return


    def transformCD(self,x,y,transwidth=None,ax=None):
        """Transform raw data into fraction unfolded per temp value, by fitting to
            a general unfolding equation that extracts baseline/slopes"""
        #fit baseline slopes and get intercepts
        A,X=Fitting.doFit(expdata=zip(x,y),model='Unfolding',conv=1e-7,noiter=100,silent=True)
        fity = X.getFitLine(x)        
        fd=X.getFitDict()
        if ax!=None:
            p=ax.plot(x,fity,'r',lw=2)
            self.drawParams(ax,fd)
        #we then use slopes and intercepts get frac unfolded at each temp        
        mn = fd['bn']; mu = fd['bd'] #slopes
        #if mu>0.01: mu = 0.01
        yn = fd['an']; yu = fd['ad'] #intercepts
        d50 = fd['d50']; m = fd['m']
        
        #try to take useful transition region of data 
        if transwidth != None:
            for i in x:
                if i>d50: 
                    mid = x.index(i)
                    break
            L=int(mid-transwidth); U=int(mid+transwidth)        
            x,y = x[L:U], y[L:U]
        
        t=[]; f=[]
        #print mu, mn
        for T,yo in zip(x,y):            
            fu = (yo-(yn+mn*T)) / ((yu+mu*T)-(yn+mn*T))
            #print fu, (yo-(yn+mn*T)), (m), mu, mn
            #if f>0:
            f.append(fu)
            t.append(T)            
        return t,f
    
    def fitVantHoff(self, E=None, d=None, xy=None, transwidth=80, invert=False,
                        show=True, figname=None):
        """Derive fraction unfolded, get K and fit to Van't Hoff.
           see http://www.jbc.org/content/277/43/40717.full
           or http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2144003/
        """
        
        if E !=None:
            ek=EkinDataset(E.getDataset(d))       
            x,y = ek.getXYActive()
        elif xy!=None:
            x,y = xy            
        if invert == True:
            y = [max(y)-i for i in y[:]] 
            
        f=plt.figure(figsize=(10,6))
        ax=f.add_subplot(121)
        p=ax.plot(x,y,'o',alpha=0.6)
        ax.set_xlabel('T')

        x,y = self.transformCD(x,y,transwidth,ax)
        
        #derive lnK vs 1/T
        t=[]; k=[]
        for T,fu in zip(x,y):
            if fu>=1 or fu<=0:
                continue 
            K = fu/(1-fu)
            klog = math.log(K)
            k.append(klog)
            t.append(1/T)
        if len(t)<2: return None, None, None
        ax1=f.add_subplot(122)
        p=ax1.plot(t,k,'x',color='black')
        ax1.set_xlabel('1/T')#(r'$1/T ($K^-1)$') 
        ax1.set_ylabel('ln K') 
        
        formatter = matplotlib.ticker.ScalarFormatter()
        formatter.set_scientific(True) 
        formatter.set_powerlimits((0,0))
        ax1.xaxis.set_major_formatter(formatter)
        for l in ax1.get_xticklabels():
            l.set_rotation(30)
            
        #fit this van't hoff plot
        A,X=Fitting.doFit(expdata=zip(t,k),model='Linear')
        fitk = X.getFitLine(t)      
        p=ax1.plot(t,fitk,'r',lw=2)
        fd=X.getFitDict()
        self.drawParams(ax1,fd)
        
        #slope is deltaH/R/1000 in kJ
        deltaH = -fd['a']*self.R/1000
        deltaS = fd['b']*self.R/1000
        f.suptitle("Method 1 - deltaH: %2.2f deltaS: %2.2f" %(deltaH,deltaS),size=18)
        if show==True:            
            self.showTkFigure(f)
            
        if figname != None:
            figname = figname.replace('.','_')
            f.savefig(figname)
            plt.close()
        if E!=None:          
            fdata = Fitting.makeFitData(X.name,vrs=X.variables)
            E.insertDataset(xydata=[t,k], newname=d+'_vanthoff',replace=True,fit=fdata)
            E.saveProject() 
        return deltaH, deltaS, ax

    def fitElwellSchellman(self,E=None, d=None, xy=None,transwidth=50,
                                invert=False,show=True,figname=None):
        """Fit entire raw data simultaneously to the three main thermodynamic
           parameters using Elwell/Schellman method"""
        if E !=None:
            ek=EkinDataset(E.getDataset(d))       
            x,y = ek.getXYActive()            
        elif xy!=None:
            x,y = xy
        else:
            return
        if invert == True:
            y = [max(y)-i for i in y[:]]
        f=plt.figure(figsize=(10,6))
        ax=f.add_subplot(121)
        ax.set_xlabel('T')
        p=ax.plot(x,y,'o',alpha=0.5)

        x,y = self.transformCD(x,y,transwidth,ax)
        
        t=[];dg=[]
        R=8.3144e-3
        for T,fu in zip(x,y):
            if fu>=1 or fu<=0:
                continue            
            K = fu/(1-fu)            
            deltaGt = -R * T * math.log(K)
            dg.append(deltaGt)
            t.append(T)

        ax1=f.add_subplot(122)
        p=ax1.plot(t,dg,'x',color='black')
        ax1.set_xlabel('T')
        ax1.set_ylabel('dG(T)')
        
        A,X=Fitting.doFit(expdata=zip(t,dg),model='elwellschellman',grad=1e-9,conv=1e-10)        
        fity = X.getFitLine(t)
        p=ax1.plot(t,fity,'r',lw=2)
        fd=X.getFitDict()
        self.drawParams(ax1,fd)
        deltaH=fd['deltaH']; deltacp=fd['deltacp']; Tm=fd['Tm']
        f.suptitle("Method 2 - deltaH: %2.2f deltaCp: %2.2f Tm: %2.2f" %(deltaH,deltacp,Tm),size=18)
        if show == True:
            self.showTkFigure(f)
        if figname != None:
            figname = figname.replace('.','_')
            f.savefig(figname)
            plt.close()            
        if E!=None:          
            fdata = Fitting.makeFitData(X.name,vrs=X.variables)
            E.insertDataset(xydata=[t,dg], newname=d+'_vanthoff2',replace=True,fit=fdata)
            E.saveProject()
        return deltaH, Tm, deltacp
    
    def breslauerMethod(self,E=None, d=None, xy=None,invert=False,
                        show=True,figname=None):
        """Finds slope of trans region and plugs this in to equation
        http://www.springerlink.com/content/r34n0201g30563u7/  """
        if E !=None:
            ek=EkinDataset(E.getDataset(d))       
            x,y = ek.getXYActive()            
        elif xy!=None:
            x,y = xy
        else:
            return
        f=plt.figure(figsize=(10,6))
        ax=f.add_subplot(111)
        ax.set_xlabel('T')
        p=ax.plot(x,y,'o',alpha=0.5)       
        A,X=Fitting.doFit(expdata=zip(x,y),model='Unfolding',conv=1e-7,noiter=40) 
        fity = X.getFitLine(x)      
        p=ax.plot(x,fity,'r',lw=2)
        fd=X.getFitDict()
        self.drawParams(ax,fd)   
        Tm = fd['d50']; m = fd['m']                 
        R = 8.3144e-3
        deltaH =  R * math.pow(Tm,2) * m
        f.suptitle("Method 4 - deltaH: %2.2f Tm: %2.2f" %(deltaH,Tm),size=18)
        if show == True:
            self.showTkFigure(f)           
        if figname != None:
            figname = figname.replace('.','_')
            f.savefig(figname)
            plt.close()
        return deltaH, Tm
        
    def fitDifferentialCurve(self, E=None, d=None, xy=None,smooth=0,
                                invert=False,show=True,figname=None):
        """Derive differential denaturation curve and fit to get deltaH
           We smooth the unfolding curve and then differentiate and finally
           fit to a 3 parameter equation.
           See http://www.ncbi.nlm.nih.gov/pubmed/10933511"""
      
        if E !=None:
            ek=EkinDataset(E.getDataset(d))       
            x,y = ek.getXYActive()
        elif xy!=None:
            x,y = xy
        else:
            return   
        if invert == True:
            y = [max(y)-i for i in y[:]] 

        leg=[]; lines=[]
        f=plt.figure(figsize=(10,6))
        ax=f.add_subplot(121)
        p=ax.plot(x,y,'x',color='black',mew=3,alpha=0.5)
        leg.append(p); lines.append('original')
        #smooth
        if smooth == 0:
            smooth=int(len(x)/15.0)
        s=self.smoothListGaussian(y,smooth)    
        p=ax.plot(x[:len(s)-1],s[:-1],lw=3)
        leg.append(p); lines.append('smoothed')
        ax.set_title("original data")
        ax.set_xlabel('T')
        ax1=f.add_subplot(122)
        #differentiate
        dx,ds = self.differentiate(x[:len(s)],s)
        #ds = [i/max(ds) for i in ds]
        ds = [i*10 for i in ds]
        cw=csv.writer(open('diffcd.csv','w'))
        for row in zip(dx,ds):
            cw.writerow(row)
        p=ax1.plot(dx,ds,'-',lw=1.5,alpha=0.7,color='black')
        leg.append(p); lines.append('differential')
        ax1.set_title("differential denaturation")
        ax1.set_xlabel('T'); ax1.set_ylabel('dsignal/dT')
        
        A,X=Fitting.doFit(expdata=zip(dx,ds),model='diffDenaturation',grad=1e-9,conv=1e-10)        
        fity = X.getFitLine(dx)
        p=ax1.plot(dx,fity,'r',lw=2)
        leg.append(p); lines.append('fit')
        t=X.getFitDict()
        self.drawParams(ax1,t)
        dHkcal=t['deltaH']/4.184
        f.suptitle('Method 3 - deltaH: %2.2f kJ/mol (%2.2f kcal) Tm: %2.2f' %(t['deltaH'],dHkcal,t['Tm']),size=18)
        ax.legend(leg,lines,loc='best',prop=FontProperties(size="smaller"))
        #f.subplots_adjust(hspace=0.8)
        if show == True:
            self.showTkFigure(f)
        if figname != None:
            figname = figname.replace('.','_')
            f.savefig(figname)
            plt.close()
        if E!=None:          
            fdata = Fitting.makeFitData(X.name,vrs=X.variables)
            E.insertDataset(xydata=[dx,ds], newname=d+'_diff',replace=True,fit=fdata)
            E.saveProject()
        return t['deltaH'],t['Tm']
    
    def differentiate(self, x,y):
        dy = numpy.diff(y,1)
        dx = x[:len(dy)]
        return dx,dy

    def smoothListGaussian(self,data,degree=5):
        """Gaussian data smoothing function"""
        #buffer data to avoid offset result
        data=list(data)        
        data = [data[0]]*(degree-1) + data + [data[-1]]*degree
        window=degree*2-1  
        weight=numpy.array([1.0]*window)  
        weightGauss=[]  
        for i in range(window):  
            i=i-degree+1  
            frac=i/float(window)  
            gauss=1/(numpy.exp((4*(frac))**2))
            weightGauss.append(gauss)  
        weight=numpy.array(weightGauss)*weight  
        smoothed=[0.0]*(len(data)-window)
        for i in range(len(smoothed)):  
            smoothed[i]=sum(numpy.array(data[i:i+window])*weight)/sum(weight)  
        return smoothed

    def invert(self,data):
        inv=[i for i in data]
        return inv
        
    def simulateCD(self,noise=1.0):
        """Simulate some CD spec data"""
        x=list(numpy.arange(290,380,0.2)); y=[]
        X=Fitting.getFitter(model='Unfolding',
                              vrs=[-16, 0.01, -11.6, 0.01, 2.7, 324])
        fity = X.getFitLine(x)
        for i in fity:
            noise=numpy.random.normal(i, 1.0/2)
            y.append(i+noise)
        cw=csv.writer(open('cd.csv','w'))
        for row in zip(x,y):
            cw.writerow(row)
        return x,y
       
    def drawParams(self,ax,d):
        ymin, ymax = ax.get_ylim()
        xmin, xmax = ax.get_xlim()
        inc=(ymax-ymin)/20
        xinc=(xmax-xmin)/20
        y=ymax-inc
        for k in d:            
            s = k+'='+str(round(d[k],3))
            ax.text(xmin+xinc,y,s,fontsize=10)
            y-=inc            
        return

    def usetex(self):
        plt.rc('text', usetex=True)
        return
   
    def doAll(self, methods=['method 1']):
        """Process all datasets in ekinprj""" 
        E=self.E
        vals={}
        from Dialogs import PEATDialog
        pb=PEATDialog(self.mainwin, option='progressbar',
                                      message='Analysing Data..')
        pb.update_progress(0)
        total = len(E.datasets); count=0
        for d in E.datasets: 
            if '_diff' in d or '_vanthoff' in d:
                continue
            vals[d]={}
            name = d
            if 'method 1' in methods:
                vals[d]['dH1'], vals[d]['dS1'], ax = self.fitVantHoff(E,d,show=False,figname=name)
            if 'method 2' in methods:
                vals[d]['dH2'], vals[d]['dTm2'], vals[d]['dCp2'] = self.fitElwellSchellman(E,d,show=False,figname=name)                                
            if 'method 3' in methods:
                vals[d]['dH3'], vals[d]['dTm3'] = self.fitDifferentialCurve(E,d,show=False,figname=name)
            count += 1
            pb.update_progress(float(count)/total*100.0)           
        pb.close()
        self.showTable(vals)
        return
    
    def showTable(self, data):
        """Show results in table"""
        from PEATDB.DictEdit import DictEditor       
        D=DictEditor(self.mainwin)
        D.loadTable(data)
        return
                
    def benchmark(self,E=None,d=None, method=1):
        """Test methods with varying paramaters, smoothing etc"""
        if E==None and self.E != None:
            E = self.E; d=self.dmenu.getcurselection()

        path='vh_benchmark'
        if not os.path.exists(path):
            os.mkdir(path)        
        dHvals=[]
        
        if method == 1:
            xlabel = 'width (K)'
            title = 'method 1: deltaH variation with trans region width fit'            
            vals=range(5,140,5)
            for w in vals:            
                dH, dS, ax = self.fitVantHoff(E,d,transwidth=w,show=False,
                                              figname=os.path.join(path,'%s_%s.png' %(d,w)))
                if dH == None: dH=0
                dHvals.append(dH)
            #take best values from middle    
            #dHvals= dHvals[5:16] 
        elif method == 2:
            xlabel = 'width (K)'
            title = 'method 2: deltaH variation with width fit'
            vals=range(5,140,5)
            for w in vals:
                dH, dcp, dTm = self.fitElwellSchellman(E,d,transwidth=w,show=False,
                                                       figname=os.path.join(path,'%s_%s.png' %(d,w)))
                dHvals.append(dH)
        elif method == 3:
            xlabel = 'smoothing degree'
            title = 'method 3: deltaH variation with degree of smoothing'
            vals=range(1,30,3)
            for s in vals:            
                dH, dTm = self.fitDifferentialCurve(E,d,smooth=s,show=False,
                                                    figname=os.path.join(path,'%s_%s.png' %(d,s)))
                dHvals.append(dH)                
        mean = numpy.mean(dHvals)
        stdev = numpy.std(dHvals)
        f=plt.figure()
        ax=f.add_subplot(111)
        ax.plot(vals, dHvals,lw=2)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('deltaH (kJ)')
        ax.set_title('mean: %2.2f stdev: %2.2f'%(mean, stdev))
        f.suptitle(title)        
        f.savefig('benchmark_%s.png' %method)
        cw=csv.writer(open('benchmark_%s.csv' %method,'w'))
        for row in zip(vals,dHvals):
            cw.writerow(row)
        return

    def benchmarkLimitedData(self, E=None,d=None, method=1):
        """test any method with varying limited data"""
        if E==None and self.E != None:
            E = self.E; d=self.dmenu.getcurselection()

        path='vh_benchmark'
        if not os.path.exists(path):
            os.mkdir(path)        
        dHvals=[]
        vals=[]
        if method == 1:
            L=range(5,140,5)
            for w in vals:
                dH, dS, ax = self.fitVantHoff(E,d,transwidth=w,show=False,
                                              figname=os.path.join(path,'%s_%s.png' %(d,w)))
        return    
    
    @classmethod 
    def plotCorrelation(self,x=None,y=None,xlabel='method1',ylabel='method2'):
        if x==None:
            data=open('compared.csv','r')
            cr=csv.reader(data)
            x=[float(r[0]) for r in cr]; data.seek(0) 
            y=[float(r[1]) for r in cr]
        f=plt.figure()
        ax=f.add_subplot(111)
     
        line = ax.scatter(x, y, marker='o',alpha=0.8)
        cl = numpy.arange(0,max(x)+50)
        ax.plot(cl, cl, 'g', alpha=0.5,lw=2)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(150,600); ax.set_ylim(150,600)        
        ax.set_title('Correlation')
        from scipy.stats import stats
        cc = str(round(pow(stats.pearsonr(x,y)[0],2),2))
        ax.text(400,180, r'$r^2= %s$' %cc, fontsize=16)
        plt.rc('text', usetex=True)        
        self.showTkFigure(f)
        return

    def showTkFigure(self, fig):
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
        fr = Toplevel()
        canvas = FigureCanvasTkAgg(fig, master=fr)
        #self.canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        mtoolbar = NavigationToolbar2TkAgg( canvas, fr )
        mtoolbar.update()
        
        return

def main():
    """Run some analysis"""

    from optparse import OptionParser
    parser = OptionParser()
    app = VantHoff()
    parser.add_option("-f", "--file", dest="file",
                        help="Open a local db")
    parser.add_option("-e", "--ekinprj", dest="ekinprj",
                        help="Open an ekin project")
    parser.add_option("-d", "--dataset", dest="dataset",
                        help="Dataset name")    
    parser.add_option("-m", "--method", dest="method", default=1, type='int',
                        help="Choose method - 1: Van't Hoff plot, 2: Differential fit")      
    parser.add_option("-b", "--benchmark", dest="benchmark", action='store_true',
                       help="Test", default=False)
    parser.add_option("-a", "--all", dest="all", action='store_true',
                       help="Do all datasets in ekinprj", default=False)
    parser.add_option("-w", "--width", dest="width", default=50, type='int',
                       help="Width of transition region to fit for method 1")
    parser.add_option("-s", "--smoothing", dest="smoothing", default=5, type='int',
                       help="Degree of smoothing to apply in method 2 (default 5)")
    parser.add_option("-i", "--invert", dest="invert", action='store_true',
                       help="Invert raw data", default=False)
    
    opts, remainder = parser.parse_args()
   
    if opts.file != None and os.path.exists(opts.file):
        app.loadDB(opts.file)        
    if opts.ekinprj != None and os.path.exists(opts.ekinprj):
        E = EkinProject()
        E.openProject(opts.ekinprj)
        d = opts.dataset
    else:
        x,y = app.simulateCD()
        E = EkinProject()
        d='cdtest'
        E.insertDataset(xydata=[x,y], newname=d)
    if opts.all == True:
        self.doAll(E, methods)             
    if opts.benchmark == True:
        app.benchmark(E,d,method=opts.method)           

        #app.plotCorrelation()
    else:        
        if opts.method == 1:
            app.fitVantHoff(E,d,transwidth=opts.width,invert=opts.invert)
        elif opts.method == 2:
            app.fitElwellSchellman(E,d,transwidth=opts.width,invert=opts.invert)
        elif opts.method == 3:
            app.fitDifferentialCurve(E,d,smooth=opts.smoothing,invert=opts.invert)
        elif opts.method == 4:
            app.breslauerMethod(E,d,invert=opts.invert)                       
            
if __name__ == '__main__':
    main()

