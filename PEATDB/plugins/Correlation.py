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
import types, os, sys, math, string, re
import math, numpy
#from scipy.stats import stats
try:
    import matplotlib
    from matplotlib.font_manager import FontProperties
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
except:
    pass
from PEATDB.TableModels import TableModel
from PEATDB.Ekin.Base import EkinDataset
import PEATDB.Ekin.Fitting as Fitting

class CorrelationAnalyser(Plugin):
    """A class for processing specific Kinetics exp data"""
    """Author: Damien Farrell"""

    capabilities = ['gui','uses_sidepane']
    requires = ['pylab','scipy']
    menuentry = 'Correlation Analysis'
    gui_methods = { 'doCorrelation':'Plot Correlation'}
           
    about = 'This plugin is designed for some correlation analysis and for 10R project'

    def __init__(self):
        self.phlist = ['6','6.5','7','7.5','8','8.5','9','10']
        return

    def main(self, parent=None):
        self.parent = parent
        if parent!=None:             
            self.DB = parent.DB
        if self.DB == None:
            print 'plugin requires DB'
            return
        self._doFrame()
        return

    def _doFrame(self):
        if self.parent == None:
            self.mainwin=Toplevel()
            self.mainwin.title('Correlation Analysis')
            self.mainwin.geometry('800x600+200+100')
        elif 'uses_sidepane' in self.capabilities:
            self.mainwin = self.parent.createChildFrame()

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
        """Call plot from gui"""
        if self.filterbymenu.getcurselection() != '-':
            filt = (self.filterbymenu.getcurselection(), self.filtervalentry.getvalue())
        else:
            filt = None
        #get x-y vals from sheet selected
        sheet=self.smenu.getcurselection()
        xcol=self.xcolsmenu.getcurselection()
        ycol=self.ycolsmenu.getcurselection()
        
        model = self.DB.getLabbookSheet(sheet)
        x,y,names,muts = model.getColumns([xcol,ycol,'name','Mutations'],
                            filterby=filt, allowempty=False) 
        if len(x)<1 or len(y)<1:
            print 'no data to plot'
            return
        labels = zip(names, muts)
        self.plotCorrelation(x,y,labels,side=TOP)
        return        
     
    def showHistogram(self, lists, labels=[],title='distr',
                      style=0,ax=None):
        if ax==None:
            f=plt.figure()
            ax=f.add_subplot(111)
        if style==0:
            n, bins, p = ax.hist(lists, 10, label=labels, histtype='bar',
                                        alpha=0.6,lw=0.2)            
        else:         
            bins=10   
            for l in lists:
                n, bins, p = ax.hist(l,bins=bins,label=labels,
                                     alpha=0.6,lw=0.2)
        ax.legend()
        ax.set_title(title)
        return ax
        
    def plotCorrelation(self, x, y, labels=None,
                          key='Mutations',
                          title='',
                          xlabel='Predicted',
                          ylabel='Experimental', ms=5,
                          err=None, axeslabels=True,
                          plotname=None, stats=True,        
                          side=LEFT, ax=None):
        """Show exp vs pred. for a set of x-y values """
               
        #check if x and y are number cols and round
        x=[round(float(i),2) for i in x]
        y=[round(float(i),2) for i in y]

        if min(x)<min(y): a=min(x)-2
        else: a=min(y)-4
        if max(x)>max(y): b=max(x)+2
        else: b=max(y)+4
       
        colors = ['b','g','r','y','m','c']        
        if ax==None:
            fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(111)
        else:       
            #fig = ax.get_figure()
            fig = None
        if axeslabels==True:    
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        legs=[]; legnames=[]
        bad=[]; good=[]
       
        if len(x) ==0:
            return
        errs = [i - j for i, j in zip(x, y)]
        line = ax.plot(x, y, 'o', mew=0, ms=ms,
                       alpha=0.6, label=plotname,
                       picker=3)  

        #draw expected correlation line with slope x
        slope=1
        #set range of axes
        ax.plot((a,b),(a,b),'g')            
        ax.axhline(y=0, color='grey'); ax.axvline(x=0,color='grey')
        ax.set_xlim(a,b); ax.set_ylim(a,b)
        ax.set_title(title)
           
        if err!=None:
            ax.plot((a,b),(a+err,b+err),'--',color='g')
            ax.plot((a,b),(a-err,b-err),'--',color='g')         
        if stats==True:
            self.addStats(ax,x,y)
        if fig!=None:
            fig.suptitle('Predicted vs Experimental')            
            from PEATDB.Actions import DBActions
            frame = DBActions.showTkFigure(fig, side=side)            
            #mh = MouseHandler(ax, self, labels, key)
            #mh.connect()
            mh = self.addMouseHandler(ax, labels, key)
        else:
            mh = frame = None
        return ax, frame, mh

    def addMouseHandler(self, ax, labels, key):
        """Add mouse event picker to plot so users can
        get info on each point, such as mutation name"""
        m = MouseHandler(ax, self, labels, key)
        m.connect()
        return m

    def addStats(self, ax, x, y):
        cc = math.pow(numpy.corrcoef(x,y)[0][1],2)
        rmse = self.rmse(x,y)
        ax.text(0.1,0.9,'$r^2=%s$' %round(cc,3),
            transform=ax.transAxes,fontsize=18)
        ax.text(0.1,0.85,'$rmse=%s$' %round(rmse,3),
            transform=ax.transAxes,fontsize=18)        
        return
    
    def rmse(self, ar1, ar2):
        """Mean squared error"""
        ar1 = numpy.asarray(ar1)
        ar2 = numpy.asarray(ar2)
        dif = ar1 - ar2
        dif *= dif
        return numpy.sqrt(dif.sum()/len(ar1))

    def getStats(self,x,y):
        x=[round(float(i),2) for i in x]
        y=[round(float(i),2) for i in y]
        errs = [i[0]-i[1] for i in zip(x,y)]
        cc = round(numpy.corrcoef(x,y)[0][1],3)        
        rmse = round(self.rmse(x,y),3)
        meanerr = round(numpy.mean(errs),3)
        return cc, rmse, meanerr

    def pearsonr(self,x,y):
        """pearson correl coeff"""
        n=len(x)
        vals=range(n)         
        #regular sums
        sumx=sum([float(x[i]) for i in vals])
        sumy=sum([float(y[i]) for i in vals])         
        #sum of the squares
        sumxSq=sum([x[i]**2.0 for i in vals])
        sumySq=sum([y[i]**2.0 for i in vals])         
        #sum of the products
        pSum=sum([x[i]*y[i] for i in vals])         
        #do pearson score
        num=pSum-(sumx*sumy/n)
        den=((sumxSq-pow(sumx,2)/n)*(sumySq-pow(sumy,2))**.5)
        if den==0: 
            return 1
        r=num/den
        return r
  
    def markOutliers(self,x,y):
        labels = model.getColumnData(columnName=labelcol,filterby=filterby)      
        
        #outliers = self.findOutliers(x,y,errs)
        '''for l in labels:
            i=labels.index(l)
            #label and list outliers
            if l in outliers:
                #c = plt.Circle((x[i], y[i]), 0.4,fill=False,alpha=0.7)
                ax.annotate(l, (x[i]+0.4, y[i]), xytext=None, textcoords='data',
                                    fontsize=8)
                #ax.add_patch(c)'''        

     
    def analyseCorrelation(self, sheet, filterby=None,
                            xcol='exp', ycol='Total', labelcol=None):
        """Analyse correlation"""
        model = self.DB.getLabbookSheet(sheet)
        x=model.getColumnData(columnName=xcol, filterby=filterby)
        y=model.getColumnData(columnName=ycol, filterby=filterby)
        x,y = zip(*self.tofloats(zip(x,y)))
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
        print 'no. of points:', len(x)
        print 'correlation co-eff:', numpy.corrcoef(x,y)[0][1]
        print 'stdev of errors:', numpy.std(errs)
        return   

    def plotNorm(self, data, title='', lw=3, ax=None):
        """Plot a normal distr given some data"""
        if ax==None:
            f=plt.figure()
            ax=f.add_subplot(111)       
        mean = numpy.mean(data)
        std = numpy.std(data)            
        ax.hist(data,bins=20,normed=1,histtype='step',lw=lw)
        #draw the normal dist    
        x=numpy.arange(mean-std*4,mean+std*4,max(data)/100)
        y=[]               
        for i in x:
            fx = 1/(math.sqrt(2*math.pi*math.pow(std,2))) * math.exp( -math.pow(i-mean,2)/(2*math.pow(std,2)))
            y.append(fx)
        ax.plot(x,y,'-',lw=lw,alpha=0.7)
        ax.set_title(title)
        return ax

    def QQplot(self, data=None, labels=None, title=None,n=3, ax=None):
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
        if ax==None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        #ax.set_xlabel('normal distr')
        #ax.set_ylabel('data')
        ax.plot(x, y, 'o', ms=2, alpha=0.7)
        ax.plot(x, x, '-', ms=2, alpha=0.7)
        if labels!=None:
            markoutliers()
        if title!=None:
            ax.set_title(title)
        #fig.savefig('qq.png')
        return outliers

    @classmethod
    def ShapiroWilk(self, a1):
        """Test departure from normality"""
        import scipy
        x1=numpy.mean(a1); n1=len(a1); s1=numpy.std(a1)
        w,p = scipy.stats.shapiro(a1)
        print 'W=%s p=%s' %(w,p)
        return round(w,3),round(p,4)

    @classmethod
    def tofloats(self, lsts):
        """Return floats of a single list or tuple of paired values"""
        def evaluate(l):
            for i in l:
                try: float(i)
                except:
                    return False
            return True        
        vals = [i for i in lsts if evaluate(i) == True]
        result=[]
        for l in vals:
            v=[]
            for i in l:
                v.append(float(i))
                result.append(v)
        return result             

    def test(self):
        x=[5,6,8,9,'t']; y=[6,'',3,7,8]
        labels=['dsad','sas','fef','xx']
        xy = self.tofloats(zip(x,y))
        x,y = zip(*xy)      

class MouseHandler:
    """Class to handle mouse click actions on plots"""

    bbox_props = dict(boxstyle="round", fc="#FFFC17", ec="0.4", alpha=0.8)
    
    def __init__(self, ax, parent, labels=None, key=None):
        self.ax = ax
        self.parent = parent
        self.press = False
        self.rect = None
        self.labels = labels
        self.key = key
        self.event = None
        self.infolabel = None
        self.circle = None
        return

    def connect(self):        
        self.cidpress = self.ax.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.ax.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.ax.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)
        self.cidpick = self.ax.figure.canvas.mpl_connect('pick_event', self.on_pick)        

    def on_pick(self, event): 
        """Handle pick event, picking on points"""

        obj = event.artist
        ind = event.ind       
        xd = obj.get_xdata()[ind[0]]
        yd = obj.get_ydata()[ind[0]]
        info = self.labels[ind[0]]        
        if type(info) is types.TupleType and not None in info:
            labels = '\n'.join(info)
        elif type(info) is types.DictType:
            pass
        labels = labels+'\n'+'pred:'+str(xd)+' exp:'+str(yd)+'\n'+'err='+str(xd-yd)
        self.infolabel = self.ax.annotate(labels, (xd+0.4, yd), xytext=None, textcoords='data',
                                    fontsize=12, bbox=self.bbox_props)
        self.circle = plt.Circle((xd, yd), 0.2,fill=False)
        self.ax.add_patch(self.circle)
        self.ax.figure.canvas.draw()
        if hasattr(self,'table'):
            ml=self.table.model
            if self.key != None:
                name = info[0]                           
                recname = ml.filterBy('name', name, op='=')[0]                
                self.table.movetoSelectedRow(recname=recname)
        return
    
    def on_press(self, event):
        """Handle mouse click"""
        self.press = True
        self.x = event.xdata; self.y = event.ydata
        self.parent.selection = None        
        if self.rect != None:
            self.rect.set_visible(False)
            #self.ax.patches = []
            if self.rect in self.ax.patches:
                self.ax.patches.remove(self.rect)
            self.ax.figure.canvas.draw()
            self.rect=None        
        return

    def on_motion(self, event):
        """Draw a selection rect"""
        if self.press == False:
            return
        x = event.xdata; y = event.ydata
        if x == None or y==None:
            return
        
        dx = x-self.x; dy=y-self.y
        if self.rect == None:
            #print 'new'
            self.rect = Rectangle((self.x,self.y),dx,dy, fc='lightblue',ec='blue',alpha=0.6)
            self.ax.add_patch(self.rect)
        else:
            self.rect.set_width(dx)
            self.rect.set_height(dy)

        #draw selection rect
        self.ax.figure.canvas.draw()
        self.parent.selection = (self.x, self.y, x, y)
        if x!=None: self.prevx=x
        if y!=None: self.prevy=y
        return

    def on_release(self, event):
        self.press = False
        if self.infolabel in self.ax.texts:
            self.ax.texts.remove(self.infolabel)
        if self.circle in self.ax.patches:
            self.ax.patches.remove(self.circle)
        return

    def disconnect(self):
        """Disconnect all the stored connection ids"""
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
        self.rect.figure.canvas.mpl_disconnect(self.cidmotion)
        self.rect.figure.canvas.mpl_disconnect(self.cidpick)
        
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

    if opts.tests == True:
        app.main()
        app.mainwin.mainloop()
        app.test()    
 
if __name__ == '__main__':
    main()

