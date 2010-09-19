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
    from matplotlib.font_manager import FontProperties
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
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
        """From gui"""
        if self.filterbymenu.getcurselection() != '-':
            filt = (self.filterbymenu.getcurselection(), self.filtervalentry.getvalue())
        else:
            filt = None
        #get x-y vals from sheet selected
        sheet=self.smenu.getcurselection()
        xcol=self.xcolsmenu.getcurselection()
        ycol=self.ycolsmenu.getcurselection()
        model = self.DB.getLabbookSheet(sheet)
        x=model.getColumnData(columnName=xcol, filterby=filt)
        y=model.getColumnData(columnName=ycol, filterby=filt)
        if len(x)<1 or len(y)<1:
            print 'no data to plot'
            return
        x,y = zip(*self.tofloats(zip(x,y)))
        labels = model.getColumnData(columnName='Mutations', filterby=filt)
        self.plotCorrelation(x,y,labels)
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
        
    def plotCorrelation(self, x, y, labels=None,                       
                          labeloutliers=False,
                          plottitle='Predicted vs Experimental',
                          xlabel='Predicted',
                          ylabel='Experimental',
                          limit=None):
        """Show exp vs pred. for a set of x-y values """
        
        colors = ['b','g','r','y','m','c']        
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        legs=[]; legnames=[]
        bad=[]; good=[]
       
        if len(x) ==0:
            return
        errs = [a - b for a, b in zip(x, y)]
        line = ax.plot(x, y, 'x', color='red', mew=2, alpha=0.7, picker=3)  

        #draw expected correlation line with slope x
        slope=1
        #set range of axes
        if max(y)>max(x): lims = min(y)-5,max(y)+5
        else: lims = min(x)-5,max(x)+5
        cx = numpy.arange(lims[0],lims[1])
        cy = [slope*i for i in cx]
        ax.plot(cx, cy, 'g', alpha=0.7)
        ax.plot(cx-4, cy, 'g', alpha=0.5,linestyle='--')
        ax.plot(cx+4, cy, 'g', alpha=0.5,linestyle='--')
        ax.axhline(y=0, color='grey'); ax.axvline(x=0,color='grey')
        #ax.set_xlim(min(x)-2,max(x)+2)
        if limit==None:
            limit=lims[1]
        ax.set_xlim(lims[0],limit)
        ax.set_ylim(lims[0],limit)           
        ax.set_title('%s points' %len(x))
        cc = str(round(pow(stats.pearsonr(x,y)[0],2),2))
        ax.text(1,16, r'$r^2= %s$' %cc, fontsize=16)
        fig.suptitle(plottitle)
        from PEATDB.Actions import DBActions
        DBActions.showTkFigure(fig)        
        m = MouseHandler(ax, self, labels)
        m.connect()
        return ax

    def addMouseHandler(self, ax):
        """Add mouse event picker to plot so users can
        get info on each point, such as mutation name"""
        m = MouseHandler(ax, self)
        m.connect()    
        return
    
    def getStats(self,x,y):
        '''finalx=[x[i] for i in range(len(labels)) if labels[i] not in outliers]
        finaly=[y[i] for i in range(len(labels)) if labels[i] not in outliers]
        cc1 = str(round(pow(stats.pearsonr(finalx,finaly)[0],2),2))
        #ax.text(1,14, r'$r^2_{final}= %s$' %cc1, fontsize=16)'''
        return
        
    def markOutliers(self,x,y):
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
      
        self.plotCorrelation(x,y)

class MouseHandler:
    """Class to handle mouse click actions on plots"""

    bbox_props = dict(boxstyle="round", fc="#FFFC17", ec="0.4", alpha=0.8)
    
    def __init__(self, ax, parent, labels=None):
        self.ax = ax
        self.parent = parent
        self.press = False
        self.rect = None
        self.labels = labels
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
        obj = event.artist
        ind = event.ind
        xd = obj.get_xdata()[ind[0]]
        yd = obj.get_ydata()[ind[1]]            
        info = self.labels[ind[0]]
        print info
        self.infolabel = self.ax.annotate(info, (xd+0.4, yd), xytext=None, textcoords='data',
                                    fontsize=12, bbox=self.bbox_props)
        self.circle = plt.Circle((xd, yd), 0.2,fill=False)
        self.ax.add_patch(self.circle)
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

