#!/usr/bin/env python
#
# DataPipeline - A data import and fitting tool
# Copyright (C) 2011 Damien Farrell
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
# Email: damien.farrell_at_ucd.ie
# Normal mail:
# Damien Farrell
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
#

import os, sys, math, random, numpy, string, types
import csv
import pickle
import numpy as np
from Tkinter import *
import tkFileDialog, tkSimpleDialog
import Pmw
from PEATDB.Ekin.Ekin_main import EkinApp, PlotPanel
from PEATDB.Ekin.Base import EkinProject, EkinDataset
from PEATDB.Ekin.Fitting import Fitting

class ModelDesignApp(Frame):
    """Simple GUI for designing new models for use in ekin fitting"""

    def __init__(self, parent=None, modelsfile=None):

        self.parent=parent
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
        else:
            self.main=Toplevel()
            self.master=self.main
        self.main.title('Model Design')
        ws = self.main.winfo_screenwidth()
        hs = self.main.winfo_screenheight()
        w = 900; h=600
        x = (ws/2)-(w/2); y = (hs/2)-(h/2)
        self.main.geometry('%dx%d+%d+%d' % (w,h,x,y))
        self.main.protocol('WM_DELETE_WINDOW',self.quit)
        self.tableformat = {'cellwidth':50, 'thefont':"Arial 10",
                            'rowheight':16, 'editable':False,
                            'rowselectedcolor':'yellow'}
        self.filename = None
        self.createModels()
        self.setupGUI()
        return

    def setupGUI(self):
        """Do GUI elements"""
        m = PanedWindow(self.main,
                           orient=HORIZONTAL,
                           sashwidth=3,
                           showhandle=True)
        m.pack(side=TOP,fill=BOTH,expand=1)
        m1 = PanedWindow(m,
                           orient=VERTICAL,
                           sashwidth=3,
                           showhandle=True)
        m.add(m1)
        f1=Frame(m1,width=200)

        self.entrywidgets = {}
        self.guessentrywidgets = {}
        fields = ['equation','varnames','description']
        heights = [100,50,80]
        for f in fields:
            i=fields.index(f)
            self.entrywidgets[f] = Pmw.ScrolledText(f1,
                labelpos = 'n',
                label_text=f,
                hull_width=400,
                hull_height=heights[i],
                usehullsize = 1)
            self.entrywidgets[f].pack(fill=BOTH,side=TOP,expand=1,padx=4)

        self.createGuessEntryWidget(f1)

        f2=Frame(m1)
        m1.add(f2)
        b=Button(f2,text='New Model',command=self.newModel)
        b.pack(side=TOP,fill=BOTH,pady=2)
        b=Button(f2,text='Load Models File',command=self.loadModelsFile)
        b.pack(side=TOP,fill=BOTH,pady=2)
        b=Button(f2,text='Save Current Models',command=self.save)
        b.pack(side=TOP,fill=BOTH,pady=2)

        self.modelselector = Pmw.ScrolledListBox(f2,
                items=self.modelsdict.keys(),
                labelpos='nw',
                label_text='Current Models:',
                listbox_height = 3,
                selectioncommand=self.loadModel,
                #dblclickcommand=self.loadModel,
                usehullsize = 1,
                hull_width = 300,
                hull_height = 100,
        )
        self.modelselector.pack()
        b=Button(f2,text='Test Model',command=self.testFitter,bg='#FFFF99')
        b.pack(side=TOP,fill=BOTH,pady=2)
        m1.add(f1)
        self.previewer = FitPreviewer(m,app=self)
        m.add(self.previewer)
        return

    def createGuessEntryWidget(self, parent):
        fr = Pmw.ScrolledFrame(parent,labelpos = 'n',
                                  label_text='guess')
        fr.pack(fill=BOTH,side=TOP,expand=1,padx=4)
        fr.configure(horizflex = 'expand')
        balloon = Pmw.Balloon(self.main)
        balloon.bind(fr, 'Models usually require initial guesses.\n'
                         'Enter a calculation or float value for any variable.')
        self.guessentry = fr
        return fr

    def createModels(self):
        """Create default models file"""
        self.modelsdict = {
            'Linear': {'Name': 'Linear',
             'description': 'straight line',
             'equation': 'a*x+b',
             'guess': {'a': 'min(y)', 'b': 'max(y)-min(y)/max(x)-min(x)'},
             'name': 'Linear',
             'varnames': 'a,b'},
            'Sigmoid': {'Name': 'Sigmoid',
             'description': 'simple sigmoid',
             'equation': 'bottom+(top-bottom)/(1+exp((tm-x)/slope))',
             'guess': {'slope': '1', 'top': 'max(y)', 'tm': '(max(x)-min(x))/2+min(x)', 'bottom': 'min(y)'},
             'name': 'Sigmoid',
             'varnames': 'tm,bottom,top,slope'},
            'Gaussian': {'Name': 'Gaussian',
             'description': 'gaussian function, a bell-shaped curve',
             'equation': 'a*exp(-(pow((x-b),2)/(pow(2*c,2))))',
             'guess': {},
             'name': 'Gaussian',
             'varnames': 'a,b,c'}
        }
        return

    def newModel(self):
        """Create new model"""
        name = tkSimpleDialog.askstring('New model','Enter a name',
                                          initialvalue='',
                                          parent=self.main)
        if not name: return
        model = { 'description': '',
             'equation': 'a*x+b',
             'guess': "",
             'name': name,
             'varnames': 'a,b'}
        self.modelsdict[name] = model
        self.updateModelSelector()
        self.modelselector.setvalue(name)
        self.loadModel()
        return

    def save(self):
        """Save current models"""
        if self.filename == None:
            self.filename = tkFileDialog.asksaveasfilename(defaultextension='.dict',
                                              initialdir=os.getcwd(),
                                              filetypes=[("dict files","*.dict"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        if not self.filename:
            return
        f=open(self.filename, 'w')
        pickle.dump(self.modelsdict, f)
        f.close()
        return

    def saveAs(self):
        """Save current models as a new file"""
        self.filename = None
        self.save()
        return

    def loadModelsFile(self):
        """Load a dict of models, the pickled format from Ekin"""
        filename=tkFileDialog.askopenfilename(defaultextension='.dict',
                                              initialdir=os.getcwd(),
                                              filetypes=[("dict files","*.dict"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        if filename:
            self.modelsdict = pickle.load(open(filename,'r'))
            self.updateModelSelector()
            self.filename = filename
        return

    def updateModelSelector(self):
        self.modelselector.setlist(self.modelsdict.keys())
        return

    def updateGuessEntryWidgets(self, model):
        """Update guess entry widget, special case as guess vals are in
        a dict"""
        for v in self.guessentrywidgets.keys()[:]:
            w = self.guessentrywidgets[v]
            del self.guessentrywidgets[v]
            w.forget()

        if not model.has_key('guess') or not type(model) is types.DictType:
            model['guess'] = {}
        fr = self.guessentry.interior()
        varnames = model['varnames'].split(',')
        for v in varnames:
            w = self.guessentrywidgets[v] = Pmw.EntryField(fr,
                labelpos = 'w',
                label_text = v,
                hull_width=400)
            if v in model['guess']:
                w.setvalue(model['guess'][v])
            w.pack(fill=BOTH,side=TOP,expand=1,pady=2)
        return

    def loadModel(self):
        """Load a model from the current dict"""

        #add confirmation dialog here

        sel = self.modelselector.getcurselection()[0]
        self.currentmodel = model = self.modelsdict[sel]
        for f in model:
            if f in self.entrywidgets:
                self.entrywidgets[f].settext(model[f])
        self.updateGuessEntryWidgets(model)
        self.previewer.replot()
        return

    def createFitter(self, model):
        """Create fitter from the current model entry values"""
        fitclass = Fitting.createClass(**model)
        x,y = self.previewer.getCurrentData()
        X = fitclass(exp_data=zip(x,y),variables=None,
                     callback=self.previewer.updateFit)
        return X

    def testFitter(self):
        """Parse the entry fields and create the fitter"""
        self.currentmodel = self.parseValues()
        X = self.createFitter(self.currentmodel)
        X.guess_start()
        X.fit(rounds=60)
        return

    def parseValues(self):
        """Parse model entry and create a new fitter"""
        model = {}
        for f in self.entrywidgets:
            model[f] = self.entrywidgets[f].get().rstrip()
        return model

    def quit(self):
        self.main.destroy()
        if not self.parent:
            sys.exit()
        return

class FitPreviewer(Frame):
    """Class to preview model fitting in plot"""
    def __init__(self, master, app=None):

        Frame.__init__(self, master)
        self.main=self.master
        self.app = app
        self.gui()
        self.sampleData()
        return

    def gui(self):
        fr = Frame(self)
        b=Button(fr,text='load csv',command=self.importCSV)
        b.pack(side=LEFT,fill=BOTH)
        b=Button(fr,text='load ekin prj',command=self.loadProject)
        b.pack(side=LEFT,fill=BOTH)
        b=Button(fr,text='update',command=self.replot)
        b.pack(side=LEFT,fill=BOTH)
        b=Button(fr,text='prev',command=self.prev)
        b.pack(side=LEFT,fill=BOTH)
        b=Button(fr,text='next',command=self.next)
        b.pack(side=LEFT,fill=BOTH)
        fr.pack(side=TOP)
        self.plotframe = PlotPanel(parent=self, side=BOTTOM, height=200)
        self.dsindex = 0
        self.opts = {'markersize':18,'fontsize':10}
        return

    def sampleData(self):
        E=self.E = EkinProject()
        x=range(20)
        y=[i+np.random.normal(.2) for i in x]
        ek = EkinDataset(xy=[x,y])
        self.E.insertDataset(ek,'testdata')
        self.plotframe.setProject(E)
        self.replot()

    def replot(self):
        dset = self.E.datasets[self.dsindex]
        self.plotframe.plotCurrent(dset,options=self.opts)

    def getCurrentData(self):
        dset = self.E.datasets[self.dsindex]
        ek = self.E.getDataset(dset)
        return ek.getxy()

    def updateFit(self, selfdiff, vrs, fitvals, c, X):
        self.plotframe.updateFit(X)
        #if self.stopfit == True:
        #    X.stop_fit=1

    def prev(self):
        if self.dsindex <= 0:
            self.dsindex = 0
        else:
            self.dsindex -= 1
        self.replot()
        return

    def next(self):
        if self.dsindex >= self.E.length-1:
            self.dsindex = self.E.length-1
        else:
            self.dsindex += 1
        self.replot()
        return

    def importCSV(self):
        """Import csv file"""
        from PEATDB.Ekin.IO import Importer
        importer = Importer(self, parent_win=self)
        importer.path = os.getcwd()
        newdata = importer.import_multiple()
        for name in newdata.keys():
            self.E.insertDataset(newdata[name], newname=name)
        self.dsindex = self.E.datasets.index(name)
        self.replot()
        return

    def loadProject(self):
        filename=tkFileDialog.askopenfilename(defaultextension='.ekinprj',
                                              initialdir=os.getcwd(),
                                              filetypes=[("Ekin files","*.ekinprj"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        if not filename: return
        self.E.openProject(filename)
        self.replot()
        return

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="file",
                        help="Raw file", metavar="FILE")

    opts, remainder = parser.parse_args()
    app = ModelDesignApp(modelsfile=opts.file)
    if opts.file != None:
        app.openRaw(opts.file)
    app.mainloop()

if __name__ == '__main__':
    main()
