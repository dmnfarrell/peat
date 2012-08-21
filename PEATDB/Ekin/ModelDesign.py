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
import csv, copy
import pickle
import numpy as np
from Tkinter import *
import tkFileDialog, tkSimpleDialog, tkMessageBox
import Pmw
from PEATDB.Ekin.Plotting import PlotPanel
from PEATDB.Ekin.Base import EkinProject, EkinDataset
import PEATDB.Ekin.Fitting as Fitting
import Ekin_images
from PEATDB.GUI_helper import *

class ModelDesignApp(Frame, GUI_help):
    """Simple GUI for designing new models for use in ekin fitting"""

    def __init__(self, parent=None):

        self.parent=parent
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
        else:
            self.main=Toplevel()
            #self.master=self.main
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
        self.modelsdict = Fitting.createModels()
        self.setupGUI()
        self.currentname = 'Linear'
        self.loadModel(self.currentname)
        return

    def setupGUI(self):
        """Do GUI elements"""
        self.createMenuBar()
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
        b=Button(f2,text='Load Models File',command=self.loadModelsFile)
        b.pack(side=TOP,fill=BOTH,pady=2)
        b=Button(f2,text='Save Current Models',command=self.save)
        b.pack(side=TOP,fill=BOTH,pady=2)
        self.filenamevar = StringVar()
        w=Label(f2, text='', textvariable=self.filenamevar, fg='blue')
        w.pack(side=TOP,fill=BOTH,pady=2)
        self.updateFileLabel()
        f3 = Frame(f2)
        f3.pack(fill=BOTH)
        self.modelselector = Pmw.ScrolledListBox(f3,
                items=self.modelsdict.keys(),
                labelpos='nw',
                label_text='Current Models:',
                listbox_height = 3,
                dblclickcommand=self.loadModel,
                listbox_selectmode = EXTENDED,
                usehullsize = 1,
                hull_width = 300,
                hull_height = 100,
        )
        self.modelselector.pack()
        b=Button(f3,text='New',command=self.newModel)
        b.pack(side=LEFT,fill=BOTH,pady=2)
        b=Button(f3,text='Delete',command=self.deleteModel)
        b.pack(side=LEFT,fill=BOTH,pady=2)
        b=Button(f3,text='Rename',command=self.renameModel)
        b.pack(side=LEFT,fill=BOTH,pady=2)
        self.currmodelvar = StringVar()
        Label(f3,text='Current Model:').pack(side=LEFT,fill=BOTH,pady=2)
        Label(f3,textvariable=self.currmodelvar,anchor=W,fg='blue').pack(side=LEFT,fill=BOTH,pady=2)
        b=Button(f2,text='Test Model',command=self.testFitter,bg='#FFFF99')
        b.pack(side=TOP,fill=BOTH,pady=2)
        b=Button(f2,text='Compare Models',command=self.findBestModel)
        b.pack(side=TOP,fill=BOTH,pady=2)
        m1.add(f1)
        self.previewer = FitPreviewer(m,app=self)
        m.add(self.previewer)
        return

    def createMenuBar(self):
        """Create the menu bar for the application"""
        self.menu=Menu(self.main)
        self.file_menu={'01Open models file':{'cmd':self.loadModelsFile},                        
                        '02Quit':{'cmd':self.quit}}
        self.file_menu=self.create_pulldown(self.menu,self.file_menu)
        self.menu.add_cascade(label='File',menu=self.file_menu['var'])
        self.help_menu={ '01Online Help':{'cmd': self.help},
                         '02About':{'cmd': self.about},}
        self.help_menu=self.create_pulldown(self.menu,self.help_menu)
        self.menu.add_cascade(label='Help',menu=self.help_menu['var'])
        self.main.config(menu=self.menu)
        return

    def updateFileLabel(self):
        """Update the file name"""
        if self.filename == None:
            self.filenamevar.set('no file loaded')
        else:
            self.filenamevar.set(self.filename)
        return

    def createGuessEntryWidget(self, parent):
        """Guess entry widget"""
        fr = Pmw.ScrolledFrame(parent,labelpos = 'n',
                                  label_text='guess')
        fr.pack(fill=BOTH,side=TOP,expand=1,padx=4)
        fr.configure(horizflex = 'expand')
        balloon = Pmw.Balloon(self.main)
        balloon.bind(fr, 'Models usually require initial guesses.\n'
                         'Enter a calculation or float value for any variable.')
        self.guessentry = fr
        return fr

    def newModel(self, name=None):
        """Create new model"""
        name = tkSimpleDialog.askstring('New model','Enter a name',
                                          initialvalue='modelname',
                                          parent=self.main)
        if not name or name == '' or name in self.modelsdict.keys():
            print 'model exists or empty'
            return
        model = { 'description': '',
                 'equation': 'a*x+b',
                 'guess': {},
                 'name': name,
                 'varnames': 'a,b'}
        self.modelsdict[name] = model
        self.updateModelSelector()
        self.modelselector.setvalue(name)
        self.loadModel()
        return

    def renameModel(self):
        """Rename a model"""
        currentname = self.modelselector.getcurselection()[0]
        name = tkSimpleDialog.askstring('Rename model','Enter a new name',
                                          initialvalue=currentname,
                                          parent=self.main)
        if not name or name == '' or name in self.modelsdict.keys():
            print 'model exists or empty'
            return
        self.modelsdict[name] = copy.deepcopy(self.modelsdict[currentname])
        del self.modelsdict[currentname]
        self.updateModelSelector()
        self.modelselector.setvalue(name)
        self.loadModel()
        return

    def deleteModel(self):
        """Delete a model"""
        names = currentname = self.modelselector.getcurselection()
        for name in names:
            del self.modelsdict[name]
        self.updateModelSelector()
        self.modelselector.setvalue(self.modelsdict.keys()[0])
        self.loadModel()
        return

    def updateModelsDict(self):
        """Send the current model to the dict"""
        name = self.currentname
        self.modelselector.setvalue(name)
        self.modelsdict[name] = self.currentmodel
        return

    def save(self):
        """Save current models"""
        if not tkMessageBox.askyesno('Overwrite Models?',
                                      'Are you sure you want to save?',
                                      parent=self.main):
            return
        if self.filename == None:
            self.filename = tkFileDialog.asksaveasfilename(defaultextension='.dict',
                                              initialdir=os.getcwd(),
                                              filetypes=[("dict files","*.dict"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        if not self.filename:
            return
        f=open(self.filename, 'w')
        self.updateModelsDict()
        pickle.dump(self.modelsdict, f)
        f.close()
        return

    def saveAs(self):
        """Save current models as a new file"""
        self.filename = None
        self.save()
        return

    def loadModelsFile(self, filename=None):
        """Load a dict of models, the pickled format from Ekin"""
        if filename==None:
            filename=tkFileDialog.askopenfilename(defaultextension='.dict',
                                              initialdir=os.getcwd(),
                                              filetypes=[("dict files","*.dict"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        if filename:
            self.modelsdict = pickle.load(open(filename,'r'))
            self.updateModelSelector()
            self.filename = filename
            self.updateFileLabel()
        return

    def updateModelSelector(self):
        modelnames = sorted(self.modelsdict.keys())
        self.modelselector.setlist(modelnames)
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
            if type(model['guess']) is types.StringType:
                model['guess'] = eval(model['guess'])
            if v in model['guess']:
                w.setvalue(model['guess'][v])
            w.pack(fill=BOTH,side=TOP,expand=1,pady=2)
        return

    def loadModel(self, sel=None):
        """Load a model from the current dict"""

        #add confirmation dialog here
        if sel == None:
            self.currentname = sel = self.modelselector.getcurselection()[0]
        self.currentmodel = model = self.modelsdict[sel]
        for f in model:
            if f in self.entrywidgets:
                self.entrywidgets[f].settext(model[f])
        self.updateGuessEntryWidgets(model)
        self.previewer.replot()
        self.currmodelvar.set(self.currentname)
        return

    def createFitter(self, model, name):
        """Create fitter from the current model entry values"""
        fitclass = Fitting.createClass(**model)
        x,y = self.previewer.getCurrentData()
        X = fitclass(exp_data=zip(x,y),variables=None,
                     callback=self.previewer.updateFit,
                     name=name)
        return X

    def testFitter(self):
        """Parse the entry fields and create the fitter"""
        self.currentmodel = self.parseValues()
        X = self.createFitter(self.currentmodel, self.currentname)
        X.guess_start()
        if self.previewer.plotfunctionvar.get() == 1:
            self.previewer.plotframe.clearData()
            self.previewer.plotframe.updateFit(X, showfitvars=True)
        else:
            X.fit(rounds=60)
            self.updateModelsDict()
            #update the global fitting models since the plotter uses them
            Fitting.updateModels(self.modelsdict)
            self.previewer.finishFit(X)
        return

    def findBestModel(self):
        """Determine best fit model using SS F-test"""
        models = self.modelselector.getcurselection()
        if len(models)<2:
            print 'select 2 or more models'
            return
        fitters= []
        for m in models:
            model = self.modelsdict[m]
            X = self.createFitter(model, m)
            fitters.append(X)
        Fitting.updateModels(self.modelsdict)
        self.previewer.findBestModel(models)
        return

    def parseValues(self):
        """Parse model entry and create a new fitter"""
        model = {}
        for f in self.entrywidgets:
            model[f] = str(self.entrywidgets[f].get().rstrip())
        model['guess']={}
        for v in self.guessentrywidgets:
            val = str(self.guessentrywidgets[v].get().rstrip())
            if val != '':
                model['guess'][v] = val
        return model

    def help(self):
        import webbrowser
        link='http://code.google.com/p/peat/wiki/ModelDesign'
        webbrowser.open(link,autoraise=1)
        return

    def about(self):
        win=Toplevel()
        win.geometry('+500+350')
        win.title('About ModelDesign')
        win.maxsize(width=400,height=400)
        '''logo = Images.logo()
        label = Label(win,image=logo)
        label.image = logo
        label.pack(fill=BOTH,padx=4,pady=4)'''
        text="""ModelDesign is an application that is used to create non-linear 
             fitting models for use in the Ekin application or elsewhere.              
             Released under GPL v3
             (C) Copyright 2012- Damien Farrell """
        text=text.replace('\t','')
        text= ' '.join(text.split())
        Label(win,text=text,wraplength=400).pack(fill=Y,side=TOP,pady=4)
        return

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
        self.E = None
        self.gui()
        self.sampleData()
        self.stop = False
        return

    def gui(self):
        fr = Frame(self)
        b=Button(fr,text='load csv',command=self.importCSV)
        b.pack(side=LEFT,fill=BOTH,padx=1)
        b=Button(fr,text='load ekin',command=self.loadProject)
        b.pack(side=LEFT,fill=BOTH,padx=1)
        self.stopbtn=Button(fr,text='stop',command=self.stopFit)#,state=DISABLED)
        self.stopbtn.pack(side=LEFT,fill=BOTH,padx=1)
        self.previmg = Ekin_images.prev()
        b=Button(fr,text='prev',image=self.previmg,command=self.prev)
        b.pack(side=LEFT,fill=BOTH)
        self.nextimg = Ekin_images.next()
        b=Button(fr,text='next',image=self.nextimg,command=self.next)
        b.pack(side=LEFT,fill=BOTH)
        fr.pack(side=TOP)
        self.dsetselector = Pmw.ComboBox(fr, entry_width=15,
                        selectioncommand = self.selectDataset)
        self.dsetselector.pack(side=LEFT,fill=BOTH,padx=1)
        fr2 = Frame(self)
        fr2.pack(side=TOP)
        self.plotfunctionvar = IntVar()
        cb=Checkbutton(fr2, text='plot function only',variable=self.plotfunctionvar)
        cb.pack(side=TOP)
        self.plotframe = PlotPanel(parent=self, side=BOTTOM, height=200)
        self.dsindex = 0
        self.opts = {'markersize':18,'fontsize':10,'showfitvars':True}
        return

    def sampleData(self):
        E =self.E = EkinProject()
        E.createSampleData()
        self.plotframe.setProject(E)
        self.datasets = sorted(self.E.datasets)
        self.replot()
        self.updateSelector()
        return

    def updateSelector(self):
        self.dsindex=0
        lst = self.datasets
        self.dsetselector.setlist(lst)
        self.dsetselector.selectitem(lst[self.dsindex])
        return

    def replot(self, dset=None):
        """Replot"""

        if dset==None:
            dset = self.datasets[self.dsindex]
        self.plotframe.plotCurrent(dset,options=self.opts)
        return

    def finishFit(self, X):
        """Call when fitting done"""
        d = self.datasets[self.dsindex]
        model = X.getResult()['model']
        fitdata = Fitting.makeFitData(model, X.variables, X.getError())
        self.E.setFitData(d, fitdata)
        self.replot()
        self.stop = False
        return

    def stopFit(self):
        self.stop = True
        return

    def getCurrentData(self):
        dset = self.datasets[self.dsindex]
        ek = self.E.getDataset(dset)
        return ek.getxy()

    def updateFit(self, selfdiff, vrs, fitvals, c, X):
        self.plotframe.updateFit(X, showfitvars=True)
        self.update_idletasks()
        self.update()
        if self.stop == True:
            X.stop_fit = 1
        return

    def prev(self):
        if self.dsindex <= 0:
            self.dsindex = 0
        else:
            self.dsindex -= 1
        self.replot()
        self.dsetselector.selectitem(self.datasets[self.dsindex])
        return

    def next(self):
        if self.dsindex >= self.E.length-1:
            self.dsindex = self.E.length-1
        else:
            self.dsindex += 1
        self.replot()
        self.dsetselector.selectitem(self.datasets[self.dsindex])
        return

    def selectDataset(self, name):
        self.dsindex = self.datasets.index(name)
        self.replot(name)
        return

    def importCSV(self):
        """Import csv file"""
        from PEATDB.Ekin.IO import Importer
        importer = Importer(self, parent_win=self)
        importer.path = os.getcwd()
        newdata = importer.import_multiple()
        for name in newdata.keys():
            self.E.insertDataset(newdata[name], newname=name)
        self.datasets = sorted(self.E.datasets)
        self.dsindex = self.datasets.index(name)
        self.replot()
        self.updateSelector()
        return

    def loadProject(self):
        filename=tkFileDialog.askopenfilename(defaultextension='.ekinprj',
                                              initialdir=os.getcwd(),
                                              filetypes=[("Ekin files","*.ekinprj"),
                                                         ("All files","*.*")],
                                              parent=self.main)
        if not filename: return
        self.E.openProject(filename)
        self.datasets = sorted(self.E.datasets)
        self.replot()
        self.updateSelector()
        return

    def findBestModel(self, models, callback=None):
        """Determine best fit model using SS F-test"""
        dset = self.datasets[self.dsindex]
        ek = self.E.getDataset(dset)
        fitdata,p = Fitting.findBestModel(ek, models=models, silent=True)
        best = fitdata['model']
        text = 'Best model is %s with p value=%2.2e' %(best,p)
        tkMessageBox.showinfo('Best fit model result',text)
        return

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="file",
                        help="Models File", metavar="FILE")

    opts, remainder = parser.parse_args()
    app = ModelDesignApp()
    if opts.file != None:
        app.loadModelsFile(opts.file)
    app.mainloop()

if __name__ == '__main__':
    main()
