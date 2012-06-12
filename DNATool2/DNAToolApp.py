#!/usr/bin/env python
#
# DNATool - A program for DNA sequence manipulation
# Copyright (C) 2012- Damien Farrell & Jens Erik Nielsen
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
# Email: farrell.damien_at_gmail.com


"""Main GUI application for DNATool"""

import os, sys, math, random, string, types
import pickle
from Tkinter import *
import Pmw
import tkFileDialog
import Dialogs
from PEATDB.GUI_helper import *
from Base import Project, Sequence, SequenceCanvas, ORFOverview
from Primer import PrimerDatabase, PrimerDBGUI, PrimerDesignGUI
from SeqIO import SequenceAnalysis, SequenceAlignmentTool, BlastInterface
import RestrictionDigest
import Images
from Prefs import Preferences
from Editor import TextEditor

class App(Frame, GUI_help):
    """Data pipe GUI for importing and fitting of raw data.
       This class uses ekin provided an automated pipeline for fitting
       raw text data files and propagating errors.
       Uses a config file to store the pipeline settings"""

    def __init__(self, parent=None, rawfile=None, conffile=None):

        self.parent=parent
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
        else:
            self.main=Toplevel()
            self.master=self.main
        self.main.title('DNATool Desktop')
        self.main.protocol('WM_DELETE_WINDOW',self.quit)

        self.defaultprefs = {'rememberwindowsize':0,'orientation':'horizontal',
                             'alignmenttool':'clustal','clustalpath':'clustalw'}
        self.defaultopts = [('rememberwindowsize','checkbutton',1,'Remember window size'),
                            ('orientation','menu',('horizontal','vertical'),
                            'Default sideframe orientation'),
                            ('alignmenttool','menu',('clustal','muscle','needle'),
                            'External alignment tool'),
                            ('clustalpath','entry','','Path to Clustal')]
        self.loadPrefs()
        self.setGeometry()
        self.P = Project()
        self.primerdb = PrimerDatabase()
        self.setupApps()
        self.setupGUI()
        self.sidepane.selectpage('PrimerDB')
        return

    def setGeometry(self):
        if self.rememberwindowsize == 1 and self.preferences.has_key('lastwindowsize'):
            lastwindowsize = self.preferences.get('lastwindowsize')
            self.winsize = lastwindowsize
            #self.w = int(lastwindowsize.split('x')[0])
        else:
            self.winsize = self.getBestGeometry()
        self.main.geometry(self.winsize)

        return

    def getBestGeometry(self):
        """Calculate optimal geometry from screen size"""
        ws = self.main.winfo_screenwidth()
        hs = self.main.winfo_screenheight()
        self.w=w=ws/1.3; h=500
        x = (ws/2)-(w/2); y = (hs/2)-(h/2)
        g = '%dx%d+%d+%d' % (w,h,x,y)
        return g

    def saveGeometry(self):
        """Save current geometry before quitting"""
        g = self.main.geometry()
        self.preferences.set('lastwindowsize', g)
        return

    def setupApps(self):
        """Creates a dict of classes so that GUI frames like primer design
           can be dynamically added without explicit methods for each"""
        self.apps = {'Primer Design': PrimerDesignGUI,
                     'ORF Overview': ORFOverview,
                     'Sequence Analysis': SequenceAnalysis,
                     'Sequence Alignment Tool': SequenceAlignmentTool,
                     'Restriction Digest Summary': RestrictionDigestSummary,
                     'Blast Interface':BlastInterface,
                     'NotePad': TextEditor}
        return

    def loadPrefs(self):
        """Load prefs from a preferences instance"""
        self.preferences = Preferences('DNATool2',{})
        temp = {}
        for key in self.defaultprefs:
            if not key in self.preferences.prefs:
                self.preferences.set(key, self.defaultprefs[key])
            temp[key] = self.preferences.get(key)
        self.applyPrefs(temp)
        return

    def savePrefs(self):
        """Save prefs to a preferences instance"""
        self.applyPrefs(self.prefsdialog.getValues())
        self.preferences = Preferences('DNATool2',{})
        for key in self.defaultprefs:
            self.preferences.set(key, self.__dict__[key])
        return

    def applyPrefs(self, prefs=None):
        """Apply prefs from a given dict"""
        if prefs == None: prefs=self.defaultprefs
        for key in prefs:
            self.__dict__[key] = prefs[key]
        return

    def setupGUI(self):
        """Do GUI elements"""
        self.visibleapps = {}
        #bottom panel
        bottom = Frame(self.main,height=50)
        bottom.pack(side=BOTTOM,fill=BOTH,pady=1)
        #status bar
        sb = self.addStatusBar(bottom)
        tb = self.addWindowToolBar(bottom)
        tb.pack(side=LEFT,anchor='e',expand=1)

        self.m = PanedWindow(self.main,
                           orient=self.orientation,
                           sashwidth=2,
                           showhandle=True)
        self.m.pack(side=TOP,fill=BOTH,pady=2,expand=1)

        sc = self.sc = SequenceCanvas(self.m,parentapp=self,
                                      width=900,height=800)
        sc.pack(side=TOP,fill=BOTH,pady=2)
        self.m.add(sc)
        self.m.paneconfigure(sc,sticky='news',min=200)

        self.createMenuBar()
        self.createToolBar()
        self.createSidePane()
        if self.orientation == 'vertical':
            self.m.paneconfigure(self.sidepane,
                                      min=200)
        self.createPrimerDBGUI()
        self.createChildFrame(name='NotePad')
        return

    def createMenuBar(self):
        """Create the menu bar for the application. """
        self.menu=Menu(self.main)
        self.file_menu={'01Open Project':{'cmd':self.openProject},
                        '02Open Sequence':{'cmd':self.openSequence},
                        '05Quit':{'cmd':self.quit}}
        self.file_menu=self.create_pulldown(self.menu,self.file_menu)
        self.menu.add_cascade(label='File',menu=self.file_menu['var'])
        self.edit_menu={'01Undo':{'cmd':self.undo},
                        '02Copy':{'cmd':self.copy},
                        '03Select All':{'cmd':self.sc.selectAll},
                        '04Configure Restriction Enzymes':{'cmd':self.restrictionEnzymesDialog}}
        self.edit_menu=self.create_pulldown(self.menu,self.edit_menu)
        self.menu.add_cascade(label='Edit',menu=self.edit_menu['var'])

        self.primer_menu={'01Primer DB':{'cmd':self.createPrimerDBGUI}}
        self.primer_menu=self.create_pulldown(self.menu,self.primer_menu)
        self.menu.add_cascade(label='Primer Tools',menu=self.primer_menu['var'])

        self.seqanal_menu={'01x':{'cmd':self.openSequence}}
        self.seqanal_menu=self.create_pulldown(self.menu,self.seqanal_menu)
        self.menu.add_cascade(label='Sequence Analysis',menu=self.seqanal_menu['var'])

        self.view_menu=Menu(self.menu)
        self.menu.add_cascade(label='Tools',menu=self.view_menu)
        self.appsvars = {}
        for i in self.apps.keys():
            self.appsvars[i] = IntVar()
            def func(args):
                def new():
                    self.toggleApps(args)
                return new
            self.view_menu.add_checkbutton(label=i, onvalue=True,
                                            offvalue=False,
                                            command=func(i),
                                            variable=self.appsvars[i])

        self.help_menu={ '01Online Help':{'cmd': self.help},
                         '02About':{'cmd': self.about},}
        self.help_menu=self.create_pulldown(self.menu,self.help_menu)
        self.menu.add_cascade(label='Help',menu=self.help_menu['var'])
        self.main.config(menu=self.menu)
        return

    def createSidePane(self, width=200):
        """Side panel for various dialogs is tabbed notebook that can hold multiple
            dialogs at once """
        self.closeSidePane()
        self.sidepane = Frame(self.m,
                              bd=1,relief=RAISED)
        self.sidepane = Pmw.NoteBook(self.m)
        self.sidepane.component('hull').configure(relief='sunken',
                                                    borderwidth=1)
        self.m.paneconfigure(self.sidepane,
                                  width=width,
                                  sticky='news')
        self.m.add(self.sidepane)
        self.sidepane.setnaturalsize()
        return self.sidepane

    def closeSidePane(self):
        """Destroy sidepine"""
        if hasattr(self, 'sidepane'):
            self.m.forget(self.sidepane)
            self.sidepane.destroy()
        return

    def createChildFrame(self, sidepane=True, width=400, name=None):
        """Create a child frame in sidepane notebook or as toplevel"""
        geometry = '600x400+500+250'
        if sidepane == True:
            cframe = self.sidepane.add(name)
            self.sidepane.tab(name).configure(font='fixed 8',anchor='w')
            self.sidepane.selectpage(name)
        else:
            cframe = Toplevel()
            if name == None:
                title = 'default'
            else:
                title = name
            cframe.title(title)
        self.__dict__[name] = cframe
        #if name is a pre-defined app we load fetch the class
        if self.apps.has_key(name):
            cls = self.apps[name]
            inst = cls(cframe, parentapp=self)
            inst.pack(fill=BOTH,expand=1)
            self.appsvars[name].set(1)
            if hasattr(inst, 'geometry'):
                geometry = inst.geometry
        if sidepane == False:
            cframe.geometry(geometry)
        #bind frame close
        def func(evt):
            if hasattr(self, name):
                del self.__dict__[name]
            if name in self.appsvars:
                self.appsvars[name].set(0)
        cframe.bind("<Destroy>", func)
        return cframe

    def closeChildFrame(self):
        """Close sidepane frame"""
        name = self.sidepane.getcurselection()
        self.sidepane.delete(name)
        return

    def detachChildFrame(self):
        """Since we cannot actually detach frames in Tkinter, we
           remove the frame in the sidepane make a new toplevel window"""
        name = self.sidepane.getcurselection()
        self.sidepane.delete(name)
        if name == 'PrimerDB':
            fr = self.createPrimerDBGUI(sidepane=False)
        else:
            fr = self.createChildFrame(sidepane=False, name=name)
        return fr

    def toggleApps(self, name):
        """Remove or add apps from view"""
        print name
        if self.appsvars[name].get() == 1:
            if not hasattr(self, name):
                fr = self.createChildFrame(sidepane=True, name=name)
        else:
            if hasattr(self, name):
                if name in self.sidepane.pagenames():
                    self.sidepane.delete(name)
                else:
                    self.__dict__[name].destroy()
                    del self.__dict__[name]
        return

    def createToolBar(self):
        """Toolbar"""
        self.toolbar = ToolBar(self.main, self)
        self.toolbar.pack(side=TOP,fill=X,pady=1,before=self.m)
        img = Images.openproject()
        hlptxt="Open Project from file"
        self.toolbar.addButton('Open Project', self.openProject, img, hlptxt)
        img = Images.saveproject()
        hlptxt="Save Project"
        self.toolbar.addButton('Save Project', self.openProject, img, hlptxt)
        img = Images.undo()
        hlptxt="Undo"
        self.toolbar.addButton('Save Project', self.openProject, img, hlptxt)
        img = Images.zoomout()
        hlptxt="Zoom Out"
        self.toolbar.addButton('Save Project', self.sc.zoomOut, img, hlptxt)
        img = Images.zoomin()
        hlptxt="Zoom In"
        self.toolbar.addButton('Save Project', self.sc.zoomIn, img, hlptxt)
        img = Images.windowprefs()
        hlptxt="Seqeuence display preferences"
        self.toolbar.addButton('Seqeuence display preferences', self.sc.showPrefs,
                         img, hlptxt)
        img = Images.prefs()
        hlptxt="Preferences"
        self.toolbar.addButton('Preferences', self.globalPrefsDialog, img, hlptxt)
        return

    def addWindowToolBar(self, parent):
        """Toolbar for changing frame views"""
        toolbar = ToolBar(parent, self)
        img = Images.closeframe()
        hlptxt="Close current sideframe"
        toolbar.addButton('tile vertical', self.closeChildFrame, img, hlptxt)
        img = Images.detachframe()
        hlptxt="Detach current sideframe as seperate window"
        toolbar.addButton('tile vertical', self.detachChildFrame, img, hlptxt)
        img = Images.tilevertical()
        hlptxt="Tile frames vertically"
        toolbar.addButton('x', self.tileVertical, img, hlptxt)
        img = Images.tilehorizontal()
        hlptxt="Tile frames horizontally"
        toolbar.addButton('tile horizontal', self.tileHorizontal, img, hlptxt)
        return toolbar

    def addStatusBar(self, parent):
        fr = Frame(parent, width=100)
        Label(fr,text='status stuff').pack()
        fr.pack(fill=X,side=LEFT,anchor='w')
        return

    def removeAll(self):
        """Remove all widgets"""
        for w in self.main.children.values():
            w.destroy()
        return

    def tileHorizontal(self):
        """Re-arrange paned widgets horizontally"""
        if self.orientation == 'horizontal':
            return
        self.removeAll()
        self.orientation = 'horizontal'
        self.setupGUI()
        return

    def tileVertical(self):
        """Re-arrange paned widgets vertically"""
        if self.orientation == 'vertical':
            return
        self.removeAll()
        self.orientation = 'vertical'
        self.setupGUI()
        return

    def openSequence(self):
        """Load a sequence"""
        seq = Sequence()
        self.sc.loadSequence(seq)
        return

    def openProject(self, filename=None):
        """Load a project"""
        if filename == None:
            filename = Dialogs.openFilename(parent=self, ext='dtp')
        if not filename:
            return
        self.P.load(filename)
        #open sequence
        seq = self.P.DNAseq
        #print self.P.used_enzymes, self.P.cut_pos
        self.sc.project = self.P
        self.sc.update()

        #update primer db
        self.primerdb = PrimerDatabase(self.P.data['primer_dict'])
        self.pdbwin.database = self.primerdb
        self.pdbwin.update()
        return

    def copy(self):
        self.sc.copySequence()
        return

    def undo(self):
        return

    def createPrimerDBGUI(self, sidepane=True):
        """We create this individually to avoid confusion"""
        if hasattr(self, 'PrimerDB'):
            self.PrimerDB.lift()
            return
        fr = self.createChildFrame(name='PrimerDB',sidepane=sidepane)
        self.pdbwin = PrimerDBGUI(fr, database=self.primerdb)
        self.pdbwin.update()
        self.pdbwin.pack(fill=BOTH,expand=1)
        #bind window closing so it's automatically recreated in the sidepane
        if sidepane == False:
            def func(evt):
                self.createPrimerDBGUI()
            self.pdbwin.bind("<Destroy>", func)
        return fr

    def applyRestrSiteSettings(self):
        if hasattr(self, 'restrsitedialog'):
            vals = self.restrsitedialog.getValues()
        self.sc.update()
        return

    def restrictionEnzymesDialog(self, sidepane=True):
        defaults = {'uniquesitesonly':1,'showonlyenzymescutmax':50,
                    'ignorecutpos':1,'minlengthrecseq':5,'excludepromiscuous':1}
        options = [('uniquesitesonly','checkbutton',1,'Show unique sites only'),
                   ('showonlyenzymescutmax','scale',(1,100),'Show only enzymes that cut max'),
                   ('ignorecutpos','checkbutton',1,'Ignore cut position when finding isoschiziomers'),
                   ('minlengthrecseq','scale',(1,10),'Min length of recognition sequence'),
                   ('excludepromiscuous','checkbutton',1,'Exclude promiscuous enzymes')]
        geometry = '400x250+500+250'
        fr = self.createChildFrame(name='Configure Restriction Enzymes',sidepane=sidepane)
        self.restrsitedialog = Dialogs.GenericDialog(fr, options, defaults)
        defaults = {'uniquesitesonly':1,'showonlyenzymescutmax':50,
                    'ignorecutpos':1,'minlengthrecseq':5,'excludepromiscuous':1}
        self.restrsitedialog.pack(fill=BOTH,expand=1)
        return

    def applySettings(self):
        self.applyPrefs(self.prefsdialog.getValues())
        #self.removeAll()
        #self.setupGUI()
        #self.globalPrefsDialog()
        return

    def globalPrefsDialog(self, sidepane=True):
        curr = {}
        for key in self.defaultprefs:
            curr[key] = self.__dict__[key]
        fr = self.createChildFrame(name='Global settings', sidepane=sidepane)
        self.prefsdialog = Dialogs.PrefsDialog(fr, self, self.defaultopts,
                                                curr, self.applySettings)
        self.prefsdialog.pack(fill=BOTH,expand=1)
        return

    def help(self):
        import webbrowser
        link='http://code.google.com/p/peat/wiki/DNAtool'
        webbrowser.open(link,autoraise=1)
        return

    def about(self):
        win=Toplevel()
        win.geometry('+500+350')
        win.title('About DNATool')
        win.maxsize(width=400,height=400)
        logo = Images.logo()
        label = Label(win,image=logo)
        label.image = logo
        label.pack(fill=BOTH,padx=4,pady=4)
        text="""DNATool App, Version 2
             is a python desktop application for DNA sequence manipulation.
             Released under GPL v3
             (C) Copyright 2012- Damien Farrell """
        text=text.replace('\t','')
        text= ' '.join(text.split())
        Label(win,text=text,wraplength=400).pack(fill=Y,side=TOP,pady=4)
        return

    def quit(self):
        if self.rememberwindowsize == 1:
            self.saveGeometry()
        self.main.destroy()
        if not self.parent:
            sys.exit()
        return

class RestrictionDigestSummary(Frame):
    """Summary of restr digest"""
    def __init__(self, parent, parentapp=None):
        self.parent = parent
        self.parentapp = parentapp
        Frame.__init__(self, parent)
        self.main=self
        self.makeGUI()
        return

    def makeGUI(self):
        self.sortvar = StringVar()
        self.sortvar.set('# of cuts')
        s = Pmw.OptionMenu(self.main,
                        labelpos = 'w',
                        label_text = 'Sort by:',
                        menubutton_textvariable = self.sortvar,
                        items = ['alphabetical','# of cuts'],
                        menubutton_width = 10)
        s.pack(side=TOP,fill=BOTH)
        self.details = Pmw.ScrolledText(self.main)
        self.details.pack(side=TOP,fill=BOTH,expand=1)
        return

class ToolBar(Frame):
    """Uses the parent instance to provide the functions"""
    def __init__(self, parent=None, parentapp=None):
        Frame.__init__(self, parent, width=600, height=30)
        self.parentframe = parent
        self.parentapp = parentapp
        return

    def addButton(self, name, callback, img=None, helptxt=None):
        if img==None:# and helptxt==None:
            b = Button(self, text=name, command=callback, height=1,
                         relief=GROOVE, overrelief=RAISED)
        else:
             b = Button(self, text=name, command=callback, width=20, height=20,
                             relief=GROOVE, image=img, overrelief=RAISED)
        b.image = img
        b.pack(side=LEFT, padx=2, pady=2, ipadx=1, ipady=1)
        #add balloon
        if helptxt!=None:
            balloon=Pmw.Balloon(self.parentframe)
            balloon.bind(b, helptxt)
        return

    def addSeperator(self):
        Label(self, text=' ').pack(side=LEFT, padx=2, pady=2)
        return

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="file",
                        help="Sequence file", metavar="FILE")
    parser.add_option("-p", "--project", dest="project",
                        help="Project file", metavar="FILE")

    opts, remainder = parser.parse_args()
    app = App()
    if opts.file != None:
        app.openSequence(opts.file)
    if opts.project != None:
        app.openProject(opts.project)
    app.mainloop()

if __name__ == '__main__':
    main()

