#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#
# Written by Damien Farrell, Feb 2010

try:
    from Plugins import Plugin
except:
    from PEATDB.Plugins import Plugin
from Tkinter import *
from PEATDB.GUI_helper import *
from PEATDB.Ekin.Titration import TitrationAnalyser
from PEATDB.Base import PDatabase
from PEATDB.Ekin.Base import EkinProject
from PEATDB.Tables import TableCanvas, ColumnHeader
from PEATDB.Ekin.Tables import EkinProjModel, EkinProjTable
from PEATDB.Ekin.Ekin_main import PlotPanel
from PEATDB.DictEdit import DictEditor
import os
import pickle

class NMRTitration(Plugin, GUI_help):
    """NMR ph titration analysis"""
    capabilities = ['gui']
    menuentry = 'NMR Titration Analysis'

    def main(self, parent=None):      
        if parent!=None:           
            self.DB = parent.DB
        self._doFrame()
        return

    def _doFrame(self):
        if self.parent == None:
            self.mainwin=Toplevel()           
        else:    
            self.mainwin=Toplevel()
        # Get platform into a variable
        import platform
        self.currplatform = platform.system()
        self.mainwin.title('NMR Titration Analysis')
        self.mainwin.geometry('800x600+200+100')
        self.createMenuBar()
        self.ekinprojects = {}      
        self.addToolBar()
        return

    # GUI stuff

    def createMenuBar(self):
        """Create the menu bar for the application. """
        self.menu=Menu(self.mainwin)

        self.file_menu={ '01Open Ekin project':{'cmd':self.loadEkinProj},            
                         '02Open DB':{'cmd':self.openDB},   
                         '03Quit':{'cmd':self.quit}}
        self.file_menu=self.create_pulldown(self.menu,self.file_menu)
        self.menu.add_cascade(label='File',menu=self.file_menu['var'])

        self.util_menu={'01Edit Dict':{'cmd': self.editDict},
                        '02Autoset Residue Names':{'cmd': self.add_ResidueNames},
                        '03Autoset Residue Numbers':{'cmd': self.add_ResidueNumbers}}
                       
        self.util_menu=self.create_pulldown(self.menu,self.util_menu)
        self.menu.add_cascade(label='Utils',menu=self.util_menu['var'])

        self.anal_menu={'01Chem. shift distrib':{'cmd': self.analysepKas},
                        '02Compare Across Nuclei':{'cmd': self.compareNuclei},
                        '03Ghost Mapping':{'cmd': self.doMapping}}
        self.anal_menu=self.create_pulldown(self.menu,self.anal_menu)
        self.menu.add_cascade(label='Analysis',menu=self.anal_menu['var'])
        self.mainwin.config(menu=self.menu)
        return
        
    def addToolBar(self):
        fr=self.toolbar=Frame(self.mainwin)
        fr.pack(side=TOP,fill=X)     
        Button(fr, text='Get Ekin prj from DB', command=self.selectfromDB,
                    relief=GROOVE, bg='#B0C4DE').pack(side=LEFT)        
        Button(fr, text='Plot Selected', command=self.plotSelected,
                    relief=GROOVE, bg='#B0C4DE').pack(side=LEFT)
        '''Button(fr, text='Fit Selected', command=self.fittSelected,
                    relief=GROOVE, bg='#B0C4DE').pack() '''
        self.dbstatus=StringVar()
        Label(self.toolbar,textvariable=self.dbstatus).pack(side=LEFT)
        return

    def updateDBStatus(self):
        if self.DB != None:
            self.dbstatus.set('db loaded with %s recs' %len(self.DB.getRecs()))
        return
    
    def selectfromDB(self):
        """Get an ekin prj from the DB"""        
        if self.DB == None:
            return        
        from PEATDB.Actions import DBActions
        fr=Toplevel()
        rbox, cbox = DBActions.getRecordsSelector(self.DB,fr)
        def loadselected():
            item = rbox.curselection()[0]
            rec = self.DB.getRecs()[int(item)]
            item = cbox.curselection()[0]
            col = self.DB.getFields()[int(item)]           
            E=self.DB[rec][col]
            self.loadEkinProj(E)
            fr.destroy()
        Button(fr, text='OK', 
                command=loadselected).grid(row=3,column=0,
                columnspan=2,sticky='news',padx=1,pady=3)       
        return
        
    def showEkinProject(self, E):
        """Show ekin prj"""
        if hasattr(self, 'pw'):
            self.pw.destroy()        
        self.pw = PanedWindow(self.mainwin,
                           orient=HORIZONTAL,
                           sashwidth=3,
                           showhandle=True,
                           opaqueresize=False)   
        self.pw.pack(side=LEFT,fill=BOTH,expand=1)  
        self.showEkinTable(self.pw, E)
        self.showPlotPanel(self.pw, E)
        return
        
    def showEkinTable(self, parent, E=None):
        """Show a list of ekin prj datasets in a table"""      
        self.currentmodel = EkinProjModel(E)
        fr=Frame(parent)
        self.pw.add(fr,minsize=200)        
        self.ekintable = EkinProjTable(fr, self.currentmodel)       
        self.ekintable.createTableFrame()
        return        
    
    def showPlotPanel(self, parent, E=None):
        """Add an ekin plot frame"""
        fr=Frame(parent)
        self.pw.add(fr)
        self.plotframe = PlotPanel(parent=fr, tools=True)
        if E != None:
            E.checkDatasets()
            self.plotframe.setProject(E)
            self.plotSelected()
        return
        
    def plotSelected(self):
        """Plot ekin """
        datasets = self.ekintable.get_selectedRecordNames()        
        self.plotframe.plotCurrent(datasets=datasets)
        return

    def add_ResidueNames(self):
        """Try to set residue names for selected datasets"""
        return

    def add_ResidueNumbers(self):
        """Try to set residue numbers from dataset names"""
        return

    def analysepKas(self):
        """Get the main pKas of all/titr group and do analysis"""
        E = self.currprj
        if E==None: return
        t = TitrationAnalyser()
        p = t.findpKas(E, titratable=False, reliable=False, minspan=0.06)       
        t.analysepKas(p)
        return

    def compareNuclei(self):
        """Compare corresponding datasets for a protein """     
        return
    
    def doMapping(self):
        return
        
    def editDict(self):
        """Edit ghost pka mapping in table"""
        D = DictEditor(self)
        return
    
    # Mainly IO methods    

    def openDB(self):
        import tkFileDialog
        filename=tkFileDialog.askopenfilename(defaultextension='.fs',
                                              filetypes=[("PEAT DB","*.fs"),
                                                        ("All files","*.*")],                                               
                                              parent=self.mainwin)  
        if filename != None:
            self.loadDB(filename)
            self.updateDBStatus()
        return    

    def loadEkinProj(self, E=None):
        """Load an ekin project file"""
        import os, types
        if E == None:
            import tkFileDialog
            filename=tkFileDialog.askopenfilename(defaultextension='.ekinprj',
                                                  filetypes=[("Ekin project","*.ekinprj"),
                                                             ("All files","*.*")],                                               
                                                  parent=self.mainwin)
            if filename != None:
                if os.path.isfile(filename):
                    fd=open(filename)
                    import pickle
                    data=pickle.load(fd)
                    E=EkinProject(data=data)               
                    self.ekinprojects[filename] = E 
                    fd.close()
            else:
                return
        self.currprj = E       
        self.showEkinProject(E)     
        return

    def addEkinProj(self, ekinproj=None, replace=1):
        """Load an ekin project file"""
        return
    
    def quit(self):
        self.mainwin.destroy()
        return

    def analyseTitDB(self, DB, col, names=None):
        """Extract titdb pKas"""
        nuclnames = {'1H NMR':'H','15N NMR':'N'}
        t = TitrationAnalyser()        
        #extract reliable pkas from selected proteins     
        #p=t.extractpKas(DB, col, names=names, titratable=False, reliable=False, minspan=0.06)
        #t.analysepKas(p)
        t.compareNuclei(DB, '15N NMR', '1H NMR', titratable=False, names=names)
        
        #ghost mapping..
        #t.mappKas(DB,col,p,names=['Protein G B1'],nucleus=nuclnames[col],calculatespans=False)       
        #t.mappKas(DB,col,p,names=['Xylanase (Bacillus subtilus)'],nucleus=nuclnames[col],calculatespans=False)
        
        return
        
    def save(self, DB, col, prot, E):
        if E==None and DB != None:
            E.saveProject(prot+'_'+col)
            DB[prot][col] = E            
        else:
            E.currentmode = 'NMR titration'
            E.saveProject()
        return
            
    def titDBUtils(self, DB=None, col=None, prot=None, a=None, E=None,
                    refit=False, addmeta=False, getexperrs=False, 
                    yuncert=None):
        """Add some meta and refit all for an ekin prj or a rec/field in db"""
        if E==None and DB != None:
            E = DB[prot][col]
            E.checkDatasets()
        t = TitrationAnalyser() 
        if refit == True:
            models = ['Linear', '1 pKa 2 Chemical shifts', 
                        '2 pKas, 3 Chemical shifts',
                        '3 pKas, 4 Chemical shifts']        
            E = t.findBest(E, models, geterrs=False)
        if addmeta == True:
            E = t.setMetaInfo(E, atom=a)
        if getexperrs == True:
            if yuncert == None:
                print 'No value for Y uncertainty!, please supply it'
                return
            print 'Using %s for y uncertainty.' %yuncert
            print
            E = t.getExpErrs(E, xuncert=0.1, yuncert=yuncert)
        self.save(DB, col, prot, E)
        #DB.commit('refit/added meta info')
        return E 
        
    def benchmarkExpErr(self, DB):
        """Test effects of model and noise on est exp error technique"""
        E = DB['HEWL']['1H NMR']
        d='D66N-HN'
        E.plotDatasets(d,filename='experrtest.png')
        print E.getFitData(d)
        ferrs = E.estimateExpUncertainty(d, runs=20, 
                xuncert=0.1, yuncert=0.03)
        print ferrs
        return

    def addpKaTables(self, DB, names, col='1H NMR'):
        """Create labbook tables for 'real' pKas for required proteins"""
        t = TitrationAnalyser()
        prots = t.getProtNames(DB)
        for name in names:
            recname = DB.getRecordName(name)
            E = DB[recname][col]
            titrresidues = t.getResidues(E, titratable=True)
            S = DB.createLabbookSheet(name+'.pKas')
            for r in titrresidues:
                d, res, resnum = r
                pKa = ''
                S.addRecord(res+resnum,pka=pKa,resname=res,
                            resnum=resnum,error='')           
            DB.saveLabbook(name+'.pKas', S)
        DB.saveLabbooktoFile('titdb.labbook')
        return
        
def main():
    """Run some analysis"""
    from optparse import OptionParser
    parser = OptionParser()
    app = NMRTitration()
    DB=None; E=None
    parser.add_option("-f", "--file", dest="file",
                        help="Open a local db")
    parser.add_option("-e", "--ekinprj", dest="ekinprj",
                        help="Open an ekin project")    
    parser.add_option("-t", "--titdb", dest="titdb", action='store_true',
                       help="titr db analysis", default=False)
    parser.add_option("-r", "--refit", dest="refit", action='store_true',
                       help="refit specific ekin data", default=False)
    parser.add_option("-u", "--getexperrs", dest="getexperrs", action='store_true',
                       help="get exp uncertainties", default=False)   
    parser.add_option("-m", "--addmeta", dest="addmeta", action='store_true',
                       help="add meta data for NMR", default=False)     
    parser.add_option("-x", "--analyse", dest="analyse", action='store_true',
                       help="extract pKas", default=False)     
    parser.add_option("-p", "--protein", dest="protein", help="protein")
    parser.add_option("-c", "--col", dest="col", help="field")
    parser.add_option("-a", "--atom", dest="atom", help="atom")
    parser.add_option("-b", "--benchmark", dest="benchmark", action='store_true',
                       help="benchmark some stuff", default=False)        
    parser.add_option("-g", "--gui", dest="gui", action='store_true',
                       help="start gui app", default=False)  
     
    opts, remainder = parser.parse_args()   
    if opts.file != None and os.path.exists(opts.file):
        app.loadDB(opts.file)
        
    if opts.gui == True:
        app.main()
        app.mainwin.mainloop()
        return
        
    #some tit db funcs    
    if opts.titdb == True:
        DB = PDatabase(server='peat.ucd.ie', username='guest',
                       password='123', project='titration_db',
                       port=8080)
        
    yuncerts = {'H':0.03,'N':0.1,'C':0.2}
    complete = ['HEWL', 'Bovine Beta-Lactoglobulin',
                'Plastocyanin (Anabaena variabilis)',
                'Plastocyanin (Phormidium)',
                'Glutaredoxin', 'CexCD (Apo)',
                'Protein G B1','Xylanase (Bacillus subtilus)']
    try:
        yuncert=yuncerts[opts.atom]
    except:
        yuncert=None
    if opts.ekinprj != None:
        E = EkinProject()
        E.openProject(opts.ekinprj)
    
    if opts.analyse == True:        
        app.analyseTitDB(DB, opts.col, complete)
    elif opts.benchmark == True:
        app.benchmarkExpErr(DB)
    else:
        app.titDBUtils(DB, opts.col, opts.protein, a=opts.atom, E=E,
                        refit=opts.refit, addmeta=opts.addmeta, 
                        getexperrs=opts.getexperrs, yuncert=yuncert)
        #app.addpKaTables(DB, complete)
        
if __name__ == '__main__':
    main()
