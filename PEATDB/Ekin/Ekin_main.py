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

import pickle, os, math
from Tkinter import*
import Pmw
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

from PEATDB.Tables import TableCanvas, ColumnHeader
from PEATDB.TableModels import TableModel
import PEATDB.Table_images as Table_images
from PEATDB.Ekin.Base import *
from PEATDB.Ekin.Tables import EkinDataModel, EkinDataTable
from PEATDB.Ekin.Fitting import Fitting
from PEATDB.Ekin.Convert import EkinConvert
from PEATDB.Ekin.Ekin_map import *
from PEATDB.Ekin.IO import *
import Ekin_images
from PEATDB.Ekin.Titration import TitrationAnalyser
from PEATDB.GUI_helper import *

class EkinApp(Frame, Ekin_map_annotate, GUI_help):

    def __init__(self, parent=None, project=None, data=None, mode=None,
                  protein=None, field=None, allowreturntoDB=None):
        "Initialize the application."
        self.parent=parent
        if not self.parent:
            Frame.__init__(self)
            self.ekin_win=self.master
        else:
            self.ekin_win=Toplevel()
            self.master=self.ekin_win
        self.ekin_win.title('Ekin - Data fitting tool')
        self.ekin_win.geometry('1000x720+200+100')
        self.ekin_win.protocol('WM_DELETE_WINDOW',self.quit)

        self.ekinfields = EkinProject.ekinfields
        self.modes = EkinProject.modes
        self.mode_definition = EkinProject.mode_definition
        self.path = os.getcwd()
       
        self.protein=protein
        self.filename=None
        self.allowreturntoDB = None
        if allowreturntoDB == 1:
            self.allowreturntoDB = 1
        if self.parent and protein:
            self.allowreturntoDB = 1
            try:
                self.protein_name=self.parent.data['DBinstance'].DB[protein]['Name']
            except:
                self.protein_name=self.parent.DB.data[protein]['name']
            #check for other open instances for the cell if called in PEAT
            if self.parent.ekin_instances.has_key(protein+field):
                if self.parent.ekin_instances[protein+field] > 1:
                    print 'EKIN ALREADY OPENED FOR THIS CELL'
                    self.allowreturntoDB = 0
        else:
            self.protein_name=protein
        self.field=field

        #create IO objects
        self.savedir = os.getcwd()
        self.importer = Importer(self)
        self.exporter = Exporter(self)
        self.importer.savedir = self.savedir
        self.exporter.savedir = self.savedir

        self.setupVars()
        self.setupGUI()

        #if ekinprojectfile provided at command line, we load it
        if project != None:
            if type(project) is EkinProject:
                self.setProject(project, mode)
            else:
                #should be changed to loadFile??
                self.loadProject(project=project)
        else:
            #just load the data
            self.loadProject(data=data, mode=mode)
            pass
        return

    def setupGUI(self):
        """Do GUI elements"""
        self.createDatasetMenu()
        self.m = PanedWindow(self.ekin_win,
                           orient=HORIZONTAL,
                           sashwidth=3,
                           showhandle=True)        
        self.m.pack(side=RIGHT, fill=BOTH,expand=1)        
        self.sidepane = Frame(self.m)
        self.m.add(self.sidepane)         
        plotpane = Frame(self.m)
        self.m.add(plotpane)
        #add tabular data view frame
        self.tableframe = DataPanel(parent=self.sidepane)
        self.tableframe.pack(side=TOP,fill=BOTH,expand=1)

        self.modespecific_frame = Frame(self.sidepane, borderwidth=2,relief=GROOVE)
        self.do_mode_specific_frame(self.modespecific_frame)
        self.modespecific_frame.pack(fill=X)

        #add plotter and fitter panels, these are independent classes
        self.plotframe = PlotPanel(plotpane, side=BOTTOM, tools=True)
        self.createFitFrame()
        #a small toolbar
        self.apptoolBar = ToolBar(plotpane, self)
        self.apptoolBar.pack(side=TOP,fill=X)       
        if self.parent and self.allowreturntoDB == 1:
            self.apptoolBar.add_button('Return Data to DB', self.return_data,
                                        help='Return current data in table to DB',
                                        side=LEFT)
        self.createMenuBar()
        self.createBindings()
        return

    def createFitFrame(self):
        self.fitframe = FitterPanel(parent=self.sidepane)
        self.fitframe.setParentApp(self)
        self.fitframe.setPlotter(self.plotframe)
        self.fitframe.pack(side=BOTTOM,fill=BOTH)
        return

    def showDatasetsTable(self):
        """Show a list of all datasets in a table"""
        if self.showdatasetstable.get() == 1:              
            from PEATDB.Ekin.Tables import EkinProjModel, EkinProjTable

            def plotselected():
                datasets = self.datasetstable.get_selectedRecordNames()        
                self.plotframe.plotCurrent(datasets=datasets, plotoption=self.overlayPlots.get())
            def createTable():
                self.currentmodel = EkinProjModel(self.E)
                self.datasetsframe = Frame(self.m)
                self.m.add(self.datasetsframe,minsize=250)
                self.datasetstable = EkinProjTable(self.datasetsframe, self.currentmodel)       
                self.datasetstable.createTableFrame()
            def refresh():
                self.currentmodel = EkinProjModel(self.E)
                self.datasetstable.redrawTable()
            createTable()
            b=Frame(self.datasetsframe)
            b.grid(row=3,column=1)
            Button(b, text='plot selected', command=plotselected).pack(side=LEFT,fill=BOTH)
            Button(b, text='refresh', command=refresh).pack(side=LEFT,fill=BOTH)
        else:
            self.m.forget(self.datasetsframe)
        return

    def setupVars(self):
       
        self.displaytabs={}
        self.currentdataset = StringVar()
        self.currenttable=None
        self.project_open=1
        if not hasattr(self,'defaultsavedir'):
            self.defaultsavedir = os.getcwd()
        self.display_multiple=IntVar()
        self.display_multiple.set(0)
        self.display_all=IntVar()
        self.display_all.set(0)
        self.no_multiplePlots=IntVar()
        self.no_multiplePlots.set(4)
        self.no_multipleCols=IntVar()
        self.no_multipleCols.set(0)
        self.overlayPlots=IntVar()
        self.overlayPlots.set(2)
        self.compare_multiple=IntVar()
        self.compare_multiple.set(0)
        self.mode_vars=[]
        self.mode_var=IntVar()
        self.mode_var.set(-1)
        self.showfitpanel=IntVar()
        self.showfitpanel.set(1)
        self.showdatasetstable=IntVar()
        self.showdatasetstable.set(0)
        
        self.kcat=StringVar()
        self.substr_conc=StringVar()
        self.enzyme_conc=StringVar()
        self.residue_type = StringVar()
        self.residue_num = StringVar()
        self.atom_type = StringVar()
        self.chain_id = StringVar()
        #pdb things for titr stuff are defined in this class

        self.residue_list = TitrationAnalyser.residue_list
        self.atom_types = TitrationAnalyser.atom_types
        return

    def createMenuBar(self):
        """Create the menu bar for the application. """
        self.menu=Menu(self.ekin_win)
        self.file_menu={ '01New':{'cmd':self.newProject},
                         '02Open':{'cmd':self.openProject},
                         '03Close':{'cmd':self.closeProject},
                         '04Save':{'cmd':self.saveProject},
                         '05Save As':{'cmd':self.save_as_project},
                         '06Append':{'cmd':self.appendProject},
                         '06Quit':{'cmd':self.quit}}

        self.file_menu=self.create_pulldown(self.menu,self.file_menu)
        self.menu.add_cascade(label='File',menu=self.file_menu['var'])

        self.data_menu={'01Open dataset':{'cmd':self.read_dataset},
                        '02Save dataset':{'cmd':self.saveDataset},
                        '03sep':{None:None},
                        '04Add New Dataset':{'cmd':self.addDataset},
                        '05Rename Dataset':{'cmd':self.renameDataset},
                        '06Copy Dataset':{'cmd':self.copy_dataset},
                        '07Delete Dataset':{'cmd':self.deleteDataset},
                        '08Delete All':{'cmd':self.deleteAll},
                        '09sep':{None:None},
                        '10Adjust Data':{'cmd':self.adjustDataDialog}}

        self.data_menu=self.create_pulldown(self.menu,self.data_menu)
        self.menu.add_cascade(label='Data',menu=self.data_menu['var'])

        self.IE_menu={  '01Import from CSV file':{'cmd':self.import_csv},
                        '02Export tab to CSV file':{'cmd':self.export_csv},
                        '03Import Sparky peak files':{'cmd':self.import_chem_shift},
                        #'04Import from chem. shift file':{'cmd':self.import_chem_shift},
                        '05Import from CcpNmr file':{'cmd':self.import_ccpnmr},
                        #'06Transfer peaks to pKaSystem':{'cmd':self.open_pKaSystem_dialog},
                        '07sep':{None:None},
                        '08Import CD temperature scan':{'cmd':self.import_CD_tempscan},
                        '09Analyse temp dependence of CD data':{'cmd':self.insert_CD_temp_datatab}
                        }

        self.IE_menu=self.create_pulldown(self.menu,self.IE_menu)
        self.menu.add_cascade(label='Import/Export',menu=self.IE_menu['var'])

        #set up mode menu
        self.mode_menu=Menu(self.menu,tearoff=0)
        self.doModeMenu()
        self.menu.add_cascade(label='Mode',menu=self.mode_menu)

        # Annotate/Map menu
        #
        self.annotate_menu={'01Edit meta data':{'cmd':self.editMetaData},
                            #'02View all meta data':{'cmd':self.viewMetaData},
                            '02Map datatab to structure':{'cmd':self.map_datatab2structure}}
        self.annotate_menu=self.create_pulldown(self.menu,self.annotate_menu)
        self.menu.add_cascade(label='Annotate/Map',menu=self.annotate_menu['var'])

        # Plot menu
        self.plot_menu=Menu(self.menu,tearoff=0)

        self.plot_menu.add_command(label='Plot options',
                                       command=self.plotframe.plotOptions)
        self.plot_menu.add_checkbutton(label='Show multiple plots',
                                       command=self.multiplePlotsOption,
                                       variable=self.display_multiple,onvalue=1,offvalue=0)
        self.plot_menu.add_checkbutton(label='Show all',
                                       command=self.redrawGraph,
                                       variable=self.display_all,onvalue=1,offvalue=0)
        self.plot_menu.add_checkbutton(label='Compare plots',
                                       command=self.compareAllDialog,
                                       variable=self.compare_multiple,onvalue=1,offvalue=0)
        self.plot_menu.add_checkbutton(label='Overlay view',
                                       command=self.redrawGraph,
                                       variable=self.overlayPlots,onvalue=3,offvalue=2)
        self.plot_menu.add('separator')
        self.plot_menu.add_command(label='Disable selected',
                                       command=self.disableSelected)
        self.plot_menu.add_command(label='Enable selected',
                                        command=lambda: self.disableSelected(1))
        self.plot_menu.add_command(label='Delete selected',
                                        command=lambda: self.deleteSelected(1))        
        self.menu.add_cascade(label='Plot',menu=self.plot_menu)

        # Fit menu
        self.fit_menu=Menu(self.menu,tearoff=0)
        self.fit_menu={ '01Find best model':{'cmd':self.fitframe.showBestModelDialog},
                        '02Fit all datasets':{'cmd':self.fitAll},
                        '03Create new fit model':{'cmd':self.createFitModel}}
        self.fit_menu=self.create_pulldown(self.menu,self.fit_menu)
        self.menu.add_cascade(label='Fit',menu=self.fit_menu['var'])

        self.view_menu=Menu(self.menu,tearoff=0)       
        self.view_menu.add_checkbutton(label='Show fit panel',command=self.updateView,
                                       variable=self.showfitpanel,onvalue=1,offvalue=0)        
        self.view_menu.add_checkbutton(label='Show datasets in table',command=self.showDatasetsTable,
                                       variable=self.showdatasetstable,onvalue=1,offvalue=0)     
        self.menu.add_cascade(label='View',menu=self.view_menu)
        
        # Help menu
        self.help_menu={'01Online Help':{'cmd':self.online_documentation},
                        '02About':{'cmd':self.about}}
        self.help_menu=self.create_pulldown(self.menu,self.help_menu)
        self.menu.add_cascade(label='Help',menu=self.help_menu['var'])
        self.ekin_win.config(menu=self.menu)

        return

    def createBindings(self):
        """Bind keys"""
        self.ekin_win.bind("<Control-n>", self.newProject)
        self.ekin_win.bind("<Control-o>", self.openProject)
        self.ekin_win.bind("<Control-s>", self.saveProject)
        #self.ekin_win.bind("<1>", lambda:self.disableSelected(0))
        self.ekin_win.bind("<Control-q>", self.quit)
        self.ekin_win.bind("<Control-i>", self.import_csv)
        self.ekin_win.bind("<Control-e>", self.export_csv)
        self.ekin_win.bind_all("<F1>", self.about)
        self.ekin_win.bind_all("<Control-r>", self.redrawGraph)
        #self.ekin_win.focus_set()
        #self.tk_focusFollowsMouse()
        return

    def do_mode_specific_frame(self, parent):
        """Create mode specific things in a frame"""
        l=Label(parent,text='      ')
        l.grid(row=0,column=0, sticky='news')
        return

    def createDatasetMenu(self):
        """Do menu for selecting the dataset and dataset operations"""

        row=0
        self.dataset_menu=Frame(self.ekin_win,
                             height=30,
                             width=300,
                             borderwidth=2)
        #self.dataset_menu.grid(row=row,column=0, columnspan=2, rowspan=1, sticky='w')
        self.dataset_menu.pack(fill=X )
        self.dataset_select_button=None


        prev = Ekin_images.prev()
        self.prevbutton = Button(self.dataset_menu,text='Previous',
                                 image=prev,
                                 #compound='left',
                                 #width=11,
                                 borderwidth=2, relief=GROOVE,
                                 command=self.prev_dataset)
        self.prevbutton.grid(row=0,column=1,sticky='nws',ipadx=2, ipady=2, padx=1,pady=1)
        self.prevbutton.image = prev
        self.prevballoon=Pmw.Balloon(self.ekin_win)
        self.prevballoon.bind(self.prevbutton, 'Previous dataset')

        next = Ekin_images.next()
        self.nextbutton = Button(self.dataset_menu,text='Next',
                                image=next,
                                borderwidth=2, relief=GROOVE,
                                command=self.next_dataset)
        self.nextbutton.image = next
        self.nextballoon=Pmw.Balloon(self.ekin_win)
        self.nextballoon.bind(self.nextbutton, 'Next dataset')
        self.nextbutton.grid(row=0,column=2,sticky='nws',ipadx=2, ipady=2, padx=1,pady=1)

        add = Ekin_images.add()
        self.addDB = Button(self.dataset_menu,text='Add dataset',
                            image=add,
                            #width=11,
                            borderwidth=2, relief=GROOVE,
                            command=self.addDataset)
        self.addDB.image = add
        self.addballoon=Pmw.Balloon(self.ekin_win)
        self.addballoon.bind(self.addDB, 'Add dataset')
        self.addDB.grid(row=row,column=3,sticky='nws',ipadx=2, ipady=2, padx=1,pady=1)

        delb = Ekin_images.delb()
        self.deleteDB = Button(self.dataset_menu,text='Delete dataset',
                            image=delb,
                            #width=11,
                            borderwidth=2, relief=GROOVE,
                            command=self.deleteDataset)
        self.deleteDB.image = delb
        self.delballoon=Pmw.Balloon(self.ekin_win)
        self.delballoon.bind(self.deleteDB, 'Delete dataset')
        self.deleteDB.grid(row=row,column=4,sticky='nws',ipadx=2, ipady=2, padx=1,pady=1)

        self.sortbyDB = Menubutton(self.dataset_menu,
                                        text='Sort by',
                                        width=12,
                                        borderwidth=2, relief=GROOVE)

        self.sortby_menu = Menu(self.sortbyDB, tearoff=0)
        self.sortbyDB['menu']=self.sortby_menu
        self.sortby_var = StringVar()
        self.sortby_var.set('name')
        self.sortby_menu.add_radiobutton(label='Name',
                           variable=self.sortby_var,
                           value='name',
                           indicatoron=1,
                           command=self.update_sortby)
        self.sortby_menu.add_radiobutton(label='First number',
                           variable=self.sortby_var,
                           value='firstnum',
                           indicatoron=1,
                           command=self.update_sortby)
        self.sortbyballoon=Pmw.Balloon(self.ekin_win)
        self.sortbyballoon.bind(self.sortbyDB, 'Sort datasets by name or 1st num')
        self.sortbyDB.grid(row=row,column=5,sticky='nws',padx=1)
        self.updateDatasetSelector()
        return

    def doModeMenu(self):
        """Put modes in menu"""
        if not hasattr(self, 'E'):
            return
        self.mode_menu.delete(0, self.mode_menu.index('last'))
        count=0        
        for mode in self.modes:
            #If a mode was set AND we are in PEAT, lock it in
            #if self.parent!=None and not self.E.currentmode or (self.E.currentmode and self.E.currentmode==mode):
            self.mode_menu.add_radiobutton(label=mode,
                               variable=self.mode_var,
                               value=count,
                               indicatoron=1,
                               command=self.updateMode)            
            count=count+1
        #print self.E.currentmode, self.E.mode
        if self.E.currentmode == None:
            self.E.currentmode = self.E.mode
        self.mode_var.set(self.modes.index(self.E.currentmode))
        return

    def updateMode(self, setmode=1):
        """Update the mode, change the selection of fitters available"""
        current = self.currentdataset.get()
        # Usually only update currentmode when menu changed
        if setmode == 1 or not self.E.currentmode:
            self.E.currentmode = self.modes[self.mode_var.get()]
            self.fitframe.setMode(self.E.currentmode)           

        if current == None:
            return
        self.clear_modespecific_frame()
        if self.E.currentmode == 'NMR titration':
            self.do_titration_stuff()
        return

    def updateView(self):
        if self.showfitpanel.get() == 1:
            self.createFitFrame()
            self.updateAll()             
        else:
            self.fitframe.destroy()
        return
        
    def updateAll(self, d=None):
        """Update all elements to reflect current dataset"""
        if len(self.E.datasets) == 0:
            return
        if d!=None: self.currentdataset.set(d)
        self.update_dataset()
        self.fitframe.setData(self.E, self.currentdataset.get())
        self.fitframe.update()
        #reset the multiple display options when someone selects a dataset
        self.compare_multiple.set(0)
        self.display_all.set(0)
        self.redrawGraph()
        self.updateModeSpecificData()
        if hasattr(self, 'M') and self.M != None:
            self.M.applyCurrent()
            self.M.updateFields(dataset=self.currentdataset.get())
        #self.showDatasetsTable()    
        return

    def update_dataset(self):
        """Update data table to reflect current dataset"""
        if len(self.E.datasets) == 0:
            return
        ekindata = self.data[self.currentdataset.get()]
        model = self.fitframe.model_type.get()
        self.tableframe.draw_Table(ekindata, model)
        self.tableframe.addBinding("<Return>", self.redrawGraph)
        #callback for table entry
        #self.tableframe.setCallback(self.redrawGraph)

        return

    def update_labels(self):
        """Reset dataset labels to reflect the new mode"""
        for dataset in self.E.datasets:
            self.data[dataset][0]['label']=self.FIT.data_labels[0]
            self.data[dataset][1]['label']=self.FIT.data_labels[1]
        return

    def updateFittedData(self, fitdata):
        """Transfer the fitted values to the fit arrays for the present tab"""
        dataset = self.currentdataset.get()
        self.E.setFitData(dataset, fitdata)
        return

    def update_sortby(self):
        """Update the sorted order for datasets in notebook"""
        self.updateDatasetSelector()
        return

    def updateDatasetSelector(self,parent=None):
        """Updates the dataset selector to reflect current datasets"""
        if not hasattr(self, 'E'):
            return
        if parent == None:
            parent = self.dataset_menu
        if self.dataset_select_button:
            self.dataset_select_button.destroy()
        if len(self.E.datasets) == 0:
            self.currentdataset.set(None)
            self.dataset_select_button = Button(parent,text='No Data')
            self.dataset_select_button.grid(row=0,column=0,columnspan=1,padx=4,sticky='nws')
            return
        self.dataset_select_button=Menubutton(parent,
                                            textvariable=self.currentdataset,
                                            width=12,
                                            borderwidth=2, relief=GROOVE,
                                            background='#488AC7',
                                            foreground='yellow',
                                            activeforeground='red')
        self.dataset_sel_menu=Menu(self.dataset_select_button,tearoff=0)
        self.dataset_select_button['menu']=self.dataset_sel_menu

        # Add the names of all datasets
        names=self.E.datasets
        if self.sortby_var.get() == 'firstnum':        
            resorted = self.sort_by_Num(names)
        else:
            names.sort()
            resorted = names
        p=0

        for keys in resorted:
            if self.sortby_var.get() == 'firstnum':
                name = keys[1]
            else:
                name = keys
            if p%30==0 and p!=0:
                self.dataset_sel_menu.add_radiobutton(label=name,
                                                  variable=self.currentdataset,
                                                  value=name,
                                                  indicatoron=0,
                                                  columnbreak=1,
                                                  command=self.updateAll)
            else:
                 self.dataset_sel_menu.add_radiobutton(label=name,
                                                  variable=self.currentdataset,
                                                  value=name,
                                                  indicatoron=0,
                                                  command=self.updateAll)

            p=p+1

        self.dataset_select_button.grid(row=0,column=0,columnspan=1,padx=4,sticky='nws')

        # Add protein name and current mode in top bar
        if self.protein!=None:
            proteinlbl=Label(parent,text=self.protein_name,fg='Blue')
            proteinlbl.grid(row=0,column=8,columnspan=1,sticky='news',padx=4,pady=2)
        if self.E.currentmode != None:
            modelbl=Label(parent,text='mode: '+self.E.currentmode,fg='Blue')
            modelbl.grid(row=0,column=6,columnspan=1,sticky='news',padx=4,pady=2)
        if self.field!=None:
            #get columnlabel
            fieldlabel = self.parent.get_field_label(self.field)
            flbl=Label(parent,text='field: '+fieldlabel,fg='Blue')
            flbl.grid(row=0,column=7,sticky='news',padx=4,pady=2)
        elif self.filename!=None:
            filenamelbl=Label(parent,text=self.filename,fg='Blue')
            filenamelbl.grid(row=0,column=8,columnspan=1,sticky='news',padx=4,pady=2)             
        return

    def sort_by_Num(self, p):
        """Sort text keys by contained numbers"""
        splitkeys={}
        import re
        r=re.compile('\D')
        for k in p:
            value=None
            for splitval in r.split(k):
                try:
                    value=int(splitval)
                except ValueError:
                    pass
            if not value:
                print 'no number in', k
            splitkeys[k]=value
        items = splitkeys.items()
        items = [(v, k) for (k, v) in items]
        items.sort()
        return items
        
    def redrawGraph(self, event=None):
        """Replot"""
        if self.display_all.get() == 1:
            total = len(self.E.datasets)
            if total >= 25:
                tkMessageBox.showwarning('Warning',
                                         'There are too many plots to display at once.\nShowing the first 25 only.',
                                         parent=self.ekin_win)
                total=25
            self.plotframe.plotCurrent(showmultiple=True, datasets=self.E.datasets[:25],
                                    plotoption=self.overlayPlots.get())

        elif self.display_multiple.get() == 1:
            #plot x multiple graphs from current dataset
            selecteddatasets=[]
            names = self.E.datasets
            names.sort()
            curr = self.currentdataset.get()
            total = len(names)
            pos = names.index(curr)
            p = self.no_multiplePlots.get()
            if p > total:
                p = total
            c = self.no_multipleCols.get()
            for i in range(pos, pos+p):
                if i < total:
                    selecteddatasets.append(names[i])
            self.plotframe.plotCurrent( datasets=selecteddatasets,
                                    plotoption=self.overlayPlots.get(), cols=c)
        elif self.compare_multiple.get()==1:
            #compare user chosen datasets
            selecteddatasets=[]
            for dataset in self.E.datasets:
                if self.displaytabs[dataset].get()==1:
                    selecteddatasets.append(dataset)
            c = int(math.ceil(math.sqrt(len(selecteddatasets))))
            self.plotframe.plotCurrent(datasets=selecteddatasets,
                                                plotoption=self.overlayPlots.get(), cols=c)
        else:           
            self.plotframe.plotCurrent(datasets=self.currentdataset.get())
        return


    def insertDataset(self, newdata, newname, fit=None, update=1):
        """Insert pre-existing data as a dataset"""

        name = newname
        if newname in self.E.datasets:
            if tkMessageBox.askyesno('Overwrite Dataset?',name+' is already present.\n'
                                      'Do you want to overwrite it?',
                                      parent=self.ekin_win):
                self.deleteDataset(name, confirm=None)                
            else:
                name = tkSimpleDialog.askstring("Dataset name?",
                                   "Enter a dataset name:",
                                    parent=self.ekin_win)
            if not name:
                return

        self.E.insertDataset(newdata, name, fit)
        self.currentdataset.set(name)
        self.updateDatasetSelector()
        if update == 1:
            self.updateAll()
        return

    def insertMultipleDatasets(self, newdata, overwrite=None):
        """Adds multiple datasets at once from a dict of ekin datasets"""
        if newdata != None:
            for name in newdata.keys():
                if name in self.E.datasets:
                    if overwrite == 1:
                        self.deleteDataset(name, confirm=None, update=False)
                        self.insertDataset(newdata[name], name, update=None)
                else:
                    self.insertDataset(newdata[name], name, update=None)
        else:
            return
        if name in self.E.datasets:
            self.currentdataset.set(name)
        else:
            self.currentdataset.set(self.E.datasets[0])
        self.updateAll()
        return

    def addDataset(self, label=None, update=1):
        """Add a new empty dataset"""

        if label==None:
            label = tkSimpleDialog.askstring("New dataset name?",
                                               "Enter dataset name:",
                                               parent=self.ekin_win)
            if not label:
                return
        if label in self.E.datasets:
            import tkMessageBox
            tkMessageBox.showwarning('Error',
                                     'Dataset name exists',
                                     parent=self.ekin_win)
            return

        self.E.addDataset(label)
        self.updateDatasetSelector()
        self.currentdataset.set(label)
        if update == 1:
            self.updateAll()
        return


    def deleteDataset(self, name=None, confirm=1, update=True):
        """Delete a dataset"""
        if name == None:
            name = self.currentdataset.get()
        if confirm == 1:
            import tkMessageBox
            if not tkMessageBox.askyesno("Deleting Dataset","Are you sure you want to delete?",
                                    parent=self.ekin_win):
                print 'cancelled'
                return
        self.E.deleteDataset(name)
        self.updateDatasetSelector()
        if update == True:
            self.prev_dataset()
            self.updateAll()
        return

    def renameDataset(self):
        """Rename a dataset"""
        currname = self.currentdataset.get()
        newname=tkSimpleDialog.askstring('Rename data series',
                                             'Enter a new name',
                                             initialvalue=currname,
                                             parent=self.ekin_win)
        if newname == None or newname == currname:
            return
        print newname    
        self.E.copyDataset(currname, newname)
        self.deleteDataset(currname, confirm=None, update=False)
        self.currentdataset.set(newname)
        self.updateAll()
        return

    def copy_dataset(self, newname=None):
        """Copy a dataset"""
        name = self.currentdataset.get()
        if newname == None:
            newname=tkSimpleDialog.askstring('Copy data series',
                                             'Enter a new name',
                                             parent=self.ekin_win)
        if newname == None:
            return
        elif newname in self.E.datasets:
            import tkMessageBox
            sorry = tkMessageBox.showinfo("Duplicate Name","Name Exists.\nChoose a different name",
                                            parent=self.ekin_win)
            return

        else:
            self.E.copyDataset(name, newname)
            self.updateDatasetSelector()
            self.updateAll()
        return

    def deleteAll(self, confirm=1):
        """Delete All"""
        if confirm == 1:
            import tkMessageBox
            if not tkMessageBox.askyesno("Deleting All","Are you sure you want\n to delete all datasets?",
                                    parent=self.ekin_win):
                return
        self.E.deleteAll()
        self.updateDatasetSelector()
        self.updateAll()
        return

    def prev_dataset(self):
        name=self.currentdataset.get()
        names=self.E.datasets
        if self.sortby_var.get() == 'firstnum':           
            resorted = self.sort_by_Num(names)
            index=0
            for nowname in resorted:
                if name==nowname[1]:
                    break
                index=index+1
            if index-1>=0:
                self.currentdataset.set(resorted[index-1][1])
                self.updateAll()
        else:
            names.sort()
            index=0
            for nowname in names:
                if name==nowname:
                    break
                index=index+1
            if index-1>=0:
                self.currentdataset.set(names[index-1])
                self.updateAll()
        return

    def next_dataset(self):
        name=self.currentdataset.get()
        names=self.E.datasets
        if self.sortby_var.get() == 'firstnum':           
            resorted = self.sort_by_Num(names)
            index=0
            for nowname in resorted:
                if name==nowname[1]:
                    break
                index=index+1
            if len(names)>index+1:
                self.currentdataset.set(resorted[index+1][1])
                self.updateAll()
        else:
            names.sort()
            index=0
            for nowname in names:
                if name==nowname:
                    break
                index=index+1
            if len(names)>index+1:
                self.currentdataset.set(names[index+1])
                self.updateAll()
        return

    def saveDataset(self,event=None):
        """Save dataset to a file"""
        import tkFileDialog, os
        filename=tkFileDialog.asksaveasfilename(defaultextension='.Ekindat',
                                                initialdir=os.getcwd(),
                                                filetypes=[("Ekin data","*.Ekindat"),
                                                           ("All files","*.*")])
        if not filename:
            return
        name=self.currentdataset.get()
        self.E.saveDataset(name, filename)
        return

    def read_dataset(self,event=None):
        """Read dataset from a file"""
        import tkFileDialog
        filename=tkFileDialog.askopenfilename(defaultextension='.Ekindat',initialdir=self.savedir,
                                              filetypes=[("Ekin data","*.Ekindat"),
                                                         ("All files","*.*")],
                                              parent=self.ekin_win)
        newdata = {}
        if os.path.isfile(filename):
            fd=open(filename)
            import pickle
            newdata=pickle.load(fd)
            fd.close()
            # Add the data tab with a filename
            showname=os.path.split(filename)[1]
            showname=showname.replace('_','-')
            
        #old format      
        if type(newdata) is types.DictType and newdata.has_key('data'):
            newdata=newdata['data']
        self.insertDataset(newdata, showname)
        return

    def multiplePlotsOption(self):
        """How many plots to display at once"""
        if self.display_multiple.get() == 1:
            self.display_all.set(0)
            '''np=tkSimpleDialog.askinteger(title='Select multiple plots',
                                             prompt='How many plots do you wish to\n view at once?',
                                             parent=self.ekin_win)'''
            mpDlg = MultipleValDialog(title='No. of plots',
                                        initialvalues=(self.no_multiplePlots.get(),self.no_multipleCols.get()),
                                        labels=('no. of plots','no. of cols'), types=('int','int'),
                                        parent=self.ekin_win)
            if mpDlg.result == True:
                np = mpDlg.results[0]
                nc = mpDlg.results[1]
            self.no_multiplePlots.set(np)
            self.no_multipleCols.set(nc)
        self.redrawGraph()
        return

    def compareAllDialog(self):
        """Choose which datasets to plot together"""
        self.display_multiple.set(0)
        self.display_all.set(0)
        if self.compare_multiple.get() == 0:
            self.redrawGraph()
            return
        if hasattr(self, 'compall_dialog') and self.compall_dialog != None:
            self.compall_dialog.deiconify()
            return

        #self.createSidePane(200)
        #self.compall_dialog = Frame(self.sidepane, width=300, height=500)
        #self.compall_dialog.pack(side=LEFT,fill=BOTH, expand=1, padx=4,pady=4)
        self.compall_dialog = Toplevel()
        self.compall_dialog.geometry('+800+200')
        self.compall_dialog.title('Compare all')

        self.datasetlist = LabelFrame(self.compall_dialog, text="Choose Datasets", padx=5, pady=5)
        self.datasetlist.grid(row=0, column=0, columnspan=2, sticky='news')
        r=0
        c=0
        i=0

        for name in self.E.datasets:
            if r>20:
                r=0
                c=c+2
            r=r+1
            Label(self.datasetlist,text=name).grid(row=r,column=c)
            self.displaytabs[name] = IntVar()
            self.displaytabs[name].set(0)
            cb=Checkbutton(self.datasetlist, onvalue=1,offvalue=0,variable=self.displaytabs[name])
            cb.grid(row=r,column=c+1,padx=2,pady=2)

        def close():
            if self.compall_dialog:
                self.compall_dialog.destroy()
                self.compall_dialog = None
                #self.resetSidePane()
            return

        b = Button(self.compall_dialog, text="Apply", command=self.redrawGraph)
        b.grid(row=1,column=1,sticky='news',padx=4,pady=4)
        c=Button(self.compall_dialog,text='Close', command=close)
        c.grid(row=1,column=0,sticky='news',padx=4,pady=4)
        return

    def adjustDataDialog(self):
        """Edit datasets globally - dialog"""
        self.adjustdialog = Toplevel()
        self.adjustdialog.title('Adjust Data')
        self.set_geometry(self.ekin_win,self.adjustdialog)
        var = DoubleVar()
        var.set(0.0)
        columnvar = IntVar()
        columnvar.set(0)
        opvar = StringVar()
        opvar.set('+')
        row=0
        column=0
        def close():
            if hasattr(self, 'adjustdialog'):
                self.adjustdialog.destroy()
        a1 = Frame(self.adjustdialog)
        a1.grid(row=row,column=0,padx=2,pady=2,sticky='news')
        Label(a1, text='Value:').pack(side=LEFT,fill=BOTH)
        Entry(a1, textvariable=var, width=10, bg='white').pack(side=LEFT,fill=BOTH)
        Label(a1, text='Column:').pack(side=LEFT,fill=BOTH)
        Entry(a1, textvariable=columnvar, width=4, bg='white').pack(side=LEFT,fill=BOTH)
        opbutton = Menubutton(a1,textvariable=opvar,
				 relief=GROOVE, width=16, bg='lightblue')
        opmenu = Menu(opbutton, tearoff=0)
        opbutton['menu'] = opmenu
        i=0
        for p in ['+','-','*','/']:
            opmenu.add_radiobutton(label=p, variable=opvar,
                                    value=p, indicatoron=1)
            i+=1
        opbutton.pack(side=LEFT,fill=BOTH)
        d = self.currentdataset.get()
        Ek = self.E.getDataset(d)
        row=1
        def apply():
            Ek.adjust(op='+',column=columnvar.get(),val=var.get())
            self.update_dataset()
            self.redrawGraph()
        Button(self.adjustdialog,text='Apply',command=apply).grid(row=row,column=0,padx=2,sticky='news')
        Button(self.adjustdialog,text='Close',command=close).grid(row=row,column=1,padx=2,sticky='news')

        return

    def clear_modespecific_frame(self):
        if hasattr(self, 'mode_widgets') and self.mode_widgets!=None:
            for widget in self.mode_widgets:
                widget.destroy()
            self.mode_widgets=None
        return

    def updateModeSpecificData(self):
        """Set specific mode data widgets when dataset changed"""
        d = self.currentdataset.get()            
        meta = self.E.getMetaData(d)
        if self.E.currentmode=='Simple enzyme kinetic' or self.E.currentmode == 'Enzyme pH-activity':
            if meta.has_key('enzyme_conc'):
                self.enzyme_conc.set(meta['enzyme_conc'])
            else:
                self.enzyme_conc.set('')
            if meta.has_key('substr_conc'):
                 self.substr_conc.set(meta['substr_conc'])
            else:
                self.substr_conc.set('')
            if meta.has_key('kcat'):
                 self.kcat.set(meta['kcat'])
            else:
                self.kcat.set('')
        elif self.E.currentmode == 'NMR titration':
            #print 'nmr'
            if meta.has_key('residue'):
                self.residue_type.set(meta['residue'])
            else:
                self.residue_type.set('')
            if meta.has_key('res_num'):
                self.residue_num.set(meta['res_num'])
            else:
                self.residue_num.set('')
            if meta.has_key('atom_type'):
                self.atom_type.set(meta['atom_type'])
            else:
                self.atom_type.set('')
            if meta.has_key('chain_id'):
                self.chain_id.set(meta['chain_id'])
            else:
                self.chain_id.set('')
        return

    def do_titration_stuff(self):
        """Add residue info for titration mode"""

        self.mode_widgets=[]
        self.chain_id_entry = Entry(self.modespecific_frame,
                                       background = 'yellow', width=2,relief=GROOVE,
                                       textvariable=self.chain_id)
        self.chain_id_entry.bind('<KeyRelease>', self.set_chain_id)
        self.mode_widgets.append(self.chain_id_entry)
        self.chain_id_entry.grid(row=0,column=0,pady=1,sticky='news')
        self.residue_menu = Pmw.OptionMenu(self.modespecific_frame,
                                labelpos = 'w',
                                label_text = 'Residue:',
                                menubutton_textvariable = self.residue_type,
                                items = self.residue_list,
                                menubutton_width = 10,
                                command = self.set_residue_type )
        self.residue_menu.grid(row=0,column=2)
        self.residue_menu.component('label').configure(background = '#ffffcc')
        self.mode_widgets.append(self.residue_menu)
        self.residue_num_entry = Entry(self.modespecific_frame,
                                       background = 'white', width=5,relief=GROOVE,
                                       textvariable=self.residue_num)
        self.residue_num_entry.bind('<KeyRelease>', self.set_residue_num)
        self.residue_num_entry.grid(row=0,column=3,pady=1,sticky='news')
        self.mode_widgets.append(self.residue_num_entry)
        self.atom_type_menu = Pmw.OptionMenu(self.modespecific_frame,
                                labelpos = 'w',
                                menubutton_textvariable = self.atom_type,
                                items = self.atom_types,
                                menubutton_width = 5,
                                command = self.set_atom_type )

        self.mode_widgets.append(self.atom_type_menu)
        self.atom_type_menu.grid(row=0,column=4)
        return

    def set_residue_type(self, event=None):
        d = self.currentdataset.get()     
        if d == None:
            return
        self.E.addMeta(d,'residue',self.residue_type.get())
        return

    def set_residue_num(self, event=None):
        d = self.currentdataset.get()        
        if d == None:
            return     
        self.E.addMeta(d,'res_num',self.residue_num.get())
        return

    def set_atom_type(self, event=None):
        d = self.currentdataset.get()       
        if d == None:
            return
        self.E.addMeta(d,'atom_type',self.atom_type.get())
        return

    def set_chain_id(self, event=None):
        d = self.currentdataset.get()       
        if d == None:
            return       
        self.E.addMeta(d,'chain_id',self.chain_id.get())
        return

    def save_comments(self,event=None):
        """Set comments for this dataset"""
        d = self.currentdataset.get()     
        if d != None:
            self.E.addMeta(d, 'comment', self.details.get(1.0, END))
        self.commentframe.destroy()
        return

    def addedit_comments(self,event=None):
        """Show window with specific comments for this dataset"""

        self.commentframe=Toplevel()
        self.commentframe.title('Add Comments')
        self.commentframe.geometry("300x300")
        self.set_geometry(self.ekin_win,self.commentframe)
        row=0
        column=0
        import tkFont
        thefont = tkFont.Font ( family="Helvetica", size=11, weight="bold" )

        yscrollbar=Scrollbar(self.commentframe,orient='vertical',width=14)
        yscrollbar.grid(row=row,column=2,sticky='news',padx=2)
        xscrollbar=Scrollbar(self.commentframe,orient='horizontal',width=14)
        xscrollbar.grid(row=row+1,column=0,columnspan=2,sticky='news')
        self.details=Text(self.commentframe,background='white',
                      foreground='black',
                      state=NORMAL,
                      exportselection=1,
                      yscrollcommand=yscrollbar.set,
                      xscrollcommand=xscrollbar.set,
                      wrap='word',
                      font=thefont,
                      bg='lightgray')
        self.details.config(width=60)
        self.details.config(width=24)
        yscrollbar.config(command=self.details.yview)
        xscrollbar.config(command=self.details.xview)
        self.details.grid(row=0,column=0,columnspan=2,sticky='NEWS',padx=2,pady=2)
        self.details.config(state=NORMAL)
        x=Button(self.commentframe,text='save & close',command=self.save_comments)
        x.grid(row=3,column=0,sticky='NEWS',padx=3,pady=3)
        self.commentframe.rowconfigure(0,weight=1)
        self.commentframe.columnconfigure(0,weight=1)

        d = self.currentdataset.get()
        meta = self.E.getMetaData(d)
        if d != None:
            if meta.has_key('comment'):
                self.details.insert(END,meta['comment'])
        return

    def editMetaData(self):
        """Edit the metadata using the MetaData class"""
        self.M = MetaData(data=self.E.__meta_data__, project=self.E)
        if len(self.M.glob) == 0:
            self.M.autoCreate(mode = self.E.mode)
        d=self.currentdataset.get()
        self.mwin = self.M.createDialog(dataset=d, parent=self.ekin_win)
        self.ekin_win.wait_window(self.mwin)        
        return

    def viewMetaData(self):
        """View summary of all meta data"""
        M = MetaData(data=self.E.__meta_data__,project=self.E)
        M.viewAll()
        return

    def fitAll(self):
        """Fit all datasets"""
    
        self.stopfit = False
        def go():
            #if findbestvar.get() == 1:
            #    pass
            models = [fitmodelvar.get()]
            tot = len(self.E.datasets); p=0
            for d in self.E.datasets:
                if self.stopfit == True:
                    return
                self.fitall_win.update_idletasks()
                if len(models)==1:
                    try:
                        fdata, X = self.E.fitDataset(d, model=models[0],silent=True,noiter=no_itervar.get(),
                                                conv=converg_diff.get())#,callback=self.fitall_win.update)
                    except:
                        pass                    
                else:
                    fdata, p = self.E.findBestModel(d, models=models)#, checkfunc=checkfunc)
                thistab=p
                p=p+1
                m=(p/float(tot))*100
                fitbar.updateProgress(newValue=m)
            self.updateAll()
            return

        def stop():
            self.stopfit=True
            self.fitall_win.destroy()
            return

        #fit all data dialog
        self.stopfitall=0
        self.fitall_win=Toplevel()
        #self.fitall_win.geometry('+300+300')
        self.set_geometry(self.ekin_win,self.fitall_win)
        self.fitall_win.title('Fit All Datasets')
        self.fitall_win.grab_set()
        self.fitall_win.transient(self.ekin_win)
        lbl0=Label(self.fitall_win,text='This will loop through all datasets and \n'
            'fit them in turn. Set a reasonable number of iterations.\n'
            'You can also change the convergence value (larger=faster).',bg='yellow')
        lbl0.grid(row=0,column=0,columnspan=2,padx=3,pady=2)
        lbl1=Label(self.fitall_win,text='Max Iterations:')
        lbl1.grid(row=1,column=0,padx=3,pady=2)
        no_itervar = IntVar(); no_itervar.set(100)
        entry1=Entry(self.fitall_win,textvariable=no_itervar,bg='white')
        entry1.grid(row=1,column=1,padx=3,pady=2)

        lbl1=Label(self.fitall_win,text='Convergence at:')
        lbl1.grid(row=2,column=0,padx=3,pady=2)
        converg_diff = DoubleVar(); converg_diff.set(1e-9)
        entry1=Entry(self.fitall_win,textvariable=converg_diff,bg='white')
        entry1.grid(row=2,column=1,padx=3,pady=2)

        #add model chooser here
        fitmodelvar = StringVar(); fitmodelvar.set('Linear')
        Label(self.fitall_win,text='Model to use:').grid(row=3,column=0,sticky='news')
        model_button = Menubutton(self.fitall_win,textvariable=fitmodelvar,relief=GROOVE,
                                    width=10)
        model_menu = Menu(model_button,tearoff=0)
        model_button['menu'] = model_menu
        modelts = EkinProject.mode_definition['General']
        for text in modelts:
            model_menu.add_radiobutton(label=text,
                            variable=fitmodelvar,
                            value=text,
                            indicatoron=1)
        model_button.grid(row=3,column=1,sticky='news',padx=2,pady=3)

        doallfitsbutton = Button(self.fitall_win, text="Go", width=20, command=go)
        doallfitsbutton.grid(row=5,column=0,padx=3,pady=2)
        cancelbutton = Button(self.fitall_win, text="Cancel", width=20, command=stop)
        cancelbutton.grid(row=5,column=1,padx=3,pady=2)

        from PEATDB.ProgressBar import ProgressBar
        Label(self.fitall_win,text='Progress:').grid(row=6,column=0,padx=2,pady=4)
        fitbar = ProgressBar(self.fitall_win,fillColor='lightblue',width=120)
        fitbar.frame.grid(row=6,column=1,padx=2,pady=4)

        return

    def disableSelected(self, status=0):
        """Disable/Enable selected points in current plot"""
        bounds = self.plotframe.selection
        name = self.currentdataset.get()
        ek = self.E.getDataset(name)
        ek.setActiveBounded(bounds,status)
        self.update_dataset()
        self.redrawGraph()
        return
        
    def deleteSelected(self, evt=None):
        """Delete selected points"""
        bounds = self.plotframe.selection
        name = self.currentdataset.get()
        ek = self.E.getDataset(name)
        ek.removeBounded(bounds)
        self.update_dataset()
        self.redrawGraph()        
        return        

    def createFitModel(self):
        from PEATDB.DictEdit import DictEditor
        app = DictEditor(self.ekin_win)
        app.loadDict(Fitting.modelsfile)
        #self.ekin_win.wait_window(app)
        #update models
        print 'model file',Fitting.modelsfile
        Fitting.presetmodels = pickle.load(open(Fitting.modelsfile,'r'))  
        return
    
    #
    # Import and export functions, most should use the Importer class to do the actual
    # importing stuff and get the returned dataset that we simply insert
    #
    def import_csv(self, event=None):
        """Will allow import of multiple or single text files"""
        self.importer.path = self.path
        newdatasets = self.importer.import_multiple()        
        self.insertMultipleDatasets(newdatasets)

        return

    def export_csv(self, event=None):
        """Will allow export of multiple or single csv files"""
        name = self.currentdataset.get()
        self.exporter.export_csv(self.data[name], self.__datatabs_fits__[name])
        return

    def import_chem_shift(self):
        """Import sparky files"""
        self.importer.open_sparky_dialog()
        return

    def import_ccpnmr(self):
        """Import CCPNmr file"""
        newdatasets = self.importer.import_ccpnmr()
        #set the mode to NMR
        self.mode_var.set(4)
        self.updateMode()
        self.insertMultipleDatasets(newdatasets)
        return

    def import_CD_tempscan(self):
        """Import CD temp scan data"""
        newdatasets = self.importer.import_CD_tempscan()
        self.mode_var.set(5)
        self.updateMode()
        self.insertMultipleDatasets(newdatasets)
        return

    def insert_CD_temp_datatab(self):

        return

    def setProject(self, E, mode='General'):
        """Set from an ekinproject object"""
        self.E = E
        self.E.checkDatasets()
        self.E.checkMeta()
        self.data = self.E.data
        
        if not hasattr(self.E, 'mode'):
            self.E.currentmode = self.E.mode = mode

        #self.M = MetaData(self.E.__meta_data__, project=self.E)
        #if len(self.M.glob) == 0:
        #    self.M.autoCreate(mode = mode)
        self.M = None
        self.project_open = 1
        self.updateDatasetSelector()
        self.plotframe.setProject(self.E)
        self.doModeMenu()
        self.updateMode(setmode=1)        
        if hasattr(self.E,'__currentdataset__'):
            self.currentdataset.set(self.E.__currentdataset__)
        elif len(self.E.datasets)>0:
            self.currentdataset.set(self.E.datasets[0])
        if hasattr(self.E, '__displaymultiple__'):
            self.display_multiple.set(self.E.__displaymultiple__)
        if hasattr(self.E, '__nomultiplePlots__'):
            self.no_multiplePlots.set(self.E.__nomultiplePlots__)
        self.updateAll()
        return

    def loadProject(self, project=None, data=None, mode=None):
        """Load ekin project and update GUI"""
 
        if mode == None:
            self.mode = mode ='General'
        if project != None and os.path.exists(project):
            self.E = EkinProject()
            self.E.openProject(project)
            self.E.checkDatasets()
            self.filename = project
            self.path = os.path.abspath(self.filename)
        elif data != None:
            self.E = EkinProject(data=data, mode=mode)
        else:
            self.E = EkinProject(mode=mode)
            self.E.addDataset('data')
        self.data = self.E.data   
        self.M = None
        if not hasattr(self.E, 'mode'):
            self.E.mode = mode

        self.project_open = 1
        self.updateDatasetSelector()
        self.plotframe.setProject(self.E)
        self.doModeMenu()
        self.updateMode(setmode=1)
        self.mode_var.set(self.modes.index(mode))
        if hasattr(self.E,'__currentdataset__'):
            self.currentdataset.set(self.E.__currentdataset__)
        elif len(self.E.datasets)>0:
            self.currentdataset.set(self.E.datasets[0])
        if hasattr(self.E, '__displaymultiple__'):
            self.display_multiple.set(self.E.__displaymultiple__)
        if hasattr(self.E, '__nomultiplePlots__'):
            self.no_multiplePlots.set(self.E.__nomultiplePlots__)

        self.updateAll()
        return


    def newProject(self, event=None):
        """Create a new project"""
        if self.project_open == 1:
            import tkMessageBox
            if tkMessageBox.askyesno("Close Project?","Close Current Project?",
                                    parent=self.ekin_win):
                self.closeProject(ask=0)
            else:
                return
        if self.project_open == 0:
            self.E.currentmode=None
            self.loadProject()
            self.updateAll()
            return
        else:
            return

    def openProject(self,filename=None):
        """Open a project"""
        import tkFileDialog, os
        if filename==None:
            # check current project first
            if self.project_open == 1:
                result = self.saveDialog()
                if result == 'Yes':
                    self.saveProject()
                    #self.closeProject(ask=0)
                #elif result == 'No':
                #    self.closeProject(ask=0)
                elif result == 'Cancel':
                    return
            filename=tkFileDialog.askopenfilename(defaultextension='.ekinprj',
                                                  initialdir=self.path,
                                                  filetypes=[("Ekin project","*.ekinprj"),
                                                             ("All files","*.*")],
                                                  parent=self.ekin_win)
        if not filename:
            return
        self.loadProject(filename)
        self.filename=filename
        self.path=os.path.abspath(filename)
        return

    def saveProject(self, event=None):
        "Save a project"""
        if hasattr(self.plotframe,'Opts') and self.plotframe.Opts != None:
            import copy
            self.E.__plotopts__ = copy.deepcopy(self.plotframe.Opts.opts)
        self.E.__currentdataset__ = self.currentdataset.get()
        self.E.__displaymultiple__ = self.display_multiple.get()
        self.E.__nomultiplePlots__ = self.no_multiplePlots.get()
        if hasattr(self, 'filename') and self.filename != None:
            self.E.saveProject(self.filename)
        else:
            self.save_as_project()
        return

    def save_as_project(self):
        """Save as project"""
        import tkFileDialog, os
        filename = tkFileDialog.asksaveasfilename(parent=self.ekin_win,
                                                defaultextension='.ekinprj',
                                                initialdir=self.defaultsavedir,
                                                filetypes=[("Ekin project","*.ekinprj"),
                                                           ("All files","*.*")])
        if not filename:
            return
        self.E.saveProject(filename)
        # Set the filename
        self.filename=filename
        return

    def appendProject(self):
        """Append another project to this one"""
        filename=tkFileDialog.askopenfilename(defaultextension='.ekinprj',
                                              initialdir=os.getcwd(),
                                              filetypes=[("Ekin project","*.ekinprj"),
                                                         ("All files","*.*")],
                                              parent=self.ekin_win)
        if filename != None:
            E1 = EkinProject()
            E1.openProject(filename)
            self.E.addProject(E1)
            self.updateDatasetSelector()
            self.updateAll()
        return

    def closeProject(self,ask=1):
        """Close Ekin project"""
        if ask == 1 and self.project_open == 1:
            import tkMessageBox
            if not tkMessageBox.askyesno("Close Project?","Close Current Project?",
                                    parent=self.ekin_win):
                return

        self.project_open = 0
        self.filename=None        
        self.newProject()
        return

    def saveDialog(self):
        """Save dialog"""
        import Pmw
        dialog = Pmw.MessageDialog(self.parent,
                    title = 'Close Project?',
                    defaultbutton = 0,
                    buttons = ('Yes', 'No', 'Cancel'),
                    message_text = 'Save Current Project before closing?')
        dialog.iconname('Simple message dialog')
        result = dialog.activate()
        return result

    def return_data(self):
        """Return to PEAT or some parent app"""
        if self.protein and self.field:
            self.E._PEATinfo = {'record':self.protein,'column':self.field}

        #remove attribs that cannot not be saved
        tempflds = ['fitline', 'ax']
        for f in tempflds:
            del self.E.__dict__[f]
        if hasattr(self.plotframe,'Opts') and self.plotframe.Opts != None:
            import copy
            self.E.__plotopts__ = copy.deepcopy(self.plotframe.Opts.opts)
        self.E.__currentdataset__ = self.currentdataset.get()
        self.E.__displaymultiple__ = self.display_multiple.get()
        self.E.__nomultiplePlots__ = self.no_multiplePlots.get()

        self.parent.Ekinprj = self.E
        self.ekin_win.destroy()
        return

    def online_documentation(self,event=None):
        """Open the online documentation"""
        import webbrowser
        link='http://enzyme.ucd.ie/main/index.php/Ekin'
        webbrowser.open(link,autoraise=1)
        return

    def about(self, event=None):
        """Display about window"""
        self.ab_win=Toplevel(self.ekin_win)
        self.ab_win.geometry('+100+350')
        self.ab_win.title('About Ekin')
        row=0

        import Ekin_images
        logo = Ekin_images.logo()
        label = Label(self.ab_win,image=logo)
        label.image = logo
        #newButton = Button(self.ab_win,image=logo)
        #label.pack()
        label.grid(row=row,column=0,sticky='news')
        row=1
        text=['Ekin Ver 3.0 ','Component of Protein Engineering and Analysis Tool',
        'A subcomponent of PEAT that lets you fit your experimental',
        'data to useful equations. Michaelis-Menten, Henderson-Hasselbalch',
        ', Bell-shaped curves determined by two pKa values',
        'Thermal denaturation, Chemical denaturation, etc.',
        'Matplotlib is required to use this application.',
        'Authors: Jens Erik Nielsen, Damien Farrell, University College Dublin',
        '(C) Copyright 2003- Jens Erik Nielsen All rights reserved']

        for line in text:
            tmp=Label(self.ab_win,text=line)
            tmp.grid(row=row,column=0,sticky='news')
            row=row+1
        return

    def get_geometry(self,widget):
        """Get the geometry of a widget
        Return width,height,xorg,yorg"""
        widget.update_idletasks()
        txt=widget.winfo_geometry()
        width=int(txt.split('x')[0])
        rest=txt.split('x')[1]
        height=int(rest.split('+')[0])
        xorg=int(rest.split('+')[1])
        yorg=int(rest.split('+')[2])
        return width,height,xorg,yorg

    def set_geometry(self,pwidget,widget):
        """Set the position of widget in the middle of pwidget"""
        w,h,x,y=self.get_geometry(pwidget)
        sw,sh,dummy,dummy2=self.get_geometry(widget)
        xoffset=int((w-sw)/2)
        yoffset=int((h-sh)/2)
        widget.geometry('+%d+%d' %(x+xoffset,y+yoffset))
        return


    def quit(self, event=None):
        """Quit application"""
        #ask to save current
        if self.project_open == 1:
            result = self.saveDialog()
            if result == 'Yes':
                self.saveProject()
            elif result == 'Cancel':
                return
        self.ekin_win.destroy()
        if not self.parent:
            sys.exit()
        return


class FitterPanel(Frame):
    """Display and manage ekin fit data - all Tk stuff here"""
    def __init__(self, parent):
        Frame.__init__(self, parent)
        self.parent=parent
        self.mode='General'
        self.createTkVars()
        self.doFrame()
        return

    def setParentApp(self, parentapp):
        self.parentapp=parentapp
        return

    def setFitData(self, data=None):
        """Set current ekin fit data"""
        self.fitdata=data
        return

    def setData(self, E, dataset):
        """Set current ekin fit data"""
        self.fitdata = E.__datatabs_fits__[dataset]
        self.data = E.data[dataset]      
        return

    def setPlotter(self, plotter):
        """Set the plotter for fit callback"""
        self.plotter = plotter
        return

    def createTkVars(self):
        """Create Tk vars for the fit panel"""
        self.model_type=StringVar()
        self.model_type.set('Linear')
        self.no_iter=IntVar()
        self.no_iter.set(300)
        self.converg_diff=DoubleVar()
        self.converg_diff.set(1e-9)
        self.guess_converg_diff=IntVar()
        self.guess_converg_diff.set(0)
        self.grad_crit=DoubleVar()
        self.grad_crit.set(1e-12)
        self.LM_damper = DoubleVar()
        self.LM_damper.set(2.0)
        self.guess_startparms=IntVar()
        self.guess_startparms.set(1)
        self.maxdiff=StringVar()
        self.maxerror=StringVar()
        self.maxdiff.set(0.3)
        self.maxerror.set(0.1)
        self.count=IntVar()
        self.count.set(0)

        self.fvars=[[IntVar(), StringVar(), DoubleVar(), DoubleVar()],
                    [IntVar(), StringVar(), DoubleVar(), DoubleVar()],
                    [IntVar(), StringVar(), DoubleVar(), DoubleVar()],
                    [IntVar(), StringVar(), DoubleVar(), DoubleVar()],
                    [IntVar(), StringVar(), DoubleVar(), DoubleVar()],
                    [IntVar(), StringVar(), DoubleVar(), DoubleVar()],
                    [IntVar(), StringVar(), DoubleVar(), DoubleVar()],
                    [IntVar(), StringVar(), DoubleVar(), DoubleVar()],
                    [IntVar(), StringVar(), DoubleVar(), DoubleVar()],
                    [IntVar(), StringVar(), DoubleVar(), DoubleVar()]
                    ]
        for f in self.fvars:
            f[0].set(1)
            f[2].set(0.0)
        return

    def doModelButton(self):
        """Mode drop down menu"""
        if hasattr(self, 'model_button'):
            self.model_button.destroy()
        self.model_button = Menubutton(self,textvariable=self.model_type,relief=RAISED,width=self.model_button_width)
        self.model_menu = Menu(self.model_button,tearoff=0)
        self.model_button['menu'] = self.model_menu

        modelts = EkinProject.mode_definition[self.mode]
        for text in modelts:
            self.model_menu.add_radiobutton(label=text,
                            variable=self.model_type,
                            value=text,
                            indicatoron=1,
                            command=lambda : self.update (reset=True))
        self.model_button_row=0
        self.model_button.grid(row=0,column=0,
                                sticky='nes',columnspan=2,
                                padx=2,pady=3)
        return

    def doFrame(self):
        """Create GUI elements for fitting frame"""
        self.model_button_width=40
        self.model_button=Menubutton(self,textvariable=self.model_type,relief=RAISED,
                                    width=self.model_button_width)
        self.model_menu=Menu(self.model_button,tearoff=0)
        self.model_button['menu']=self.model_menu
        row=1
        self.doModelButton()
        self.optionsButton()
        paramsFrame = Frame(self,bg='lightblue')
        paramsFrame.grid(row=row,column=0,columnspan=3)

        text=['Include','Name','Start value','','Fitted value']
        col = 0
        for t in text:
            Label(paramsFrame,text=t,font='times12bold').grid(row=0,column=col,sticky='news')
            col = col+1

        # Make space for parameters
        r=1
        for i in range(7):
            Checkbutton(paramsFrame,variable=self.fvars[i][0]).grid(row=r,column=0,sticky='news')
            Label(paramsFrame,textvariable=self.fvars[i][1]).grid(row=r,column=1,sticky='news')
            Entry(paramsFrame,width=8,justify='center',text=self.fvars[i][2],
                        bg='#FFFF99').grid(row=r,column=2,sticky='news')
            def transferVal(var=i):
                """Set the starting value to the fitted value"""
                if self.fvars[var][3].get() != '-':
                    self.fvars[var][2].set(self.fvars[var][3].get())
                return
            Button(paramsFrame,text='<--',command=transferVal,relief='groove').grid(row=r,column=3)
            Label(paramsFrame,textvariable=self.fvars[i][3]).grid(row=r,column=4,sticky='news')
            r=r+1

        # Button for calculating fit
        row=row+1
        self.sqdiff=StringVar()
        self.sqdiff.set('No fit')
        Label(self,text='Sum of sq diffs').grid(row=row,column=1,padx=2,pady=1,sticky='news')
        Label(self,textvariable=self.sqdiff).grid(row=row,column=2,padx=2,pady=1,sticky='news')
        row=row+1
        Button(self,text='Calculate best fit',command=self.doFit,relief='groove',
                    bg='lightblue').grid(row=row,column=0,padx=2,pady=2,sticky='news')

        # Button for just updating the fit
        row=row+1
        b=Button(self,text='Preview fit', command=self.previewFit,relief='groove',
                         bg='lightblue')
        b.grid(row=row,column=0,padx=2,pady=2,sticky='news',columnspan=1)
        Pmw.Balloon(self).bind(b, 'Update curve to reflect manual start parameters')

        # Button for estimation of experimental uncertainty
        b=Button(self,text='Estimate exp. uncertainty', command=self.showExpUncertaintyDialog,
                    relief='groove',bg='lightblue')
        b.grid(row=row,column=1,padx=2,pady=2,sticky='news',columnspan=2)
        Pmw.Balloon(self).bind(b, 'Get an error estimate for the fit params')

        # Radio button to select fitting method
        self.fitting_method = IntVar()
        row=row+1
 
        n=Frame(self)
        Label(n,text='Rounds:').pack(side=LEFT)
        Entry(n,textvariable=self.no_iter,width=6).pack(side=LEFT)
        n.grid(row=row,column=0,padx=2,pady=2,sticky='news')

        Label(self,text='LM damper start value').grid(row=row,column=1)
        Entry(self,width=8,justify='center',
              textvariable=self.LM_damper).grid(row=row,column=2,sticky='ns')

        lbl1=Label(self,text='Conv. criteria (error/gradient):')
        lbl1.grid(row=row+1,column=0,padx=3,pady=2)
        guessconvchk = Checkbutton(self,onvalue=1,offvalue=0,variable=self.guess_converg_diff)
        guessconvchk.grid(row=row+1,column=3,padx=1,pady=2)
        entry1=Entry(self,textvariable=self.converg_diff,width=8,justify='center',bg='white')
        entry1.grid(row=row+1,column=1,sticky='ns',padx=3,pady=2)

        entry2=Entry(self,textvariable=self.grad_crit,width=8,justify='center',bg='white')
        entry2.grid(row=row+1,column=2,sticky='ns',padx=3,pady=2)
        return

    def optionsButton(self):
        """Options menu"""
        self.optsbutton = Menubutton(self,text='Options',relief=RAISED,width=8,bg='#FFFF78')
        self.optsmenu = Menu(self.optsbutton,tearoff=0)
        self.optsbutton['menu']=self.optsmenu
        self.optsbutton.grid(row=0,column=2,
                                sticky='nes',columnspan=1,
                                padx=2,pady=3)

        self.optsmenu.add_checkbutton(label='Guess start values',
                                       variable=self.guess_startparms,onvalue=1,offvalue=0)
        self.optsmenu.add_checkbutton(label='Guess error crit.',
                                       variable=self.guess_converg_diff,onvalue=1,offvalue=0)
        # Send current fit to comments
        self.markasfitted=IntVar()
        self.markasfitted.set(0)
        self.optsmenu.add_checkbutton(label='Mark as Fitted',onvalue=1,offvalue=0,
                                variable=self.markasfitted, command=self.send_fit_to_comments)
        self.optsmenu.add_command(label='Find best model',
                                       command=self.showBestModelDialog)
        return

    def showStopFitFrame(self, parentframe):
        """Show the fit status frame while fitting"""
        if parentframe == None:
            self.count_win=Frame()
            self.count_win.pack()
        else:
            self.count_win=Frame(parentframe)
            self.count_win.grid(row=3,column=1,columnspan=2,sticky='news')
        Label(self.count_win,text='Iteration #').pack(fill=BOTH,side=LEFT)
        Label(self.count_win,textvariable=self.count).pack(fill=BOTH,side=LEFT)
        b=Button(self.count_win,text='Stop fit',command=self.stopFit,
                    relief=GROOVE,bg='#E55451',highlightbackground='#E55451')
        b.pack(fill=BOTH,side=LEFT)

        return

    def stopFit(self):
        """Stop the fit"""
        self.count_win.destroy()
        self.stopfit=True
        return

    def doFit(self):
        """Fit the currently loaded dataset and update plotter"""
        self.stopfit=False
        currfitdata = self.fitdata
        model = self.model_type.get()
        if self.guess_converg_diff.get() == 1:
            conv = None
        else:
            conv = self.converg_diff.get()
        grad = self.grad_crit.get()
        noiter = self.no_iter.get()
        guess = self.guess_startparms.get()
        damper = self.LM_damper.get()

        def updatenow(diff, vrs, fitvals, c, X):
            """callback for fit update"""
            fdata=Fitting.makeFitData(model, vrs, diff) #get ekin fitdata from vars            
            self.update(fitdata=fdata, startvars=False, X=X)
            self.update_idletasks()
            self.plotter.updateFit(X)
            self.count.set(c)
            self.count_win.update()
            if self.stopfit == True:
                X.stop_fit=1
            return

        self.plotter.plotCurrent()  #update plot first, just in case
        self.showStopFitFrame(self)

        changevars=[]
        startvals=[]
        for i in self.fvars:
            changevars.append(i[0].get())
            startvals.append(i[2].get())

        fitresult, X = Fitting.doFit(self.data, model=model,
                                    noiter=noiter, conv=conv, grad=grad, LM_damper=damper,
                                    guess=guess,
                                    silent=False, callback=updatenow,
                                    changevars=changevars,
                                    startvalues=startvals)

        if fitresult == None:
            self.stopFit()
            return         
        
        #send new fit to current ekin dataset
        self.parentapp.updateFittedData(self.fitdata)
        #update fitting frame when finished fit?
        self.stopFit()
        return

    def findBestModel(self, doall=False, callback=None):
        """Use fitting class to determine best fit model using f-testing"""
        models=[]
        E = self.parentapp.E
        for m in self.tempmodelvars:
            if self.tempmodelvars[m].get() == 1:
                models.append(m)
        if len(models)==0:
            return None        
        if doall == False:
            d = self.parentapp.currentdataset.get()
            res, p = E.findBestModel(d, models=models, silent=True)
            callback(res, d)
        else:             
            for d in E.datasets:
                if self.stopvar.get()==1:
                    break
                res, p = E.findBestModel(d, models=models, silent=True)
                callback(res, d)            
        return 

    def showBestModelDialog(self):
        """Find the best model using f-test criteria"""

        self.stopvar=IntVar(); self.stopvar.set(0)
        models = EkinProject.mode_definition[self.mode]
        self.fbm_win = Toplevel(self.parent.master)
        self.parentapp.set_geometry(self.parentapp.ekin_win, self.fbm_win)
        self.fbm_win.title('Find best model')

        frame1 = Frame(self.fbm_win)
        frame1.grid(row=0,column=0)
        info = '''Select models to include in the test, models\n must have increasing degrees of freedom\n and must be nested in some meaningful sense\n for this test to be useful'''
        Label(frame1,text=info,bg='yellow').grid(row=0,column=0,columnspan=3)
        self.tempmodelvars={}
        r=1; i=0
        for l in ['name','params',''] :
            Label(frame1,text=l).grid(row=r,column=i)
            i+=1
        r=2
        for m in models:
            self.tempmodelvars[m]=IntVar(); self.tempmodelvars[m].set(0)
            Label(frame1,text=m).grid(row=r,column=0)
            X = Fitting.getFitter(model=m)            
            n = len(X.variables)
            Label(frame1,text=n).grid(row=r,column=1)
            Checkbutton(frame1,variable=self.tempmodelvars[m],
                                onvalue=1,offvalue=0).grid(row=r,column=2)
            r+=1

        self.fbresults=Pmw.ScrolledText(self.fbm_win,
                                labelpos = 'n',
                                usehullsize = 1,
                                hull_width = 400,
                                hull_height = 400)
        self.fbresults.grid(row=0,column=1,sticky='news')
        # buttons
        frame3 = Frame(self.fbm_win)
        frame3.grid(row=1,column=0,columnspan=2)
        
        def update(res, d):
            print 'res',res
            if res == None:
                self.fbresults.insert(END, 'No models selected?\n')
            else:
                self.fbresults.insert(END, '%s: %s' %(d,res['model']))
                self.fbresults.insert(END, '\n')
            self.fbresults.yview("moveto",1.0)   
            #we dont update atm because of memory leak    
            #self.parentapp.updateAll(d)
            self.update_idletasks()
            frame1.update()
            
        def go():
            self.stopvar.set(0)
            self.findBestModel(doall=doallvar.get(), callback=update)             

        def stop():
            self.stopvar.set(1)
            
        Button(frame3, text='Go', relief=GROOVE, bg='#B0C4DE',
               command=go).pack(side=LEFT)
        def close():
            self.fbm_win.destroy()
        Button(frame3,text='Close', relief=GROOVE, bg='#B0C4DE',
               command=close).pack(side=LEFT)
        
        Button(frame3,text='Stop', relief=GROOVE, bg='red',
               command=stop).pack(side=LEFT)
        doallvar = IntVar(); doallvar.set(0)
        Label(frame3,text='do all').pack(side=LEFT)
        Checkbutton(frame3,variable=doallvar,
                            onvalue=1,offvalue=0).pack(side=LEFT)
        return

    def update(self, reset=False, fitdata=None, data=None,
                   startvars=True, X=None):
        """Reset the current model and/or set new data, if provided.
           reset=True resets all vars, ie. a new model has been selected"""

        if reset == True:
            model = self.model_type.get()
            #resets current fitdata to generic start values for that model            
            self.fitdata = Fitting.makeFitData(model=model)

        if fitdata != None:
            self.fitdata = fitdata
        if data != None:
            self.setData(data)
        self.updateVars(self.fitdata, reset=reset, X=X)
        if startvars == True and self.fitdata!=None:
            self.setStartVars()
        return

    def setMode(self, mode):
        self.mode = mode
        self.doModelButton()
        return

    def updateVars(self, fitdata, reset=True, X=None):
        """Update with current fitdata, only done when we load the dataset"""
        if fitdata != None and fitdata.has_key('model'):
            model = fitdata['model']            
            self.model_type.set(model)
            if X==None: X = Fitting.getFitter(model=model)
            if fitdata.has_key('error'):
                self.sqdiff.set('%9.4e' %float(fitdata['error']))
            else:
                self.sqdiff.set('N/A')
            for i in range(len(self.fvars)):
                if fitdata.has_key(i):
                    #in case fit param is stored as string
                    if isinstance(fitdata[i],str):
                        fitdata[i]=float(fitdata[i])
                    val=fitdata[i]
                    if reset:
                        self.fvars[i][0].set(1)
                    self.fvars[i][3].set('%9.4e' %val)
                else:
                    self.fvars[i][0].set(0)
                    self.fvars[i][3].set('-')
            #param labels
            for i in range(len(self.fvars)):
                try:
                    val = X.varnames[i]
                    self.fvars[i][1].set(val)
                except:
                    self.fvars[i][1].set('')
        else:
            #no fit data passed, so we set to default vals
            for i in range(len(self.fvars)):
                self.fvars[i][3].set('-')
        return

    def setStartVars(self, usefitted=False):
        """Update start vars only"""

        if usefitted == True:
            for i in range(len(self.fvars)):
                self.fvars[i][2].set(self.fvars[i][3].get())

        else:
            #otherwise get default start vars
            model = self.model_type.get()
            X = Fitting.getFitter(model)
            vrs = X.variables
            #print vrs
            for i in range(len(self.fvars)):
                try:
                    self.fvars[i][2].set(vrs[i])
                except:
                    self.fvars[i][2].set(0)
        return

    def showExpUncertaintyDialog(self):
        """Show dialog for doing Exp Uncertainty"""
        data = self.data; fitdata = self.fitdata
        model = fitdata['model']
        X = Fitting.getFitter(model)        
        vrs = X.variables
        names = X.getVarNames()
        self.eeu_win = Toplevel(self.parent.master)
        self.parentapp.set_geometry(self.parentapp.ekin_win, self.eeu_win)
        self.eeu_win.title('Estimate experimental uncertainty')

        frame1 = Frame(self.eeu_win)
        frame1.pack()
        self.eeu_x_uncertainty = DoubleVar()
        self.eeu_x_uncertainty.set(0.01)
        self.eeu_y_uncertainty = DoubleVar()
        self.eeu_y_uncertainty.set(0.01)
        self.eeu_no_runs = IntVar()
        self.eeu_no_runs.set(10)
        row=1; col=1

        Label(frame1,text='Uncertainty of x-values:').grid(row=row,column=col);  col+=1
        Entry(frame1,textvariable=self.eeu_x_uncertainty,bg='white').grid(row=row,column=col); col=1
        row+=1
        Label(frame1,text='Uncertainty of y-values:').grid(row=row,column=col); col+=1
        Entry(frame1,textvariable=self.eeu_y_uncertainty,bg='white').grid(row=row,column=col); col=1
        row+=1
        Label(frame1,text='Number of runs:').grid(row=row,column=col); col+=1
        Entry(frame1,textvariable=self.eeu_no_runs,bg='white').grid(row=row,column=col); col=1
        row+=1

        frame2 = Frame(self.eeu_win)
        frame2.pack()
        no_runs = IntVar()
        no_runs.set(1)
        self.current_run = IntVar()
        self.current_run.set(0)
        row=0; col=0
        modellabel = Label(frame2,text='Model: '+model, font="bold").grid(row=row,column=col)
        row+=1; col=1
        self.eeu_averages={}; self.eeu_stdevs={}
        Label(frame2,text='Average',foreground='blue').grid(row=row,column=col,columnspan=1); col+=1
        Label(frame2,text='Standard dev.',foreground='blue').grid(row=row,column=col,columnspan=1)
        row+=1; col=0
        #fill in average and standard deviations for each var in model
        for n in names:
            self.eeu_averages[n] = DoubleVar()
            self.eeu_averages[n].set(0.0)
            self.eeu_stdevs[n] = DoubleVar()
            self.eeu_stdevs[n].set(0.0)
            Label(frame2, text='%s'%n ).grid(row=row,column=0,sticky='w',columnspan=3)
            Label(frame2,textvariable=self.eeu_averages[n]).grid(row=row,column=1,sticky='w')
            Label(frame2,textvariable=self.eeu_stdevs[n]).grid(row=row,column=2,sticky='w')
            row+=1

        # Errors and standard deviations
        self.eeu_sumsqaverages = DoubleVar()
        self.eeu_sumsqstdevs = DoubleVar()
        col=0
        Label(frame2,text='Sum of sq errors:').grid(row=row,column=col,sticky='w',columnspan=3)
        col+=1
        Label(frame2,textvariable=self.eeu_sumsqaverages).grid(row=row,column=col,sticky='w',)
        col+=1
        Label(frame2,textvariable=self.eeu_sumsqstdevs).grid(row=row,column=col,sticky='w')
        row+=1;col=0
        Label(frame2,text='Run number ').grid(row=row,column=col); col+=1
        Label(frame2,textvariable=self.current_run).grid(row=row,column=col); col+=1
        Label(frame2,text=' of ').grid(row=row,column=col); col+=1
        Label(frame2,textvariable=self.eeu_no_runs).grid(row=row,column=col); col=1
        row+=1;

        # stop and quit buttons
        frame3 = Frame(self.eeu_win)
        frame3.pack()
        row=1; col=1
        Button(frame3, text='Start', relief=GROOVE, bg='#B0C4DE',
               command=self.startestimateUncertainty).grid(row=row,column=col)
        col+=1
        Button(frame3,text='Stop', relief=GROOVE, bg='#B0C4DE',
               command=self.stopestimateUncertainty).grid(row=row,column=col)
        col+=1
        def close():
            self.eeu_win.destroy()
        Button(frame3,text='Close', relief=GROOVE, bg='#B0C4DE',
               command=close).grid(row=row,column=col)

        return

    def stopestimateUncertainty(self):

        return

    def startestimateUncertainty(self):
        """start estimate uncertainty runs"""
        Fitting.estimateExpUncertainty(self.data, self.fitdata,
                                        xuncert=self.eeu_x_uncertainty.get(),
                                        yuncert=self.eeu_y_uncertainty.get(),
                                        runs=self.eeu_no_runs.get(),
                                        callback=self.updateEstimateUncertainty)
        return

    def updateEstimateUncertainty(self, r, stats):
        """callback for estimate uncertainty GUI"""
        self.current_run.set(r+1)
        sumsqaverages=0; sumsqstdevs=0
        for v in stats:
            self.eeu_averages[v].set('%5.4f' %stats[v][0])
            self.eeu_stdevs[v].set('%6.3e' %stats[v][1])
            sumsqaverages += pow(stats[v][0],2)
            sumsqstdevs += pow(stats[v][1],2)

        self.eeu_sumsqaverages.set('%5.4f' %math.sqrt(sumsqaverages))
        self.eeu_sumsqstdevs.set('%5.4f' %math.sqrt(sumsqstdevs))
        if hasattr(self,'eeu_win'):
            self.eeu_win.update_idletasks()
            self.eeu_win.update()
        return

    def previewFit(self):
        """Just update fit curve to reflect current paramenters in start vars"""
        vrs = self.getStartVars()
        model = self.model_type.get()        
        x,y = self.data.getxy()
        X=Fitting.getFitter(model, vrs=vrs, expdata=zip(x,y))
        self.plotter.updateFit(X)
        return

    def getStartVars(self):
        """Get start vars"""
        vrs=[]
        for i in range(len(self.fvars)):
            vrs.append(self.fvars[i][2].get())
        return vrs

    def send_fit_to_comments(self):
        return
    
class PlotPanel(Frame):
    """Plotter for ekin data using pylab - resuable"""
    def __init__(self, parent=None, width=400, height=100, side=TOP, tools=False):
        Frame.__init__(self, parent=None, width=400, height=100)

        self.parent = parent
        self.setupCanvas(side, tools)
        self.normalise = IntVar()
        self.normalise.set(0)
        self.show_legend = IntVar()
        self.show_legend.set(0)
        self.viewoption = IntVar()
        self.viewoption.set(0)
        self.options = None
        self.Opts = Options(redraw=self.plotCurrent) 
        self.selection = None
        self.m = None
        self.plotoption = 2
        self.datasets = None
        return

    def setupCanvas(self, side, tools):
        plt.close()
        from matplotlib import figure
        self.fig = figure.Figure(figsize=(5,4), dpi=80)       
        # create a tk.DrawingArea
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.parent)
        #self.canvas.show()
        self.canvas.get_tk_widget().pack(side=side, fill=BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=side, fill=BOTH, expand=1)
        mtoolbar = NavigationToolbar2TkAgg( self.canvas, self.parent )
        mtoolbar.update()
        if tools == True:
            self.addToolBar()
        return
    
    def addToolBar(self):
        """Extra toolbar with save button etc"""        
        t=Frame(self.parent)
        t.pack(fill=X)
        dpivar=IntVar()
        dpivar.set(300)
        Label(t, text='dpi:').pack(side=LEFT)
        Entry(t, textvariable=dpivar).pack(side=LEFT)
        b=Button(t, text='Save',width=5,command=lambda:self.saveFigure(dpivar.get()))
        b.pack(side=LEFT,padx=2,pady=2)
        b=Button(t, text='Options',command=self.plotOptions)
        b.pack(side=LEFT,padx=2,pady=2)
        b=Button(t, text='Annotate',command=self.annotate)
        b.pack(side=LEFT,padx=2,pady=2)        
        return
        
    def setProject(self, E):
        self.E = E
        if hasattr(self.E, '__plotopts__') and self.E.__plotopts__ != None:
            self.Opts = Options(redraw=self.plotCurrent, opts=self.E.__plotopts__)
            self.options = self.E.__plotopts__
        if hasattr(self.E, '__currentdataset__'):
            self.setCurrent(self.E.__currentdataset__)
        return

    def setCurrent(self, dataset):
        self.datasets = dataset
        return

    def plotCurrent(self, datasets=None, cols=0,
                    plotoption=None, options=None):
        """Plot current datapoints and fits"""
        
        if options==None:
            options = self.Opts.opts

        #remember for next time if we plot multiple - hacky  
        if datasets == None:
            datasets = self.datasets
        else:
            self.datasets = datasets
        if plotoption == None:
            plotoption = self.plotoption
        else:
            self.plotoption = plotoption
        if datasets == None: return

        #Note: TkAgg backend causes memory leak with repeated calls to plot
        self.fig.clear()
        #plt.close(self.fig)
        
        if type(datasets) is types.ListType and len(datasets) > 1:
            self.ax = self.E.plotDatasets(datasets, figure=self.fig, plotoption=plotoption, cols=cols,                               
                                **options)
        else:
            self.ax = self.E.plotDatasets(datasets, figure=self.fig, **options)
            self.addSelectionHandler()
        self.canvas.draw()
        #self.ax.hold(False)
        return

    def plotData(self, ekindata, fitdata):
        """Just plot supplied ekin data"""
        self.canvas.draw()
        return

    def updateFit(self, X):
        """Just update fit line"""
        self.E.updateFit(X)
        self.canvas.draw()
        return

    def updatePoints(self):
        """Update current points"""
        print 'update'
        self.E.updatePlot()
        self.canvas.draw()
        return

    def plotOptions(self):
        """Create options for plotting using Pylab class"""
        from Pylab import Options        
        if self.Opts == None:            
            self.Opts = Options(redraw=self.plotCurrent)            
        self.Opts.plotSetup()
        return

    def addSelectionHandler(self):
        """Add selection handler"""
        from Pylab import MouseMonitor
        if self.m != None:
            try:
                self.m.disconnect()
            except:
                pass
        self.m = MouseMonitor(self.ax, self)
        self.m.connect()
        return

    def saveFigure(self, dpi):
        """Save current figure"""
        import tkFileDialog, os
        filename=tkFileDialog.asksaveasfilename(parent=self.parent,
                                                defaultextension='.png',
                                                initialdir=os.getcwd(),
                                                filetypes=[("png","*.png"),
                                                           ("jpg","*.jpg"),
                                                           ("tiff","*.tif"),
                                                           ("All files","*.*")])        
        if not filename:
            return                 
        self.fig.savefig(filename, dpi=dpi)
        return    

    def annotate(self):
        """Basic adding of annotation to a plot"""

        axes = self.fig.get_axes()
        def addObject(objtype='text',text='test',x=1,y=1):
            
            bbox_props = dict(boxstyle="round", fc="#FFFC17", ec="0.4", alpha=0.8)
            if objtype == 'text':
                axes[0].text(x,y, text, ha="center", va="center", size=12, bbox=bbox_props)
            elif objtype == 'arrow':
                axes[0].annotate("",
                        xy=(x,y), xycoords='data',
                        xytext=(x,y), textcoords='data',
                        arrowprops=dict(arrowstyle="->", facecolor='black',
                                        connectionstyle="arc3"),
                        )
            self.canvas.draw()

        main = Toplevel()
        main.geometry('100x200+100+100')
        x=DoubleVar(); y=DoubleVar()
        txt=StringVar()
        objtype=StringVar()

        x=Pmw.EntryField(main,labelpos = 'w',
                label_text = 'x:',
                value = '0')
        x.pack()
        y=Pmw.EntryField(main,labelpos = 'w',
                label_text = 'y:',
                value = '0')
        y.pack()
        txt=Pmw.EntryField(main,labelpos = 'w',
                label_text = 'text:',
                value = '0')
        txt.pack()        
        Pmw.OptionMenu(main,labelpos = 'w',
                            label_text = 'type:',
                            menubutton_textvariable = objtype,
                            items = ['text','arrow'],
                            menubutton_width = 8).pack()
        Button(main,text='Add',command=lambda: addObject(objtype=objtype.get(),
                                                    text=txt.getvalue(),
                                                    x=x.getvalue(),y=y.getvalue())).pack()
        return

    
class DataPanel(Frame):
    """Display ekin data"""
    def __init__(self, parent):
        Frame.__init__(self, parent)
        #self.parent = parent
        self.create_tablebuttons()
        self.main = Frame(self, relief=RAISED)
        self.main.pack(side=RIGHT,fill=BOTH, expand=1)
        return

    def draw_Table(self, ekindata=None, fitmodel=None):
        """Draw a table for ekin datapoints"""
        self.currenttable = None
        self.fitmodel = fitmodel
        self.ekindata = ekindata

        self.currentmodel = EkinDataModel(ekindata)
        self.currenttable = EkinDataTable(self.main, model=self.currentmodel,
                                             width=250, height=220) 
        self.currenttable.loadPrefs()
        self.currenttable.createTableFrame()

        if self.currentmodel.getRowCount() == 0:
            self.currenttable.autoAdd_Rows(10)

        return

    def create_tablebuttons(self):
        def add_button(parent, name, callback, img=None, helptxt=None):
            if img==None and helptxt==None:
                b = Button(parent, text=name, command=callback,
                             relief=GROOVE)
            else:
                b = Button(parent, text=name, command=callback,
                                 relief=GROOVE,
                                 image=img)
            b.image = img
            b.pack(ipadx=1, ipady=1, padx=1,pady=1)
            if helptxt!=None:
                balloon=Pmw.Balloon(parent)
                balloon.bind(b, helptxt)
            return

        #add buttons
        self.buttonframe=Frame(self)
        #self.buttonframe.grid(row=0,column=0,sticky='ns')
        self.buttonframe.pack(side=LEFT)
        parent= self.buttonframe
        img = Table_images.add_row()
        hlptxt="Add a new row/record"
        add_button(parent,'Add record', self.add_Row, img, hlptxt)
        #img = Table_images.add_col()
        #hlptxt="Add a new column"
        #add_button(parent,'Add col', self.add_Column, img, hlptxt)
        img = Table_images.del_row()
        hlptxt="Delete selected row/record"
        add_button(parent,'Delete record', self.delete_Row, img, hlptxt)
        #img = Table_images.del_col()
        #hlptxt="Delete selected column"
        #add_button(parent,'Delete col', self.delete_Column, img, hlptxt)
        return

    def add_Row(self, name=None):
        """get next number and add row"""
        self.currenttable.add_Row()
        return

    def delete_Row(self):
        """get next number and add row"""
        self.currenttable.delete_Row()
        return

    def setCallback(self, func):
        self.currentmodel.setCallback(func)

    def addBinding(self, key, func):
        self.currenttable.bind_all(key, func)
        return

class ToolBar(Frame):
    """Uses the parent instance to provide the functions"""
    def __init__(self, parent=None, parentapp=None):
        Frame.__init__(self, parent, width=600, height=40)
        self.parentframe = parent
        self.parentapp = parentapp
        #add buttons
        img = Table_images.new_proj()
        #self.add_button('Plot Options', self.parentapp.plotframe.plotOptions)
        self.add_button('Specify Exp Conditions', self.parentapp.editMetaData,help='Specify experimental conditions')
        self.add_button('Show Comments', self.parentapp.addedit_comments, help='Show comments')
        #self.add_button('Close Window', self.parentapp.quit)

        return

    def add_button(self, name, callback, img=None, help=None, side=LEFT):
        """Add button"""
        if img==None:
            b = Button(self, text=name, command=callback,
                            relief=GROOVE,
                            bg='#B0C4DE')
        else:
            b = Button(self, text=name, command=callback,
                             relief=GROOVE,
                             image=img)
        b.image = img
        b.pack(side=side, fill=BOTH, padx=2, pady=2)
        if help != None:
            balloon=Pmw.Balloon(self.parentframe)
            balloon.bind(b,help)
        return

class MultipleValDialog(tkSimpleDialog.Dialog):
    """Simple dialog to get multiple values"""

    def __init__(self, parent, title=None, initialvalues=None, labels=None, types=None):
        if labels != None and types != NoneType:
            self.initialvalues = initialvalues
            self.labels = labels
            self.types = types
        tkSimpleDialog.Dialog.__init__(self, parent, title)

    def body(self, master):

        r=0
        self.vrs=[];self.entries=[]
        for i in range(len(self.labels)):
            Label(master, text=self.labels[i]).grid(row=r, column=0)
            self.vrs.append(IntVar())
            self.vrs[i].set(self.initialvalues[i])
            #self.vrs.set('')
            self.entries.append(Entry(master, textvariable=self.vrs[i], bg='white'))
            self.entries[i].grid(row=r, column=1,padx=2,pady=2,sticky='news')
            r+=1

        return self.entries[0] # initial focus

    def apply(self):
        self.result = True
        self.results = []
        for i in range(len(self.labels)):
            self.results.append(self.vrs[i].get())
        #print self.vrs
        return


def main():
    "Run the application"
    import sys
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="ekinfile",
                        help="Open an ekin proj", metavar="FILE")
    parser.add_option("-a", "--show all", dest="ekinfile",
                        help="Show all plots", metavar="BOOLEAN")
    opts, remainder = parser.parse_args()
    
    if opts.ekinfile != None:
        print opts.ekinfile
        app=EkinApp(project=opts.ekinfile)
        app.mainloop()
    else:
        app=EkinApp()
        app.mainloop()

    return

if __name__ == '__main__':
    main()

