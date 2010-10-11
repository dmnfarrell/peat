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

from Tkinter import *
import Pmw
import os
import tkSimpleDialog, tkFileDialog, tkMessageBox
from datetime import datetime
try:
    import matplotlib
    matplotlib.use('TkAgg')
    from matplotlib.font_manager import FontProperties
    import matplotlib.pyplot as plt   
except:
    pass
from Base import PDatabase
from PEATDB.Record import PEATRecord
from PEATTables import PEATTableModel
from PEATTables import PEATTable
from GUI_helper import *
from PEATDB.Dialogs import *
import PEAT_images as PEAT_images
import Table_images as Table_images
from Actions import DBActions
from Extfile import FileHandler
from Prefs import Preferences
from Plugins import Plugin

class App(Frame, GUI_help):
    """
    Main PEAT App for new ZODB-based, DB class is in Base.py
    Author: Damien Farrell, December 2009
    """
    def __init__(self, parent=None, DB=None):
        """Initialize the application."""
        self.parent=parent
        #If there is data to be loaded, show the dialog first
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
            self.peatinfo=None
        else:
            self.main=Toplevel()
        self.title = 'PEAT v2 alpha'

        self.MAIN_width=900; self.MAIN_height=600
        self.main.geometry('900x600+200+100')

        self.preferences=Preferences('PEAT',{})
        self.server = 'localhost'
        self.port = 8080
        self.username = self.getUserID()
        self.password = ''
        self.project = 'test'   #project name on server
        self.backend = 'mysql'

        self.filename = None       
        self.DB = None

        self.getPrefs()
        self.setupVars()
        self.setupGUI()
        self.doBindings()
        self.discoverPlugins()

        self.ekin_instances = {}
        self.DNAtool_state = None

        #after loading the database instance is self.DB
        if DB != None:
            self.loadDB(DB)

        self.setTitle()
        return

    def getUserID(self):
        """Try to fill in the user name"""
        try:
            import getpass
            return getpass.getuser()
        except:
            return ''

    def getPrefs(self):
        """Try to get some saved prefs"""
        items = {'username':'string','password':'string',
                 'blobdir':'string',
                 'memcachesize':'string',
                 'recent':'list','showDialogsinSidePane':'boolean',
                 'thumbsize':'string',
                 'molgraphApplication':'string','molgraphAppPath':'string'}
        for i in items:
            if not self.preferences.has_key(i):
                if items[i] == 'string':
                    self.preferences.set(i,'')
                elif items[i] == 'list':
                    self.preferences.set(i,[])
                elif items[i] == 'boolean':
                    self.preferences.set(i,True)
            self.__dict__[i] = self.preferences.get(i)

        return

    def discoverPlugins(self):
        """Discover plugins"""
        from Plugins import init_plugin_system, find_plugins, get_plugins_by_capability
        #peatpath = os.path.dirname(sys.modules['PEATDB'].__file__)
        peatpath = os.path.split(__file__)[0]
        homedir = os.path.expanduser("~")
        paths = [peatpath, os.getcwd(), homedir]
        pluginpaths = [os.path.join(p, 'plugins') for p in paths]
        print pluginpaths
        try:
            failed = init_plugin_system(pluginpaths)
        except:            
            return
        plgmenu = self.Pluginmenu['var']
        for plugin in get_plugins_by_capability('gui'):
            def func(p, **kwargs):
                def new():
                   p.main(**kwargs)
                return new
            plgmenu.add_command(label=plugin.menuentry,
                                command=func(plugin, parent=self))
        if len(failed)>0:
            for f in failed:
                self.recordEvent('failed to load plugin %s, error: %s' %(f[0],f[1]))
        return

    def updatePlugins(self):
        """Update plugins"""             
        self.Pluginmenu['var'].delete(2, self.Pluginmenu['var'].index(END))            
        self.discoverPlugins()
        return
        
    def setupVars(self):
        """Set some tk vars"""
        self.yasara=False
        self.pageviewvar = BooleanVar()
        self.pageviewvar.set(False)
        self.displaycachevar = BooleanVar()
        self.displaycachevar.set(False)
        return

    def setupGUI(self):
        """Do GUI elements"""
        self.addOptions(self.main)
        self.createMenuBar()
        self.canvas_x_size=self.MAIN_width-100
        self.canvas_y_size=self.MAIN_height-100
        self.canvas_x=1900
        self.canvas_y=20000
        self.canvas_border_x=250
        self.canvas_border_y=85

        self.masterframe = PanedWindow(self.main,
                                       orient=HORIZONTAL,
                                       sashwidth=3,
                                       showhandle=True)
        self.masterframe.grid(row=1,column=0,columnspan=2,rowspan=1,
                              sticky='news',
                              pady=2,ipady=2)

        self.createToolBar()
        self.createStatusPane()
        self.createInfoPane()

        self.main.rowconfigure(1,weight=1)
        self.main.columnconfigure(0,weight=1)
        self.cv = self.showBlankCanvas(self.main,row=1)
        self.welcomeLogo()

        return

    def createToolBar(self):
        self.toolbar = ToolBar(self.main, self)
        self.toolbar.grid(row=0,column=0,columnspan=1,rowspan=1,
                              sticky='news',
                              pady=2,ipady=2)

        img = PEAT_images.folder_database()
        hlptxt="Open DB from file"
        self.toolbar.add_button('Open local DB', self.openLocal, img, hlptxt)
        img = PEAT_images.connect()
        hlptxt="Connect to DB"
        self.toolbar.add_button('Connect to DB', self.connectDialog, img, hlptxt)
        img = PEAT_images.database_save()
        hlptxt="Save current changes to DB"
        self.toolbar.add_button('Save Changes', self.saveDB, img, hlptxt)
        img = PEAT_images.table_refresh()
        hlptxt="Refresh table"
        self.toolbar.add_button('Refresh', self.updateView, img, hlptxt)
        img = PEAT_images.table_row_insert()
        hlptxt="Add a new row/record"
        self.toolbar.add_button('Add record', self.addRecordDialog, img, hlptxt)
        img = PEAT_images.table_row_delete()
        hlptxt="Delete selected row/record"
        self.toolbar.add_button('Delete record', self.deleteRecord, img, hlptxt)
        img = PEAT_images.table_add()
        hlptxt="Add a new column"
        self.toolbar.add_button('Add col', self.addFieldDialog, img, hlptxt)
        img = PEAT_images.table_delete()
        hlptxt="Delete selected column"
        self.toolbar.add_button('Delete col', self.deleteField, img, hlptxt)
        img = PEAT_images.arrow_undo()
        hlptxt="Undo current changes (since last save)"
        self.toolbar.add_button('Undo', self.undo, img, hlptxt)
        img = PEAT_images.magnifier()
        hlptxt="Search"
        self.toolbar.add_button('Filter Records', self.showFilteringBar, img, hlptxt)
        #img=PEAT_images.DNAtool()
        img = None
        hlptxt="DNAtool: sequence design"
        self.toolbar.add_button('DNAtool', self.startDNAtool, img, hlptxt)
        hlptxt="Ekin: curve fitting"
        self.toolbar.add_button('Ekin', self.startEkin, img, hlptxt)
        hlptxt="Labbook: tabulated shared data"
        self.toolbar.add_button('Labbook', lambda: self.startLabbook(protein='ALL'), img, hlptxt)
        self.toolbar.add_seperator()
        return

    def createStatusPane(self):
        """Small status pane for db specific info"""
        statusframe = Frame(self.main)
        statusframe.grid(row=0,column=1,columnspan=1,rowspan=1,
                              sticky='news', padx=2)
        self.sizevar=StringVar()
        self.numrecsvar=StringVar()
        self.storagetypevar=StringVar()
        self.updateStatusPane()
        Label(statusframe, text='Size:',fg='gray50').pack(side=LEFT)
        l=Label(statusframe, textvariable=self.sizevar, fg='darkblue')
        l.pack(side=LEFT)
        Label(statusframe, text='Recs:',fg='gray50').pack(side=LEFT)
        l=Label(statusframe, textvariable=self.numrecsvar, fg='darkblue')
        l.pack(side=LEFT)
        Label(statusframe, text='Storage:',fg='gray50').pack(side=LEFT)
        l=Label(statusframe, textvariable=self.storagetypevar, fg='darkblue')
        l.pack(side=LEFT)
        return

    def updateStatusPane(self):
        """Update status pane"""
        if self.DB == None:
            self.sizevar.set(0)
            self.numrecsvar.set(0)
            self.storagetypevar.set('')
        else:
            self.sizevar.set(str(self.DB.getSize())+' KB')
            self.numrecsvar.set(self.DB.length())
            self.storagetypevar.set(self.DB.storagetype)
        return

    def createInfoPane(self):
        """Make info panel for logging etc"""
        self.eventlog = Pmw.ScrolledText(self.main,
                                          labelpos = 'nw',
                                          label_text='Messages',
                                          usehullsize = 1,
                                          hull_width=self.MAIN_width,
                                          hull_height=80)
        self.eventlog.component('text').configure(background='#FFFAF0')
        self.eventlog.grid(row=3,column=0,columnspan=2,
                              sticky='news', pady=2,ipady=2)
        self.recordEvent('Welcome to PEAT')
        return

    def createSidePane(self, width=20):
        """Side panel for various dialogs"""
        if hasattr(self, 'sidepane'):
            self.masterframe.forget(self.sidepane)
            self.sidepane.destroy()
        self.sidepane = Frame(self.masterframe, height=self.MAIN_height,
                              bd=1,relief=RAISED)
        self.masterframe.paneconfigure(self.sidepane,
                                       width=width,
                                       sticky='news')
        self.masterframe.add(self.sidepane)
        return self.sidepane

    def resetSidePane(self, width=20):
        self.createSidePane(width=width)
        if hasattr(self,'tableframe'):
            self.masterframe.paneconfigure(self.sidepane, before=self.tableframe)
        return

    def createMenuBar(self):
        """Create the menu bar for the application. """
        self.menu=Menu(self.main)
        self.file_menu={ '01Connect to DB':{'cmd':self.connectDialog},
                         '02User Setup':{'cmd':self.loginSetup},
                         '03Open local DB':{'cmd':self.openLocal},
                         '04New DB':{'cmd':self.newDB},
                         '05Close DB':{'cmd':self.closeDB},
                         '06Save Changes':{'cmd':self.saveDB},
                         '07Save a Copy As':{'cmd':self.saveAs},
                         '08Create DB on Server':{'cmd':self.createServerProject},
                         '09Undo Current Changes':{'cmd':self.undo},
                         '10Undo Previous Commits':{'cmd':self.showUndoLog},
                         '11Quit':{'cmd':self.quit}}
        self.file_menu=self.create_pulldown(self.menu,self.file_menu)
        self.menu.add_cascade(label='Database',menu=self.file_menu['var'])
        self.addPrevProjects(self.file_menu['var'])

        self.rec_menu={ '01Add Protein':{'cmd':self.addRecordDialog},
                        '02Delete Protein':{'cmd':self.deleteRecord},
                        '03Add Field':{'cmd':self.addFieldDialog},
                        '04Delete Field':{'cmd':self.deleteField},
                        '05Show/Hide Fields':{'cmd':self.hideFieldsDialog},
                        '06Add/Change PDB File':{'cmd':self.addPDBFile},
                        '07Add Mutant':{'cmd':self.addMutantDialog},
                        '08sep':'',
                        '09Find in table':{'cmd':self.createSearchBar},
                        '10Filter Records':{'cmd':self.showFilteringBar},
                        '11Advanced Search':{'cmd':self.advancedSearch}}
        self.rec_menu=self.create_pulldown(self.menu,self.rec_menu)
        self.menu.add_cascade(label='Records',menu=self.rec_menu['var'])

        self.IO_menu={'01Import from Text/CSV':{'cmd':self.importCSV},
                      #'02Import Mutants':{'cmd':self.importMutants},
                      #'02Export to csv file':{'cmd':self.exportCSV},
                      '03Import External/Binary files':{'cmd':self.importFileset},
                      '04Export Binary Files':{'cmd':self.exportExtFiles}}
        self.IO_menu=self.create_pulldown(self.menu,self.IO_menu)
        self.menu.add_cascade(label='Import/Export',menu=self.IO_menu['var'])

        self.view_menu={ '01Show Change Log':{'cmd':self.showUndoLog},
                         '02Show Changed Items':{'cmd':self.showChanged},
                         '03Show Binary files':{'cmd':self.showExtFiles},
                         '04New Table View':{'cmd':self.createTableView}}
        self.view_menu=self.create_pulldown(self.menu,self.view_menu)
        self.view_menu['var'].add_checkbutton(label="Use Page View", variable=self.pageviewvar,
                                       command=self.pageView)
        self.menu.add_cascade(label='View',menu=self.view_menu['var'])

        self.settings_menu={ '01General':{'cmd':self.showSettings},
                             '02DB Specific Settings':{'cmd':self.editDBInfo},
                             '03Table Prefs':{'cmd':self.showtablePrefs},
                             '04sep':'',
                             '05Inspect DB Meta':{'cmd':self.inspectDBMeta},
                             '06Pack Database':{'cmd':self.packDB},
                             '07Testing':{'cmd':self.testsDialog}, }
        self.settings_menu=self.create_pulldown(self.menu,self.settings_menu)
        self.menu.add_cascade(label='Settings',menu=self.settings_menu['var'])

        self.tools_menu={'01Compare Records':{'cmd':self.startComparator},
                     '02Set Reference Protein':{'cmd':self.setRefProtein},
                     '03Remodel Mutant Structures':{'cmd':self.remodelPEATModels},
                     '04Fetch PDB':{'cmd':self.fetchPDB},
                     '05Update Mutation Codes':{'cmd':self.updateMutations},
                     '06Update Mutant Sequences':{'cmd':self.updateAASequences}}                     
        self.tools_menu=self.create_pulldown(self.menu,self.tools_menu)
        self.menu.add_cascade(label='Tools',menu=self.tools_menu['var'])        
        
        self.labbookmenu = Menu(self.tools_menu['var'], tearoff=1)
        self.tools_menu['var'].add_cascade(label='Labbooks', menu=self.labbookmenu)
       
        self.Pluginmenu={'01Update Plugins':{'cmd':self.updatePlugins},
                         '02sep':'',}
        self.Pluginmenu=self.create_pulldown(self.menu,self.Pluginmenu)
        self.menu.add_cascade(label='Plugins',menu=self.Pluginmenu['var'])

        self.help_menu=Menu(self.menu,tearoff=0)
        self.help_menu.add_command(label='Online documentation',command=self.onlineDocumentation)
        self.help_menu.add_command(label='Report Bug',command=self.gotoBugzilla)
        self.help_menu.add_command(label="About PEATDB",command=self.aboutPEAT)
        self.menu.add_cascade(label='Help',menu=self.help_menu)

        self.main.config(menu=self.menu)
        return

    def updateMenus(self, load=False, close=False):
        """Do whatever changes to menus required wehen we load and close"""
        if load==True:
            self.updateLabbooksMenu()
        return

    def doBindings(self):
        #self.master.bind("<Configure>", self.resize)
        self.master.bind("<Control R>", self.updateTable)
        self.main.protocol('WM_DELETE_WINDOW', self.quit)

        return

    def createSearchBar(self, event=None):
        """Add a find entry box"""
        if self.DB == None:
            return
        frame = Frame(self.main)
        row=0
        def close():
            frame.destroy()
        self.findtext=StringVar()
        self.findbox=Entry(frame,textvariable=self.findtext,width=30,bg='white')
        self.findbox.grid(row=row,column=1,sticky='news',columnspan=2,padx=2,pady=2)
        self.findbox.bind('<Return>',self.dofindText)
        Label(frame,text='Find:').grid(row=row,column=0,sticky='ew')
        self.findagainbutton=Button(frame,text='Find Again', command=self.dofindAgain)
        self.findagainbutton.grid(row=row,column=3,sticky='news',padx=2,pady=2)
        self.cbutton=Button(frame,text='Close', command=close)
        self.cbutton.grid(row=row,column=4,sticky='news',padx=2,pady=2)
        frame.grid(row=2,column=0)
        return

    def showFilteringBar(self):
        """Add filter bar"""
        if self.DB == None:
            return
        #supply any non ekin fields
        fields = []
        model = self.table.getModel()
        for f in model.columnNames:            
            if model.columntypes[f] in PEATRecord.ekintypes:
                pass
            else:
                fields.append(f)        
        self.filterframe = self.table.createFilteringBar(self.main, fields=fields)
        self.filterframe.grid(row=2,column=0)
        return
    
    def addOptions(self, root):
        root.option_add('*Button*padx', 4)
        root.option_add('*Button*pady', 4)
        return

    def addPrevProjects(self, parentmenu):
        """Add previous projects from prefs file to the DB menu"""
        import string, types
        self.recentmenu = submenu = Menu(parentmenu, tearoff=1)
        for s in self.preferences.get('recent'):
            if s and type(s) is types.TupleType:
                server,port,project,backend=s
                text=server+':'+str(port)+':'+project
                def func(**kwargs):
                    def new():
                       self.connect(**kwargs)
                    return new
                submenu.add_command(label=text,
                                    command=func(server=server,port=port,
                                                 project=project,backend=backend))
            elif type(s) is types.StringType and os.path.exists(s):
                if len(s)>30: text = s[:10]+'..'+s[20:]
                else: text = s
                def func(arg):
                    def new():
                        self.openLocal(arg)
                    return new
                submenu.add_command(label=text, command=func(s))
        parentmenu.add_separator()
        parentmenu.add_command(label='Clear Recent', command=self.clearPrevProjects)
        parentmenu.add_cascade(label='Recent Projects', menu=submenu)
        return

    def clearPrevProjects(self):
        self.preferences.set('recent', [])
        self.recentmenu.delete(0, self.recentmenu.index('last'))
        return
    
    def connect(self, server=None, port=None, project=None, backend=None):
        """Connect with no dialog"""

        if server != None:
            self.server = server; self.project = project
            self.backend = backend

            DB = PDatabase(server=server, port=port,
                      username=self.username,
                      password=self.password,
                      project=project,
                      backend=self.backend,
                      blob_dir=self.preferences.get('blobdir'))
            if DB.connection == None:
                tkMessageBox.showwarning("Connection Error",
                         'Could not connect to %s:%s.\n'
                         'Error message returned: %s'
                          %(project,server,DB.errormessage))

                return
            self.closeDB()
            self.recordEvent('Connected to '+server)
            self.loadDB(DB)

        return

    def connectDialog(self, event=None):
        """Get a server/port, storage and open DB"""
        backends = PDatabase.backends        
        mpDlg = MultipleValDialog(title='Connect to DB',
                                    initialvalues=(self.username, self.password, self.server,
                                                   self.port, self.project,
                                                   backends),
                                    labels=('username','password','server','port','project','backend'),
                                    types=('string','password','string','int','string','list'),
                                    parent=self.main)

        if mpDlg.result == True:
            self.username = mpDlg.results[0]
            self.password = mpDlg.results[1]
            self.server = mpDlg.results[2]
            self.port = mpDlg.results[3]
            self.project = mpDlg.results[4]
            self.backend = mpDlg.results[5]
        else:
            return
        if self.server != None:
            DB = PDatabase(server=self.server, port=self.port,
                      username=self.username,
                      password=self.password,
                      project=self.project,
                      backend=self.backend)

            if DB.connection == None:
                tkMessageBox.showwarning("Connection Error",
                         'Could not connect to %s:%s.\n'
                         'Error message returned: %s'
                          %(self.project,self.server,DB.errormessage))
                self.recordEvent('Connection failure to '+self.server)
                return
            #We can try to close current db now
            answer = self.closeDB()
            if answer == 'cancel':
                return
            self.recordEvent('Connected to '+self.server)
            self.loadDB(DB)
            self.addtoRecent((self.server,self.port,self.project,self.backend))

        return

    def addtoRecent(self, value):
        """"Store loaded project info in prefs"""
        prev = self.preferences.get('recent')
        if not value in prev:
            prev.append(value)
        self.preferences.set('recent', prev)
        return

    def loadDB(self, DB):
        """Load peat db and display it in a table"""
        #createtableframe should be optimised for reading large dbs
        self.DB = DB
        self.removeBlankCanvas()
        self.createSidePane()
        self.tablemodel = PEATTableModel(DB)
        self.tableframe = Frame(self.masterframe, bd=1,relief=RAISED)
        self.table = PEATTable(self.tableframe, self.tablemodel, parentapp=self)
        self.table.loadPrefs(self.preferences)
        self.table.createTableFrame()
        if self.preferences.get('thumbsize') == '':
            self.preferences.set('thumbsize',200)
        self.table.thumbsize = int(self.preferences.get('thumbsize'))
        #we are using a panedwindow
        self.masterframe.paneconfigure(self.tableframe,
                                       width=850,
                                       sticky='news')
        self.masterframe.add(self.tableframe)
        self.setTitle()
        self.recordEvent('Loaded db ok')
        self.updateStatusPane()
        self.updateMenus(load=True)
        return

    def setTitle(self):
        """set title of window"""
        if self.DB != None:
            if self.filename != None:
                self.main.title(self.title+': '+self.filename)
            else:
                self.main.title(self.title+': '+self.server+' '+str(self.project))
        else:
            self.main.title(self.title)
        return

    def newDB(self):
        """Create a new local DB"""
        if self.DB != None:
            answer=self.closeDB()
            if answer == 'cancel':
                return
        filename=tkFileDialog.asksaveasfilename(defaultextension='.fs',
                                   initialdir=os.getcwd(),
                                   filetypes=[("zodb fs","*.fs"),("All files","*.*")])
        if filename:
            if os.path.exists(filename):
                for i in ['.lock','.index','']:
                    os.remove(filename+i)
            DB = PDatabase(local=filename)
            self.loadDB(DB)
        return

    def openLocal(self, filename=None):
        """Open local DB"""
        print 'local file'
        print filename
        if filename == None:
            filename=tkFileDialog.askopenfilename(defaultextension='.fs',
                                       initialdir=os.getcwd(),
                                       filetypes=[("zodb fs","*.fs"),("All files","*.*")])
        if filename:
            if self.DB != None:
                answer = self.closeDB()
                if answer == 'cancel':
                    return
            self.filename = filename
            DB = PDatabase(local=filename)
            self.loadDB(DB)
            self.addtoRecent(filename)
        return

    def saveDB(self, event=None):
        """Save current changes, conflicts are handled here"""

        if self.DB.isChanged() == False:
            return
        #cache changed objects in case of conflict
        import copy
        #saved = copy.deepcopy(self.DB.getChanged())
        
        if self.username == '':
            tkMessageBox.showwarning("Set a user name",
                                     'You need to set a user name before saving.\n'
                                     'So that PEAT can log who made the change.')
            return

        def handleconflict(result):
            '''handle this inside db class??'''
            dialog.deactivate()
            if result == 'Discard':
                self.DB.abort()
                #self.DB.connection.sync()
                self.updateTable()
            elif result == 'Overwrite':
                #first cache changed objects, abort
                #and update then we try to apply users changes
                self.DB.abort()
                #self.DB.connection.sync()
                #reapply changes..
                for s in saved:
                    print s
                #result = self.DB.commit(user=self.username, note=comment)
                #if result == False:
                    #failed.. give up

            return

        askcomment = self.preferences.get('promptforlogcomments')
        comment = None
        if self.DB != None:
            if askcomment == True:
                comment = tkSimpleDialog.askstring('Log comment',
                                             'Enter a comment for the log',
                                             initialvalue='minor',
                                             parent=self.main)
            if comment == None:
                return
            result = self.DB.commit(user=self.username, note=comment)
            if result == False:
                #a conflict occured, we give user a choice
                dialog = Pmw.Dialog(self.main,
                          buttons = ('Discard My Changes', 'Overwrite', 'Cancel'),
                          defaultbutton = 'OK',
                          title = 'Conflict Error!',
                          command = handleconflict)
                w = Label(dialog.interior(),
                      text = 'A conflict occured, someone else has changed \n'
                             'the same record/field since you last synced the DB.\n'
                             'You can discard your changes OR overwrite with \n'
                             'your changes. Pressing cancel will take no action.',
                      background = '#FFFAF0',
                      width=50,
                      pady=10)
                w.pack(expand=1, fill=BOTH, padx=2, pady=2)
                dialog.activate(geometry='centerscreenalways')
        else:
            #update any synced changes if using remote zeo
            #self.DB.connection.sync()
            self.updateTable()
        self.updateStatusPane()
        self.updateChangelog()
        return

    def saveAs(self):
        """Save a local copy of the DB, if the DB is remote, this acts as a local
           backup, but may take some time to get the data from the server, need to
           warn user before we start"""
        if self.DB == None:
            return
        size = round(self.DB.getSize()/1024,1)
        if size > 30:
            answer = tkMessageBox.askyesno("DB Save Warning",
                                       'This Database is %s MB in size.\n'
                                       'Are you sure you want to save locally?' %size)
            if answer == False:
                return

        import copy
        filename=tkFileDialog.asksaveasfilename(defaultextension='.fs',
                                   initialdir=os.getcwd(),
                                   filetypes=[("zodb fs","*.fs"),("All files","*.*")])
        if filename:
            if os.path.exists(filename):
                for i in ['.lock','.index','']:
                    os.remove(filename+i)
            newDB = PDatabase(local=filename)
            #show progress dialog
            from Dialogs import PEATDialog
            pb=PEATDialog(self.master, option='progressbar',
                                          message='Getting records..')
            pb.update_progress(0)
            total = len(self.DB.getRecs())
            count=0
            for r in self.DB.getRecs():
                rec = copy.deepcopy(self.DB.data[r])
                newDB.add(r, record=rec)
                count+=1
                pb.update_progress(float(count)/float(total)*100.0)
            pb.close()
            for m in self.DB.meta.keys():
                newDB.meta[m] = copy.deepcopy(self.DB.meta[m])
            #newDB.cellcache = copy.deepcopy(self.DB.cellcache)
            newDB.commit()
            newDB.close()
            self.recordEvent('Saved current DB as %s' %filename)
        return

    def closeDB(self, event=None):
        """Close the open DB"""
        if self.DB == None:
            return
        if self.DB.isChanged() == True:
            answer = tkMessageBox.askquestion("Save recent changes?",
                                              "Save pending changes before closing?",
                                              icon=tkMessageBox.QUESTION,
                                              type=tkMessageBox.YESNOCANCEL)            
            if str(answer) == 'cancel':                
                return 'cancel'
            elif str(answer) == 'yes':
                self.saveDB()

        self.DB.close()
        self.DB = None
        self.removeTable()        
        self.recordEvent('Closed db')
        self.filename = None
        self.setTitle()
        self.cv = self.showBlankCanvas(self.main,row=1)
        self.welcomeLogo()
        self.updateStatusPane()
        return True

    def packDB(self):
        if self.DB == None:
            return
        days = tkSimpleDialog.askinteger('Pack DB',
                                 'Reduces the size of the DB by removing old\n'
                                 'copies. This step cannot be reversed!\nNo. of days back to keep:',
                                 initialvalue=7,
                                 parent=self.master)

        if days:
            self.DB.pack(days=days)
            self.updateStatusPane()
            self.updateChangelog()
        return

    def addRecordDialog(self, event=None):
        """Add a new protein to our database"""
        if self.DB == None:
            return
        self.added=False

        def closewin(event=None):
            enterwin.destroy()
            return

        def add():
            self.added = self.addRecord(new_name.get())
            if self.added == False:
                tkMessageBox.showwarning("Key present",
                                         "A record with this key is already in the DB")
            else:
                enterwin.destroy()

            return

        import pdb
        enterwin=Toplevel()
        enterwin.title('Choose data')

        top=self.master.winfo_toplevel()
        rootx=top.winfo_rootx()
        rooty=top.winfo_rooty()
        enterwin.geometry('+%d+%d' %(rootx+100,rooty+100))
        enterwin.transient(self.master)
        enterwin.grab_set()

        row=1
        l=Label(enterwin,text='What type of data will you add?')
        l.grid(row=row,column=0,sticky='news',columnspan=2)
        options=[('DNA sequence','DNA'),
                 ('Protein sequence','Prot'),
                 ('PDB file','PDB'),
                 ('Name only','Name')]
        enter_choice=StringVar()
        enter_choice.set('Name')
        row=row+1
        for text,mod in options:
            b=Radiobutton(enterwin,text=text,variable=enter_choice,value=mod)
            b.grid(row=row,column=0,sticky='w')
            row=row+1

        l=Label(enterwin, text="Enter name")
        l.grid(row=row,column=0, sticky='news')
        new_name=StringVar()
        e=Entry(enterwin,textvariable=new_name)
        e.grid(row=row,column=1,sticky='news')

        row=row+1
        b=Button(enterwin,text='Cancel',command=closewin)
        b.grid(row=row,column=1,sticky='news')
        b2=Button(enterwin,text='Enter',command=add)
        b2.grid(row=row,column=0,sticky='news')

        e.focus_set()
        self.master.wait_window(enterwin)
        if self.added == True:
            self.addRecordAction(new_name.get(), action=enter_choice.get())
        return

    def addRecordAction(self, name, action=None):
        """Additional actions when adding a new rec"""
        if action == None:
            return
        if action == 'DNA':
            tkMessageBox.showinfo("Add DNA sequence",
                                  'I will now start DNAtool. '
                                  'Load your DNA sequence, select the open '
                                  'reading frame and press "Return to DataBase" to continue')
            done = False
            while not done:
                self.startDNAtool()
                DNAdata=self.DNAtool_state
                if DNAdata.has_key('ORF_selected'):
                    done=1
                else:
                    if not tkMessageBox.askokcancel("ORF not selected",
                                             "You did not select an ORF.\nReturn to DNAtool?"):
                        return

        elif action == 'Prot':
            DBActions.addProteinSeq(self.DB, name)

        elif action == 'PDB':
            DBActions.addPDBFile(self.DB, name)

        return

    def addRecord(self, name):
        """Add a record"""
        result = self.DB.add(name)
        if result == False:
            return result

        self.updateTable()
        self.table.setSelectedRow(self.tablemodel.getRecordIndex(name))
        self.updateTable()
        return True

    def deleteRecord(self, event=None):
        if getattr(self,'table',None):
            self.sel_prots = self.table.get_selectedRecordNames()

        answer = tkMessageBox.askyesno('Delete records?',
                                 'This will delete all information on all'
                                 'selected records\n\nDo you want to continue?.',
                                 parent=self.master)
        if answer:
            for protein in self.sel_prots:
                self.DB.delete(protein)
            self.sel_prots = None
            self.updateTable()

        return

    def addMutantDialog(self):
        """Dialog for adding a mutant"""
        from Mutate import PointMutate
        p = PointMutate()
        name = self.table.get_selectedRecordNames()[0]
        mframe = self.createChildFrame()
        p.addMutant(self.DB, name, parentframe=mframe)
        return

    def addFieldDialog(self):
        """Get a description of the field and add it"""
        if self.DB == None:
            return

        def close(event=None):
            dataf_spec.destroy()
            return

        def add():
            self.DB.addField(name=name_var.var.get(), fieldtype=type_var.var.get())
            dataf_spec.destroy()
            self.updateTable()
            return

        dataf_spec=Toplevel()
        dataf_spec.title('Specify type of data field')
        top=self.master.winfo_toplevel()
        rootx=top.winfo_rootx()
        rooty=top.winfo_rooty()
        dataf_spec.geometry('+%d+%d' %(rootx+100,rooty+100))
        dataf_spec.transient(self.master)
        dataf_spec.grab_set() # Grab all input

        # Get the name
        exclude=['Ekintype']
        self.column_types=['text']+self.table.peatactions.keys()+self.table.ekin_actions.keys()
        self.column_types.remove('Ekintype')
        name_var=entry_field(dataf_spec,row=1,column=0,name='Name',ftype='text')
        type_var=entry_field(dataf_spec,
                                  row=2,
                                  column=0,
                                  name='Data field type',
                                  ftype='menu',
                                  items=self.column_types,
                                  default_item='text')
        self.default_var=entry_field(dataf_spec,row=3,column=0,name='Default value',ftype='text')
        Button(dataf_spec,
               command=add,
               text='Add data field',fg='red').grid(row=5,column=0,columnspan=1,
                                                    padx=2,pady=2,sticky='news')
        Button(dataf_spec,
               command=close,
               text='Cancel',fg='green').grid(row=5,column=1,columnspan=1,
                                              padx=2,pady=2,sticky='news')


    def addField(self):
        """Add a new data field to the database"""
        return

    def deleteField(self):
        """Delete selected field"""
        col = self.table.getSelectedColumn()
        colname = self.tablemodel.getColumnName(col)
        if colname in self.DB.meta.staticfields:
            tkMessageBox.showinfo('Cannot delete',
                                   "This field can't be removed.\n"
                                   'You may mark it as hidden.',
                                   parent=self.master)
            return
        ans =  tkMessageBox.askyesno("Delete",  "Delete This Column?")
        if ans:
            #if self.DB.length() > 20:
            from Dialogs import PEATDialog
            pb=PEATDialog(self.master, option='progressbar',
                                          message='Deleting this field for all records..')
            pb.update_progress(0)            
            callback = pb.update_progress             
            self.DB.deleteField(colname, callback=callback)
            self.DB.data._p_jar.cacheGC()
            self.updateTable()
            self.updateStatusPane()
            if pb:
               pb.close()
        return

    def undo(self):
        """Undo all changes since last commit to DB"""
        if self.DB == None or self.DB.isChanged() == False:
            return
        self.DB.abort()
        self.updateTable()
        return

    def createChildFrame(self, width=400, title=''):
        if self.showDialogsinSidePane == True:
            self.resetSidePane(width=width)
            cframe = Frame(master=self.sidepane)
            cframe.pack(fill=BOTH,expand=1)
        else:
            cframe = Toplevel()
            cframe.geometry('+100+450')
            cframe.title(title)
        return cframe

    def showChanged(self):
        """Show list of changed records"""
        if self.DB == None:
            return

        cframe = self.createChildFrame()  #we should make sidepane a class

        def close():
            cframe.destroy()
            self.resetSidePane(width=20)
            return
        def update():
            items = self.DB.getChanged()
            try:
                self.chgeditemswin.destroy()
            except:
                pass
            self.chgeditemswin = Pmw.ScrolledText(cframe,
                                           labelpos = 'nw',
                                           label_text='Changed Items Since last Save',
                                           usehullsize = 1,
                                           hull_width=self.MAIN_width,
                                           hull_height=200)
            self.chgeditemswin.pack(fill=BOTH, side=BOTTOM, padx=2,pady=2)
            for i in items:
                self.chgeditemswin.insert(END, i+'\n')
            return
        Button(cframe,text='Update',command=update).pack(fill=BOTH, side=BOTTOM, padx=2,pady=2)
        Button(cframe,text='Close',command=close).pack(fill=BOTH, side=BOTTOM, padx=2,pady=2)
        update()
        return

    def showUndoLog(self):
        """Show db commit log"""
        if self.DB == None:
            return
        from datetime import datetime
        self.changelogframe = self.createChildFrame()
        cframe = self.changelogframe
        def close():
            cframe.destroy()
            self.changelogframe = None
            self.resetSidePane(width=20)
            return
        def errormessage():
            tkMessageBox.showwarning("Undo Error",
                  'Could not undo this transaction.\n'
                  'Later commits probably altered the same data.')
        def tryUndo():
            i = self.ltable.get_selectedRecordNames()[0]
            undoid = self.ulogs[i]['id']
            comment = tkSimpleDialog.askstring('Log comment',
                                               'Enter a comment for the log',
                                                initialvalue='',
                                                parent=self.main)
            self.DB.undo(id=undoid, user=self.username,
                        note=comment,
                        callback=errormessage)
            self.updateTable()
            self.recordEvent('Undid transaction %s' %undoid)
            self.updateChangelog()

        if self.DB.supportsUndo():
            Button(cframe,text='Update',command=self.updateChangelog).grid(row=3,column=0,columnspan=3,
                                                         sticky='news')
        Button(cframe,text='Undo Selected',command=tryUndo).grid(row=4,column=0,columnspan=3,
                                                           sticky='news')
        Button(cframe,text='Close',command=close).grid(row=5,column=0,columnspan=3,
                                                           sticky='news')
        w='Use the undo feature carefully.\nYou should be aware of what data changes\n'
        w+='are being undone. Provide a comment\nwith the undo to indicate what you reverted.'
        Label(cframe,text=w,bg='lightyellow').grid(row=6,column=0,columnspan=3,sticky='news')
        self.updateChangelog()
        return

    def updateChangelog(self):

        if not hasattr(self, 'changelogframe') or self.changelogframe == None:
            return
        try:
            self.changelog.destroy()
        except:
            pass
        #show logs in a table
        from Tables import TableCanvas
        self.ulogs={}; l=self.DB.undoLog()
        i=0
        for r in range(len(l)-1,-1,-1):
            l[r]['time'] = self.formatTime(l[r]['time'])
            l[r]['user_name'] = l[r]['user_name'].strip('/ ')
            self.ulogs[i] = l[r]; i+=1

        self.ltable = TableCanvas(self.changelogframe, newdict=self.ulogs, namefield='rev',
                                  cellwidth=50, cellbackgr='#E3F6CE',
                                  thefont="Arial 10",rowheight=16, editable=False,
                                  rowselectedcolor='yellow',reverseorder=1)
        self.ltable.createTableFrame()
        return

    def formatTime(self, s):
        x = datetime.fromtimestamp(s).strftime("%d-%m-%Y %H:%M:%S")
        return x

    def updateView(self):
        """Update the table if no changes pending"""
        if self.DB==None:
            return
        if self.DB.isChanged() == True:
            tkMessageBox.showinfo('Pending changes',
                                   'You have pending changes and should save these first.',
                                   parent=self.master)
            return
        self.DB.connection.sync()
        self.updateTable()
        self.updateStatusPane()
        self.updateChangelog()
        return

    def get_field_label(self, field):
        """Return the field label"""
        model = self.table.getModel()
        return model.columnlabels[field]

    def updateTable(self, protein=None, field_name=None):
        """To be called whenever a change to DB needs to be reflected in
            table. """
        if self.DB == None:
            return
        DB=self.DB

        if hasattr(self, 'table') and self.table != None:
            model = self.tablemodel            
            #Update the model to reflect changes in DB and redraw the table
            if self.table.rows < 10000:
                sortcol = self.table.sortcol
            else:
                sortcol = None
            model.update_reclist(DB.data.keys(), self.table.sortcol)
            model.update_columnNames()
            model.update_colors()
            self.table.redrawTable()
        return

    def removeTable(self):
        if hasattr(self.table, 'filterframe') and self.filterframe != None:
            self.table.filterframe.destroy()
        self.table = None
        try:
            if self.tableframe != None:
                self.tableframe.destroy()
                self.tableframe = None
        except:
            print 'no tableframe'
        return

    def createTableView(self):
        """Create a new table view"""
        if self.DB == None or not hasattr(self, 'table'):
            return
        DB=self.DB
        newframe = Toplevel()         
        tableframe = Frame(newframe, bd=1,relief=RAISED)
        tableframe.pack(fill=BOTH,expand=1)
        table = PEATTable(tableframe, self.tablemodel, parentapp=self)
        table.loadPrefs(self.preferences)
        table.createTableFrame()
        if self.preferences.get('thumbsize') == '':
            self.preferences.set('thumbsize',200)
        table.thumbsize = int(self.preferences.get('thumbsize'))                           
        return
    
    def pageView(self, event=None):
        if self.pageviewvar.get() == False:
            self.table.paging = 0
        else:
            self.table.paging = 1
        self.table.redrawTable()
        return

    def resize(self,event):
        """Update the scrollbars"""
        if event.widget==self.master:
            if event.width>100 and event.height>100:
                self.masterframe.configure(hull_width=event.width,
                                                 hull_height=event.height)
                try:
                    self.tableframe.configure(width=event.width-50)
                    self.tableframe.configure(height=event.height-100)
                except:
                    pass
                try:
                    self.cv.configure(width=event.width-50)
                    self.cv.configure(height=event.height-100)
                    self.cv.configure(scrollregion=(0,0,self.canvas_x,self.canvas_y))
                except:
                    pass
        return

    def showBlankCanvas(self, frame, row):
        """Create a blank canvas for introduction"""
        width=900
        height=750
        cv=Canvas(frame,
                       width=self.canvas_x_size-self.canvas_border_x,
                       height=self.canvas_y-self.canvas_border_y,
                       scrollregion=(0,0,self.canvas_x_size,self.canvas_y),
                       bd=0, relief=GROOVE, bg='#9999CC',
                       highlightthickness=0)
        cv.grid(row=row,column=0,columnspan=10,rowspan=1,
                    sticky='news',
                    pady=3,ipady=2)
        #frame.add(cv)
        cv.configure(width=width-50)
        cv.configure(height=height-300)
        return cv

    def removeBlankCanvas(self):
        self.cv.destroy()
        return

    def welcomeLogo(self):
        """Show the welcome logo and text"""

        self.removeLogo()
        self.logo=None
        import tkFont
        try:
            logo = PEAT_images.logo_large_mono()
            self.cv.create_image(self.canvas_x_size/2,self.canvas_y_size/2,image=logo)
            self.cv.image=logo
        except:
            import tkFont
            font=tkFont.Font(family='Arial',size=38)
            self.cv.create_text(self.canvas_x_size/2,self.canvas_y_size/2-50,
                                                  text='P E A T',
                                                  font=font,anchor='n',
                                                  fill='black',
                                                  tag='logo')
        text=['Welcome to Protein Engineering and Analysis Tool  (PEAT)',
              'Authors: Jens Erik Nielsen, Damien Farrell and Michael Johnston',
              'Copyright Jens Erik Nielsen, University College Dublin 2003-, All rights reserved']
        self.y=20
        self.x=380
        ifont=tkFont.Font(family='Arial',size=12)
        for line in text:
            self.cv.create_text(self.x,self.y,text=line,fill='white',font=ifont, tag='logo')
            self.y=self.y+15

        # Remove the logo after some time
        #self.wait=self.after(4000, self.show_startinstructions)
        return

    def removeLogo(self,event=None):
        """Remove the logo and welcome message."""
        if self.cv == None:
            return
        self.cv.delete('logo')
        self.cv.image=None
        self.master.unbind('<KeyPress>')
        self.master.unbind('<Button-1>')
        try:
            self.after_cancel(self.wait)
        except:
            pass
        return

    def onlineDocumentation(self,event=None):
        """Open the online documentation"""
        import webbrowser
        link='http://enzyme.ucd.ie/main/index.php/PEAT_DB'
        webbrowser.open(link,autoraise=1)
        return

    def gotoBugzilla(self, event=None):
         """Open bugzilla site"""
         import webbrowser
         link='http://peat.ucd.ie/bugzilla/'
         webbrowser.open(link,autoraise=1)
         return

    def aboutPEAT(self):
        self.ab_win=Toplevel()
        self.ab_win.geometry('+100+350')
        self.ab_win.title('About PEAT')

        logo = PEAT_images.logo()
        label = Label(self.ab_win,image=logo)
        label.image = logo
        label.grid(row=0,column=0,sticky='news',padx=4,pady=4)

        text=['P E A T ','Protein Engineering and Analysis Tool',
        'Version 2 ','A database tool for analysing the effect of point mutations',
        'on the catalytic and structural characteristics of proteins and enzymes',
        'Authors: Jens Erik Nielsen, Damien Farrell and Chresten Soendergaard',
        'University College Dublin','(C) Copyright 2003- Jens Erik Nielsen All rights reserved']
        row=1
        for line in text:
            tmp=Label(self.ab_win,text=line)
            tmp.grid(row=row,column=0,sticky='news',padx=4)
            row=row+1
        return

    def recordEvent(self,message='Timestamp'):
        """Display a message in the eventlog"""
        import time
        tid=time.strftime('%H:%M:%S %d/%m',time.localtime(time.time()))
        self.eventlog.insert('0.0','%s: %s\n' %(tid,message))
        self.master.update_idletasks()
        return

    def findValue(self):
        self.currenttable.findValue()
        return

    def dofindText(self, event=None):
        """Find the text in the table"""
        if not hasattr(self,'currenttable'):
            return
        import string
        if string.strip(self.findtext.get())=='':
            return
        searchstring=self.findtext.get()
        if self.currenttable!=None:
            self.currenttable.findValue(searchstring)
        return

    def dofindAgain(self, event=None):
        """Find again"""
        if not hasattr(self,'currenttable'):
            return
        searchstring=self.findtext.get()
        if self.currenttable!=None:
            self.currenttable.findValue(searchstring, findagain=1)
        return

    def loginSetup(self):
        """Login setttings for authentication to ZEO"""

        mpDlg = MultipleValDialog(title='Login settings',
                                    initialvalues=(self.username,
                                                   self.password),
                                    labels=('user name','password'),
                                    types=('string','password'),
                                    parent=self.main)

        if mpDlg.result == True:
            self.username = mpDlg.results[0]
            self.password = mpDlg.results[1]
            self.preferences.set('username',self.username)
            self.preferences.set('password',self.password)
        else:
            return
        return

    def showSettings(self):
        """Settings dialog"""
        import Prefs
        if self.showDialogsinSidePane == True and self.DB!=None:
            self.resetSidePane(width=320)
            X=Prefs.preferences_dialog(parent=self,parentframe=self.sidepane,
                                       subset='PEAT',callback=self.getPrefs)
        else:
            X=Prefs.preferences_dialog(parent=self,subset='PEAT',callback=self.getPrefs)
        return

    def hideFieldsDialog(self):
        """Allow fields to be hidden"""
        fr=Toplevel() 
        fr.geometry('300x300+300+200')
        fr.title('Show/hide fields')
        fr.grab_set()
        fr.transient()
        userfields = self.DB.meta.userfields
        def apply():
            show=list(checkbuttons.getcurselection())            
            for f in userfields:
                if f not in show:                    
                    self.DB.showField(f, False)
                else:
                    self.DB.showField(f)                
            self.updateTable()
            return
        def close():
            fr.destroy()

        Button(fr, text='Apply', command=apply).pack(fill=BOTH,expand=1,padx=2, pady=2)
        Button(fr, text='Close', command=close).pack(fill=BOTH,expand=1,padx=2, pady=2)                    
        checkbuttons = Pmw.RadioSelect(fr,
                buttontype = 'checkbutton',
                orient = 'vertical',
                labelpos = 'n',                
                label_text = 'Deselect fields to hide them:')          
        checkbuttons.pack(fill=BOTH,expand=1)        
        for f in userfields:            
            checkbuttons.add(f)
            if userfields[f]['show'] == True:
                checkbuttons.invoke(f)
        checkbuttons.pack(fill=X, padx=5, pady=5)
        return
        
    def quit(self,event=None):
        """Close DB and quit"""
        answer = self.closeDB()
        if answer == 'cancel':
            return
        self.main.destroy()
        self.preferences.save_prefs()
        if not self.parent:
            sys.exit()
        return

    def showtablePrefs(self):
        if self.table:
            self.table.showtablePrefs()
        return

#
#move these to an actions class..?
#calls to updateTable could be replaced with an updateCell function
#

    def importCSV(self):
        """Import CSV files dialog"""
        if self.DB == None:
            return
        from PEATDB.IO import Importer
        IM = Importer()
        IM.getLines()
        IM.showimportDialog()
        self.wait_window(IM.ask_csv)
        #also get field for unique key
        vals = IM.data[0].keys()
        mpDlg = MultipleValDialog(title='Provide a Key Field',
                                    initialvalues=([vals]),
                                    labels=(['key field:']),
                                    types=(['list']),
                                    parent=self.main)
        if mpDlg.result == True:
            key = mpDlg.results[0]
        else:
            key = 'name'    
        self.DB.importDict(importdata=IM.data, namefield=key)
        self.updateTable()
        return

    def importMutants(self):
        """Specialised import method for adding mutant data from a csv
        file. Requires a wt structure to create mutant AA sequences"""
        
        def getpdb():
            return
        def doimport():
            #self.DB.importDict(importdata=IM.data, namefield='name')
            #use mutations field to derive aa seq..
            return
                
        fr=Frame()
        win=Toplevel()
        Label(win, text='This dialog allows you to import a set of mutant data\n'                                
                        'along with some associated experimental values. You will\n'
                        'need to provide the the structure for creating an AA\n'
                        'sequence for each mutant. Mutation codes should be of\n' 
                        'the form chain:residuenumber:code',bg='#ccFFFF').pack(fill=BOTH,expand=1)
        Button(win,text='Set a wt PDB',command=getpdb).pack(fill=BOTH,expand=1)
        fr=Frame(win); fr.pack(fill=X,expand=1)
        self.useref = IntVar()
        Label(fr,text='Use current reference protein').pack(side=LEFT)
        Checkbutton(fr,variable=self.useref).pack(side=LEFT)       
        self.set_centered_geometry(self.main,win)
        Button(win,text='Continue',command=doimport).pack(fill=BOTH,expand=1)
        return        
    
    def importFileset(self):
        """Use filehandler class to import a set of external files"""
        if self.DB == None:
            return
        fh = FileHandler(parent=self)
        if self.showDialogsinSidePane == True and self.DB!=None:
            self.resetSidePane(width=300)
            f=self.sidepane
        else:
            f=None
        fh.importFileset(DB=self.DB, parentframe=f, callback=self.updateTable)
        #self.updateTable()
        return

    def exportExtFiles(self):
        """Get all stored ext files in DB and save to folder"""
        fh = FileHandler(parent=self)
        fh.exportExtFiles(DB=self.DB)
        return

    def showExtFiles(self):
        """Display all stored ext files in DB"""
        fh = FileHandler(parent=self)
        fh.displayExtFiles(DB=self.DB)
        return

    def inspectDBMeta(self):
        """View DB meta like userfields, mainly for debugging"""
        if self.DB == None:
            return
        #We reuse the ekin meta data class here.. should make a base class
        #and inherit ekin and this one from it..
        from Meta import PEATMetaData
        self.PM = PEATMetaData(data=self.DB.meta)
        self.mwin = self.PM.createDialog(parent=self.main)
        self.main.wait_window(self.mwin)
        return

    def editDBInfo(self):
        """Edit info field of Meta for DB"""
        if self.DB == None:
            return
        from Meta import PEATMetaData

        self.PM = PEATMetaData(data=self.DB.meta.info)
        self.mwin = self.PM.createDialog(parent=self.main, label='DB Specific Info')
        self.main.wait_window(self.mwin)
        return

    def open_link(self, protein, field_name):
        """Open a hyperlink"""
        D = self.DB.data
        try:
            link = D[protein_name][field_name]['link']
            print link
            import webbrowser
            webbrowser.open(link,autoraise=1)
        except:
            pass
        return

    def edit_link(self, protein, field_name):
        """Edit a hyperlink"""
        D = self.DB.data
        if D[protein].has_key(field_name):
            data=D[protein][field_name]
        else:
            data=None

        self.editlinkframe=Toplevel()
        # set the position of the window
        top=self.winfo_toplevel()
        rootx=top.winfo_rootx()
        rooty=top.winfo_rooty()
        self.editlinkframe.geometry('+%d+%d' %(rootx+200,rooty+200))

        self.editlinkframe.title('Edit hyperlink')
        try:
            currtext = data['text']
            currlink = data['link']
        except:
            currtext = ''
            currlink = ''
        self.newtextvar=StringVar()
        self.newtextvar.set(currtext)
        self.newlinkvar=StringVar()
        self.newlinkvar.set(currlink)
        def getdata():
            self.linkdata={}
            self.linkdata['text'] = self.newtextvar.get()
            self.linkdata['link'] = self.newlinkvar.get()
            self.editlinkframe.destroy()
            return
        def close():
            self.editlinkframe.destroy()
            return
        Label(self.editlinkframe,text='Text:').grid(row=0,column=0,padx=2,pady=2)
        Entry(self.editlinkframe,textvariable=self.newtextvar,width=40).grid(row=0,column=1,padx=2,pady=2)
        Label(self.editlinkframe,text='Link:').grid(row=1,column=0,padx=2,pady=2)
        Entry(self.editlinkframe,textvariable=self.newlinkvar,width=40).grid(row=1,column=1,padx=2,pady=2)
        Button(self.editlinkframe,text='OK',command=getdata).grid(row=2,column=0,padx=2,pady=2,sticky='NEWS')
        Button(self.editlinkframe,text='Cancel',command=close).grid(row=2,column=1,padx=2,pady=2,sticky='NEWS')
        self.editlinkframe.columnconfigure(1,weight=1)
        self.editlinkframe.focus_set()
        self.master.wait_window(self.editlinkframe)
        if not getattr(self,'linkdata',None):
            self.linkdata=None
        if self.linkdata:
           # Got the record and column - add the data
            self.DB.data[protein][field_name] = self.linkdata
            delattr(self,'linkdata')
            #self.updateTable(protein, field_name)
            self.table.redrawCell(recname=protein, colname=field_name)
        return

    def edit_notes(self, protein, field_name):
        """edit notes in a text frame"""
        D = self.DB.data
        if D[protein].has_key(field_name):
            data=D[protein][field_name]
        else:
            data=None

        info={}
        info['record']=protein
        info['column']=field_name
        import textFrame
        self.notesframe=textFrame.textFrame(parent=self,title=protein+' notes',parent_info=info)
        if data!=None:
            self.notesframe.load_text(data['text'])
        self.master.wait_window(self.notesframe.frame)

        # did we get any data back?
        if not getattr(self,'notesdata',None):
            self.notesdata=None
        if self.notesdata:
            print self.notesdata
            if self.notesdata.has_key('parent_info'):
                notes_info=self.notesdata['parent_info']
                protein=notes_info['record']
                column=notes_info['column']
                newnotes=self.notesdata['text']
                del self.notesdata['parent_info']
            # Got the record and column - add the data

            self.DB.data[protein][field_name]  = self.notesdata
            delattr(self,'notesdata')
            #self.updateTable(protein, field_name)
            self.table.redrawCell(recname=protein, colname=field_name)
        return

    def display_file(self, protein, field_name):
        """Uses a new method of file storage with blobs
           We just get the blob from the filerecord object"""
        D = self.DB.data
        if D[protein].has_key(field_name):
            filerec=D[protein][field_name]
        else:
            return
        fh = FileHandler(self)
        fh.displayFile(filerec)
        return

    def add_file(self,protein,field_name):
        """Uses a new method of file storage with blobs
           We just get the blob from the filerecord object"""
        D = self.DB.data
        if D[protein].has_key(field_name):
            filerec=D[protein][field_name]
        else:
            filerec=None
        fh = FileHandler(self)
        newfile = fh.addFileDialog(filerec)
        #self.main.wait_window(win)
        if newfile != None:
            self.DB.addBlob(protein, field_name, fh.newfile)
            #self.updateTable(protein, field_name)
            self.table.redrawCell(recname=protein, colname=field_name)
        return

    def save_file(self,protein,field_name):
        D = self.DB.data
        if D[protein].has_key(field_name):
            filerec=D[protein][field_name]
        else:
            return
        import shutil
        myblob = filerec.blob
        ext = filerec.ext
        if myblob != None:
            f = myblob.open("r")
            name = filerec.name
            pth=tkFileDialog.asksaveasfilename(defaultextension=ext,
                                   initialfile=name,
                                   initialdir=os.getcwd(),
                                   filetypes=[("All files","*.*")])
            shutil.copyfile(f.name, pth)
        return

    def addPDBFile(self):
        """Add pdb file"""
        sel_prot = self.table.get_selectedRecordNames()[0]
        DBActions.addPDBFile(self.DB, sel_prot)
        self.updateTable()
        return

    def save_structure(self,protein,field_name,filename=None):
        """Save a PDB file from a Structure cell"""
        pdblines,X=self.DB.getStructure(protein,field_name)
        if not pdblines:
            return
        if filename==None:
            import tkFileDialog, os
            savedir=self.preferences.get('datadir')
            filename=tkFileDialog.asksaveasfilename(defaultextension='.pdb',
                                                    initialdir=savedir,
                                                    filetypes=[("PDB files","*.pdb"),
                                                               ("All files","*.*")])
        if filename:
            if filename[-4:]!='.pdb':
                filename=filename+'.pdb'
            fd=open(filename,'w')
            for line in pdblines:
                fd.write(line)
            fd.close()
        return

    def viewStructureText(self, protein, field_name):
        """Show PDB text"""
        pdblines,X=self.DB.getStructure(protein,field_name)
        if not pdblines:
            return
        import textFrame
        tf = textFrame.textFrame(parent=self.main,title=protein+' Structure')
        tf.load_text(pdblines)
        return        
    
    def displayStructure(self, protein, field_name):
        """Display a structure"""
        app = self.preferences.get('molgraphApplication')
        path = self.preferences.get('molgraphAppPath')
        DBActions.DB = self.DB
        DBActions.displayStructure(protein, field_name, app, path)
        return

    #these 2 should be peatrecord methods

    def get_protein_sequence(self, protein_name):
        """Get the amino acid sequence for this protein"""
        D = self.DB.data
        if D[protein_name].has_key('aaseq'):
            aaseq=D[protein_name]['aaseq']
        else:
            aaseq=None
        return aaseq

    def format_protein_sequence(self,sequence):
        """ Returns a string containing the one-letter sequence"""
        if not sequence:
            return ''
        seq=''
        self.aa={'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H',
                'ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',
                'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

        for resnum,restyp in sequence:
            if self.aa.has_key(restyp):
                seq=seq+self.aa[restyp]
            else:
                seq=seq+'?'
        return seq

    def display_sequence(self,protein,field_name):
        """Display the sequence of the current protein regardless of which column is chosen"""
        D = self.DB.data
        if D.has_key(protein):
            sequence=self.get_protein_sequence(protein)
            sequence=self.format_protein_sequence(sequence)
            """"Get font"""
            import tkFont
            font=tkFont.Font(family='Courier',size=14,weight='bold')
            seqwin=Toplevel()
            seqwin.geometry('600x400+100+100')
            seqwin.title('Amino acid sequence')
            ST=Text(seqwin,font=font)
            ST.grid(row=0,column=0)
            text='Protein sequence\n\n'
            text=text+'             10         20         30         40\n'
            text=text+'              |          |          |          |'
            count=0
            for s in sequence:
                if count%40==0:
                    text=text+'\n'
                    text=text+'%4d ' %count
                text=text+s
                count=count+1
                if count%10==0:
                    text=text+' '
            text=text+'\n'
            text=text+'\nLength: %d \n' %len(sequence)
            ST.insert(END,text)
        return

    def startDNAtool(self, data=None):
        """Start DNAtool - this function is far too long and should be simplified"""
        self.DNAtool_instance = None
        DB = self.DB
        D = DB.data
        M = DB.meta  #DNAtool/primers stuff is now stored in DB.meta

        #change this later
        self.DBfiledir = os.getcwd()

        if not self.DNAtool_instance:

            # Do we have data to pass?
            data_passed = {}
            if self.DNAtool_state:
                data_passed = self.DNAtool_state

            # Load the primer dictionary into the dictionary that are being passed
            data_passed['PEAT_DBrecord']='None'
            data_passed['primer_dict']={}
            if M:
                if M.has_key('DNAtool_primers'):
                    data_passed['primer_dict'] = M['DNAtool_primers'].copy()
            else:
                data_passed['PEAT_DBrecord'] = 'No database'

            # Delete any data that we only transfer if a protein is selected
            if data_passed.has_key('DNAseq'):
                del data_passed['DNAseq']
            if data_passed.has_key('ORF_selected'):
                del data_passed['ORF_selected']

            #Get the DNA sequence and ORF info from the currently selected protein
            sel_prot = self.table.get_selectedRecordNames()[0]

            # Now pass the data
            if sel_prot:
                data_passed['PEAT_DBrecord'] = sel_prot
                import types
                
                if D[sel_prot].has_key('DNAseq') and D[sel_prot].DNAseq != None:
                    data_passed['DNAseq'] = D[sel_prot].DNAseq

                    # ORF info
                    data_passed['ORF_selected']={}
                    keys=['start','aastart_number','frame','stop','length','aaseq','aaseq3']
                    for key in keys:
                        if D[sel_prot]['ORF_selected'].has_key(key):
                            data_passed['ORF_selected'][key] = D[sel_prot]['ORF_selected'][key]

                    # If we didn't find any data, then don't copy anything
                    if len(data_passed['ORF_selected'].keys())==0:
                        del data_passed['ORF_selected']

            # Now start DNAtool
            import DNAtool.DNAtool as DNAtool
            self.DNAtool_instance=DNAtool.MainWindow(self, data_passed)
            self.master.wait_window(self.DNAtool_instance.master) # Wait for DNAtool to finish
            self.DNAtool_instance=None

            # We came back from DNAtool.
            if self.DNAtool_state:
                if self.DNAtool_state['DNAseq_status']=='MODIFY RECORD':
                    # Modify the present record
                    self.DNAtool_state['DNAseq_status']='OK'

                elif self.DNAtool_state['DNAseq_status']=='NEW SEQ':
                    # Added this in as a fix to modify record when no seq previously
                    if sel_prot:
                        DNAdata = self.DNAtool_state
                        ok = DB.addDNAseq(sel_prot, DNAdata['DNAseq'])
                        data = DNAdata['ORF_selected']
                        ok = DB.addProperty(sel_prot, 'ORF_selected', data)
                       
                        if not ok:
                            raise 'Something went wrong when adding ORF info'

                        #Add the protein sequence
                        ok=DB.addProtseq(sel_prot, DNAdata['ORF_selected']['aaseq3'],
                                                    DNAdata['ORF_selected']['aastart_number'])
                        if not ok:
                            raise 'Something went wrong when adding the protein sequence'

                #this deals with changing the record in-place only when ORF is changed
                elif self.DNAtool_state['DNAseq_status']=='MODIFY ORF':
                    DNAdata=self.DNAtool_state
                    #Now do a check to see if frame has been changed
                    if D[sel_prot].has_key('ORF_selected'):
                        currframe=D[sel_prot]['ORF_selected']['frame']
                        newframe=self.DNAtool_state['ORF_selected']['frame']

                        if currframe!=newframe:
                            import tkMessageBox
                            tkMessageBox.showerror(title='Frame Changed',
                                message='Frame has been changed. You cannot save this in the same record!',
                                parent=self.master)
                            pass
                        else:
                            D[sel_prot]['ORF_selected']=DNAdata['ORF_selected']

                elif self.DNAtool_state['DNAseq_status']=='ADD NEW RECORD':

                    # Add a new protein record
                    #print "DNAtool_state['DNAseq_mutations']", self.DNAtool_state['DNAseq_mutations']
                    self.DNAtool_state['DNAseq_status']='OK'
                    initialvalue=self.DB.get_next_freename_newprot()
                    if self.DNAtool_state.has_key('DNAseq_mutations'):
                        mutations=self.DNAtool_state['DNAseq_mutations']
                        import DNAtool.mutation
                        mutant_text=''
                        if sel_prot:
                            mutant_text=sel_prot+'+'
                        for mut in mutations:
                            old=DNAtool.mutation.three_to_one[mut.split(':')[1]]
                            number=mut.split(':')[2]
                            new=DNAtool.mutation.three_to_one[mut.split(':')[-1]]
                            mutant_text=mutant_text+'%s%s%s+' %(old,number,new)
                        initialvalue=mutant_text[:-1]

                    # Get the name for the mutation
                    ok=None
                    while not ok:
                        import tkSimpleDialog
                        record_name=tkSimpleDialog.askstring('Record name',
                                                'Give a name for the new record',
                        initialvalue=initialvalue,parent=self.master)
                        if not record_name:
                            return

                        # Store the data
                        ok=DB.newmolecule(record_name)
                        if not ok:
                            self.warning('Name already in use',
                                         '%s is already in use' %record_name)

                    # Got the name, now just fill in the fields
                    DNAdata=self.DNAtool_state
                    ok=self.DB.addDNAseq(record_name, DNAdata['DNAseq'])
                    data=DNAdata['ORF_selected']
                    ok=self.DB[record_name]['ORF_selected'] = data
                    if not ok:
                        raise 'Something went wrong when adding ORF info'

                    # Add the protein sequence
                    ok=self.DB.addProtseq(record_name, DNAdata['ORF_selected']['aaseq3'],
                                                   DNAdata['ORF_selected']['aastart_number'])
                    if not ok:
                        raise 'Something went wrong when addding the protein sequence'                    
                    
                    #add a code for the mutation, this can be parsed to create
                    #a PEATSA mutationset later
                    ref = self.DB.meta.refprotein      
                    #DBActions.checkMutation(self.DB, record_name, ref)
                    self.updateTable(record_name, 'Mutations')                    
        return

    def startEkin(self, protein=None, field_name=None, mode=None):
        """Start ekin from peat cell"""
        import copy
        from Ekin.Ekin_main import EkinApp

        D = self.DB.data
        #Get the data of the cell
        if protein and field_name:
            if D[protein].has_key(field_name):
                E = copy.deepcopy(D[protein][field_name])
                #E = D[protein][field_name]
            else:
                E = None
            #add to open instances
            if not self.ekin_instances.has_key(protein+field_name):
                self.ekin_instances[protein+field_name]=1
            else:
                self.ekin_instances[protein+field_name] += 1
        else:
            E=None
            protein=None
            field=None

        EK = EkinApp(parent=self, project=E, protein=protein, field=field_name, mode=mode)
        self.Ekinprj=None
        self.master.wait_window(EK.ekin_win)

        # Ekin is done - did we get any data back?
        if protein and field_name:
            if hasattr(self,'ekin_instances') and protein != None and field_name!= None:
                if self.ekin_instances[protein+field_name] >1:
                    self.ekin_instances[protein+field_name] -= 1
                else:
                    del self.ekin_instances[protein+field_name]

        if not getattr(self, 'Ekinprj', None):
            self.Ekinprj = None
        if self.Ekinprj:
            # Yes, we have data - which protein/column should we associate it with?
            PEAT_info = self.Ekinprj._PEATinfo
            protein = PEAT_info['record']
            column = PEAT_info['column']
            # Got the record and column - simply add the data
            del self.Ekinprj._PEATinfo
            D[protein][field_name] = self.Ekinprj
            delattr(self,'Ekinprj')
            #self.updateTable(protein, field_name)
            self.table.redrawCell(recname=protein, colname=field_name)
        return

    def startLabbook(self, protein=None, field_name=None):
        """Edit a labbook style sub-table"""

        D = self.DB.data
        tabledata = None
        info = None
        #protein is ALL means we are using the main labbook
        if protein == 'ALL':
            tabledata = self.DB.meta['labbook']
            if tabledata == {}: tabledata=None
        elif hasattr(self.DB.meta, 'labbooks') and protein in self.DB.meta.labbooks:
            tabledata = self.DB.meta.labbooks[protein]               
        elif protein and field_name:
            if D[protein].has_key(field_name):
                tabledata = D[protein][field_name]
            info={}
            info['record']=protein
            info['column']=field_name

        from Labbook import LabbookApp
        if tabledata != None and isinstance(tabledata, dict):
            labbook = LabbookApp(parent=self, peatinfo=info, data=tabledata)
        else:
            labbook = LabbookApp(parent=self, peatinfo=info)
        self.master.wait_window(labbook.labbook_win)

        # did we get any data back?
        if not hasattr(self, 'templabbookdata'):
            self.templabbookdata=None
        if not self.templabbookdata is None:
            if protein == 'ALL':
                self.DB.meta['labbook'] = self.templabbookdata
            elif hasattr(self.DB.meta, 'labbooks') and protein in self.DB.meta.labbooks:
                self.DB.meta.labbooks[protein] = self.templabbookdata                
                self.updateLabbooksMenu()
            else:
                D[protein][field_name] = self.templabbookdata
            delattr(self,'templabbookdata')
            #self.updateTable(protein, field_name)
            self.table.redrawCell(recname=protein, colname=field_name)           
        return

    def updateLabbooksMenu(self ):
        """Add existing laboooks to tools menu"""
        import string, types
       
        self.labbookmenu.delete(0, self.labbookmenu.index('last'))
        self.labbookmenu.add_command(label='New Labbook', command=self.addLabbook)
        self.labbookmenu.add_command(label='Manage Labbooks', command=self.manageLabbooks)       

        if not hasattr(self.DB.meta,'labbooks'):
            return        
        self.labbookmenu.add_separator()
        for name in self.DB.meta.labbooks.keys():
            def func(**kwargs):
                def new():
                   self.startLabbook(**kwargs)
                return new
            self.labbookmenu.add_command(label=name, command=func(protein=name) )              
        return

    def addLabbook(self):
        """Create an extra labbook"""
        name = tkSimpleDialog.askstring('Name of Labbook',
                                             'Enter a name',
                                             initialvalue='',
                                             parent=self.main)
        if name:
            self.DB.addLabbook(name)
            self.startLabbook(name)            
        return

    def manageLabbooks(self):
        L=LabbookManager(self.main, self.DB)
        return

    def labbookSheetsSelector(self, frame):
        """Show list of labbook sheets"""        
        names = self.DB.meta.labbook.keys() 
        if hasattr(self.DB.meta, 'labbooks'):
            for l in self.DB.meta.labbooks:
                sheets = [l+'_'+ s for s in self.DB.meta.labbooks[l].keys()]
                names.extend(sheets)
            
        labbooklist = Pmw.ScrolledListBox(frame,              
                labelpos='nw',
                label_text='Tables',
                listbox_height = 8)
        labbooklist.pack(fill=Y,expand=1)
        labbooklist.setlist(names)
        return labbooklist           
            
    def alanineScan(self):
        return

    def setRefProtein(self):
        """Set a reference protein for making point mutation models"""
        if self.DB == None:
            return
        cframe = self.createChildFrame(200)
        Label(cframe,text='Ref Prot:').pack()
        o = self.createRecsOptsMenu(cframe)
        o.pack()
        return

    def updateMutations(self):
        """Check sequences and create mutation names"""
        if not hasattr(self.DB.meta, 'refprotein'):
            return
        ref = self.DB.meta.refprotein
        for p in self.DB.getRecs():
            DBActions.checkMutation(self.DB, p, ref)
        self.updateTable()
        return
    
    def updateAASequences(self):
        """Update the AA sequences for mutants, usually applied to 
        those mutants that don't have a DNA sequence, so at any time we
        can change their AA seq if the ref prot is altered
        requires: reference protein to be set"""
        if not hasattr(self.DB.meta, 'refprotein') or self.DB.meta.refprotein == None:
            tkMessageBox.showinfo('No ref protein',
                                  'Set a reference (wt) protein first')
            return 
        sel = self.table.get_selectedRecordNames()    
        DBActions.setSequencesfromMutationCodes(self.DB, selected=sel)    
        return
        
    def createRecsOptsMenu(self, parent, callback=None):
        """Get an option menu with a list of the records"""
        def setref(evt=None):
            prot = var.get()
            if self.DB[prot].Structure == None:              
                tkMessageBox.showinfo("Warning",
                                     'There is no structure for this record! '
                                      'You should add one.')                
            self.DB.meta.refprotein = prot            
        if not hasattr(self.DB.meta, 'refprotein'):
            self.DB.meta.refprotein = None

        var=StringVar()
        var.set(self.DB.meta.refprotein)
        o = Menubutton(parent, textvariable=var,
                            borderwidth=2,relief=GROOVE,
                            bg='#ccFFFF',
                            width=12)
        m=Menu(o,tearoff=0)
        o['menu']=m
        p=0
        for name in self.DB.getRecs():
            if p%30==0 and p!=0:
                colbrk=1
            else:
                colbrk=0
            m.add_radiobutton(label=name,                                                 
                              value=name,
                              variable=var,
                              columnbreak=colbrk,
                              command=setref)
            p+=1
        return o
    
    def remodelPEATModels(self):
        """Remodel current mutant sequences using ref protein"""
        if self.DB == None:
            return

        # Delete the temp models
        if getattr(self.DB,'3Dmodels', None):
            self.DB.tmp_3Dmodels={}
       
        c = self.createChildFrame()
        from Dialogs import PEATDialog

        Label(c, text='Construct Homology Models').pack(side=TOP,pady=4)
        o = self.createRecsOptsMenu(c)
        o.pack()
        doselectedvar=IntVar(); doselectedvar.set(0)
        c1=Frame(c); c1.pack(side=TOP,pady=4)
        Label(c1, text='Do selected only').pack(side=LEFT)
        Checkbutton(c1, variable=doselectedvar).pack(side=LEFT)
        pb=PEATDialog.createProgressBar(c)
        pb.frame.pack(padx=5,pady=5)
        def go():
            if doselectedvar.get() == True:
                sel = self.table.get_selectedRecordNames()
            else:
                sel = None
            self.stdout2Log()
            DBActions.checkModels(DB=self.DB,callback=pb.updateProgress,
                                    selected=sel)
            self.log2Stdout()
            self.updateTable()            
        Button(c,text='Go',command=go).pack()
        self.log = self.createLogWin(c)
        #log.delete(1.0,END)
        self.log.pack(fill=BOTH,padx=4,pady=4)
        return

    def createLogWin(self, parent):
        log = Pmw.ScrolledText(parent,
            labelpos = 'n',
            label_text='Log',
            usehullsize = 1,
            hull_width = 800,
            hull_height = 500,
            text_wrap='word')
        return log        

    def stdout2Log(self):
        """Redirect stdout to app control"""
        sys.stdout = self
        sys.stderr = self
        return
       
    def log2Stdout(self):
        """return to stdout"""
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        return
      
    def write(self, txt):
        """Handle stdout if required"""        
        self.log.appendtext(txt)
        self.log.update_idletasks()
        return

    def flush(self):
        pass
        return

    def fetchPDB(self):
        """Allow user to fetch a pdb file from net and store to record"""
        mpDlg = MultipleValDialog(title='Fetch PDB',
                                    initialvalues=('', self.DB.getRecs()),
                                    labels=('pdb id','record'),
                                    types=('string','list'),
                                    parent=self.main)

        if mpDlg.result == True:
            pdbid = mpDlg.results[0]
            recname = mpDlg.results[1]
            
        pdb = DBActions.fetchPDB(pdbid)       
        DBActions.addPDBFile(self.DB, recname, pdbdata=pdb)
        self.updateTable()
        return
    
    def startComparator(self):
        import sys
        from Compare import Comparator
        if self.showDialogsinSidePane == True:
            self.resetSidePane(width=500)
            C = Comparator(self, parentframe=self.sidepane, DB=self.DB)
        else:
            C = Comparator(self, parent=self, DB=self.DB)
        return

    def testsDialog(self):
        """Allow users to set up a test DB"""
        import zodbtest
        def createDB():
            mpDlg = MultipleValDialog(title='Create Test DB',
                                        initialvalues=('testdb.fs',
                                                       1000),
                                        labels=('filename','no. recs'),
                                        types=('string','int'),
                                        parent=self.main)

            if mpDlg.result == True:
                filename = mpDlg.results[0]
                norecs = mpDlg.results[1]
            else:
                return
            tstDB = zodbtest.createdb(filename, norecs)
            self.closeDB()
            self.loadDB(tstDB)
        def close():
            frame.destroy()
            self.resetSidePane(width=20)
            return

        if self.showDialogsinSidePane == True and self.DB!=None:
            self.resetSidePane(width=300) #clears it
            frame = Frame(master=self.sidepane)
            frame.pack(fill=BOTH)
        else:
            frame = Toplevel()
            frame.geometry('+100+450')
            frame.title('Tests')
        Button(frame,text='Create Test DB',command=createDB).pack(fill=BOTH,side=TOP)
        Button(frame,text='Close',command=close).pack(fill=BOTH,side=TOP)
        return

    def advancedSearch(self):
        """Search DB with current search class, which should be replaceable by
           another search class if needed"""
        if self.DB==None:
            return
        self.resetSidePane(width=400)
        def close():
            frame.destroy()
            self.resetSidePane(width=20)
            return
        
        from Search import  searchDialog
        S = searchDialog(self.sidepane, self.DB)        
        S.pack()        
        return

    def createServerProject(self):
        """Allow users with admin passwd to create a project
           on the remote MySql DB. Relstorage only."""        
        import MySQLdb as mysql      
        mpDlg = MultipleValDialog(title='Create DB on Server',
                                    initialvalues=(self.server, self.port,
                                                   self.project,'root','',''),
                                    labels=('server','port','project',
                                            'user','password','access list'),
                                    types=('string','int','string',
                                           'string','password','string'),
                                    parent=self.main)
        if not mpDlg.result:
            return
        server = mpDlg.results[0]
        port = mpDlg.results[1]
        dbname = mpDlg.results[2]
        user = mpDlg.results[3]
        passwd = mpDlg.results[4]
        access = mpDlg.results[5]
        try:
            db = mysql.connect(user=user, host=server,
                               passwd=passwd, port=port)
            c = db.cursor()    
            cmd = "create database " + dbname + ";"
            c.execute(cmd)
            cmd = "grant all privileges on " + dbname + ".* to" + access + "@%;"
            c.execute(cmd)
        except mysql.OperationalError, e:
            tkMessageBox.showinfo('Error',e)
            return
        self.connect(self, server=self.server, port=self.port,
                project=dbname, backend=relstorage)
        return        
    
    def set_geometry(self,pwidget,widget):
        """Set the position of widget in the middle of pwidget"""
        w,h,x,y=get_geometry(pwidget)
        sw,sh,dummy,dummy2=get_geometry(widget)
        xoffset=int((w-sw)/2)
        yoffset=int((h-sh)/2)
        widget.geometry('+%d+%d' %(x+xoffset,y+yoffset))
        return

#project class should handle which windows open and when to warn users
#that they should be careful in syncing when using sub apps.

class Project(object):
    """Handle an open project inside the application, also we can use this
       to save a project state and reload it later.."""
    def __init__(self, parent):
        self.parent=parent
        return

    def doRecentMenu(self, preferences):
        """Returns a menu of recent projects"""
        preferences.get('recent')
        menu = Menu()
        return menu

#
# Utility functions - move
#
def get_geometry(widget):
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


# move these

class LabbookManager(Toplevel):
    def __init__(self, parent, DB):
        Toplevel.__init__(self, parent)
        #self.transient(parent)
        self.DB = DB
        self.title('Manage Labbooks')
        self.parent = parent        
        self.geometry('400x300+300+200') 
        self.body = Frame(self)
        self.initial_focus = self.body
        self.body.pack(fill=BOTH,expand=1,padx=5, pady=5)        
        #self.grab_set()
        if not self.initial_focus:
            self.initial_focus = self
        self.protocol("WM_DELETE_WINDOW", self.cancel)
        self.updateList()
        self.labbooklist.pack(side=LEFT,fill=BOTH,expand=1)
        
        bf=Frame(self.body); bf.pack(side=RIGHT)
        Button(bf,text='New',command=self.new).pack(fill=X)
        Button(bf,text='Open',command=self.open).pack(fill=X)
        Button(bf,text='Remove',command=self.remove).pack(fill=X)
        Button(bf,text='Close',command=self.cancel).pack(fill=X)
        return

    def updateList(self):
        if hasattr(self.DB.meta,'labbooks'):
            names = self.DB.meta.labbooks.keys()
        else:
            names = []
        self.labbooklist = Pmw.ScrolledListBox(self.body,                
                labelpos='nw',
                label_text='Labbooks',
                listbox_height = 6,
                dblclickcommand=self.open)
        self.labbooklist.setlist(names)
        return

    def new(self):
        self.parent.addLabbook()
        return
    
    def open(self):
        name = self.labbooklist.curselection()
        self.parent.startLabbook(protein=name)
        return

    def remove(self):
        name = self.labbooklist.curselection()
        if not hasattr(self.DB.meta,'labbooks'):
            return
        if name in self.DB.meta.labbooks.keys():
            self.DB.meta.labbooks.remove(name)
        self.updateList()
        return
    
    def cancel(self):
        self.parent.focus_set()
        self.destroy()      
        return


class ToolBar(Frame):
    """Uses the parent instance to provide the functions"""
    def __init__(self, parent=None, parentapp=None):
        Frame.__init__(self, parent, width=600, height=18)
        self.parentframe = parent
        self.parentapp = parentapp

        return

    def add_button(self, name, callback, img=None, helptxt=None):
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

    def add_seperator(self):
        Label(self, text=' ').pack(side=LEFT, padx=2, pady=2)
        return

def main():
    "Run the application."

    import sys, os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-s", "--server", dest="server",
                        help="Open a remote db")
    parser.add_option("-o", "--port", dest="port",
                        help="Use port other than 3306")
    parser.add_option("-u", "--username", dest="username",
                        help="Access with a username")
    parser.add_option("-d", "--passwd", dest="passwd",
                        help="Use password")
    parser.add_option("-p", "--project", dest="project",
                        help="Specify project", default=1)
    parser.add_option("-f", "--file", dest="file",
                        help="Open a local db")
    opts, remainder = parser.parse_args()
    if opts.server != None:
        DB = PDatabase(server=opts.server, port=opts.port,
                  username=opts.username, password=opts.passwd,
                  project=opts.project)
        app=App(DB=DB)
    elif opts.file != None:
        if not os.path.exists(opts.file):
            print 'File does not exist.'
            DB=None

        else:
            try:
                DB = PDatabase(local=opts.file)
            except:
                print 'Could not open the file, it may be in use'
                return
        app=App(DB=DB)
        app.addtoRecent(os.path.abspath(opts.file))
    else:
        app=App()
    app.mainloop()
    return

if __name__ == '__main__':
    main()


