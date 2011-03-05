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
import tkFileDialog, tkMessageBox, tkSimpleDialog
import Pmw
import re
import os
import time
import math
import pickle

from Tables_IO import TableImporter, TableExporter
from Prefs import Preferences
from GUI_helper import *

class LabbookApp(Frame, GUI_help):
    """
    Labbook app
    """
    appname = 'Labbook client'

    def __init__(self, parent=None, peatinfo=None, title=None,
                     data=None, datafile=None):
        "Initialize the application."
        self.parent=parent

        #If there is data to be loaded, show the dialog first
        if self.parent == None:
            Frame.__init__(self)
            self.labbook_win=self.master
            self.peatinfo=None
        else:
            self.labbook_win=Toplevel()
            self.master = self.labbook_win
            self.peatinfo=peatinfo      #reference to peat protein/field
            self.DB = self.parent.DB    #ref to peat

        if title != None:
            self.title = 'subtable_' + title
        else:
            self.title = 'Labbook Tool'

        self.ID='Labbook'
        # Get platform into a variable
        import platform
        self.currplatform=platform.system()
        if not hasattr(self,'defaultsavedir'):
            self.defaultsavedir = os.getcwd()

        self.preferences=Preferences('Labbook',{'check_for_update':1})
        #self.loadprefs()
        if self.peatinfo:
            protein = self.peatinfo['record']
            field_name = self.peatinfo['column']
            self.labbook_win.title('Labbook: '+protein+'_'+field_name)
        else:
            self.labbook_win.title(self.title)
        self.labbook_win.geometry('+200+100')
        self.x_size=1000
        self.y_size=500
        self.createMenuBar()
        self.apptoolBar = ToolBar(self.labbook_win, self)
        self.apptoolBar.pack(fill=BOTH, expand=NO)
        #add find bar
        #self.createSearchBar()

        if data != None:
            self.data = data
            self.new_project(data)
        elif datafile != None:
            self.open_project(filename=datafile)
        else:
            self.new_project()
        #add return to main PEAT db button
        if self.parent:
            Button(self.labbook_win,text="Return to Database",command=self.return_data).pack()
        self.labbook_win.protocol('WM_DELETE_WINDOW',self.quit)
        self.createBindings()
        self.clipboard = None
        return

    def createBindings(self):
        """Bind keys"""
        self.labbook_win.bind("<Control-n>", self.new_project)
        self.labbook_win.bind("<Control-o>", self.open_project)
        self.labbook_win.bind("<Control-s>", self.save_project)
        self.labbook_win.bind("<Control-q>", self.quit)
        self.labbook_win.bind("<Control-i>", self.import_csv)
        self.labbook_win.bind("<Control-e>", self.export_csv)

        return

    def createMenuBar(self):
        """Create the menu bar for the application. """
        self.menu=Menu(self.labbook_win)
        self.proj_menu={'01New':{'cmd':self.new_project},
                        '02Open':{'cmd':self.open_project},
                        '03Close':{'cmd':self.close_project},
                        '04Save':{'cmd':self.save_project},
                        '05Save As':{'cmd':self.save_as_project},
                        '06Preferences..':{'cmd':self.showPrefsDialog},
                        '08Quit':{'cmd':self.quit}}
        if self.parent:
            self.proj_menu['08Return to Database']={'cmd':self.return_data}
        self.proj_menu=self.create_pulldown(self.menu,self.proj_menu)
        self.menu.add_cascade(label='Project',menu=self.proj_menu['var'])
        
        self.edit_menu={'01Copy':{'cmd':self.copy},
                        '02Paste':{'cmd':self.paste},
                        '03Copy Columns': {'cmd':self.copyColumns},
                        '04Paste Columns': {'cmd':self.pasteColumns}}
        self.edit_menu=self.create_pulldown(self.menu,self.edit_menu)
        self.menu.add_cascade(label='Edit',menu=self.edit_menu['var'])
        
        self.records_menu={'01Add Row':{'cmd':self.add_Row},
                         '02Delete Row':{'cmd':self.delete_Row},
                         '03Add Column':{'cmd':self.add_Column},
                         '04Delete Column':{'cmd':self.delete_Column},
                         '05Auto Add Rows':{'cmd':self.autoAdd_Rows},
                         '06Auto Add Columns':{'cmd':self.autoAdd_Columns},
                         '07Find Cell':{'cmd':self.createSearchBar},
                         '08Filter':{'cmd':self.showFilteringBar} 
                         }
        self.records_menu=self.create_pulldown(self.menu,self.records_menu)
        self.menu.add_cascade(label='Records',menu=self.records_menu['var'])

        self.sheet_menu={'01Add Sheet':{'cmd':self.add_Sheet},
                         '02Remove Sheet':{'cmd':self.delete_Sheet},
                         '03Copy Sheet':{'cmd':self.copy_Sheet},
                         '04Rename Sheet':{'cmd':self.rename_Sheet},
                         '05Merge Sheet':{'cmd':self.merge_Sheet},
                         }
        self.sheet_menu=self.create_pulldown(self.menu,self.sheet_menu)
        self.menu.add_cascade(label='Sheet',menu=self.sheet_menu['var'])

        self.IO_menu={'01Import from csv file':{'cmd':self.import_csv},
                      '02Export to csv file':{'cmd':self.export_csv},
                      '03Import external fileset':{'cmd':self.import_fileset},
                      }

        self.IO_menu=self.create_pulldown(self.menu,self.IO_menu)
        self.menu.add_cascade(label='Import/Export',menu=self.IO_menu['var'])


        self.view_menu = Menu(self.menu, tearoff=0)
        self.view_menu.add_radiobutton(label="Normal View", state=ACTIVE,command=self.normal_view)
        self.view_menu.add_radiobutton(label="Page View", command=self.page_view)
        self.view_menu.add_command(label="Show All", command=self.showAll)  
        self.menu.add_cascade(label='View',menu=self.view_menu)

        #
        # Help menu
        #
        self.help_menu={'01Online Help':{'cmd':self.online_documentation},
                        '02About':{'cmd':self.about_Labbook}}
        self.help_menu=self.create_pulldown(self.menu,self.help_menu)
        self.menu.add_cascade(label='Help',menu=self.help_menu['var'])

        self.labbook_win.config(menu=self.menu)

        return

    def createSearchBar(self, event=None):
        """Add a find entry box"""
        frame = Frame(self.labbook_win)
        row=0
        def close():
            frame.destroy()
        self.findtext=StringVar()
        self.findbox=Entry(frame,textvariable=self.findtext,width=30,bg='white')
        self.findbox.grid(row=row,column=1,sticky='news',columnspan=2,padx=2,pady=2)
        self.findbox.bind('<Return>',self.do_find_text)
        Label(frame,text='Find:').grid(row=row,column=0,sticky='ew')
        self.findagainbutton=Button(frame,text='Find Again', command=self.do_find_again)
        self.findagainbutton.grid(row=row,column=3,sticky='news',padx=2,pady=2)
        self.cbutton=Button(frame,text='Close', command=close)
        self.cbutton.grid(row=row,column=4,sticky='news',padx=2,pady=2)
        frame.pack(fill=BOTH, expand=NO)
        return

    def showFilteringBar(self):
        s = self.notebook.getcurselection()
        page = self.notebook.page(s)
        if not hasattr(self.currenttable, 'filterframe') or self.currenttable.filterframe == None:
            frame = self.currenttable.createFilteringBar(page)          
            frame.grid(row=5,column=0,columnspan=3)
        return

    def showAll(self):
        """Show all recs"""
        self.currenttable.showAll()
        return
    
    def loadprefs(self):
        """Setup default prefs file if any of the keys are not present"""
        defaultprefs = {'textsize':14,
                         'windowwidth': 800 ,'windowheight':600}
        for prop in defaultprefs.keys():
            try:
                self.preferences.get(prop)
            except:               
                self.preferences.set(prop, defaultprefs[prop])
        print self.preferences.get('textsize')
        return

    def showPrefsDialog(self):
        self.prefswindow = self.currenttable.showtablePrefs()

        return


    def new_project(self, data=None):
        """Create a new table, with model and add the frame"""
        if hasattr(self,'currenttable'):
            self.notebook.destroy()
            self.currenttable.destroy()

        #Create the sheets dict
        self.sheets = {}
        self.notebook = Pmw.NoteBook(self.labbook_win, raisecommand=self.setcurrenttable)
        self.notebook.pack(fill='both', expand=1, padx=4, pady=4)
        if data !=None:
            for s in data.keys():
                sdata = data[s]
                #try:
                self.add_Sheet(s ,sdata)
                #except:
                #    print 'skipping'
        else:
            #do the table adding stuff for the initial sheet
            self.add_Sheet('sheet1')
        self.notebook.setnaturalsize()
        self.setcurrenttable()
        return

    def open_project(self, event=None, filename=None):
        import os
        if filename == None:
            import tkFileDialog
            filename=tkFileDialog.askopenfilename(defaultextension='.labbook',
                                                      initialdir=os.getcwd(),
                                                      filetypes=[("Pickle file","*.labbook"),
                                                                 ("All files","*.*")],
                                                      parent=self.labbook_win)
        if os.path.isfile(filename):
            fd=open(filename)
            import pickle
            data=pickle.load(fd)
            fd.close()
        self.new_project(data)
        self.filename=filename
        return


    def save_project(self, event=None):
        if not hasattr(self, 'filename'):
            self.save_as_project()
        elif self.filename == None:
            self.save_as_project()
        else:
            self.do_save_project(self.filename)

        return

    def save_as_project(self):
        """Save as a new filename"""
        import tkFileDialog, os
        filename=tkFileDialog.asksaveasfilename(parent=self.labbook_win,
                                                defaultextension='.labbook',
                                                initialdir=self.defaultsavedir,
                                                filetypes=[("Labbook project","*.labbook"),
                                                           ("All files","*.*")])
        if not filename:
            print 'Returning'
            return
        self.filename=filename
        self.do_save_project(self.filename)
        return

    def do_save_project(self, filename):
        """Get model dicts and write all to pickle file"""
        data={}
        for s in self.sheets.keys():
            currtable = self.sheets[s]
            model = currtable.getModel()
            data[s] = model.getData()
  
        fd=open(filename,'w')
        import pickle
        pickle.dump(data,fd)
        fd.close()
        return

    def close_project(self):
        if hasattr(self,'currenttable'):
            #self.notebook.destroy()
            self.currenttable.destroy()
        return

    def import_csv(self):
        importer = TableImporter()
        #just use the dialog to load and import the file
        importdialog = importer.import_Dialog(self.labbook_win)
        self.labbook_win.wait_window(importdialog)
        modeldata = importer.modeldata
        self.add_Sheet(sheetdata = modeldata)
        return

    def export_csv(self):
        exporter = TableExporter()
        exporter.ExportTableData(self.currenttable)
        return

    def import_fileset(self):
        """Import a series of external files in a folder"""
        if self.parent == None:
            import tkMessageBox
            tkMessageBox.showwarning("Not available", "You can't use this feature outside PEAT.")
            return

        from Extfile import FileHandler
        fh = FileHandler(parent=self)
        fh.importFileset(DB=self.DB,
                         callback=self.currenttable.redrawTable)

        return

    def add_Sheet(self, sheetname=None, sheetdata=None, model=None):
        """Add a new sheet - handles all the table creation stuff"""
        def checksheet_name(name):
            if name == '':
                tkMessageBox.showwarning("Whoops", "Name should not be blank.")
                return 0
            if self.sheets.has_key(name):
                tkMessageBox.showwarning("Name exists", "Sheet name already exists!")
                return 0
        noshts = len(self.notebook.pagenames())
        if sheetname == None:
            sheetname = tkSimpleDialog.askstring("New sheet name?", "Enter sheet name:",
                                                initialvalue='sheet'+str(noshts+1),
                                                parent=self.labbook_win)
        checksheet_name(sheetname)
        page = self.notebook.add(sheetname)
        #Create the model if none provided
        if model == None:
            if sheetdata != None:
                model = LabbookTableModel(sheetdata)
            else:
                model = LabbookTableModel(rows=10,columns=5)

        #create the table: we pass the parent instance of peat
        #and peat row/col if present
        self.currenttable = LabbookTable(parent=page, model=model, sheetname=sheetname,
                                        peat=self.parent, peatinfo=self.peatinfo)

        #Load preferences into table
        self.currenttable.loadPrefs(self.preferences)
        #This handles all the canvas and header in the frame passed to constructor
        self.currenttable.createTableFrame()
        #add the table to the sheet dict
        self.sheets[sheetname] = self.currenttable
        self.saved = 0
        return sheetname

    def delete_Sheet(self, s=None):
        """Delete a sheet"""
        if s == None:
            s = self.notebook.getcurselection()
        print s
        self.notebook.delete(s)
        del self.sheets[s]
        return

    def copy_Sheet(self, newname=None):
        """Copy a sheet"""
        newdata = self.currenttable.getModel().getData().copy()
        if newname==None:
            self.add_Sheet(None, newdata)
        else:
            self.add_Sheet(newname, newdata)
        return

    def rename_Sheet(self):
        """Rename a sheet"""
        s = self.notebook.getcurselection()
        newname = tkSimpleDialog.askstring("New sheet name?", "Enter new sheet name:",
                                                initialvalue=s)
        if newname == None:
            return
        self.copy_Sheet(newname)
        self.delete_Sheet()
        return

    def merge_Sheet(self):
        """Merge 2 sheets"""
        if len(self.sheets)<2:
            return
        def domerge():
            name = shs.getvalue()[0]
            newmodel = self.sheets[name].model
            curr.model.merge(newmodel)
            self.currenttable.redrawTable()
            return
        
        s = self.notebook.getcurselection()
        curr = self.sheets[s]
        fr=Toplevel()
        fr.title('Select sheet to merge')
        self.set_centered_geometry(self.labbook_win,fr)
        shs = self.sheetsSelector(fr)
        Button(fr,text='Merge',command=domerge).pack(side=TOP)
        return
        
    def setcurrenttable(self, event=None):
        """Set the currenttable so that menu items work with visible sheet"""
        try:
            s = self.notebook.getcurselection()
            self.currenttable = self.sheets[s]
        except:
            pass
        return

    def add_Row(self):
        """Add a new row"""
        self.currenttable.add_Row()
        self.saved = 0
        return

    def delete_Row(self):
        """Delete currently selected row"""
        self.currenttable.delete_Row()
        self.saved = 0
        return

    def add_Column(self):
        """Add a new column"""
        self.currenttable.add_Column()
        self.saved = 0
        return

    def delete_Column(self):
        """Delete currently selected column in table"""
        self.currenttable.delete_Column()
        self.saved = 0
        return

    def autoAdd_Rows(self):
        """Auto add x rows"""
        self.currenttable.autoAdd_Rows()
        self.saved = 0
        return

    def autoAdd_Columns(self):
        """Auto add x rows"""
        self.currenttable.autoAdd_Columns()
        self.saved = 0
        return

    def copy(self):
        """Copy current selection to internal clipboard"""
        T = self.currenttable
        self.clipboard = T.copyCell()        
        return
    
    def paste(self):
        print self.clipboard        
        return
    
    def copyColumns(self):
        """Copy columns"""       
        self.currenttable.getSelectedColumn()
        self.clipboard = self.currenttable.copyColumns()        
        return

    def pasteColumns(self):
        """Paste columns"""        
        self.currenttable.pasteColumns(self.clipboard)
        return
        
    def findValue(self):
        self.currenttable.findValue()
        return

    def do_find_text(self, event=None):
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

    def do_find_again(self, event=None):
        """Find again"""
        if not hasattr(self,'currenttable'):
            return
        searchstring=self.findtext.get()
        if self.currenttable!=None:
            self.currenttable.findValue(searchstring, findagain=1)
        return

    def plot(self, event=None):
        self.currenttable.plot_Selected()
        return

    def plotSetup(self, event=None):
        self.currenttable.plotSetup()
        return

    def normal_view(self,event=None):
        self.currenttable.paging_Off()
        return

    def page_view(self,event=None):
        self.currenttable.paging = 1
        self.currenttable.redrawTable()
        return

    def return_data(self,event=None):
        """Return the data to PEAT"""
        returndata={}
        for s in self.sheets.keys():
            currtable = self.sheets[s]
            model = currtable.getModel()
            returndata[s] = model.getData()

        self.parent.templabbookdata=returndata.copy()
        self.quit()
        return

    def sheetsSelector(self, frame):
        """Show list of sheets"""        
        names = self.sheets.keys()            
        labbooklist = Pmw.ScrolledListBox(frame,              
                labelpos='nw',
                label_text='Sheets',
                listbox_height = 8)
        labbooklist.pack(fill=BOTH,expand=1)
        labbooklist.setlist(names)
        return labbooklist
        
    def about_Labbook(self):
        self.ab_win=Toplevel()
        self.ab_win.geometry('+100+350')
        self.ab_win.title('About Labbook')

        import PEAT_images
        logo = PEAT_images.labbook_logo()
        label = Label(self.ab_win,image=logo)
        label.image = logo
        label.grid(row=0,column=0,sticky='news',padx=4,pady=4)

        text=['Labbook ','Part of Protein Engineering and Analysis Tool',
                         'A spreadsheet-like table for PEAT.']
        row=1
        for line in text:
            tmp=Label(self.ab_win,text=line)
            tmp.grid(row=row,column=0,sticky='news',padx=4)
            row=row+1
        return

    def online_documentation(self,event=None):
        """Open the online documentation"""
        import webbrowser
        link='http://enzyme.ucd.ie/PEAT/'
        webbrowser.open(link,autoraise=1)
        return

    def quit(self, event=None):
        self.labbook_win.destroy()

        return

class ToolBar(Frame):
    """Uses the parent instance to provide the functions"""
    def __init__(self, parent=None, parentapp=None):
        Frame.__init__(self, parent, width=600, height=40)
        import Table_images
        self.parentframe = parent
        self.parentapp = parentapp
        #add buttons
        img = Table_images.new_proj()
        hlptxt="Create a new project"
        self.add_button('New Project', self.parentapp.new_project, img, hlptxt)
        img = Table_images.open_proj()
        hlptxt="Open a saved project"
        self.add_button('Open Project', self.parentapp.open_project, img, hlptxt)
        img = Table_images.save_proj()
        hlptxt="Save current project"
        self.add_button('Save Project', self.parentapp.save_project, img, hlptxt)
        img = Table_images.add_row()
        hlptxt="Add a new row/record"
        self.add_button('Add record', self.parentapp.add_Row, img, hlptxt)
        img = Table_images.add_col()
        hlptxt="Add a new column"
        self.add_button('Add col', self.parentapp.add_Column, img, hlptxt)
        img = Table_images.del_row()
        hlptxt="Delete selected row/record"
        self.add_button('Delete record', self.parentapp.delete_Row, img, hlptxt)
        img = Table_images.del_col()
        hlptxt="Delete selected column"
        self.add_button('Delete col', self.parentapp.delete_Column, img, hlptxt)
        img = Table_images.search()
        hlptxt="Filter records"
        self.add_button('Filter', self.parentapp.showFilteringBar, img, hlptxt)
        img = Table_images.plot()
        hlptxt="Plot selected cells using pylab"
        self.add_button('Plot', self.parentapp.plot, img, hlptxt)
        img = Table_images.plotprefs()
        hlptxt="Plot options"
        self.add_button('Plot Prefs', self.parentapp.plotSetup, img, hlptxt)
        img = Table_images.prefs()
        hlptxt="Table preferences"
        self.add_button('Table Prefs', self.parentapp.showPrefsDialog, img, hlptxt)        
        return

    def add_button(self, name, callback, img=None, helptxt=None):
        if img==None and helptxt==None:
            b = Button(self, text=name, command=callback,
                         relief=GROOVE)
        else:
            b = Button(self, text=name, command=callback,
                             relief=GROOVE,
                             image=img)
        b.image = img
        b.pack(side=LEFT, padx=2, pady=2, ipadx=3, ipady=3)
        #add balloon
        if helptxt!=None:
            balloon=Pmw.Balloon(self.parentframe)
            balloon.bind(b, helptxt)

        return

from Tables import TableCanvas, ColumnHeader

class LabbookTable(TableCanvas):
    """Sub-class of Tablecanvas, with some changes in behaviour to make
       a more spreadsheet-like interface for labbook"""
    def __init__(self, parent=None, model=None, sheetname=None, peat=None, peatinfo=None):
        TableCanvas.__init__(self, parent, model)
        if sheetname != None:
            self.sheetname=sheetname
        if parent != None:
            self.parent = parent
        else:
            self.parent = None
        #inside peat?
        if peat != None:
            self.peatapp = peat
            self.peatinfo = peatinfo 
            self.DB = peat.DB
        else:
            self.peatapp = None
            self.peatinfo = None
            self.DB = None
        self.cellbackgr = '#FFFAF0'
        self.entrybackgr = 'white'
        self.selectedcolor = 'yellow'
        self.rowselectedcolor = '#B0E0E6'
        self.multipleselectioncolor = '#ECD672'
        self.fileicons = {'txt': 'txticon','img':'imageicon','pdf':'pdficon',
                        'other':'blankicon','zip':'zipicon','doc':'docicon',
                        'xls':'ssheeticon'}
        self.columnactions['File'] = {"Add File": 'addFile' }# , "View File" : 'displayFile' }
        self.columnactions['Table'] = {"Edit Table": 'edit_Table' }
        self.columnactions['Ekin'] = {"Open Ekin": 'open_Ekin' }

        return

    def draw_file_tooltip(self, row, col):
        """Draw a tooltip for image and file previews"""
        import os
        h=self.height
        x1,y1,x2,y2 = self.getCellCoords(row,col)
        w=x2-x1
        absrow = self.get_AbsoluteRow(row)
        text = self.model.getValueAt(absrow,col)
        if text == '' or text == None:
            return

        #actual blob data is stored in parent db blob attr
        cell = self.model.getCellRecord(absrow,col)
        key = cell['filename']
        rec = self.DB.blobs[key]
        myblob = rec.blob

        if myblob != None:
            f = myblob.open("r")
            filename = f.name
        else:
            filename=None

        extension = text.split('.')[-1]
        if extension in ['jpg','jpeg','tiff','tif','gif','ps','png','ico']:
            extension = 'img'
        if extension in ['gz','bz','jar']:
            extension = 'zip'
        if not self.fileicons.has_key(extension):
            extension = 'other'

        #Check if its an image and try to make thumbnail.. needs PIL
        if extension == 'img' and filename != None:
            self.pic = self.get_image_thumb(filename)
            isthumb=1
        #if can't get thumb, just create icon for the file type
        if self.pic == None:
            import File_images
            func = getattr(File_images, self.fileicons[extension])
            self.pic = func()
            isthumb=0

        imgobj = self.create_image(x1+w/1.5,y2+2, image=self.pic, anchor='nw')
        if isthumb == 1:
            box = self.bbox(imgobj)
            xa=box[0]-1
            ya=box[1]-1
            xb=box[2]+1
            yb=box[3]+1
            self.create_rectangle(xa,ya,xb,yb,tag='tooltip',fill='black')
            self.lift(imgobj)

        #now draw text
        import tkFont
        sfont = tkFont.Font (family='Arial', size=12, weight='bold')
        obj=self.create_text(x1+w/1.5,y2,text=text,
                                anchor='sw',
                                font=sfont,tag='tooltip')
        box = self.bbox(obj)
        x1=box[0]
        y1=box[1]
        x2=box[2]+1
        y2=box[3]+1
        self.create_rectangle(x1+1,y1+1,x2+1,y2+1,tag='tooltip',fill='black')
        self.create_rectangle(x1,y1,x2,y2,tag='tooltip',fill='lightyellow')
        self.lift(obj)
        return

    def get_image_thumb(self, infile):
        """Create an image thumbnail for img types"""
        #print 'getting thumbnail......'
        try:
            from PIL import Image, ImageTk
        except:
            return None
        import os
        size = 200,200
        file, ext = os.path.splitext(infile)
        im = Image.open(infile)
        im.thumbnail(size, Image.ANTIALIAS)
        tkim = ImageTk.PhotoImage(im)
        #im.save(file + ".thumbnail", "GIF")
        return tkim

    def handle_double_click(self, event):
        """Do double click stuff. Selected row/cols will already have
           been set with single click binding"""
        row = self.get_row_clicked(event)
        col = self.get_col_clicked(event)
        absrow = self.get_AbsoluteRow(row)
        model=self.getModel()
        coltype = model.getColumnType(col)
        cellvalue = self.model.getCellRecord(absrow, col)
        if Formula.isFormula(cellvalue):
            self.formula_Dialog(row, col, cellvalue)

        #if coltype == 'text' or coltype == 'number':
        #    self.draw_cellentry(self.currentrow, self.currentcol)
        elif coltype == 'File':
            try:
                self.view_File(self.currentrow, self.currentcol)
            except:
                pass
        elif coltype == 'Table':
            self.edit_Table(self.currentrow, self.currentcol)
        elif coltype == 'Ekin':
            self.open_Ekin(self.currentrow, self.currentcol)
        return

    def handle_motion(self, event):
        """Handle mouse motion overridden to handle file type tooltips"""
        self.delete('tooltip')
        self.pic = None
        row = self.get_row_clicked(event)
        col = self.get_col_clicked(event)
        if col == None:
            return
        if self.check_PageView(row) == 1:
            return
        coltype = self.model.getColumnType(col)
        if 0 <= row < self.rows and 0 <= col < self.cols:
            if coltype == 'File':
                self.draw_file_tooltip(row, col)
            else:
                self.draw_tooltip(row, col)
        return

    def addBlob(self, row, col, filename):
        """Adds an external file as a blob to the DB"""
        from Record import FileRecord
        f=os.path.split(filename)[1]
        model = self.getModel()
        blob = self.getBlob(filename)
        f = FileRecord(name=filename, blob=blob)
        rec = {'text': f.name, 'filename': filename}
        model.setValueAt(rec, row, col)
        #we add the blob to the parent DB blob struct rather than to table model
        self.DB.blobs[filename] = f
        #print f, blob, model.getCellRecord(row,col)
        return

    def getBlob(self, filename):
        """Create a blob from a file"""
        from ZODB.blob import Blob
        myblob = Blob()
        b=myblob.open('w')
        o=open(filename)
        data = o.read()
        b.write(data)
        b.close()
        return myblob

    def addFile(self, row, col):
        """Uses a new method of file storage with blobs
           We just get the blob from the filerecord object"""
        from Extfile import FileHandler
        model = self.getModel()
        celldata = model.getCellRecord(row, col)
        fh = FileHandler(self)
        newfile = fh.addFileDialog()

        print newfile
        if newfile != None:
            self.addBlob(row, col, fh.newfile)
            self.redrawTable()
            print model.getValueAt(row, col)
        return

    def displayFile(self, row, col):
        """Uses a new method of file storage with blobs
           We just get the blob from the filerecord object"""
        from Extfile import FileHandler
        model = self.getModel()
        celldata = model.getCellRecord(row, col)
        if celldata != None:
            filerec=celldata
        else:
            return
        fh = FileHandler(self)
        fh.displayFile(filerec)
        return

    def edit_Table(self, row, col):
        """Edit a nested table.
        Done by opening up another instance of labbook and
        returning the data to the parent"""

        #get cell data and pass to labbook as a new table
        model = self.getModel()
        tabledata = model.getCellRecord(row,col)
        refname = str(row)+'_'+str(col)
        if self.peatinfo != None:
            refname =  self.peatinfo['record']+'_'+self.peatinfo['column']+'_'+refname
        from Labbook import LabbookApp
        if tabledata:
            self.sublabbook = LabbookApp(parent=self, title=refname, data=tabledata)
        else:
            self.sublabbook = LabbookApp(parent=self, title=refname)
        self.master.wait_window(self.sublabbook.labbook_win)

        # did we get any data back?
        if not hasattr(self,'templabbookdata'):
            self.templabbookdata=None
        else:
            model.setValueAt(self.templabbookdata, row, col)
            self.redrawTable()
            self.templabbookdata=None

        return

    def open_Ekin(self, row, col):
        """Edit Ekin data, can be used outside PEAT, but does not contain protein and
           field references"""
        model = self.getModel()
        cell_data = model.getCellRecord(row,col)
        if not cell_data:
            ed = EkinModesDialog(title="Select Mode",
                                  parent=self.parentframe)
            mode = ed.result
            print 'MODE',mode
            E = None
        else:
            #mode=cell_data['mode']
            #ekindata=cell_data['ekindata']
            E = cell_data
        from PEATDB.Ekin.Ekin_main import EkinApp
        EK=EkinApp(parent=self, project=E, allowreturntoDB=1)
        self.Ekinprj=None
        self.master.wait_window(EK.ekin_win)
        if not getattr(self, 'Ekinprj', None):
            self.Ekinprj = None
        else:
            #newdata={}
            #newdata['mode']=mode
            #newdata['ekindata']=self.Ekindata
            model.setValueAt(self.Ekinprj, row, col)
            self.Ekinprj=None
            self.redrawTable()
        return

class EkinModesDialog(tkSimpleDialog.Dialog):
    def __init__(self, parent, title=None):
        self.ekinmodes = ['General', 'Simple enzyme kinetic',
                    'Enzyme pH-activity', 'pH-stability',
                    'NMR titration', 'Protein stability',
                    'Amyloid formation']
        tkSimpleDialog.Dialog.__init__(self, parent, title)

    def body(self, master):

        import Pmw
        self.themode=StringVar()
        self.mode_choice = Pmw.OptionMenu(
            parent = master,labelpos = 'w',
            label_text = 'Enter mode:',
            menubutton_textvariable = self.themode,
            items = self.ekinmodes,
            initialitem = 'General',
            menubutton_width = 16)
        self.mode_choice.pack(pady=5)
        return self.mode_choice

    def apply(self):
        self.result = self.themode.get()
        return

from TableModels import TableModel
from TableFormula import Formula
from types import *

class LabbookTableModel(TableModel):
    """A model for managing the data in a the Labbook Table.
       Subclassed from TableModel."""

    def __init__(self, newdict=None, rows=None, columns=None):
        TableModel.__init__(self, newdict, rows, columns)
        #default types for this model
        self.defaulttypes = ['text', 'number', 'File', 'Table', 'Ekin']
        return

    def getRecordAttributeAtColumn(self, rowIndex=None, columnIndex=None,
                                        recName=None, columnName=None):
         """Get the attribute of the record at the specified column index.
            This determines what will be displayed in the cell"""

         value = None   # Holds the value we are going to return
         if columnName != None and recName != None:
             if not self.data[recName].has_key(columnName):
                 return ''
             cell = self.data[recName][columnName]
         else:
             cell = self.getCellRecord(rowIndex, columnIndex)
             columnName = self.getColumnName(columnIndex)

         if cell == None:
            return ''
         # Set the value based on the data record field

         coltype = self.columntypes[columnName]
         if Formula.isFormula(cell) == True:  #change this to e.g. cell.isFormula() ?
             value = self.doFormula(cell)
             return value
         #if not type(cell) is DictType:
         if coltype == 'text' or coltype == 'Text':
             value = cell
         elif coltype == 'number':
             value = str(cell)
         elif coltype == 'File':
             value = cell['text']
         elif coltype == 'Table':
             value = str(len(cell.keys())) + ' sheets'
         elif coltype == 'Ekin':
             value = self.ekin_show_length(cell)
         else:
             value=''
         if value==None:
             value=''

         return value

    def setValueAt(self, value, rowIndex, columnIndex, **kwargs):
        """Changed the dictionary when cell is updated by user - overridden"""
        name = self.reclist[rowIndex]
        colname = self.getColumnName(columnIndex)
        coltype = self.columntypes[colname]
        if coltype == 'number':
            try:
                if value == '': #need this to allow deletion of values
                    self.data[name][colname] = ''
                else:
                    self.data[name][colname] = float(value)
            except:
                pass
        else:
            self.data[name][colname] = value
        return

    def ekin_show_length(self,data):
        """Show only the number of data points/records, taken from DB_actions"""

        return '%s recs' %data.length

    def get_external_filenames(self):
        """Returns all external filenames stored"""
        #find cols of type File
        extfilenames= []
        filecols = []
        for key in self.columntypes.keys():
            if self.columntypes[key]=='File':
                filecols.append(key)
        if len(filecols) == 0:
            return None
        else:
            for colname in filecols:
                for recname in self.reclist:
                    if self.data[recname].has_key(colname):
                        if self.data[recname][colname].has_key('location'):
                            loc = self.data[recname][colname]['location']
                            if loc != None and loc != '':                           
                                extfilenames.append(self.data[recname][colname]['location'])

        return extfilenames

class Labbook(object):
    """Labbook object mainly used from other scripts"""
    def __init__(self, data=None):
        self.data=data
        self.filename=None
        if data != None:
            self.sheets = data.keys()
        return
       
    def load(self, filename):
        """Load from file"""
        self.data = pickle.load(open(filename,'r'))
        self.filename = filename
        return
        
    def save(self, filename=None):
        """Save labbook to a file"""
        if filename == None: filename=self.filename
        fd=open(filename,'w')        
        pickle.dump(self.data,fd)
        fd.close()        
        return

    def getSheet(self, name):
         """Get a sheet - returns a tablemodel"""
         return TableModel(self.data[name])
         
    def saveSheet(self, name, model=None, data=None):
        """Save sheet using model or dict"""
        if model != None:
            self.data[name] = model.getData()
        return
        
# Main function, run when invoked as a stand-alone Python program.

def main():
    "Run the application."

    import sys, os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="tablefile",
                        help="Open a table file", metavar="FILE")
    opts, remainder = parser.parse_args()
    if opts.tablefile != None:
        app=LabbookApp(datafile=opts.tablefile)
    else:
        app=LabbookApp()
    app.mainloop()
    return

if __name__ == '__main__':
    main()
