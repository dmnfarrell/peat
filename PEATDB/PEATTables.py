#
# This file is part Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#
#
# Written by D Farrell, Jan 2008
# Updated Dec 2009 for new ZODB format
#

from Tkinter import *
import tkSimpleDialog, tkFileDialog, tkMessageBox
from Tables import TableCanvas
from TableModels import TableModel
from TableFormula import Formula
from PEATDB.Ekin.Base import EkinProject
from Base import zDatabase
from Record import PEATRecord


class PEATTable(TableCanvas):
    """Sub-class of Tablecanvas, with some additions for PEAT main table
       We mostly override the methods that deal with the extra PEAT fields"""
    def __init__(self, parent=None, model=None, parentapp=None):
        if parentapp != None:
            self.parentapp = parentapp
        TableCanvas.__init__(self, parent, model)
        self.sortcol = 0
        self.grid_color = '#8996ab'       
        self.rowselectedcolor = '#f6fcd1'

        self.peatactions = {'PDB' : {"Save PDB file":'save_structure',
                                    "View Structure": 'displayStructure',
                                    "View PDB Text": 'viewStructureText'},
                    'Notes' : {"Edit notes": 'edit_notes' },
                    'Link' : {"Edit Link" : 'edit_link', "Open link" : 'open_link'},
                    'File' : {"View File" : 'display_file', "Add File" : 'add_file',
                            "Save File" : 'save_file'},
                    'Table' : {"Edit Table" : 'startLabbook' },              
                    'Ekintype' : {"Open Ekin" : 'startEkin' },
                    'Sequence' : {"Show sequence" : 'display_sequence'},           
                    }

        self.ekin_actions={'NMR titration':'NMR titration',
                              'Protein Stability':'Protein Stability',
                              'pH Activity profile':'Enzyme pH-activity',
                              'Michaelis-Menten kinetics':'Simple enzyme kinetic',
                              'pH stability profile':'pH-stability',
                              'Amyloid formation':'Amyloid formation',
                               'General':'General'}

        self.fileicons = {'txt': 'txticon','img':'imageicon','pdf':'pdficon',
                        'other':'blankicon','zip':'zipicon','doc':'docicon',
                        'xls':'ssheeticon'}
        self.thumbsize = 200
        return

    def do_bindings(self):
        """Bind keys and mouse clicks - overriden"""
        self.bind("<Button-1>",self.handle_left_click)
        self.bind("<Double-Button-1>",self.handle_double_click)
        self.bind("<Control-Button-1>", self.handle_left_ctrl_click)
        self.bind("<Shift-Button-1>", self.handle_left_shift_click)

        self.bind("<ButtonRelease-1>", self.handle_left_release)
        if self.ostyp=='mac':
            #For mac we bind Shift, left-click to right click
            self.bind("<Button-2>", self.handle_right_click)
            self.bind('<Shift-Button-1>',self.handle_right_click)
        else:
            self.bind("<Button-3>", self.handle_right_click)

        self.bind('<B1-Motion>', self.handle_mouse_drag)
        self.bind('<Motion>', self.handle_motion)
        self.bind("<Tab>", self.gotonextCell)

        self.bind("<Control-x>", self.delete_Row)
        self.bind("<Control-n>", self.add_Row)
        self.bind("<Delete>", self.delete_Cells)
        #bind to parentapp instead of parentframe for PEAT
        self.parentapp.master.bind("<Right>", self.handle_arrow_keys)
        self.parentapp.master.bind("<Left>", self.handle_arrow_keys)
        self.parentapp.master.bind("<Up>", self.handle_arrow_keys)
        self.parentapp.master.bind("<Down>", self.handle_arrow_keys)
        self.parentapp.master.bind("<KP_8>", self.handle_arrow_keys)

        self.parentframe.master.bind("<Return>", self.handle_arrow_keys)
        self.parentframe.master.bind("<Tab>", self.handle_arrow_keys)
        if 'windows' in self.platform:
            self.bind("<MouseWheel>", self.mouse_wheel)
        self.bind('<Button-4>', self.mouse_wheel)
        self.bind('<Button-5>', self.mouse_wheel)
        return

    def add_Row(self, recname=None):
        """Add a new rec in PEAT, need special dialog"""
        def checkrow_name(rowname):
            if rowname == '':
                tkMessageBox.showwarning("Whoops", "Name should not be blank.")
            if self.getModel().data.has_key(rowname):
                tkMessageBox.showwarning("Name exists", "Record already exists!")

        if recname == None:
            recname = tkSimpleDialog.askstring("New rec name?", "Enter rec name:")
        checkrow_name(recname)
        self.model.addRow(recname)
        self.setSelectedRow(self.model.getRecordIndex(recname))
        self.redrawTable()
        return

    def add_Column(self, newname=None, fieldtype=None, defaultval=None):
        """Add a new column"""
        self.parentapp.addFieldDialog()
        return

    def delete_Column(self):
        """Delete currently selected column"""
        col = self.getSelectedColumn()
        colname = self.model.getColumnName(col)
        #if colname in self.getModel().static_fields.keys():
        if colname == 'name':
            tkMessageBox.showwarning("Static field", "This field can't be deleted.")
            return
        n =  tkMessageBox.askyesno("Delete",  "Delete This Column?")
        if n:
            self.model.deleteColumn(col)            
            self.redrawTable()
        return

    def handle_double_click(self, event):
        """Do double click stuff. Selected row/cols will already have
           been set with single click binding"""
        row = self.get_row_clicked(event)
        col = self.get_col_clicked(event)
        absrow = self.get_AbsoluteRow(row)
        if col == None:
            return
        model = self.getModel()
        cellvalue = model.getCellRecord(absrow, col)
        coltype = model.getColumnType(col)
        colname = model.getColumnName(col)
        protein = model.getRecName(absrow)

        def do_peatcommand(fieldtype):
            functions = self.peatactions[fieldtype]
            kys = functions.keys()
            default = kys[0]
            func = getattr(self.parentapp, functions[default])
            if fieldtype == 'Ekintype':
                ekinmode = self.ekin_actions[coltype]
                func(protein=protein, field_name=colname, mode=ekinmode)
            else:
                func(protein=protein, field_name=colname)
            return
        if Formula.isFormula(cellvalue):
            self.formula_Dialog(row, col, cellvalue)
        elif coltype in self.model.ekintypes:
            do_peatcommand('Ekintype')
        elif coltype.lower() == 'text':
            self.draw_cellentry(row, col)
        elif self.peatactions.has_key(coltype):
            do_peatcommand(coltype)
        else:
            raise Exception('Unknown field')
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

    def popupMenu(self, event, rows=None, cols=None, outside=None):
        """Add right click behaviour for canvas"""
        if outside == None:
            row = self.get_row_clicked(event)
            col = self.get_col_clicked(event)
            coltype = self.model.getColumnType(col)
        popupmenu = Menu(self, tearoff = 0)
        def popupFocusOut(event):
            popupmenu.unpost()

        if outside == 1:
            pass
        else:
            def add_peatcommand(fieldtype):
                """Adds all commands per column type based on self.peatactions"""
                model = self.getModel()
                colname = model.getColumnName(col)
                protein = model.getRecName(row)
                functions = self.peatactions[fieldtype]

                for f in functions.keys():
                    func = getattr(self.parentapp, functions[f])
                    if f == 'Open Ekin':
                        ekinmode = self.ekin_actions[coltype]
                        def thisfunc(self=self,function=func,protein=protein,field_name=colname):
                            return function(protein=protein,field_name=field_name,mode=ekinmode)
                    else:
                        def thisfunc(self=self,function=func,protein=protein,field_name=colname):
                            return function(protein=protein,field_name=field_name)
                    popupmenu.add_command(label=f, command=thisfunc)
                return

            def add_commands(fieldtype):
                """Add from columnactions if present"""
                functions = self.columnactions[fieldtype]
                for f in functions.keys():
                    func = getattr(self, functions[f])
                    def thisfunc(self=self,function=func,row=row,col=col):
                        return function(row=row,col=col)
                return

            def add_defaultcommands():
                """now add general actions for all cells"""
                order = [ "Set Fill Color","Set Text Color","Edit","Copy","Paste","Fill Down","Fill Right", "Clear Data",
                         "View Record", "Select All",
                         "Plot Selected","Plot Options"]
                defaultactions={"Edit" : lambda: self.draw_cellentry(row,col),
                                "Set Fill Color" : lambda : self.setcellColor(rows,cols,key='bg'),
                                "Set Text Color" : lambda : self.setcellColor(rows,cols,key='fg'),
                                "Copy" : lambda : self.copyCell(rows, cols),
                                "Paste" : lambda : self.pasteCell(rows, cols),                                
                                "Fill Down" : lambda : self.fill_down(rows, cols),
                                "Fill Right" : lambda : self.fill_across(cols, rows),
                                "View Record" : lambda : self.getRecordInfo(row),
                                "Clear Data" : lambda : self.delete_Cells(rows, cols),
                                "Select All" : self.select_All,
                                "Plot Selected" : self.plot_Selected,
                                "Plot Options" : self.plotSetup }
                for action in order:
                    if action == 'Fill Down' and (rows == None or len(rows) <= 1):
                        continue
                    if action == 'Fill Right' and (cols == None or len(cols) <= 1):
                        continue
                    else:
                        popupmenu.add_command(label=action, command=defaultactions[action])
                return

            if self.columnactions.has_key(coltype):
                add_commands(coltype)
            if self.peatactions.has_key(coltype):
                add_peatcommand(coltype)
            elif coltype in self.model.ekintypes:
                add_peatcommand('Ekintype')
            add_peatcommand('Sequence')
            add_defaultcommands()

        popupmenu.add_command(label="Show Prefs", command=lambda : self.showtablePrefs())
        popupmenu.add_command(label="Export Table", command= self.exportTable)
        popupmenu.add_command(label="Filter Recs", command= self.showFilteringBar)
        popupmenu.bind("<FocusOut>", popupFocusOut)
        popupmenu.focus_set()
        popupmenu.post(event.x_root, event.y_root)
        return popupmenu

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

        fopen = None
        rec = self.model.getCellRecord(absrow,col)
        myblob = rec.blob

        extension = text.split('.')[-1]
        if extension in ['jpg','jpeg','tiff','tif','gif','ps','png','ico']:
            extension = 'img'
        if extension in ['gz','bz','jar']:
            extension = 'zip'
        if not self.fileicons.has_key(extension):
            extension = 'other'

        #Check if its an image and try to make thumbnail.. needs PIL
        if extension == 'img':
            if myblob != None:
                fopen = myblob.open("r")
                filename = fopen.name
            else:
                filename=None
            if filename != None:
                try:
                    self.pic = self.get_image_thumb(filename)
                    isthumb=1
                except:
                    fopen.close()
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
        if fopen != None:  fopen.close()
        return

    def get_image_thumb(self, infile):
        """Create an image thumbnail for img types"""
        #print 'getting thumbnail......'
        try:
            from PIL import Image, ImageTk
        except:
            print 'you need PIL'
            return None
        import os
        size = self.thumbsize, self.thumbsize
        file, ext = os.path.splitext(infile)        
        im = Image.open(infile)
        im.thumbnail(size, Image.ANTIALIAS)
        tkim = ImageTk.PhotoImage(im)
        return tkim

#
# Model class
#

class PEATTableModel(TableModel):
    """A model for managing the data in a the PEAT Table.
       Subclassed from TableModel.
       The actual data will always be passed from PEAT, so there must
       always be a DB instance passed to this constructor.
       The data structure is an OOBTree instance"""

    def __init__(self, DB=None, parentapp=None):
        """Constructor overridden"""
        import copy
        self.parentapp=parentapp
        self.DB = DB
        self.data = DB.data    #holds the data
        self.meta = DB.meta    #holds table formatting

        self.colors = {}
        self.colors['fg']={}
        self.colors['bg']={}
        self.columnlabels={}
        self.columnwidths={} 
        #list of editable column types
        self.editable={}
        #columntypes from PEAT
        self.ekintypes = ['General', 'NMR titration', 'Protein Stability',
                      'pH Activity profile', 'Michaelis-Menten kinetics',
                      'pH stability profile', 'Amyloid formation']
        self.ekinfields = EkinProject.ekinfields + ['mode']

        if self.data != None:           
            self.userfields = self.meta.userfields          
            self.columnNames = self.meta.staticfields.keys()
            self.columnOrder = None
            for key in self.userfields:
                if not self.userfields[key].has_key('show'):
                    self.userfields[key]['show'] = True
                if self.userfields[key]['show'] == True:
                    self.columnNames.append(key)

            #add or get colors dict
            self.update_colors()

            #finally set up column types dict
            self.columntypes=copy.deepcopy(self.meta.staticfields)
            for key in self.userfields:
               self.columntypes[key]=self.userfields[key]['field_type']               

        #Store the sorted list of names
        self.reclist = list(self.data.keys())

        self.filteredrecs = None
        #sort by name
        #self.setSortOrder(0)
        #restore col order
        if self.columnOrder:
            self.columnNames=[]
            for i in self.columnOrder.keys():
                self.columnNames.append(self.columnOrder[i] )
                i=i+1
        self.defaulttypes = ['text', 'number','Structure', 'Mutations']
        return

    def update_reclist(self, keys, sortcol=None):
        """Update the reclist when DB has been changed externally,
           ie. not from the model directly"""
        self.reclist = self.data.keys()
        if sortcol!=None:        
            self.setSortOrder(sortcol)
        return

    def update_columnNames(self):
        """Update the columns from userfields when data has been changed
           externally"""
        self.userfields = self.meta.userfields
        currentkeys = self.userfields.keys()
        currentkeys.extend(self.meta.staticfields.keys())
        #add any keys not present in new userfields
        for key in self.userfields:            
            if not key in self.columnNames:
                self.columnNames.append(key)
                self.columnlabels[key]=key                
                self.columntypes[key]=self.userfields[key]['field_type']
            if self.userfields[key]['show'] == False:
                self.columnNames.remove(key)                
        #remove keys if not in userfields anymore
        for key in self.columnNames:    
            if not key in currentkeys:
                self.columnNames.remove(key)
                del self.columnlabels[key]
                del self.columntypes[key]
        return

    def update_colors(self):
        if self.meta.has_key('table') and self.meta['table'].has_key('colors'):
            self.colors = self.meta['table']['colors']
            self.columnlabels = self.meta['table']['columnlabels']
            if self.meta['table'].has_key('columnorder'):
                self.columnOrder = self.meta['table']['columnorder']
        else:
            self.meta['table']={}
            self.meta['table']['colors'] = self.colors
            for colname in self.columnNames:
                self.columnlabels[colname]=colname
            self.meta['table']['columnlabels'] = self.columnlabels

    def save_columnOrder(self):
        """Save the column order, done before peat db saved"""
        #save current col order
        data['columnorder']={}
        i=0
        for name in self.columnNames:
            data['columnorder'][i] = name
            i=i+1
        return

    def getData(self):
        """Return the current data for saving, overrides base method"""
        import copy
        data = copy.deepcopy(self.data)
        data['colors'] = self.colors
        data['columnlabels'] = self.columnlabels
        #save current col order
        self.save_columnOrder()
        return data

    def setValueAt(self, value, rowIndex, columnIndex):
        """Changed the dictionary when cell is updated by user"""
        name = self.reclist[rowIndex]
        colname = self.getColumnName(columnIndex)
        coltype = self.columntypes[colname]       
        if coltype == 'number':
            try:
                if value == '': 
                    self.data[name][colname] = ''
                else:
                    self.data[name][colname] = float(value)
            except:
                pass
        else:
            self.data[name][colname] = value  
        return
    
    def deleteCellRecord(self, rowIndex, columnIndex):
        """Remove the cell data at this row/column"""
        colname = self.getColumnName(columnIndex)
        coltype = self.columntypes[colname]
        name = self.reclist[rowIndex]       
        if self.data[name].has_key(colname):           
            self.data[name].delete(colname)            
        return

    def getRecordAttributeAtColumn(self, rowIndex=None, columnIndex=None,
                                         recName=None, columnName=None):
        """Get the attribute of the record at the specified row/column index.
            This determines what will be displayed in the cell"""

        value = None
        if columnIndex == None and columnName != None:
            colname = columnName
        else:
            colname = self.getColumnName(columnIndex)
        coltype = self.columntypes[colname]
        if rowIndex == None and recName != None:
            name = recName
        else:
            name = self.getRecName(rowIndex)

        p = self.data[name]
        #display info is now handled in the record class
        #try:           
        value = p.getDisplayAttribute(colname)
        #except:
        #    value = ''

        return value

    def getKeys(self, name, colname):
        """Get Userdict keys"""

        if not self.data[name].has_key(colname):
            return None

        if isinstance(self.data[name], dict) or isinstance(self.data[name], PEATRecord):
            numkeys = len(self.data[name][colname])
        else:
            numkeys = 0
        return numkeys

    def ekin_show_length(self, data):
        """Show only the number of data points/records, taken from DB_actions"""
        dkeys=data.keys()
        count=0
        for key in dkeys:
            if len(key)>2:
                if key[:2]!='__':
                    count=count+1
            else:
                count=count+1
        dkeys.sort()
        dkeys_txt=''
        for key in dkeys:
            dkeys_txt=dkeys_txt+'%s\n' %key

        return '%s recs' %count

    def setSortOrder(self, columnIndex, reverse=0):
        """Changes the order that records are sorted in, which will
           be reflected in the table upon redrawing - overridden"""
        super(PEATTableModel, self).setSortOrder(columnIndex, reverse)
        '''for key in self.nodisplay:
            if key in self.reclist:
                self.reclist.remove(key)'''
        return

    def addRow(self, name=None, rectype=None):
        """Add a row"""
        index = self.getRowCount() + 1
        if name == None:
            name=str(index)
        if rectype == None:
            rectype='text'
        self.data[name]={}
        self.data[name]['name'] = name
        self.data[name]['field_type'] = rectype
        self.data[name]['Structure'] = None
        self.reclist = self.data.keys()
        
        return

    def addColumn(self, colname=None, coltype='text'):
        """Add a column"""
       
        self.DB.addField(colname, coltype)
        self.update_columnNames()
        return

    def deleteColumn(self, columnIndex):
        """delete a column"""
        colname = self.getColumnName(columnIndex)
        print colname
        self.columnNames.remove(colname)
        del self.columnlabels[colname]
        del self.columntypes[colname]
        self.DB.deleteField(colname)

        if hasattr(self, 'sortcolumnIndex') and columnIndex == self.sortcolumnIndex:
            self.setSortOrder()
        print 'column deleted'
        print 'new cols:', self.columnNames
        return

    def display_mutations(self,protein):
        """Display the mutations for this record as compared to the referece protein"""
        if self.parentapp == None:
            return 'NA'
        refprot=self.parentapp.refprotein.get()
        if refprot is None or not self.data.has_key(refprot):
            return 'NA'

        # Get the mutations
        DBi=self.parentapp.data['DBinstance']
        is_parent,operations=DBi.isparent(protein,refprot)
        operations=DBi.shorten_operations(operations)
        if is_parent:
            import string
            mut_text=string.join(operations,'+')
        else:
            mut_text='NA'
        return mut_text

    def resetcolors(self):
        """Remove all color formatting - overriden """
        self.colors={}
        self.colors['fg']={}
        self.colors['bg']={}
        self.meta['table']['colors'] = self.colors
        self.meta._p_changed = 1
        return

    def setColorAt(self, rowIndex, columnIndex, color, key='bg'):
        """Override to set dirty bit"""
        TableModel.setColorAt(self, rowIndex, columnIndex, color, key=key)
        self.meta._p_changed = 1
        return

    def simpleCopy(self, include=None):
        """Return a simple copy of the PEAT table model, with only
           text and int fields and not using the OOBTree structure"""
        M=TableModel()
        fields = self.DB.getSimpleFields()
        if include!=None:
            fields.extend(include)
        for rec in self.reclist:
            data={}
            for f in fields:
                if self.data[rec].has_key(f):
                   data[f] = self.data[rec][f]
            M.addRecord(rec,**data)        
        return M
