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
import math, sys, os, types
import pickle
import csv
import random, string
from PEATDB.Tables import TableCanvas, ColumnHeader
from PEATDB.TableModels import TableModel
from PEATDB.GUI_helper import *


class DictModel(TableModel):
    """Table model for standard dictionary assuming a dict of subdicts"""
    def __init__(self, dictdata=None):
        self.dict = dictdata
        TableModel.__init__(self)
        if dictdata == None:
            return        
        for key in dictdata:
            #assumes D is dict with some text fields
            D = dictdata[key]           
            if D.has_key('name'):
                del D['name']
            self.addRecord(name=key, **D)
            
        return

    def importDict(self, importdata, namefield='name'):
        """Import list of dicts, each one is a record"""
        if len(importdata) == 0:
            return        
        fields = importdata[0].keys()
        if not namefield in fields or not 'name' in fields:
            print 'no such field for keyname field'
            namefield=None
        for D in importdata:
            if namefield==None:
                name = ''.join(random.choice(string.digits) for x in range(8))
            else:    
                name = D[namefield]          
                del D[namefield]
            self.addRecord(name=name,**D)
        return
    
class DictEditor(Frame, GUI_help):
    """Simple app to allow edit/save dicts in table"""

    def __init__(self, parent=None):
        self.parent = parent
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
        else:
            self.main=Toplevel()
            
        self.main.title('Dict Editor')
        self.main.geometry('600x600+200+200')        
        bf=Frame(self.main)
        bf.pack(side=RIGHT)
        Button(bf,text='Load Pickle',command=self.loadDict).pack(side=TOP,fill=X)
        Button(bf,text='Save',command=self.saveDict).pack(side=TOP,fill=X)
        Button(bf,text='Import CSV',command=self.loadCSV).pack(side=TOP,fill=X)        
        Button(bf,text='Save as CSV',command=self.saveCSV).pack(side=TOP,fill=X)
        Button(bf,text='Add Row',command=self.addRow).pack(side=TOP,fill=X)
        Button(bf,text='Quit',command=self.quit).pack(side=TOP,fill=X)
        self.tframe=Frame(self.main)
        self.tframe.pack(side=LEFT,expand=1,fill=BOTH)
        self.filename = None
        self.table=None
        return

    def loadDict(self, filename=None):
        if filename==None:
            import tkFileDialog
            filename=tkFileDialog.askopenfilename(defaultextension='.pickle',
                                                  filetypes=[("pickle","*.pickle"),
                                                            ("All files","*.*")],
                                                  parent=self.main)
        if filename != None:
            self.dictdata = pickle.load(open(filename,'r'))
            self.loadTable(self.dictdata)
            self.filename = filename
        return

    def saveDict(self):
        """Pickle dict from table contents"""
        if self.table == None:
            return
        if self.filename == None:
            import tkFileDialog
            self.filename=tkFileDialog.asksaveasfilename(defaultextension='.pickle',
                                                  filetypes=[("pickle","*.pickle"),
                                                            ("All files","*.*")],
                                                  parent=self.main)            
        pickle.dump(self.model.data,open(self.filename,'w'))
        return

    def loadCSV(self):
        """Load a csv file as a dict and display"""
        from PEATDB.IO import Importer
        IM = Importer()
        IM.getLines()
        IM.showimportDialog()        
        self.wait_window(IM.ask_csv)    
        model = DictModel()
        model.importDict(IM.data)
        self.filename=None
        self.loadTable(model=model)
        return
    
    def saveCSV(self):
        """Export csv from table contents"""
        if self.filename == None:
            return
        self.table.exportTable()
        return
    
    def loadTable(self, dictdata=None, model=None):
        if model != None:
            self.model = model
        else:
            self.model = DictModel(dictdata)
        
        self.model.setSortOrder(0)
        self.table = TableCanvas(self.tframe, self.model)
        self.table.createTableFrame()
        return

    def addRow(self):
        """add a row"""
        if self.table == None:
            return
        self.table.add_Row()    
        return
        
    def quit(self):
        self.main.destroy()
        return

def main():
    """Open editor"""
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="file",
                        help="Open a pickle file")
    opts, remainder = parser.parse_args()
    app = DictEditor()
    if opts.file != None and os.path.exists(opts.file):
        app.loadDict(opts.file)
    app.mainloop()

if __name__ == '__main__':
    main()
