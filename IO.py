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

import os, types
from Tkinter import *
import Pmw
import string
import csv
import tkFileDialog, tkMessageBox, tkSimpleDialog

class Importer:
    """Generic Import helper class that can also be used outside PEAT"""
    def __init__(self, DB=None, parent=None):
        self.DB = DB
        self.parent = parent
        self.data = None
        return
   
    def doImport(self, event=None):
        """This does the import using the current dialog options
           and shows a preview. Only works within GUI context."""

        lines = self.lines
        data = []
        #set up start row        
        try:
            self.startnum = int(self.linenum.get())
        except:
            self.startnum = 0
        if self.startnum < 0 or self.startnum > len(lines):
            return -1
        #set up cols
        self.colnames = []
        s = self.separator.get()
        if s == '':
            s=' '

        line = lines[self.startnum].replace('#','')

        if self.column_row.get() == 0:
            if s == ' ':
                #splits by any non-single whitespace
                self.colnames = string.split(line)
            else:
                self.colnames=string.split(string.strip(line), s)
            self.do_checkboxes()
        #data in rows
        else:
            if self.xinfirstcol.get() == 1:   #use first row as x-column
              
                self.colnames = ['x','y']
            else:
                #x and y cols are in every row, so use manual entry of no cols
                print self.numcolsinrow.get()
                self.colnames=range(self.numcolsinrow.get())
            print 'self.colnames', self.colnames
            self.do_checkboxes()

        self.previewtext.component('columnheader').delete("1.0",END)
        headerline=''
        for c in range(len(self.colnames)):
            headerline = headerline + '%-7s   ' %self.colnames[c]
        self.previewtext.component('columnheader').insert('0.0', headerline)
        if self.column_row.get() == 0:
            #data in cols
            self.numrows = len(lines)
            self.previewtext.delete("1.0",END)
            for row in range(self.startnum, self.numrows):
                line = string.strip(lines[row])
                rowdata={}
                if self.ignorecomments.get() == 1 and line.startswith('#'):
                    continue
                else:
                    line = line.replace('#','')
                textline = ''
                if s == ' ':
                    items=string.split(line)
                else:
                    items=string.split(line, s)
                for c in range(len(items)):
                    textline = textline + '%-7s   ' %items[c]
                    colname = self.colnames[c]
                    rowdata[colname] = items[c]
                self.previewtext.insert(END, textline+ '\n')
                data.append(rowdata)
        else:
            #data in rows
            line = string.strip(lines[self.startnum])
            if s == 'tab':
                s='\t'
            xdata = string.split(line,s)
            self.numrows = len(xdata)
            print 'xdata',xdata
            print self.numrows
            datasets={}
            for row in range(self.startnum+1, len(lines)):
                line = string.strip(lines[row])
                if self.ignorecomments.get() == 1 and line.startswith('#'):
                    continue
                else:
                    line = line.replace('#','')
                items = string.split(line,s)
                name=items[0]
                data.remove(name)
                if len(data) >0:
                    datasets[name] = items
         
            self.previewtext.delete("1.0",END)
            
            #  Sort the datasets by name            
            datanames=datasets.keys()
            datanames.sort()
            for name in datanames:
                self.previewtext.insert(END, 'Dataset: ' + name + '\n')
                textline = ''
                for i in range(self.numrows):
                    try:
                        textline = str(xdata[i]) + '%+10s   ' %datasets[name][i]
                    except:
                        textline = ' '
                    self.previewtext.insert(END, textline + '\n')
                self.previewtext.insert(END, '---------------' + '\n')
            #keep ref to these 2 for import
            self.inrows_datasets = datasets
            self.inrows_xdata = xdata
        print data      

        self.data = data
        return 1
        
    def getLines(self, filename=None):
        """Get file into list"""
        if filename==None:
            filename=tkFileDialog.askopenfilename(defaultextension='.csv',
                                       initialdir=os.getcwd(),
                                       filetypes=[("csv","*.csv"),("All files","*.*")])
        if os.path.isfile(filename):
            fd=open(filename)
            lines = fd.readlines()
            fd.close()
        self.lines = lines
        return lines

    def do_checkboxes(self):
        #add check boxes for columns
        self.colchecksframe = Frame(self.ask_csv, borderwidth=2, relief=GROOVE)
        self.colschosen = []
        r=0
        c=0
        for i in range(len(self.colnames)):
            var = IntVar()
            var.set(1)
            chk = Checkbutton(self.colchecksframe, text=str(self.colnames[i])[:8],
                                variable=var)
            if c>=8:
                r=r+1
                c=0
            chk.grid(column=c,row=r)
            self.colschosen.append(var)
            c=c+1
        self.colchecksframe.grid(row=4,column=2,columnspan=3,sticky='news', padx=2,pady=2)
        #print map((lambda var: var.get()), self.colschosen)
        return
        
    def showimportDialog(self):
        """Dialog for allowing user to preview import file and make adjustments"""

        lines = self.lines
        self.ask_csv=Toplevel()
        self.ask_csv.title('Import Text/CSV')
       
        #add other options
        self.column_row=IntVar()
        self.column_row.set(0)
        self.multiple_tabs=IntVar()
        self.multiple_tabs.set(0)
        self.ignorecomments = IntVar()
        self.ignorecomments.set(1)
        self.separator=StringVar()
        self.separator.set(',')
        self.secondcolpos=IntVar() 
        self.secondcolpos.set(0)
        self.balloon = Pmw.Balloon(self.ask_csv)

        #add scrolled text areas
        self.text= Pmw.ScrolledText(self.ask_csv, labelpos = 'n',
                                    label_text='File Contents',
                                    usehullsize=1,
                                    hull_width = 300,
                                    hull_height = 300,
                                    text_padx = 4,
                                    text_pady = 4)
        self.text.grid(row=3,column=0,rowspan=2,columnspan=2,sticky='news',padx=2,pady=2)
        self.text.component('text').configure(background = '#ffffcc')
        self.previewtext= Pmw.ScrolledText(self.ask_csv, labelpos = 'n',
                                    label_text='Preview',
                                    columnheader = 1,
                                    rowheader = 1,
                                    rowcolumnheader = 1,
                                    rowheader_width = 3,
                                    Header_foreground = 'blue',
                                    rowcolumnheader_width = 3,
                                    usehullsize=1,
                                    text_wrap='none',
                                    hull_width = 300,
                                    hull_height = 300)
        self.previewtext.grid(row=3,column=2,columnspan=3,sticky='news',padx=2,pady=2)
        self.previewtext.component('text').configure(background = '#ffffcc')
        self.ask_csv.columnconfigure(2,weight=1)
        self.ask_csv.rowconfigure(3, weight=1)


        #show file contents in first pane
        def show_file():
            self.text.delete("1.0",END)
            import string
            count=1
            for line in lines[:]:
                self.text.insert(END,'LINE%d : %s\n' %(count,string.strip(line[:60])))
                count=count+1
            return

        f1=Frame(self.ask_csv, borderwidth=2, relief=GROOVE)
        Radiobutton(f1,text='Data in columns',variable=self.column_row,
                                    command=self.doImport,value=0).pack()
        Radiobutton(f1,text='Data in rows',variable=self.column_row,
                                    command=self.doImport,value=1).pack()
        f1.grid(column=0,row=0,rowspan=2,padx=4,pady=4)
        #add a frame for row options
        inrowsframe=Frame(self.ask_csv, borderwidth=2, relief=GROOVE)
        inrowsframe.grid(column=2,row=0,rowspan=3,padx=2,pady=2)
        Label(inrowsframe, text='Data in Rows Options', fg='blue').pack()
        self.rowdatastart = IntVar()
        self.rowdatastart.set(0)
        Label(inrowsframe, text='Data starts at col:').pack()
        rowdatastartEntry = Entry(inrowsframe, textvariable=self.rowdatastart, bg='white')
        rowdatastartEntry.pack()
        self.numcolsinrow = IntVar()
        self.numcolsinrow.set(0)
        Label(inrowsframe, text='No. of cols:').pack()
        numcolsinrowEntry = Entry(inrowsframe, textvariable=self.numcolsinrow, bg='white')
        numcolsinrowEntry.pack()

        self.xinfirstcol = IntVar()
        self.xinfirstcol.set(0)
        Checkbutton(inrowsframe,text='First row is x-data',
                    variable=self.xinfirstcol, command=self.doImport).pack()

        genoptsframe = Frame(self.ask_csv,relief=GROOVE)
        genoptsframe.grid(column=1,row=0,rowspan=2,padx=2,pady=2)

        sepEntry = Pmw.EntryField(genoptsframe, labelpos = 'w',
                                        label_text = 'Separator:',
                                        entry_textvariable = self.separator,
                                        validate = None, command = self.doImport)
        sepEntry.configure(entry_background = 'white')
        sepEntry.pack(fill=X,padx=2,pady=2)
        sepEntry.component('entry').bind('<KeyPress>', self.doImport)
        sepEntry.component('entry').bind('<KeyRelease>', self.doImport)
        Checkbutton(genoptsframe,text='Ignore comments',variable=self.ignorecomments,
                    command=self.doImport).pack(fill=X,padx=2,pady=2)

        show_file()

        # Variables
        self.import_start=0
        self.name_column=0
        self.names_in_sheet=IntVar()
        self.names_in_sheet.set(1)

        # Counter field
        self.linenum = Pmw.Counter(genoptsframe,
                                   increment=1,
                                   labelpos = 'w',
                                   label_text = 'Start at line',
                                   entryfield_value = "0",
                                   datatype ='numeric',
                                   entryfield_validate = self.doImport,
                                   entryfield_command = self.doImport)
        self.linenum.pack(fill=X,padx=2,pady=2)

        def cancel_import():
            #self.doImport=None
            self.ask_csv.destroy()
            return 0

        def selectNone():
            #deselect all checkboxes
            for c in self.colschosen:
                c.set(0)
            return

        imp=Button(self.ask_csv,text='Import',command=self.ask_csv.destroy,
                                              bg='lightblue',
                                              borderwidth=2, relief=GROOVE)
        imp.grid(row=6,column=2,columnspan=1,sticky='news',padx=2,pady=2)
        allb=Button(self.ask_csv,text='Select None',command=selectNone,
                                               bg='lightblue',
                                              borderwidth=2, relief=GROOVE)
        allb.grid(row=6,column=3,columnspan=1,sticky='news',padx=2,pady=2)
        canc=Button(self.ask_csv,text='Cancel',command=cancel_import,
                                               bg='lightblue',
                                              borderwidth=2, relief=GROOVE)
        canc.grid(row=6,column=4,columnspan=1,sticky='news',padx=2,pady=2)
        #self.doImport()
        return


