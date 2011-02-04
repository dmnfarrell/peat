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

import os, sys
from Tkinter import *
import Pmw
import csv
import tkFileDialog, tkMessageBox, tkSimpleDialog
from Dataset import EkinDataset

class Importer:
    """
    Import helper class for Ekin app, redone by Damien, June 08
    This class is designed to be interactively run from the GUI.
    For the standard non GUI methods, see EkinProject class
    """

    def __init__(self, parent,parent_win=None):
        self.parent = parent
        if not parent_win:
            self.ekin_win = parent.ekin_win
        else:
            self.ekin_win = parent_win 
        self.filename = None
        self.path = os.path.expanduser("~")
        return

    def get_lines(self, filename):
        """Get file into list"""
        if os.path.isfile(filename):
            fd=open(filename)
            lines = fd.readlines()
            fd.close()
        return lines

    def import_csv_dialog(self):
        """Dialog for allowing user to preview import file and make adjustments"""
        import string
        lines = self.lines
        self.do_import=True
        self.ask_csv=Toplevel()
        self.ask_csv.title('Import Text Options')
        try:
            self.parent.set_geometry(self.ekin_win, self.ask_csv)
        except:
            pass
        #add other options
        self.column_row=IntVar()
        self.column_row.set(0)
        self.multiple_tabs=IntVar()
        self.multiple_tabs.set(0)
        self.ignorecomments = IntVar()
        self.ignorecomments.set(1)
        self.separator=StringVar()
        self.separator.set(',')
        self.secondcolpos=IntVar()   #manually set pos of 2nd col
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

        def do_checkboxes():
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

        #show file contents in first pane
        def show_file():
            self.text.delete("1.0",END)
            import string
            count=1
            for line in lines[:]:
                self.text.insert(END,'LINE%d : %s\n' %(count,string.strip(line[:60])))
                count=count+1
            return

        #Func for doing preview in text area
        #This does most of the work in defining the import using the current options
        def do_preview(event=None):
            #set up start row
            try:
                self.startnum = int(self.linenum.get())
            except:
                self.startnum = 0
            try:
                self.endnum = int(self.endnum.get())
            except:
                self.endnum = 0
            if self.startnum < 0 or self.startnum > len(lines):
                return -1
            if self.endnum < 0 or self.endnum > self.startnum:
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
                do_checkboxes()
            #data in rows
            else:
                if self.xinfirstcol.get() == 1:   #use first row as x-column                   
                    self.colnames = ['x','y']
                else:
                    #x and y cols are in every row, so use manual entry of no cols
                    print self.numcolsinrow.get()
                    self.colnames=range(self.numcolsinrow.get())
                print 'self.colnames', self.colnames
                do_checkboxes()

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
                    if self.ignorecomments.get() == 1 and line.startswith('#'):
                        continue
                    else:
                        line = line.replace('#','')
                    textline = ''
                    if s == ' ':
                        data=string.split(line)
                    else:
                        data=string.split(line, s)
                    for c in range(len(data)):
                        textline = textline + '%-7s   ' %data[c]
                    self.previewtext.insert(END, textline+ '\n')
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
                    data=string.split(line,s)
                    name=data[0]
                    data.remove(name)
                    if len(data) >0:
                        datasets[name] = data                
                self.previewtext.delete("1.0",END)
                
                # Sort the datasets by name                
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
            return 1

        f1=Frame(self.ask_csv, borderwidth=2, relief=GROOVE)
        Radiobutton(f1,text='Data in columns',variable=self.column_row,
                                    command=do_preview,value=0).pack()
        Radiobutton(f1,text='Data in rows',variable=self.column_row,
                                    command=do_preview,value=1).pack()
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
                    variable=self.xinfirstcol, command=do_preview).pack()

        genoptsframe = Frame(self.ask_csv,relief=GROOVE)
        genoptsframe.grid(column=1,row=0,rowspan=2,padx=2,pady=2)

        sepEntry = Pmw.EntryField(genoptsframe, labelpos = 'w',
                                        label_text = 'Separator:',
                                        entry_textvariable = self.separator,
                                        validate = None, command = do_preview)
        sepEntry.configure(entry_background = 'white')
        sepEntry.pack(fill=X,padx=2,pady=2)
        sepEntry.component('entry').bind('<KeyPress>', do_preview)
        sepEntry.component('entry').bind('<KeyRelease>', do_preview)
        Checkbutton(genoptsframe,text='Ignore comments',variable=self.ignorecomments,
                    command=do_preview).pack(fill=X,padx=2,pady=2)

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
                                   entryfield_validate = do_preview,
                                   entryfield_command = do_preview)
        self.linenum.pack(fill=X,padx=2,pady=2)
        self.endnum = Pmw.Counter(genoptsframe,
                                   increment=1,
                                   labelpos = 'w',
                                   label_text = 'End at line',
                                   entryfield_value = "0",
                                   datatype ='numeric',
                                   entryfield_validate = do_preview,
                                   entryfield_command = do_preview)
        self.endnum.pack(fill=X,padx=2,pady=2)        
        impmultbutton = Checkbutton(genoptsframe, text='Import into multiple datatabs',variable=self.multiple_tabs)
        impmultbutton.pack(fill=X,padx=2,pady=2)
        self.balloon.bind(impmultbutton, 'Import each column into a separate dataset \nUses the fist row as dataset names')

        def cancel_import():
            self.do_import=None
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
        do_preview()
        return

    def do_import_data(self, filename=None):
        """Returns a list of new ekin datasets this functions requires the use of
           the vars from the dialog, so is not independently callable"""
        import string

        if not self.do_import:
            return
        s=self.separator.get()

        # Get a root name for the datatabs if data in rows

        if self.column_row.get() == 1:
            import tkSimpleDialog
            showname=tkSimpleDialog.askstring('Tab name',
                                       'Root name for datasets?',
                                       parent=self.ekin_win)
        datasets = {}

        # Data in rows?
        if self.column_row.get()==1:
            #All the work is done in dopreview,
            #here we just really need to create the ekin datasets
            x = [i for i in self.inrows_xdata]            
            for name in self.inrows_datasets.keys():
                y = [i for i in self.inrows_datasets[name]]
		datasets[name]=EkinDataset(xy=[x,y])
        else:
           #Data in Cols, do all files in file list
           for filename in self.file_list:               
               showname=os.path.split(filename)[1]
               showname=showname.replace('_','-')               
               showname=os.path.splitext(showname)[0]
               lines = self.get_lines(filename)
               importdata = []
               x=[]
               y=[]
               self.numrows = len(lines)
               for row in range(self.startnum, self.numrows):
                    line = string.strip(lines[row])
                    if lines[row].startswith('#'):
                        continue
                    if s == ' ':
                        data=string.split(line)
                    else:
                        data=string.split(line, s)

                    #now get chosen cols using checkboxes, may be multiple datasets created
                    temp=[]
                    for c in range(len(data)):
                        if c < len(self.colschosen) and self.colschosen[c].get() == 1:
                            temp.append(data[c])
                    importdata.append(temp)
               
               if self.multiple_tabs.get() == 1:
                   #we want to put all the cols into new datasets                   
                   for i in range(1,len(importdata)):
                       try:
                           x.append(float(importdata[i][0]))
                       except:
                           x.append(None)
                   cols = importdata[0]
                   for j in range(1,len(cols)): #iterate over each col
                       y=[]
                       for i in range(1,len(importdata)):
                           try:
                               y.append(float(importdata[i][j]))
                           except:
                               y.append(None)
                       name = cols[j]                      
		       datasets[name] = EkinDataset(xy=[x,y])

               else:
                   #otherwise just take the 1st 2 and ignore the rest
                   for i in range(len(importdata)):
                       try:
                           xval=float(importdata[i][0])
                           yval=float(importdata[i][1])
                           x.append(xval)
                           y.append(yval)
                       except:
                           pass              
		   datasets[showname] = EkinDataset(xy=[x,y])
        #print 'datasets', datasets
        return datasets

    def import_multiple(self,event=None, default_extension='.txt'):
        """Import a single or multiple files"""

        # a list with the current sequences loaded for display
        self.file_list=[]
        # populate file list here
        self.file_list=self.open_filelist(default_extension=default_extension)
        #if self.filename == None:
        if len(self.file_list) > 0:
            self.filename = self.file_list[0]
        else:
            return

        #get all lines from file and put in variable
        self.lines = self.get_lines(self.filename)
        #get options for import and then do import
        status = self.import_csv_dialog()
        self.ekin_win.wait_window(self.ask_csv)
        if status == 0:
            return None
        datasets = self.do_import_data()
        return datasets


    def import_ccpnmr(self, event=None):
        """Import ccpnmr data"""
        #load from file
        import tkFileDialog, os
        importfile=tkFileDialog.askopenfilename(initialdir=self.savedir,title='Select file',
                                                filetypes=[("txt files","*.txt"),
                                                           ("All files","*.*")],
                                                parent=self.ekin_win)
        if not importfile:
            return

        fd=open(importfile)
        lines = fd.readlines()
        fd.close()
        length=len(lines)
        if length>1:
            chmshftcol = tkSimpleDialog.askstring("Which data?",
                                               "Found "+ str(length)+" lines of data.\nWhich chem shift do you want\n H or N?" ,
                                               parent=self.ekin_win)
        if chmshftcol == None:
            return
        datasets={}
        import re
        s = re.compile('\s')
        count=0
        ph=[]
        for line in lines:
            count=count+1

            #split into list and remove all whitespace elements
            line=line.replace('\n','')
            row=s.split(line)
            for r in row[:]:
                if r =='' or r =='!':
                    row.remove(r)
            #print row
            #read ph values, first row
            if count==1:
                for r in row:
                    ph.append(r)
                print 'phLIST', ph
                continue
            #get residue name and number
            name=row[2]
            num=row[1]
            resname=str(num)+name
            print 'name=',name
            #get x and y into lists - pH and Chem shft
            itr=range(5,len(row))

            N=[]
            H=[]
            for i in itr:
                if i%2==0:
                    H.append(row[i])
                else:
                    N.append(row[i])

            print 'N',N
            print 'H',H
            #convert to ekin format and create the data in ekin
            if resname in self.parent.datasets:
                continue
            if len(ph) != len(H):
                continue
            dp=0
            x=[]
            y=[]
            for i in range(len(ph)):

                x.append(ph[i])
                if chmshftcol == 'N':
                    y.append(N[i])
                elif chmshftcol == 'H':
                    y.append(H[i])
                dp=dp+1

            #datasets[resname] = EkinConvert.xy2ekin([x,y])
	    datasets[resname]=EkinDataset(xy=[x,y])

        return datasets


    def open_sparky_dialog(self,event=None):
        """Dialog to select sparky dir"""

        #open dir first
        import tkFileDialog, os
        sparkydir=tkFileDialog.askdirectory(initialdir=self.savedir,title='Select directory containing Sparky peak files'
                                            ,mustexist=1,parent=self.ekin_win)
        if not sparkydir:
            return
        self.sparky_options=Toplevel()
        self.sparky_options.title('Sparky import options')
        self.sparky_options.geometry('+0+200')
        self.assigncol_var = StringVar()
        self.assigncol_var.set('N')
        self.pattern_var = StringVar()
        self.fval = IntVar()
        self.fval.set(5)

        def read_all():
            if self.pattern_var.get() is None:
                return
            else:
                pattern = self.pattern_var.get()

            if self.assigncol_var.get is None:
                self.assigncol = 'N'
            import sparky
            self.peaks=sparky.read_all(self.ekin_win,sparkydir,pattern)
            #show import window if all ok
            self.import_sparky()
        def cancel():
            self.sparky_options.destroy()
            return 0

        self.assigncol_menu = Pmw.OptionMenu(self.sparky_options,
                labelpos = 'w',
                label_text = 'Assigned col:',
                menubutton_textvariable = self.assigncol_var,
                items = ['H', 'N', 'both'],
                menubutton_width = 10)

        self.assigncol_menu.grid(row=0,column=0,columnspan=2,sticky='nes',padx=2,pady=4)
        Label(self.sparky_options,text='Pattern:').grid(row=1,column=0,sticky='news',padx=3,pady=4)
        Entry(self.sparky_options,textvariable=self.pattern_var).grid(row=1,column=1,sticky='news',padx=3,pady=4)
        Label(self.sparky_options,text='F value:').grid(row=2,column=0,sticky='news',padx=3,pady=4)
        Entry(self.sparky_options,textvariable=self.fval).grid(row=2,column=1,sticky='news',padx=3,pady=4)
        Button(self.sparky_options,text='Read Files',command=read_all).grid(row=4,column=0,sticky='news',padx=3,pady=4)
        Button(self.sparky_options,text='Close',command=cancel).grid(row=4,column=1,sticky='news',padx=3,pady=4)
        return


    def import_sparky(self,event=None):
        """Dialog showing selected directory containing Sparky peak files"""

        #
        # Show what we got
        #
        self.overwrite_all = IntVar()
        self.overwrite_all.set(0)
        self.removecurrent = IntVar()
        self.removecurrent.set(0)
        self.assigncol = self.assigncol_var.get()

        self.sparky_win=Toplevel()
        self.sparky_win.title('Titration curves read')
        self.sparky_win.geometry('+0+200')
        Label(self.sparky_win,text='I read titration curves for the following peaks').grid(row=0,column=0,columnspan=3)
        Label(self.sparky_win,text='Select one to analyse').grid(row=1,column=0,columnspan=3)
        #
        # List all of the groups
        #
        pHs=self.peaks.keys()
        pHs.sort()
        #
        # Loop over all pH values to get all the residues
        #
        residues=[]
        for pH in pHs:
            ress=self.peaks[pH].keys()
            for res in ress:
                if not res in residues:
                    residues.append(res)
        residues.sort()
        self.peak_residues=residues
        count=0
        self.peak=IntVar()
        self.peak.set(0)
        yscrollbar=Scrollbar(self.sparky_win,orient='vertical',width=10)
        yscrollbar.grid(row=2,column=3,sticky='nws')
        #
        self.peaks_read=Listbox(self.sparky_win,bg='white',
                                fg='black',
                                height=30,
                                width=20,
                                yscrollcommand=yscrollbar.set,
                                selectmode=EXTENDED)
        self.peaks_read.grid(row=2,column=0,columnspan=3,sticky='news',padx=3,pady=4)
        yscrollbar.config(command=self.peaks_read.yview)
        self.sparky_win.grid_rowconfigure(2, weight=1)
        #
        for residue in residues:
            self.peaks_read.insert(END,residue)

        overwritecheck = Checkbutton(self.sparky_win, text="Overwrite All", variable=self.overwrite_all)
        overwritecheck.grid(row=3,column=0,columnspan=1,padx=1,pady=1)
        removecurrentcheck = Checkbutton(self.sparky_win, text="Remove Current", variable=self.removecurrent)
        removecurrentcheck.grid(row=3,column=1,columnspan=2,padx=1,pady=1)

        #
        # Buttons for selecting
        #
        Button(self.sparky_win,text='Do Import',command=self.copy_peak).grid(row=4,column=0,sticky='news',padx=3,pady=4)
        #Button(self.sparky_win,text='Copy to pKaTool',command=self.copy_peaks_2pKaTool).grid(row=4,column=1,sticky='news',padx=3,pady=4)
        Button(self.sparky_win,text='Done',command=self.sparky_win.destroy).grid(row=4,column=2,sticky='news',padx=3,pady=4)

        #self.peaks_read.bind("<ButtonRelease-1>",self.copy_peak)
        return

    #
    # Changed this so that it asks first if u want to overwrite the current datasets already there.
    #

    def copy_peak(self,event=None):
        """Copy the selected peak to the data"""
        selection=self.peaks_read.curselection()
        residues=[]
        data={}
        f = self.fval.get()
        for sel_number in selection:
            residue=self.peak_residues[int(sel_number)]
            #
            # Construct the titration curve
            #
            data[residue]={}
            pHs=self.peaks.keys()
            pHs.sort()
            for pH in pHs:
                if self.peaks[pH].has_key(residue):
                    data[residue][pH]=self.peaks[pH][residue]

            #print data[residue]
            # Get combined shift if user chose 'both'
            if self.assigncol == 'both':
                import sparky
                Hlist=[]
                Nlist=[]
                for pH in data[residue].keys(): #create lists for H and N
                    Hlist.append(data[residue][pH]['H'])
                    Nlist.append(data[residue][pH]['N'])
                href=min(Hlist)
                nref=min(Nlist)
                for pH in data[residue].keys():
                    h=data[residue][pH]['H']
                    n=data[residue][pH]['N']
                    data[residue][pH]['both']=sparky.calc_delta_shift(h, n, href, nref ,f)

            #print data[residue]
        #
        # Delete all the other tabs first if needed
        #
        if self.removecurrent.get()==1:
            self.parent.delete_all()

        #
        # Insert the data if it's not already there, but ask to overwrite if it is!
        #
        datasets = {}
        for peak in data.keys():
            num_vars=2
            x=[]
            y=[]
            count=0
            for pH in data[peak].keys():
                x.append(pH)
                y.append(data[peak][pH][self.assigncol])
                count=count+1
            
            #datasets[peak] = EkinConvert.xy2ekin([x,y])
	    datasets[peak]=EkinDataset(xy=[x,y])

        #add to ekin from here - easier for now...
        self.parent.mode_var.set(4)
        #self.parent.update_mode()
        #self.parent.insert_multiple_datasets(datasets, overwrite=self.overwrite_all.get())
        for d in datasets:
            self.parent.insertDataset(datasets[d],d)

        return

    def import_chem_shift(self,event=None):
        """Import kinetic data from a nmr file"""

        import tkFileDialog
        filename=tkFileDialog.askopenfilename(defaultextension='.txt',initialdir=self.savedir,
                                              filetypes=[("NMR file","*.txt"),
                                                         ("All files","*.*")],
                                              parent=self.ekin_win)
        if os.path.isfile(filename):
            fd=open(filename)

            line = fd.readline()
            line = line.strip('\t\n')
            residues = line.split()
            values = []
            while len(line)>2:
                line = fd.readline()
                line = line.strip('\t\n')
                lvalues = line.split()
                values.append(lvalues)
            fd.close()
            data = {}

            for r in range(1,len(residues)):
                data[residues[r]]={}

            for i in range(1,len(values)):
                for j in range(1,len(values[i])):
                    #print values[i][0]
                    pH = values[i][0]
                    res = residues[j]
                    print 'data[%s][%s] = %s'%(res,pH,values[i][j])
                    data[residues[j]][pH]=values[i][j]

            #
            # Insert the data if it's not already there
            #
            tabnames=self.nb.pagenames()
            for residue in data.keys():
                dps=len(data[residue].keys())
                self.reprint_parameters()
                if not residue in tabnames:
                    newtab_num=self.add_datalist(residue,data_points=dps)
                    count=0
                    for pH in data[residue].keys():
                        self.data[newtab_num][0][count]['var'].set(pH)
                        self.data[newtab_num][1][count]['var'].set(data[residue][pH])
                        count=count+1
            #self.redraw_graph(fitting=1)
        return


    def open_filelist(self,default_extension='.txt'):
        """Open multiple filenames list dialog"""
        import tkFileDialog, os
        filelist=tkFileDialog.askopenfilenames(defaultextension=default_extension,
                                                initialdir=self.path,
                                                filetypes=[("csv files","*.csv"),
                                                           ("csvx files","*.csvx"),
                                                            ("txt files","*.txt"),
                                                           ("All files","*.*")],
                                                parent=self.ekin_win)
        return filelist



    def import_CD_tempscan(self):
        """Import a temperature scan from a Jasco CD spec"""
        import tkFileDialog, os
        filename=tkFileDialog.askopenfilename(defaultextension='.txt',initialdir=os.getcwd(),
                                              filetypes=[("Jasco txt","*.txt"),
                                                         ("All files","*.*")],
                                              parent=self.ekin_win)
        if not filename:
            return
        #
        # If we got a filename then read the stuff
        #
        import os
        if not os.path.isfile(filename):
            import tkMessageBox
            tkMessageBox.showwarning('File not found',
                                     'I could not find %s' %filename,
                                     parent=self.ekin_win)
            return
        #
        # Open and read file
        #
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        #
        # Parse file
        #
        # Format is:
        # 00 TITLE	0.6mg/ml hewl sigma
        # 01 DATA TYPE
        # 02 ORIGIN	JASCO
        # 03 OWNER
        # 04 DATE	06/10/17
        # 05 TIME	17:59:52
        # 06 SPECTROMETER/DATA SYSTEM	J-810
        # 07 DELTAT	0
        # 08 TUNITS	Temperature [C]
        # 09 FIRSTT	25
        # 10 LASTT	35
        # 11 NPOINTST	3
        # 12 DELTAX	-0.5
        # 13 XUNITS	Wavelength [nm]
        # 14 YUNITS	CD[mdeg]
        # 15 YUNITS2	HT[V]
        # 16 FIRSTX	350
        # 17 LASTX	250
        # 18 NPOINTS	201
        # 19 XYDATA1	25	30	35
	#
        # For a file with scans of three temperatures
        #
        data={}
        data['TITLE']=lines[0]
        data['ORIGIN']=lines[1]
        data['DATE']=lines[4]+lines[5]
        data['SPEC']=lines[6]
        #
        # Get the range of temps
        #
        t_temps=lines[19].split()[1:]
        temps=[]
        for t in t_temps:
            temps.append(float(t))
        #
        # Get the number of wavelengths
        #
        lambda_points=int(lines[18].split()[-1])
        #
        # Read the data
        #
        raw_data={}
        for x in range(lambda_points):
            line=lines[x+20]
            line=line.split()
            wavlen=float(line[0])
            count=0
            for temp in temps:
                count=count+1
                mdeg=float(line[count])
                if not raw_data.has_key(temp):
                    raw_data[temp]=[]
                raw_data[temp].append([wavlen,mdeg])
        #
        # Insert the tabs
        #
        datasets = {}
        temp1=temps[0]

        for temp in temps:
            name = 'CD(T'+str(temp)+')'
            x=[]
            y=[]
            count=0
            for wavlen,mdeg in raw_data[temp]:
                x.append(wavlen)
                y.append(mdeg)
                count=count+1

            #datasets[name] = EkinConvert.xy2ekin([x,y])
	    datasets[resname]=EkinDataset(xy=[x,y])

        return datasets

    def handleUTF(self, val):
        """handle special chars"""
        val = unicode(val, "utf-8",errors='replace')
        val = val.encode('ascii','ignore')
        try:
            val = float(val)
        except:
            val = None
        return val

class Exporter:

    def __init__(self, parent):
        self.parent = parent
        self.ekin_win = parent.ekin_win
        self.filename = None
        return

    def export_csv(self, data, fitdata=None):
        """Export the current tab to a csv file"""

        import tkFileDialog
        filename=tkFileDialog.asksaveasfilename(defaultextension='.csv',initialdir=self.savedir,
                                              filetypes=[("CSV file","*.csv"),
                                                         ("All files","*.*")],
                                              parent=self.ekin_win)
        if not filename:
            return
        labels=data[0]['label'],data[1]['label']
        #x,y,a = EkinConvert.ekin2xy(data)
	x,y = data.getxy()
        #print filename
        try:
            import csv
            csvwriter = csv.writer(open(filename, "wb"))
            #csvwriter = csv.writer(sys.stdout)
        except:
            tkMessageBox.showwarning("Save file",
                            "Cannot open this file for writing\n(%s)" % filename)
            return
        csvwriter.writerow(labels)
        for i in range(len(x)):
            csvwriter.writerow([x[i],y[i]])

        # Ask if we should export the fit too
        export_fit=tkMessageBox.askyesno("Export fit",
                                "Export fitted values?",parent=self.ekin_win)
        if export_fit and fitdata != None:
            f=[]
            for i in fitdata.keys():
                f.append(fitdata[i])
            csvwriter.writerow(f)
        return

