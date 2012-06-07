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

"""Classes for primer handling"""

import os, pickle
from Tkinter import *
import tkFileDialog
import Dialogs 
import Pmw
import fileinput

class PrimerDatabase(object):
    def __init__(self, data=None):
        if data==None:
            self.data = {}
        else:
            self.data = data
        return

    def getPrimerDetails(self, name):
        info = self.data[name]        
        text = name+'\n'
        for i in info:
            text+=i+': '+str(info[i])+'\n'
        return text

    def load(self, filename):
        """load a db"""
        fd=open(filename,'rb')
        self.data = pickle.load(fd)
        fd.close()       
        self.filename = filename
        return

    def save(self):
        fd=open(self.filename,'w')
        pickle.dump(fd, self.data)        
        return

    def add(self, primer):
        return

    def __repr__(self):
        return str(self.data.keys())

class PrimerDBGUI(Frame):
    """Class to provide a frame for text editing"""

    def __init__(self, parent, database=None):
        self.parent=parent      
        Frame.__init__(self, parent)
        self.main=self
        if database == None:
            self.database = PrimerDatabase()
        else:
            self.database = database
        self.makeGUI()
        return
        
    def makeGUI(self):
        self.menuBar()
        self.m = PanedWindow(self.main,
                           orient=HORIZONTAL,
                           sashwidth=2,
                           showhandle=True)
        self.m.grid(row=1,column=0,sticky='NEWS',pady=2)
        # Create listbox for the results   
        self.primerlist = Pmw.ScrolledListBox(self.m,
                        items=(),
                        labelpos='nw',                  
                        listbox_height = 10,
                        selectioncommand=self.selectPrimer)     
        '''self.primerlist = Listbox(self.m,bg='white',
                             fg='black',
                             height=10,width=20,
                             selectmode=EXTENDED)'''
        self.m.add(self.primerlist)      
        #Detailed results window
        self.details = Pmw.ScrolledText(self.m,
                labelpos = 'n',
                label_text='Details',
                usehullsize = 1,         
                hull_height = 250,
                text_wrap='word')
        self.m.add(self.details)
        self.m.paneconfigure(self.primerlist,sticky='news')
        fr=Frame(self.main)
        fr.grid(row=2,column=0,sticky='NEWS',pady=2)
        self.doButtons(fr)
        self.main.columnconfigure(0,weight=1)
        self.main.rowconfigure(1,weight=1)
        return
        
    def menuBar(self):
        menuBar = Pmw.MenuBar(self.main,                            
                            hull_borderwidth = 1,
                            hull_width=300)
        menuBar.grid(row=0,column=0,sticky='NEWS')

        menuBar.addmenu('File', 'primer database IO')
        menuBar.addmenuitem('File', 'command', 'Load primer database',
                    command = self.loadPrimerDB,
                    label = 'Load primer DB')               
        menuBar.addmenuitem('File', 'command', 'Save primer database',
                command = self.savePrimerDB,
                label = 'Save primer DB')
        menuBar.addmenuitem('File', 'command', 'Add primers from file',
                command = self.addPrimers,
                label = 'Add primers')
        menuBar.addmenuitem('File', 'command', 'Save to csv file',
                command = self.saveAsCSV,
                label = 'Save to csv file')
        menuBar.addmenuitem('File', 'command', 'Add primers from csv file',
                command = self.addPrimersfromCSV,
                label = 'Add primers from csv file')
        return

    def doButtons(self, frame):
        """Buttons"""
        buttons=[['Add',self.newPrimer],
                 ['Edit',self.newPrimer],
                 ['Rename',self.newPrimer],
                 ['Apply',self.applyPrimer],
                 ['Find',self.doSearch]]     

        for button,command in buttons:
            x=Button(frame,text=button,command=command)
            x.pack(side=LEFT,fill=BOTH,padx=2,pady=2)
        return

    def newPrimer(self):
        return

    def applyPrimer(self):
        return

    def doSearch(self):
        return
    
    def selectPrimer(self):
        name = self.primerlist.getcurselection()[0]
        data = self.database.data[name]
        text = self.database.getPrimerDetails(name)
        self.details.settext(text)
        return

    def update(self):
        """Update with current primer db"""
        primers = self.database.data
        self.primerlist.setlist(primers.keys())
        self.details.clear()
        return

    def loadPrimerDB(self):
        """Load a new primer database"""
        filename = Dialogs.openFilename(parent=self,ext='primerDB')
        if not filename:
            return
        self.database.load(filename)
        #print self.database
        self.update()
        return

    def savePrimerDB(self):
        return
        
    def addPrimers(self):
        return

    def saveAsCSV(self):
        return
        
    def addPrimersfromCSV(self):
        return   
        
class PrimerDesignGUI(Frame):

    def __init__(self, parent, parentapp):
        self.parent=parent   
        self.parentapp = parentapp   
        Frame.__init__(self, parent)
        self.main=self
        self.makeGUI()
        return
        
    def makeGUI(self):
        return

