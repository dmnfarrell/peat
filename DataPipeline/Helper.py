#!/usr/bin/env python
#
# DataPipeline - A data import and fitting tool
# Copyright (C) 2011 Damien Farrell
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
# Email: damien.farrell_at_ucd.ie
# Normal mail:
# Damien Farrell
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
# 

from Tkinter import *
import Images
import Pmw
from Base import Pipeline
from Editor import TextEditor

"""Implements a simple wizard in the application to allows users to
select and appropriate Importer and some other base settings"""

class HelperDialog(Frame):
    
    def __init__(self, parent=None):
        self.parent=parent         
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
        else:
            self.main=Toplevel()
            self.master=self.main
        self.main.title('Configuration Helper')
        self.main.geometry('680x500+300+200')
        self.formatPanel()
        self.filename = 'mysettings.conf'
        self.p = Pipeline()
        self.p.createConfig(self.filename)        
        
        if self.parent:
            self.parent.wait_window(self.main)
        return
        
    def formatPanel(self):
        """Showing format options"""
        panel1 = Frame(self.main)
        panel1.pack(fill=BOTH,side=BOTTOM)
        self.img1 = Images.formatsDiagram()
        label1 = Label(self.main,image=self.img1)
        label1.pack(fill=BOTH,padx=4,pady=4,ipadx=2,ipady=2)        
        self.formatsbox = Pmw.ScrolledListBox(panel1,
                        items=(range(1,9)),
                        labelpos='n',
                        label_text='Format Option:',
                        listbox_height = 4,                                        
                        usehullsize = 1,
                        hull_width = 100,
                        hull_height = 100 )
        self.formatsbox.grid(row=0,column=0,padx=4,pady=4)
        
        showconfigbtn = Button(panel1, text='Show Config',command=self.showConfig)
        showconfigbtn.grid(row=0,column=3,padx=4,pady=4)        
        writeconfigbtn = Button(panel1, text='Write Config',command=self.writeConfig)
        writeconfigbtn.grid(row=0,column=4,padx=4,pady=4)
        
        return        
      
    def addButton(self, parent, img, row, col):
        button = Button(parent,image=img)
        button.grid(row=row,column=col,padx=4,pady=4,ipadx=2,ipady=2,sticky='ew')
        return button
    
    def writeConfig(self):
        """Write out a config using current settings"""
        
        return
 
    def showConfig(self):
        """Show current config"""        
        t = TextEditor(parent=self,title=self.filename)
        t.openFile(self.filename)        
        return
     
