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
import tkSimpleDialog, tkFileDialog, tkMessageBox
import Pmw

class TopLevelModalDialog(Toplevel):
    def __init__(self, parent, width=300, height=100):
        Toplevel.__init__(self, parent)
        self.transient(parent)
        self.title('')
        self.parent = parent
        x=(parent.winfo_rootx()+parent.winfo_width()/2)-width/2
        y=(parent.winfo_rooty()+parent.winfo_height()/2)-height/2        
        self.geometry('%sx%s+%s+%s' %(width,height,x,y)) 
        self.body = Frame(self)
        self.initial_focus = self.body
        self.body.pack(fill=BOTH,expand=1,padx=5, pady=5)        
        self.grab_set()
        if not self.initial_focus:
            self.initial_focus = self
        self.protocol("WM_DELETE_WINDOW", self.close)        
        return

    def close(self):
        self.parent.focus_set()
        self.destroy()      
        return

class ProgressDialog(TopLevelModalDialog):
    def __init__(self, parent, message='Working'):
        TopLevelModalDialog.__init__(self, parent)
        self.title(message)
        progrlbl = Label(self.body,text='Progress:')
        progrlbl.pack(fill=BOTH,padx=2,pady=4)
        import ProgressBar
        self.bar = ProgressBar.ProgressBar(self.body)
        self.bar.frame.pack(fill=Y,padx=2,pady=4)
        return


