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

try:
    from Plugins import Plugin
except:
    from PEATDB.Plugins import Plugin
import os, types, copy, pickle
from Tkinter import *
import tkFileDialog
import Pmw
import PEATSA.WebApp.Data
import PEATSA.WebApp.UtilityFunctions
import PEATSA.Core as Core
from PEATDB.Dialogs import MultipleValDialog
from PEATDB.Actions import DBActions
from PEATDB.TableModels import TableModel
from PEATDB.Tables import TableCanvas
import tkMessageBox, tkSimpleDialog

class SummaryPlugin(Plugin):
    """Plugin for PEAT that enables several projects to be sumamrised together"""
    capabilities = ['gui']
    requires = ['PEATSA']
    menuentry = 'Summary Plugin'
    gui_methods = {'help':'Help',
                   'quit':'Close Window'}
    buttonorder = ['createJobDialog','fetchJob','editConfigFile','help','quit']
    about = 'This plugin allows several projects to be sumamrised together'    
    
    def main(self, parent=None, DB=None):       
        if parent==None:
            if DB!=None:
                self.DB = DB
            else:
                return 
        else:
            self.parent = parent
            self.DB = parent.DB
            self.parentframe = None            
            self._doFrame()
        return

    def _doFrame(self):
        self.mainwin = self.parent.createChildFrame(width=460,title='PEATSA Plugin')        
        methods = self._getmethods()
        methods = [m for m in methods if m[0] in self.gui_methods.keys()]        
        l=Label(self.mainwin, text='PEATSA Interface')
        l.pack(side=TOP,fill=BOTH)

        self.log = self.createLogWin(self.mainwin)
        self.log.pack(side=TOP,fill=BOTH,expand=1)         
        self.stdout2Log()
        return
    

    
