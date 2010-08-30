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

from Plugins import Plugin
from Tkinter import *
from GUI_helper import *
from PEATDB.Yasara import YasaraControl

class YasaraPlugin(Plugin, GUI_help):
    """Yasara plugin for PEAT App"""
    capabilities = ['gui','uses_sidepane']
    requires = ['yasaramodule']
    menuentry = 'Yasara Control'
    gui_methods = []
    about = 'Handles an open yasara instance'    

    def main(self, parent):      
        if parent==None:
            return
        self.parent = parent
        self.DB = parent.DB
        self.parentframe = None
        
        try:
            self._doImport()
        except:
            import tkMessageBox
            tkMessageBox.showwarning('Failed',
                                     'No yasara specified in settings!')
            return
        self._doFrame() 
        return

    def _doImport(self):
        """Same as DBActions"""
        import sys,os
        self.preferences=self.parent.preferences       
        yasaradir = self.preferences.get('molgraphAppPath')
        if not os.path.isdir(yasaradir):
            yasaradir=os.path.split(yasaradir)[0]
        dirname=os.path.split(yasaradir)[1]
        if dirname.lower()=='yasara.app':
            yasaradir=os.path.join(yasaradir,'yasara')
       
        sys.path.append(os.path.join(yasaradir,'pym'))
        sys.path.append(os.path.join(yasaradir,'plg'))
        import yasaramodule as yasara
        self.yasara = yasara
        return
    
    def _doFrame(self):
        if 'uses_sidepane' in self.capabilities:
            self.mainwin = self.parent.createChildFrame()    
        else:
            self.mainwin=Toplevel()
            self.mainwin.title('Yasara control')
            self.mainwin.geometry('800x600+200+100')            

        methods = self._getmethods()
        methods = [m for m in methods if m[0] in self.gui_methods]       
        l=Label(self.mainwin, text='Yasara plugin')
        l.pack(side=TOP,fill=BOTH)
        yc = YasaraControl(self.mainwin, self.yasara)
        yc.pack()
        return
      
    def quit(self):
        self.mainwin.destroy()
        return
    
