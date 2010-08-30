#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#
# Written by Damien Farrell, Feb 2010

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
    
