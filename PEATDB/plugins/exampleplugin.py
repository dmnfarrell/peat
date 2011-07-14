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

class ExamplePlugin(Plugin, GUI_help):
    """Template GUI plugin for PEAT App"""
    capabilities = ['gui','uses_sidepane']
    requires = ['']
    menuentry = 'Test Plugin'
    gui_methods = {}
    about = 'This plugin is a template'    

    def main(self, parent):      
        if parent==None:
            return
        self.parent = parent
        self.DB = parent.DB
        self.parentframe = None
        self._doFrame()
        return   

    def _doFrame(self):
        if 'uses_sidepane' in self.capabilities:
            self.mainwin = self.parent.createChildFrame()    
        else:
            self.mainwin=Toplevel()
            self.mainwin.title('A PEAT Plugin')
            self.mainwin.geometry('800x600+200+100')            
             
        self.ID='Base GUI Plugin'
        #self._createMenuBar()
        methods = self._getmethods()
        methods = [m for m in methods if m[0] in self.gui_methods]
        #self._createButtons(methods)
        l=Label(self.mainwin, text='This is a template plugin')
        l.pack(side=TOP,fill=BOTH)
        self.mainwin.bind("<Destroy>", self.quit)
        return
  
    def _createButtons(self, methods):
        """Dynamically create buttons for supplied methods"""
        for m in methods:           
            b=Button(self.mainwin,text=m[0],command=m[1])
            b.pack( side=TOP,fill=BOTH)            
        return
    
    def _createMenuBar(self):
        """Create the menu bar for the application. """
        self.menu=Menu(self.mainwin)
        self.file_menu={ '01Quit':{'cmd':self.quit}}
        self.file_menu=self.create_pulldown(self.menu,self.file_menu)
        self.menu.add_cascade(label='File',menu=self.file_menu['var'])
        
        self.mainwin.config(menu=self.menu)
        return
   
    def loadDB(self):
        """Load a DB"""
        from PEATDB.Base import PDatabase
        if local != None:
            self.DB = PDatabase(local=local)
        return     
        
    def quit(self, evt=None):
        """Override this to handle pane closing"""
        return
    
    
def main():
    """Run some analysis"""
    from optparse import OptionParser
    parser = OptionParser()
    app = ExamplePlugin()
    parser.add_option("-f", "--file", dest="file",
                        help="Open a local db")
    if opts.file != None and os.path.exists(opts.file):
        app.loadDB(opts.file)
        

if __name__ == '__main__':
    main()       
