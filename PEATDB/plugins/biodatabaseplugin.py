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
import Pmw
from GUI_helper import *
from Bio import Entrez

class BioDBQuery(Plugin, GUI_help):
    """Bio database query plugin for PEAT"""
    capabilities = ['gui','uses_sidepane']
    requires = ['biopython']
    menuentry = 'Bio DB Query'
    gui_methods = {'quit':'Quit'}
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
        self.ID='Bio DB Query Plugin'
        if 'uses_sidepane' in self.capabilities:
            self.mainwin = self.parent.createChildFrame()
        else:
            self.mainwin=Toplevel()
            self.mainwin.title(self.ID)
            self.mainwin.geometry('800x600+200+100')

        methods = self._getmethods()
        methods = [m for m in methods if m[0] in self.gui_methods]
        self._createButtons(methods)
        self.queryvar=StringVar()
        Entry(self.mainwin,textvariable=self.queryvar, width=20,bg='white').pack(side=TOP,fill=BOTH,pady=4)
        b=Button(self.mainwin,text='search',command=self.doSearch)
        b.pack(side=TOP,fill=BOTH,pady=2)
        self.dbvar=StringVar()
        self.dbvar.set('pubmed')
        opts = Pmw.OptionMenu(self.mainwin,
                labelpos = 'w',
                label_text = 'Database:',
                menubutton_textvariable = self.dbvar,
                items = ['pubmed', 'ncbisearch', 'protein', 'nucleotide',
                         'structure', 'genome', 'books',
                          'gene', 'genomeprj'],
                menubutton_width = 10)
        opts.pack(side=TOP,fill=BOTH,pady=2)
        self.results = Pmw.ScrolledText(self.mainwin,
                labelpos = 'n',
                label_text='Search results',
                usehullsize = 1,
                hull_width = 400,
                hull_height = 500,
                text_wrap='word')
        self.results.pack(side=TOP,fill=BOTH,padx=4,pady=4)
        return

    def _createButtons(self, methods):
        """Dynamically create buttons for supplied methods, which is a tuple
            of (method name, label)"""     
        for m in methods:           
            b=Button(self.mainwin,text=self.gui_methods[m[0]],command=m[1])
            b.pack(side=BOTTOM,fill=BOTH)
        return
    
    def doSearch(self, query=None):
        import types
        fields = ['Title','AuthorList','DOI']
        structfields = ['PdbDescr', 'ExpMethod', 'EC',
                        'OrganismList', 'ProteinChainCount', 'PdbDepositDate',
                        'Resolution', 'Id']

        if query == None:
            query = self.queryvar.get()
        db = self.dbvar.get()
        if db == 'structure': fields = structfields
        self.results.delete(1.0,END)
        Entrez.email = "peat_support@ucd.ie"
        handle = Entrez.esearch(db=db,term=query)
        record = Entrez.read(handle)
        idlist = record["IdList"]

        for i in idlist:
            handle = Entrez.esummary(db=db, id=i)
            rec = Entrez.read(handle)
            for f in fields:
                try:
                    if type(rec[0][f]) is Entrez.Parser.ListElement:
                        for a in rec[0][f]:
                            self.results.insert(END,a+', ')
                    else:
                        self.results.insert(END,rec[0][f]+'\n')
                except:
                    pass
            self.results.insert(END,'\n'+'------------------------'+'\n')

        return

    def quit(self):
        self.mainwin.destroy()
        return

