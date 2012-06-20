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

import os, sys, random
import types
import numpy
from Base import PDatabase
from Record import PEATRecord
from Tkinter import *
import Pmw

class searchHandler(object):
    """Base search class for peat sub objects, we inherit from this"""
    def __init__(self, DB, fields=None):
        self.DB = DB
        if fields == None:
            self.fields = DB.getFields()
        else:
            self.fields = fields
        self.attributes = [] #possible attributes to search on
        self.data = DB.data
        self.reclist = DB.getRecs()

        self.funcs = {'contains':Operators.contains,'=':Operators.equals,
                   '>':Operators.greaterthan,'<':Operators.lessthan,
                   'starts with':Operators.startswith,
                   'ends with':Operators.endswith,
                   'has length':Operators.haslength}
        return

    def searchfunc(self, col, attr, val, op):
        """Must be overriden, this code is an example only
           col: field in DB
           attr: name of attribute to search
           val: value to search for
           op: relational operator to apply e.g. '>'
           returns a list of records that are True for the search
        """
        func = self.funcs[op]
        for rec in self.reclist:
            if not self.data[rec].has_key(col): continue
            if func(val,data[rec][col]):
                names.append(rec)
        return names

    def doSearch(self, filters=None):
        """If we need to search directly from handler"""
        if filters==None:
            F = [self.getFilter]
        else:
            F = filters

        SF = [self.searchfunc for i in F]
        names = Base.Filter(F, SF)
        self.names = names
        return names

    def getWidget(self, parent, callback=None):
        """Return a frame with search widgets specific to this handler"""
        #for now each handler creates one widget by default
        if len(self.fields)==0:
            self.frame = None
        else:
            self.frame = Frame(master=parent)
            self.widget = searchWidget(self.frame, self.fields, self.attributes, callback)
            self.widget.pack()
        return self.frame

    def getFilter(self):
        return self.widget.getFilter()

class simpleSearch(searchHandler):
    """Example for simple string and float searches"""
    def __init__(self, DB, fields=None):
        searchHandler.__init__(self, DB, fields)
        if fields == None:
            self.fields = DB.getTextFields()
        return

    def searchfunc(self, col, attr, val, op):
        """Simple search for strings, attr is always the string itself,
           so it can be passed as empty"""
        floatops = ['=','>','<']
        func = self.funcs[op]
        data = self.data
        names=[]
        for rec in self.reclist:
            if data[rec].has_key(col):
                #try to do numerical comparisons if required
                if op in floatops:
                    try:
                        item = float(data[rec][col])
                        v=float(val)
                        if func(v, item) == True:
                            names.append(rec)
                        continue
                    except:
                        pass
                if col == 'name':
                    item = rec
                else:
                    item = str(data[rec][col])
                if func(val, item):
                    names.append(rec)
        return names

class dictSearch(searchHandler):
    """Standard python dictionary searching"""
    def __init__(self, DB, fields=None):
        searchHandler.__init__(self, DB, fields)
        if fields == None:
            self.fields = DB.getDictFields()
        self.attributes = ['name','link','sequence','label']
        return

    def searchfunc(self, col, attr, val, op):
        """attr is a subkey in the dict that can be recursively found """
        floatops = ['=','>','<']
        func = self.funcs[op]
        data = self.data
        names=[]
        subdata = {}
        print attr
        for rec in self.reclist:
            if data[rec].has_key(col):
                item = self.data[rec][col]
                print type(item)
                #if not type(item) is types.DictType():
                #    continue
                #do recursive search if sub dicts
                if attr in item:
                    if func(val, item[attr]):
                        names.append(rec)

        return names

class ekinSearch(searchHandler):
    """Ekin search handler"""
    def __init__(self, DB, fields=None):
        searchHandler.__init__(self, DB, fields)
        if fields == None:
            self.fields = DB.getEkinFields()
        self.attributes = ['dataset', 'model']
        return

    def searchfunc(self, col, attr, val, op):
        """Search an ekin project"""
        names = []
        datasets = []
        func = self.funcs[op]

        for rec in self.reclist:
            if not self.data[rec].has_key(col): continue
            E = self.data[rec][col]
            print rec, E
            if attr == 'dataset':
                for d in E.datasets:
                    if func(val, d):
                        names.append(rec)
                        datasets.append(d)
            elif attr == 'model':
                for d in E.datasets:
                    F = E.getFitData(d)
                    if F == {}: continue
                    model = F['model']
                    if func(val, model):
                        names.append(rec)
                        datasets.append(d)
            else:
                #we try metadata for search attribute
                for d in E.datasets:
                    M = E.getMetaData(d)
                    if attr in M.keys() and M[attr] != None:
                        if func(val, M[attr]):
                            names.append(rec)
                            datasets.append(d)

        self.datasets = datasets
        return names

    def getWidget(self, parent, callback=None):
       """Add some stuff to widget for ekin handler"""
       searchHandler.getWidget(self, parent, callback)
       return self.frame


class structureSearch(searchHandler):
    """Structure search handler"""
    def __init__(self, DB, fields=None):
        searchHandler.__init__(self, DB, fields)
        if fields == None:
            self.fields = ['Structure']
        self.attributes = ['length']
        return

    def searchfunc(self, col, attr, val, op):
        """Search a structure"""
        names = []
        func = self.funcs[op]
        for rec in self.reclist:
            if not self.data[rec].has_key(col): continue

            if self.data[rec].hasStructure() == 'available':
                struct = self.data[rec][col]
            else:
                continue
            #do something with the structure here..
            print struct, attr
            #if func(val, struct):
            #    names.append(rec)
        return names

class sequenceSearch(searchHandler):
    """Sequence search handler"""
    def __init__(self, DB, fields=None):
        searchHandler.__init__(self, DB, fields)
        if fields == None:
            self.fields = ['aaseq']
        self.attributes = ['length']
        return

class searchWidget(Frame):
    """Class providing filter widgets to retrieve and return search requests.
       This class is called by a search handler to provide specific widgets and
       can be inherited or extended"""
    operators = ['contains','=','>','<','starts with',
                 'ends with','has length']
    booleanops = ['AND','OR','NOT']

    def __init__(self, parent, fields, attributes, callback=None):
        Frame.__init__(self, parent)
        self.parent=parent
        self.filtercol=StringVar()
        if 'name' in fields:
            initial = 'name'
        else:
            initial = fields[0]
        self.attributes = attributes
        self.callback = callback

        row=0
        self.booleanop=StringVar()
        booleanopmenu = Pmw.OptionMenu(self,
                menubutton_textvariable = self.booleanop,
                items = self.booleanops,
                initialitem = 'AND',
                menubutton_width = 6)
        booleanopmenu.grid(row=row,column=0,sticky='news',padx=2,pady=2)
        filtercolmenu = Pmw.OptionMenu(self,
                labelpos = 'w',
                label_text = 'Column:',
                menubutton_textvariable = self.filtercol,
                items = fields,
                initialitem = initial,
                menubutton_width = 10)
        filtercolmenu.grid(row=row,column=1,sticky='news',padx=2,pady=2)
        cbutton=Button(self,text='-', command=self.close)
        cbutton.grid(row=row,column=2,sticky='news',padx=2,pady=2)

        row=1; col=0
        self.attribute = StringVar()
        if len(self.attributes)>0:
            #if we have attributes of the object defined, add entry
            self.attribute = StringVar()
            attrmenu = Pmw.OptionMenu(self,
                    menubutton_textvariable = self.attribute,
                    items = self.attributes,
                    initialitem = self.attributes[0],
                    menubutton_width = 6)
            attrmenu.grid(row=row,column=col,sticky='news',padx=2,pady=2)
            attrentry = Entry(self,textvariable=self.attribute,width=6)
            attrentry.grid(row=row,column=col+1,sticky='news',padx=2,pady=2)
            col=2
        self.operator=StringVar()
        operatormenu = Pmw.OptionMenu(self,
                menubutton_textvariable = self.operator,
                items = self.operators,
                initialitem = 'contains',
                menubutton_width = 8)
        operatormenu.grid(row=row,column=col,sticky='news',padx=2,pady=2)
        self.filtercolvalue=StringVar()
        valsbox=Entry(self,textvariable=self.filtercolvalue,width=20,bg='white')
        valsbox.grid(row=row,column=col+1,sticky='news',padx=2,pady=2)
        valsbox.bind("<Return>", self.callback)
        return

    def close(self):
        """Destroy and remove from parent"""
        self.parent.filters.remove(self)
        self.destroy()
        return

    def getFilter(self):
        """Get filter values for this instance"""
        col = self.filtercol.get()
        val = self.filtercolvalue.get()
        op = self.operator.get()
        attr = self.attribute.get()
        booleanop = self.booleanop.get()
        return col, attr, val, op, booleanop

class searchDialog(Frame):
    """Global dialog to handle searches in PEAT using registered
       search handlers. The multiple filters are then put together"""
    def __init__(self, parent, DB):
        Frame.__init__(self, parent)
        self.DB = DB
        #add handlers we have defined here
        self.handlers = ['simple', 'dict', 'ekin', 'structure']
        self.searchbars = []
        self.currenthandlers = []
        self.doButtons()
        self.doframe()
        return

    def doframe(self):

        rf = LabelFrame(self,text='Results:')
        self.resultsvar = StringVar()
        Label(rf, textvariable=self.resultsvar,fg='blue').pack(fill=BOTH,pady=4)
        rf.pack(fill=BOTH,expand=1,padx=2,pady=4)
        scr = Pmw.ScrolledFrame(self,
                labelpos = 'n',
                label_text = 'Filters:',
                usehullsize = 1,
                hull_width = 600,
                hull_height = 600)
        scr.component('borderframe').configure(relief='groove')
        scr.pack(fill=BOTH,expand=1)
        self.sframe = scr.interior()
        return

    def doButtons(self):
        fr=Frame(self)
        self.handler = StringVar()
        handlermenu = Pmw.OptionMenu(fr,
                menubutton_textvariable = self.handler,
                items = self.handlers,
                labelpos = 'w',
                label_text = 'Add Filter:',
                command=self.addSearch,
                initialitem = 'ekin',
                menubutton_width = 8)
        handlermenu.pack(side=LEFT)
        Button(fr, text='Clear', command=self.clearSearches).pack(side=LEFT,padx=2,pady=2)
        Button(fr, text='Close', command=self.close).pack(side=LEFT,padx=2,pady=2)
        Button(fr, text='Search', bg='lightblue',
                command=self.doSearch).pack(side=LEFT,padx=2,pady=2)
        fr.pack(side=TOP,fill=BOTH, expand=1)
        return

    def addSearch(self, evt=None):
        """Show filter dialog for current handler"""
        if self.handler.get() == 'ekin':
            S=ekinSearch(self.DB)
        elif  self.handler.get() == 'simple':
            S=simpleSearch(self.DB)
        elif  self.handler.get() == 'dict':
            S=dictSearch(self.DB)
        elif  self.handler.get() == 'structure':
            S=structureSearch(self.DB)
        frame=S.getWidget(parent=self.sframe, callback=self.doSearch)
        if frame == None:
            return
        frame.pack(anchor='nw',expand=1)
        self.currenthandlers.append(S)
        self.searchbars.append(frame)
        return

    def doSearch(self, evt=None):
        """Put searches together for all search handlers present"""
        if len(self.currenthandlers)==0:
            return
        F=[]
        SF=[]
        for s in self.currenthandlers:
            F.append(s.getFilter())
            SF.append(s.searchfunc)
        names = Base.Filter(F, SF)
        self.updateResults(names)
        print names
        return

    def clearSearches(self):
        for fr in self.searchbars:
            fr.destroy()
        return

    def updateResults(self, recs):
        self.resultsvar.set('found %s records' %len(recs))

    def close(self):
        """Destroy and remove from parent"""
        self.destroy()
        return

class Base(object):
    """Base methods for searching"""
    def __init__(self):
        return

    @classmethod
    def Filter(self, filters=None, searchfuncs=None):
        """Do multiple search terms"""
        sets = []
        for f,sfunc in zip(filters, searchfuncs):
            col, attr, val, op, boolean = f
            names = sfunc(col, attr, val, op)
            sets.append((set(names), boolean))
        names = sets[0][0]

        for s in sets[1:]:
            b=s[1]
            if b == 'AND':
                names = names & s[0]
            elif b == 'OR':
                names = names | s[0]
            elif b == 'NOT':
                names = names - s[0]
        names = list(names)
        return names

class Operators(object):

    def __init__(self):
        return

    @classmethod
    def contains(self,v1,v2):
        if v1 in v2:
            return True
    @classmethod
    def equals(self,v1,v2):
        if v1==v2:
            return True
    @classmethod
    def greaterthan(self,v1,v2):
        if v2>v1:
            return True
        return False
    @classmethod
    def lessthan(self,v1,v2):
        if v2<v1:
            return True
        return False
    @classmethod
    def startswith(self,v1,v2):
        if v2.startswith(v1):
            return True
    @classmethod
    def endswith(self,v1,v2):
        if v2.endswith(v1):
            return True
    @classmethod
    def haslength(self,v1,v2):
        if len(v2)>v1:
            return True

def test():
    DB = PDatabase(local='titdb.fs')

    #s = simpleSearch(DB)
    #names = s.doSearch(filters=[('stab', '3', '>', 'AND'),('choice', 'a', '=', 'AND')])
    e=ekinSearch(DB)
    names = e.doSearch(filters=[('1H NMR', 'residue', 'ASP', 'contains', 'AND')])
                               # ('pKas', 'model', '1 pKa', 'contains', 'AND'),
                                #('N+H pKas', 'model', '1 pKa', 'contains', 'AND')])
    print names

if __name__ == '__main__':

    test()
