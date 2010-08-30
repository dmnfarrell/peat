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

from Tkinter import *
import Pmw
from types import *

class FilterFrame(Frame):
        
    def __init__(self, parent, fields, callback=None, closecallback=None):
        """Create a filtering gui frame
        Callback must be some method that can accept tuples of filter
        parameters connected by boolean operators  """
        Frame.__init__(self, parent)
        self.parent = parent
        self.callback = callback
        self.closecallback = closecallback
        self.fields = fields
        self.filters = []        
        self.addFilterBar()
        addbutton=Button(self,text='Go', command=self.callback,bg='lightblue')
        addbutton.grid(row=0,column=0,sticky='news',padx=2,pady=2)         
        addbutton=Button(self,text='+Add Filter', command=self.addFilterBar)
        addbutton.grid(row=0,column=1,sticky='news',padx=2,pady=2)       
        cbutton=Button(self,text='Close', command=self.close)
        cbutton.grid(row=0,column=2,sticky='news',padx=2,pady=2)
        self.resultsvar=IntVar()
        Label(self,text='found:').grid(row=0,column=3,sticky='nes')
        Label(self,textvariable=self.resultsvar).grid(row=0,column=4,sticky='nws',padx=2,pady=2)
        return 
           
    def addFilterBar(self):
        """Add filter"""        
        index = len(self.filters)
        f=FilterBar(self, index, self.fields)
        self.filters.append(f)
        f.grid(row=index+1,column=0,columnspan=5,sticky='news',padx=2,pady=2)      
        return
    
    def close(self):
        """Close frame and do stuff in parent app if needed"""
        self.closecallback()
        self.destroy()

    def doFiltering(self, searchfunc, filters=None):       
        """Filter recs by several filters using user provided search function.
        Provides a list of tuples with filter values"""
        F=[]       
        for f in self.filters:
            F.append(f.getFilter())
        print F     
        sets = []
        for f in F:            
            col, val, op, boolean = f         
            names = searchfunc(col, val, op)
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
        self.updateResults(len(names))
        return names
        
    def updateResults(self, i):
        self.resultsvar.set(i)
        return
        
class FilterBar(Frame):
    """Class providing filter widgets"""
    operators = ['contains','=','>','<','starts with',
                 'ends with','has length']
    booleanops = ['AND','OR','NOT']                 
    def __init__(self, parent, index, fields):
        Frame.__init__(self, parent)
        self.parent=parent
        self.index = index
        self.filtercol=StringVar()
        if 'name' in fields:
            initial = 'name'
        else:
            initial = fields[0]        
        filtercolmenu = Pmw.OptionMenu(self,
                labelpos = 'w',
                label_text = 'Column:',
                menubutton_textvariable = self.filtercol,
                items = fields,
                initialitem = initial,
                menubutton_width = 10)
        filtercolmenu.grid(row=0,column=1,sticky='news',padx=2,pady=2)
        self.operator=StringVar()
        operatormenu = Pmw.OptionMenu(self,                
                menubutton_textvariable = self.operator,
                items = self.operators,
                initialitem = 'contains',
                menubutton_width = 8)
        operatormenu.grid(row=0,column=2,sticky='news',padx=2,pady=2)        
        self.filtercolvalue=StringVar()        
        valsbox=Entry(self,textvariable=self.filtercolvalue,width=20,bg='white')
        valsbox.grid(row=0,column=3,sticky='news',padx=2,pady=2)
        valsbox.bind("<Return>", self.parent.callback)
        self.booleanop=StringVar()
        booleanopmenu = Pmw.OptionMenu(self,                
                menubutton_textvariable = self.booleanop,
                items = self.booleanops,
                initialitem = 'AND',
                menubutton_width = 6)
        booleanopmenu.grid(row=0,column=0,sticky='news',padx=2,pady=2)        
        cbutton=Button(self,text='-', command=self.close)
        cbutton.grid(row=0,column=5,sticky='news',padx=2,pady=2)
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
        booleanop = self.booleanop.get()
        return col, val, op, booleanop
    
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
