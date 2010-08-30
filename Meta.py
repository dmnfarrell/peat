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

import types
import Ekin.Utils as Utils
from Tkinter import *
import tkFileDialog, tkMessageBox, tkSimpleDialog
import Pmw

class MetaData(object):
    """Base Meta Data class - allows us to do simple/generic data management on dicts
       e.g for manipulating meta data by providing access views to the dict
       We should sublclass this for usage in peat. Ekin has it's own.
    """

    def __init__(self, data=None):
        if data != None:
            self.data = data
        else:
            self.initialise()
        return

    def initialise(self):
        """Create base empty data"""
        self.data = {}
        return

    def add(self, fieldname, data):
        """Add a new field, if no key name is provided, the field
           is global to all keys"""
        d=self.data
        if fieldname==None or data==None:
            print 'error, provide field name and data'
            return
        if not self.data.has_key(fieldname):
            d[fieldname] = {}
        d[fieldname] = data
        return

    def remove(self, field):
        """Remove a field for all datsets or just supplied name"""
        del self.data[key][field]
        return

    def rename(self, old, new):
        """rename a field"""
        if key!=None and self.data[key].has_key(old):
            temp = copy.deepcopy(self.data[old])
            self.data[new] = temp
            del self.data[old]

        return

    def get(self, fieldname=None):
        """Get a field"""
        d=self.data
        if fieldname==None:
            return d[key]
        else:
            return self.data[fieldname]
        return

    def getFields(self):
        fields = []
        for i in self.data:
            fields.append(i)
        print fields
        return fields

    def returnData(self):
        return self.data

    def printData(self):
        print self.data
        return

    def __repr__(self):
        """Overload string repr of class"""
        d = len(self.data)-1
        return "Meta data with %s keys" %d

    def prettyPrint(self):
        """prints the entire dict"""
        x=Utils.prettyPrint(self.data)
        print x
        return


    # --- Following methods are used from a graphical app and require tkinter

    def createDialog(self, showglob=True, parent=None, label=None):
        """Create a dialog for inputing fields, include all current global
           fields and ones local to the supplied key"""
        d=self.data
        self.metawin=Toplevel(parent)
        metawin = self.metawin
        metawin.title('Edit/View meta data')
        metawin.geometry('+200+200')
        self.keylabel=StringVar()
        self.keylabel.set(label)
        Label(metawin, font='Arial 16', fg='blue',
                    textvariable=self.keylabel).grid(row=0,column=0)


        self.updateFields()
        bf = Frame(metawin); bf.grid(row=0,column=2,rowspan=3)
        addfieldbtn=Button(bf,text='Add field',command=self.addFieldDialog)
        addfieldbtn.pack(fill=BOTH,expand=1)
        removefieldbtn=Button(bf,text='Remove field',command=self.removeFieldDialog)
        removefieldbtn.pack(fill=BOTH,expand=1)
        okbtn=Button(bf,text='Apply',command=self.applyCurrent)
        okbtn.pack(fill=BOTH,expand=1)
        Closebtn=Button(bf,text='Close',command=self.close)
        Closebtn.pack(fill=BOTH,expand=1)
        metawin.rowconfigure(0,weight=1)
        metawin.columnconfigure(0,weight=1)
        return metawin

    def updateFields(self):
        """Do field entry widgets for standard meta data items of the form
           data, attrib. Other fields are not displayed but are allowed """
        gvars=self.gvars={}
        cvars=self.cvars={}
        tvars=self.tvars={}
        d=self.data
        metawin = self.metawin
        glbals = self.data
        availabletypes = [types.StringType, types.FloatType, types.IntType]

        def addwidget(frame, fname, fdata, r):
            var=None
            tvar=IntVar(); tvar.set(0)
            if type(fdata) in availabletypes:
                if type(fdata) == types.IntType:
                    var=IntVar()
                else:
                    var=StringVar()
                var.set(fdata)
                Label(frame, text=fname).grid(row=r,column=0,sticky='news',padx=4,pady=3)
                Entry(frame, textvariable=var, bg='white').grid(row=r,column=1,
                                            sticky='news',padx=4,pady=3)
            else:
                Label(frame, text=fname).grid(row=r,column=0,sticky='news',padx=4,pady=3)
                st = Pmw.ScrolledText(frame, borderframe = 1,
                                        labelpos = 'n', label_text=str(type(fdata)),
                                        text_wrap='none')
                st.grid(row=r,column=1,sticky='news',padx=4,pady=3)
                p = Utils.prettyPrint(fdata) #more readable
                st.setvalue(p)
                st.component('text').configure(height=6)
                st.component('text').configure(state=DISABLED)

            return var, tvar

        #self.keylabel.set(key)

        gf = LabelFrame(metawin, text='fields')
        gf.grid(row=1,column=0,sticky='news',padx=4,pady=3)
        gf.columnconfigure(1,weight=1)
        r=0
        for g in glbals.keys():
            var, tvars[g] = addwidget(gf, g, glbals[g], r)
            if var != None:
                gvars[g]=var
            r+=1

        return

    def applyCurrent(self, key=None):
        """Apply current tkinter var data to field"""
        d=self.data
        for g in self.gvars:
            self.data[g] = self.gvars[g].get()

        return

    def close(self):
        self.metawin.destroy()

    def addFieldDialog(self):
        """Add a new field"""
        sd = NewFieldDialog(title="New field",
                              parent=self.metawin)
        if sd.result == None:
            return
        else:
            fieldtype = sd.result[0]
            newname = sd.result[1]
            data = sd.result[2]
            attrib = sd.result[3]
        if newname != None:
            if newname in self.getFields():
                tkMessageBox.showwarning("Name exists",
                                         "Name already exists!",
                                         parent=self.metawin)
            else:
                if fieldtype == 'string' or fieldtype == 'int':
                    self.add(newname, data)
                #else:
                #    self.add(newname, data)
        self.updateFields()
        return

    def removeFieldDialog(self):
        """Remove any checked field"""
        selected=[]
        #get list of checked fields
        for t in self.tvars:
            if self.tvars[t].get() == 1:
                selected.append(t)
        if tkMessageBox.askyesno("Delete?","Do you want to remove these fields?",
                                    parent=self.metawin):
            for s in selected:
                self.remove(s, self.currkey)
            self.updateFields(self.currkey)
        return


    def viewAll(self, parent=None):
        """Display all meta data in a table"""
        from PEATDB.Tables import TableCanvas
        from PEATDB.TableModels import TableModel

        tp = Toplevel(master=parent)
        tframe = Frame(tp)
        tframe.pack()
        m = TableModel()
        m.importDict(self.data)
        table = TableCanvas(tframe, model=m)
        table.createTableFrame()

        return


class NewFieldDialog(tkSimpleDialog.Dialog):
    """Simple dialog to get data for new fields"""

    def __init__(self, parent, title=None):
        self.items=['string','int','dict']
        tkSimpleDialog.Dialog.__init__(self, parent, title)

    def body(self, master):

        Label(master, text="Field type:").grid(row=0)
        Label(master, text="Name:").grid(row=1)
        Label(master, text="Data:").grid(row=2)

        self.v1=StringVar()
        self.v1.set('string')
        self.b1 = Menubutton(master,textvariable=self.v1,relief=RAISED)
        self.menu=Menu(self.b1,tearoff=0)
        self.b1['menu']=self.menu

        for option in self.items:
            self.menu.add_radiobutton(label=option,
                                          variable=self.v1,
                                          value=option,
                                          indicatoron=1)
        self.e2 = Entry(master);self.e3 = Entry(master);self.e4 = Entry(master)
        self.b1.grid(row=0, column=1,padx=2,pady=2,sticky='news')
        self.e2.grid(row=1, column=1,padx=2,pady=2,sticky='news')
        self.e3.grid(row=2, column=1,padx=2,pady=2,sticky='news')

        return self.b1 # initial focus

    def apply(self):
        first = self.v1.get()
        second = self.e2.get()
        third = self.e3.get()
        fourth = self.e4.get()
        self.result = first, second, third, fourth
        return

class PEATMetaData(MetaData):
    def __init__(self, data=None):
        MetaData.__init__(self, data=data)
        return

