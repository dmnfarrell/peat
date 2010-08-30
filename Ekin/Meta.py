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
import Utils
from Tkinter import *
import tkFileDialog, tkMessageBox, tkSimpleDialog
import Pmw

class MetaData(object):
    """Ekin Meta Data class - allows us to manage meta data stored in ekin
       (or some similar app), by providing access methods to the metadata etc.
    """
    modetemplate = ['General',
                'Simple enzyme kinetic',
                'Enzyme pH-activity',
                'pH-stability',
                'NMR titration',
                'Protein Stability',
                'Amyloid formation']

    def __init__(self, data=None, project=None):

        if data != None:
            self.data = data
        else:
            self.initialise()
        if not self.data.has_key('glob'):
            self.data['glob'] = {}
        self.glob = self.data['glob']
        if project != None:
            self.project = project
        return

    def initialise(self):
        """Create base empty meta data"""
        self.data = {}
        self.data['glob'] = {}
        return

    def add(self, fieldname, data, attrib='', dataset=None):
        """Add a new field, if no dataset name is provided, the field
           is global to all datasets"""
        d=self.data
        if fieldname==None or data==None:
            print 'error, provide field name and data'
            return
        if dataset==None:
            if not d.has_key('glob'):
                self.glob = {}
            self.glob[fieldname]={}
            self.glob[fieldname]['data'] = data
            self.glob[fieldname]['attrib'] = attrib
        else:
            if not d.has_key(dataset):
                d[dataset] = {}
            d[dataset][fieldname] = {}
            d[dataset][fieldname]['data'] = data
            d[dataset][fieldname]['attrib'] = attrib
        return

    def remove(self, field, dataset=None):
        """Remove a field for all datsets or just supplied name"""
        if field in self.glob:
            del self.glob[field]
        if dataset==None:
            for d in self.data:
                if field in self.data[d]:
                    #print d,field
                    del self.data[d][field]
        else:
            del self.data[dataset][field]
        return

    def rename(self, old, new, dataset=None):
        """rename a field"""
        if dataset!=None and self.data[dataset].has_key(old):
            temp = copy.deepcopy(self.data[old])
            self.data[new] = temp
            del self.data[old]
        elif field in self.glob:
            pass

        return

    def setAttrib(self, field, attrib, dataset=None):
        """set an attribute"""
        if field in self.glob:
            self.glob[field]['attrib']=attrib
        elif dataset!=None and self.data[dataset].has_key(field):
            self.data[dataset][field]['attrib']=attrib

        return

    def get(self, fieldname=None, dataset=None):
        """Get a field"""
        d=self.data
        if dataset!=None:
            if not d.has_key(dataset):
                print 'no such dataset'
                return
            else:
                return d[dataset][fieldname]['data']
        else:
            if fieldname==None:
                return d[dataset]
            elif fieldname not in self.glob:
                print 'no such field'
                return
            else:
                return self.glob[fieldname]['data']
        return

    def getFields(self):
        fields = self.glob.keys()
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
        g = len(self.glob)
        return "Meta data with %s datasets and %s global fields" %(d, g)

    def prettyPrint(self):
        """prints the entire dict"""
        x=Utils.prettyPrint(self.data)
        print x
        return

    #
    # --- Following methods are used from a graphical app and require tkinter
    #

    def createDialog(self, dataset=None, showglob=True, parent=None):
        """Create a dialog for inputing fields, include all current global
           fields and ones local to the supplied dataset"""
        d=self.data
        
        self.currdataset = dataset
        self.metawin=Toplevel(parent)
        metawin = self.metawin
        metawin.title('Edit/View meta data')
        metawin.geometry('+200+200')
        self.datasetlabel=StringVar()
        self.datasetlabel.set(dataset)
        Label(metawin, font='Arial 16', fg='blue',
                    textvariable=self.datasetlabel).grid(row=0,column=0)

        if not d.has_key('glob'):
            d['glob'] = {}
            self.glob = d['glob']

        self.updateFields(dataset)
        bf = Frame(metawin); bf.grid(row=0,column=2,rowspan=3)
        addfieldbtn=Button(bf,text='Add field',command=self.addFieldDialog)
        addfieldbtn.pack(fill=BOTH,expand=1)
        removefieldbtn=Button(bf,text='Remove field',command=self.removeFieldDialog)
        removefieldbtn.pack(fill=BOTH,expand=1)
        editattribbtn=Button(bf,text='Edit attrib',command=self.editAttribDialog)
        editattribbtn.pack(fill=BOTH,expand=1)
        okbtn=Button(bf,text='Apply',command=self.applyCurrent)
        okbtn.pack(fill=BOTH,expand=1)
        Closebtn=Button(bf,text='Close',command=self.close)
        Closebtn.pack(fill=BOTH,expand=1)
        metawin.rowconfigure(0,weight=1)
        metawin.columnconfigure(0,weight=1)
        return metawin

    def updateFields(self, dataset=None):
        """Do field entry widgets for standard meta data items of the form
           data, attrib. Other fields are not displayed but are allowed """
        gvars=self.gvars={}
        cvars=self.cvars={}
        tvars=self.tvars={}
        d=self.data
        metawin = self.metawin
        glbals = self.glob
        self.currdataset = dataset
        def addwidget(frame, fname, fdata, fattrib, r):
            var=None
            tvar=IntVar(); tvar.set(0)
            if type(fdata) is types.StringType or type(fdata) is types.FloatType:
                var=StringVar()
                var.set(fdata)
                Label(frame, text=fname).grid(row=r,column=0,sticky='news',padx=4,pady=3)
                Entry(frame, textvariable=var, bg='white').grid(row=r,column=1,
                                            sticky='news',padx=4,pady=3)
            else:
                Label(frame, text=c).grid(row=r,column=0,sticky='news',padx=4,pady=3)
                st = Pmw.ScrolledText(frame, borderframe = 1,
                                        labelpos = 'n', label_text=str(type(fdata)),
                                        text_wrap='none')
                st.grid(row=r,column=1,sticky='news',padx=4,pady=3)
                p = Utils.prettyPrint(fdata) #more readable
                st.setvalue(p)
                st.component('text').configure(height=6)
                st.component('text').configure(state=DISABLED)
            #attributes
            Label(frame, text=fattrib).grid(row=r,column=2,
                                            sticky='news',padx=4,pady=3)
            Checkbutton(frame, variable=tvar).grid(row=r,column=3,sticky='news',padx=2,pady=2)


            return var, tvar

        self.datasetlabel.set(dataset)

        gf = LabelFrame(metawin, text='global')
        gf.grid(row=1,column=0,sticky='news',padx=4,pady=3)
        gf.columnconfigure(1,weight=1)
        r=0
        for g in glbals:
            var, tvars[g] = addwidget(gf, g, glbals[g]['data'], glbals[g]['attrib'], r)
            if var != None:
                gvars[g]=var
            r+=1
        if dataset!=None and d.has_key(dataset):
            currdata=d[dataset]
            cf = LabelFrame(metawin, text='dataset specific')
            cf.grid(row=2,column=0,sticky='news',padx=4,pady=3)
            cf.columnconfigure(1,weight=1)
            r=0
            for c in currdata:
                try:
                    itemdata = currdata[c]['data']
                    var, tvars[c] = addwidget(cf, c, itemdata, currdata[c]['attrib'], r)
                except:
                    itemdata = currdata[c]
                    var, tvars[c] = addwidget(cf, c, itemdata, None, r)
                if var != None:
                    cvars[c]=var
                r+=1
        return

    def applyCurrent(self, dataset=None):
        """Apply current tkinter var data to field"""
        d=self.data
        if dataset==None:
            dataset=self.currdataset
        glbals = self.glob

        for g in self.gvars:
            self.glob[g]['data'] = self.gvars[g].get()
        for c in self.cvars:
            d[dataset][c]['data'] = self.cvars[c].get()
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
                if fieldtype == 'global':
                    self.add(newname, data, attrib, None)
                elif fieldtype == 'all datsets':
                    for d in self.project.datasets:
                        self.add(newname, data, attrib, d)
                else:
                    self.add(newname, data, attrib, self.currdataset)
        self.updateFields(self.currdataset)
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
                self.remove(s, self.currdataset)
            self.updateFields(self.currdataset)
        return

    def editAttribDialog(self):
        """Edit attribute of field(s)"""
        selected=[]
        d=self.currdataset
        for t in self.tvars:
            if self.tvars[t].get() == 1:
                selected.append(t)
        val = tkSimpleDialog.askstring("New attribute value",
                                            "Enter new value:",
                                             parent=self.metawin)
        for s in selected:
            self.setAttrib(s, val, d)
        self.updateFields(d)
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

    def convertfromOld(self, E=None):
        """Convert from old meta data format """
        if E == None:
            E = self.project
        old = E.__Exp_Meta_Dat__
        #print old
        for f in old['Operator']:
            self.add(f, old['Operator'][f])
        for i in old['Meta_Data']:
            #print i
            f = i['Name']
            try:
                units=i['units']
            except:
                units=None
            #print f
            self.add(f, i['value'], units)

        return

    def autoCreate(self, mode):
        import time
        """Auto create some meta data specific to the ekin mode"""
        if mode=='General' or mode=='NMR titration':
            self.data['glob'] = {  'Buffer': {   'attrib': 'mM', 'data': 'KCl'},
                        'Comment': {   'attrib': None, 'data': ''},
                        'D2O concentration': {   'attrib': '%', 'data': ''},
                        'Experimental protocol': {   'attrib': None, 'data': ''},
                        'Literature Reference': {   'attrib': None, 'data': ''},
                        'Protein concentration': {   'attrib': 'mM', 'data': ''},
                        'Sample volume': {   'attrib': 'mL', 'data': ''},
                        'Temperature': {   'attrib': 'C', 'data': '35'},
                        'date': {   'attrib': '', 'data': time.asctime()},
                        'operator': {   'attrib': '', 'data': ''},
                        'pH': {   'attrib': 'pHU', 'data': ''},
                       'status': {   'attrib': '', 'data': ''}}
        elif mode=='Simple enzyme kinetic':
            self.data['glob'] = {
                        'Comment': {   'attrib': None, 'data': ''},
                        'Experimental protocol': {   'attrib': None, 'data': ''},
                        'Literature Reference': {   'attrib': None, 'data': ''},
                        'Protein concentration': {   'attrib': 'mM', 'data': ''},
                        'Sample volume': {   'attrib': 'mL', 'data': 'x'},
                        'Temperature': {   'attrib': 're', 'data': '35'},
                        'date': {   'attrib': '', 'data': time.asctime()},
                        'operator': {   'attrib': '', 'data': ''},
                        'pH': {   'attrib': 'pHU', 'data': ''},
                        'status': {   'attrib': '', 'data': ''}}

        elif mode=='Enzyme pH-activity' or mode=='Protein Stability':
            self.data['glob'] = {
                        'Comment': {   'attrib': None, 'data': ''},
                        'Enzyme concentration': {   'attrib': 'mM', 'data': ''},
                        'Experimental protocol': {   'attrib': None, 'data': ''},
                        'Literature Reference': {   'attrib': None, 'data': ''},
                        'Temperature': {   'attrib': 'C', 'data': '35'},
                        'date': {   'attrib': '', 'data': time.asctime()},
                        'operator': {   'attrib': '', 'data': ''},
                        'pH': {   'attrib': 'pHU', 'data': ''},
                        'status': {   'attrib': '', 'data': ''}}

        else:
            print 'no presets for this mode'
        return


class NewFieldDialog(tkSimpleDialog.Dialog):
    """Simple dialog to get data for new fields"""

    def __init__(self, parent, title=None):
        self.items=['global','current dataset','all datsets']
        tkSimpleDialog.Dialog.__init__(self, parent, title)

    def body(self, master):

        Label(master, text="Field type:").grid(row=0)
        Label(master, text="Name:").grid(row=1)
        Label(master, text="Data:").grid(row=2)
        Label(master, text="Attrib:").grid(row=3)

        self.v1=StringVar()
        self.v1.set('global')
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
        self.e4.grid(row=3, column=1,padx=2,pady=2,sticky='news')
        return self.b1 # initial focus

    def apply(self):
        first = self.v1.get()
        second = self.e2.get()
        third = self.e3.get()
        fourth = self.e4.get()
        self.result = first, second, third, fourth
        return

