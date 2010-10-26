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

import math, sys, os, types
from Base import EkinProject, EkinDataset
import IO
#from Convert import EkinConvert
from Fitting import Fitting
from Meta import MetaData
from PEATDB.Tables import TableCanvas, ColumnHeader
from PEATDB.TableModels import TableModel

            
class EkinDataModel(TableModel):
    """A model for managing the ekin datapoints"""
    def __init__(self, dataset=None, fitmodel=None):

        """create a model dict from ekin dataset and allow simple add/remove rows"""        
        ek = self.ek = dataset             
        modeldata = self.createData(ek)
        TableModel.__init__(self, modeldata)        
        return

    def createData(self, ek):
        """Populate table model from ekin dataset"""
        modeldata={}
        label0,label1 = ek.labels
        label2='xerr'
        label3='yerr'
        labels = self.labels = [label0,label1,label2,label3]        
        modeldata['columnnames']=labels
        modeldata['columntypes']={}; modeldata['columnorder']={}
        modeldata['columnlabels']={}
        i=0
        for l in labels:
            modeldata['columntypes'][l] = 'number'
            modeldata['columnorder'][i]=l
            i+=1
        for colname in modeldata['columnnames']:
            modeldata['columnlabels'][colname]=colname
        
        for i in range(ek.length()): 
            x=ek.x[i]; y=ek.y[i]
            modeldata[i] = {label0:x,label1:y}
            modeldata[i][label2] = ek.xerr(i)
            modeldata[i][label3] = ek.yerr(i)
        
        return modeldata
        
    def checkactivepoint(self, row):
        """Need to set color in a row """       
        i = self.reclist[row]        
        return self.ek.active[i]

    def setValueAt(self, value, rowIndex, columnIndex):
        """Overridden for ekin updates"""
        TableModel.setValueAt(self,value,rowIndex,columnIndex)
        
        ek = self.ek
        i = self.reclist[rowIndex]
        col=columnIndex
        if value=='':
            value=None
        else:
            value=float(value)
        if col==0:
            ek.x[i]=value
        elif col==1:
            ek.y[i]=value
        else:
            ek.errors[col-2][i] = value
          
        if hasattr(self, 'callback'):
            self.callback()
        return

    def setCallback(self, func):
        self.callback = func
        return

    def addRow(self, name):
        """Add row to end"""
        ek=self.ek   
        i = ek.length()
        TableModel.addRow(self, i)
        ek.add()
        return

    def deleteRow(self, row, update=True):
        """Delete row from the dataset and update table"""        
        i=self.reclist[row]
        self.ek.remove(i)
        if update==True:
            modeldata = self.createData(self.ek)
            self.setupModel(modeldata) 
        return

    def deleteRows(self, rowlist=None):  
        """Delete multiple or all rows"""
        if rowlist == None:
            rowlist = range(len(self.reclist))        
        self.ek.removeMultiple(rowlist)
        modeldata = self.createData(self.ek)
        self.setupModel(modeldata)      
        return
        

class EkinDataTable(TableCanvas):
    """Sub-class of Tablecanvas for display of datasets"""
    def __init__(self, parent=None, model=None, width=250,height=100):

        TableCanvas.__init__(self, parent, model, width=width,height=height)
        self.cellwidth=50
        self.rowheight=18
        self.thefont = "Arial 10"
        self.cellbackgr = 'white'
        self.autoresizecols = 1
        self.parent = parent
        return

    def createTableFrame(self):
        TableCanvas.createTableFrame(self)
        self.drawInactive()
        return

    def add_Row(self, rowname=None):
        """Add a new row"""
        self.model.addRow(rowname)
        self.redrawTable()
        return

    def drawInactive(self):
        """Need to set color for inactive points"""
        for row in self.rowrange:
            if self.model.checkactivepoint(row) == 0:
                self.setcellColor(row, newColor='gray90', key='bg', redraw=False)
        self.redrawTable()
        return

    def handle_point_activate(self, row, active=1):
        """Need to set color in a row if the corresponding point clicked in plotter"""
        if active == 0:
            self.setcellColor(row, newColor='yellow', key='bg')
        else:
            self.setcellColor(row, newColor=self.cellbackgr, key='bg')
        return

class EkinProjModel(TableModel):
    def __init__(self, E):
        self.E = E        
        TableModel.__init__(self) 
        self.addColumn('name'); self.addColumn('model')
        for d in E.datasets:            
            fdata = E.getFitData(d)            
            M = E.getMetaData(d)
            if E.currentmode == 'NMR titration':
                try:
                    edata = E.getDataset(d)
                    M['residue'] = edata['residue']
                    M['resnum'] = edata['res_num']
                    M['atom'] = edata['atom_type']
                except:
                    pass                
            if fdata != {} and fdata!= None:
                model = fdata['model']
            else:
                model=''
            self.addRecord(name=d, **M)
                
        return
    
class EkinProjTable(TableCanvas):
    """Sub-class of Tablecanvas """
    def __init__(self, parent=None, model=None, width=300,height=100):
        TableCanvas.__init__(self, parent, model, width=width,height=height)
        self.parent = parent
        return

