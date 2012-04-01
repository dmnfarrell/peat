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

from TableFormula import Formula
from types import *
import copy

class TableModel(object):
    """A base model for managing the data in a TableCanvas class"""

    keywords = {'columnnames':'columnNames', 'columntypes':'columntypes',
               'columnlabels':'columnlabels', 'columnorder':'columnOrder',
               'colors':'colors'}

    def __init__(self, newdict=None, rows=None, columns=None):
        """Constructor"""
        self.initialiseFields()
        self.setupModel(newdict, rows, columns)
        #finally set default sort order as first col
        #self.setSortOrder()
        return

    def setupModel(self, newdict, rows=None, columns=None):
        """Create table model"""
        if newdict != None:
            self.data = copy.deepcopy(newdict)
            for k in self.keywords:
                if self.data.has_key(k):
                    self.__dict__[self.keywords[k]] = self.data[k]
                    del self.data[k]
            #read in the record list order
            if self.data.has_key('reclist'):
                temp = self.data['reclist']
                del self.data['reclist']
                self.reclist = temp
            else:
                self.reclist = self.data.keys()
        else:
            #just make a new empty model
            self.createEmptyModel()
       
        if not set(self.reclist) == set(self.data.keys()):
            print 'reclist does not match data keys'
        #restore last column order
        if hasattr(self, 'columnOrder') and self.columnOrder != None:
            self.columnNames=[]
            for i in self.columnOrder.keys():
                self.columnNames.append(self.columnOrder[i])
                i=i+1
        self.defaulttypes = ['text', 'number']
        #setup default display for column types
        self.default_display = {'text' : 'showstring',
                                'number' : 'numtostring'}
        #add rows and cols if they are given in the constructor
        if newdict == None:
            if rows != None:
                self.auto_AddRows(rows)
            if columns != None:
                self.auto_AddColumns(columns)
        self.filteredrecs = None          
        return
        
    def initialiseFields(self):
        """Create base fields, some of which are not saved"""
        self.data = None    # holds the table dict
        self.colors = {}    # holds cell colors
        self.colors['fg']={}
        self.colors['bg']={}
        #default types
        self.defaulttypes = ['text', 'number']
        #list of editable column types
        self.editable={}
        self.nodisplay = []
        self.columnwidths={}  #used to store col widths, not held in saved data
        return

    def createEmptyModel(self):
        """Create the basic empty model dict"""
        self.data = {}
        # Define the starting column names and locations in the table.
        self.columnNames = []
        self.columntypes = {}
        self.columnOrder = None
        #record column labels for use in a table header
        self.columnlabels={}
        for colname in self.columnNames:
            self.columnlabels[colname]=colname
        self.reclist = self.data.keys()
        return

    def importDict(self, newdata, namefield=None):
        """Try to create a table model from some arbitrary dict"""
        import types
        if namefield == None:
            namefield = 'name'
        #get cols from sub data keys
        colnames = []
        colnames.append(namefield)
        
        for k in newdata:
            fields = newdata[k].keys()
            for f in fields:
                if not f in colnames:
                    colnames.append(f)

        for c in colnames:
            self.addColumn(c)

        #add the data
        #print colnames
        for k in newdata:
            self.addRow(k)
            for c in colnames:
                if c == namefield:
                    self.data[k][namefield] = str(k)
                else:
                    self.data[k][c] = str(newdata[k][c])

        return
    
    def getDefaultTypes(self):
        """Get possible field types for this table model"""
        return self.defaulttypes

    def getData(self):
        """Return the current data for saving"""
        import copy
        data = copy.deepcopy(self.data)
        data['colors'] = self.colors
        data['columnnames'] = self.columnNames
        #we keep original record order
        data['reclist'] = self.reclist
        #record current col order
        data['columnorder']={}
        i=0
        for name in self.columnNames:
            data['columnorder'][i] = name
            i=i+1
        data['columntypes'] = self.columntypes
        data['columnlabels'] = self.columnlabels
        return data

    def getAllCells(self):
        """Return a dict of the form rowname: list of cell contents
          Useful for a simple table export for example"""
        records={}
        for row in range(len(self.reclist)):
            recdata=[]
            for col in range(len(self.columnNames)):
                recdata.append(self.getValueAt(row,col))
            records[row]=recdata
        return records

    def getColCells(self, colIndex):
        """Get the viewable contents of a col into a list"""
        collist = []
        if self.getColumnType(colIndex) == 'Link':
            return ['xxxxxx']
        else:        
            for row in range(len(self.reclist)):
                v = self.getValueAt(row, colIndex)
                collist.append(v)
        return collist

    def getlongestEntry(self, columnIndex):
        """Get the longest cell entry in the col"""
        collist = self.getColCells(columnIndex)        
        maxw=5
        for c in collist:
            try:                
                w = len(str(c))
            except UnicodeEncodeError:
                pass
            if w > maxw:
                maxw = w
        #print 'longest width', maxw
        return maxw

    def getRecordAtRow(self, rowIndex):
        """Get the entire record at the specifed row."""
        name = self.getRecName(rowIndex)
        record = self.data[name]
        return record

    def getCellRecord(self, rowIndex, columnIndex):
        """Get the data held in this row and column"""
        value = None
        colname = self.getColumnName(columnIndex)
        coltype = self.columntypes[colname]
        name = self.getRecName(rowIndex)    
        if self.data[name].has_key(colname):
            celldata=self.data[name][colname]
        else:
            celldata=None
        return celldata

    def deleteCellRecord(self, rowIndex, columnIndex):
        """Remove the cell data at this row/column"""
        colname = self.getColumnName(columnIndex)
        coltype = self.columntypes[colname]
        name = self.getRecName(rowIndex)
        if self.data[name].has_key(colname):
            del self.data[name][colname]
        return

    def getRecName(self, rowIndex):
        """Get record name from row number"""
        if len(self.reclist)==0:
            return None
        if self.filteredrecs != None:
            name = self.filteredrecs[rowIndex]
        else:
            name = self.reclist[rowIndex]
        return name

    def setRecName(self, newname, rowIndex):
        """Set the record name to another value - requires re-setting in all
           dicts that this rec is referenced"""
        if len(self.reclist)==0:
            return None
        currname = self.getRecName(rowIndex)
        self.reclist[rowIndex] = newname
        import copy
        temp = copy.deepcopy(self.data[currname])
        self.data[newname] = temp
        self.data[newname]['Name'] = newname
        del self.data[currname]
        for key in ['bg', 'fg']:
            if self.colors[key].has_key(currname):
                temp = copy.deepcopy(self.colors[key][currname])
                self.colors[key][newname] = temp
                del self.colors[key][currname]
        print 'renamed'
        #would also need to resolve all refs to this rec in formulas here!

        return

    def getRecordAttributeAtColumn(self, rowIndex=None, columnIndex=None,
                                        recName=None, columnName=None):
         """Get the attribute of the record at the specified column index.
            This determines what will be displayed in the cell"""

         value = None         
         if columnName != None and recName != None:
             if not self.data[recName].has_key(columnName):
                 return ''
             cell = self.data[recName][columnName]
         else:
             cell = self.getCellRecord(rowIndex, columnIndex)
             columnName = self.getColumnName(columnIndex)
         if cell == None:
             cell=''
         # Set the value based on the data record field
         coltype = self.columntypes[columnName]
         if Formula.isFormula(cell) == True:
             value = self.doFormula(cell)
             return value
           
         if not type(cell) is DictType:             
             if coltype == 'text' or coltype == 'Text':
                 value = cell
             elif coltype == 'number':
                 value = str(cell)
             else:
                 value = 'other'
         if value==None:
             value=''

         return value

    def getRecordIndex(self, recname):
        rowIndex = self.reclist.index(recname)
        return rowIndex

    def setSortOrder(self, columnIndex=0, reverse=0):
        """Changes the order that records are sorted in, which will
           be reflected in the table upon redrawing"""
        
        self.sortcolumnIndex = columnIndex    
        sortkey = self.getColumnName(columnIndex)
        recnames = self.reclist
        
        self.reclist = self.createSortMap(self.reclist, sortkey, reverse)
        if self.filteredrecs != None:
            self.filteredrecs = self.createSortMap(self.filteredrecs, sortkey, reverse)
        return

    def createSortMap(self, names, sortkey, reverse=0):
        """Create a sort mapping for given list"""
        import operator
        recdata = []
        for rec in names:            
            recdata.append(self.getRecordAttributeAtColumn(recName=rec, columnName=sortkey))
        #try create list of floats if col has numbers only
        try:            
            recdata = self.toFloats(recdata)
        except:
            pass
        smap = zip(names, recdata)
        #sort the mapping by the second key
        smap = sorted(smap, key=operator.itemgetter(1), reverse=reverse)        
        #now sort the main reclist by the mapping order
        sortmap = map(operator.itemgetter(0), smap)
        return sortmap
    
    def toFloats(self, l):
        x=[]
        for i in l:
            if i == '':
                x.append(0.0)
            else:    
                x.append(float(i))
        return x
    
    def getSortIndex(self):
        """Return the current sort order index"""
        if self.sortcolumnIndex:
            return self.sortcolumnIndex
        else:
            return 0

    def moveColumn(self, oldcolumnIndex, newcolumnIndex):
        """Changes the order of columns"""
        self.oldnames = self.columnNames
        self.columnNames=[]
        #self.columnOrder=[]

        print oldcolumnIndex, newcolumnIndex
        #write out a new column names list - tedious
        moved = self.oldnames[oldcolumnIndex]
        del self.oldnames[oldcolumnIndex]
        print self.oldnames
        i=0
        for c in self.oldnames:
            if i==newcolumnIndex:
                self.columnNames.append(moved)
                #self.columnOrder.append(moved)
            self.columnNames.append(c)
            #self.columnOrder.append(c)
            i=i+1
        #if new col is at end just append
        if moved not in self.columnNames:
            self.columnNames.append(moved)
            #self.columnOrder.append(moved)

        print self.columnNames
  
        return

    def addRow(self, name=None, **kwargs):
        """Add a row"""
        if self.data.has_key(name) or name in self.reclist:
            print 'name already present!!'
            return
        self.data[name]={}
        if name != None:
            self.data[name]['Name'] = name
        else:
            self.data[name]['Name'] = ''
        self.reclist.append(name)       

        return

    def deleteRow(self, rowIndex, update=True):
        """Delete a row"""
        name = self.getRecName(rowIndex)
        del self.data[name]
        if update==True:
            self.reclist = self.data.keys()

        return

    def deleteRows(self, rowlist=None):
        """Delete multiple or all rows"""
        if rowlist == None:
            rowlist = range(len(self.reclist))
        #print 'deleting' , rowlist
        #print 'reclist', self.reclist
        for row in rowlist:
            self.deleteRow(row, update=False)
        self.reclist = self.data.keys()

        return

    def addColumn(self, colname=None, coltype=None):
        """Add a column"""
        index = self.getColumnCount()+ 1
        if colname == None:
            colname=str(index)
        if colname in self.columnNames:
            #print 'name is present!'
            return
        self.columnNames.append(colname)
        self.columnlabels[colname] = colname
        if coltype == None:
            self.columntypes[colname]='text'
        else:
            self.columntypes[colname]=coltype

        return

    def deleteColumn(self, columnIndex):
        """delete a column"""
        colname = self.getColumnName(columnIndex)
        print colname
        self.columnNames.remove(colname)
        del self.columnlabels[colname]
        del self.columntypes[colname]
        #remove this field from every record
        for recname in self.reclist:
            if self.data[recname].has_key(colname):
                del self.data[recname][colname]

        if hasattr(self, 'sortcolumnIndex') and columnIndex == self.sortcolumnIndex:
            self.setSortOrder()
        print 'column deleted'
        print 'new cols:', self.columnNames
        return

    def deleteMultipleColumns(self, cols=None):
        """Remove all cols or list provided"""

        while self.getColumnCount() > 0:
            self.deleteColumn(0)
        return

    def auto_AddRows(self, numrows=None):
        """Automatically add x number of records"""
        import string
        alphabet = string.lowercase[:26]
        rows = self.getRowCount()

        if rows <= 25:
            i=rows
            j=0
        else:
            i=int(rows%25)
            j=int(round(rows/25,1))
        #print i, j
        for x in range(numrows):
            if i >= len(alphabet):
                i=0
                j=j+1
            name = alphabet[i]+str(j)
            if name in self.reclist:
                pass
            else:
                self.addRow(name)
            i=i+1
            #print self.reclist
        return

    def auto_AddColumns(self, numcols=None):
        """Automatically add x number of cols"""
        import string
        alphabet = string.lowercase[:26]
        currcols=self.getColumnCount()
        #find where to start
        start = currcols + 1
        end = currcols + numcols + 1
        new = []
        for n in range(start, end):
            new.append(str(n))
        #check if any of these colnames present
        common = set(new) & set(self.columnNames)
        extra = len(common)
        end = end + extra
        for x in range(start, end):
            self.addColumn(str(x))
        return

    def relabel_Column(self, columnIndex, newname):
        """Change the column label - can be used in a table header"""
        colname = self.getColumnName(columnIndex)
        self.columnlabels[colname]=newname
        return

    def getColumnType(self, columnIndex):
        """Get the column type"""
        colname = self.getColumnName(columnIndex)
        coltype = self.columntypes[colname]
        return coltype

    def getColumnCount(self):
         """Returns the number of columns in the data model."""
         return len(self.columnNames)

    def getColumnName(self, columnIndex):
         """Returns the name of the given column by columnIndex."""        
         return self.columnNames[columnIndex]

    def getColumnLabel(self, columnIndex):
        """Returns the label for this column"""
        colname = self.getColumnName(columnIndex)
        return self.columnlabels[colname]

    def getColumnIndex(self, columnName):
        """Returns the column index for this column"""
        colindex = self.columnNames.index(columnName)
        return colindex

    def getColumnData(self, columnIndex=None, columnName=None,
                      filterby=None):
        """Return the data in a list for this col,
            filterby is a tuple of key, value that allows to filter by rec attribute
        """
        import types
        if columnIndex != None:
            columnName = self.getColumnName(columnIndex)
        coldata = []
        if filterby != None and filterby[0] in self.columnNames:
            f, vals = filterby
            if type(vals) is not types.ListType:
                vals = [vals]
            for r in self.reclist:
                if self.data[r].has_key(f) and self.data[r][f] in vals:
                    if not self.data[r].has_key(columnName):
                        coldata.append('')
                    else:
                        coldata.append(self.data[r][columnName])
        else:
            for r in self.reclist:                                  
                if not self.data[r].has_key(columnName):
                    coldata.append('')
                else:
                    coldata.append(self.data[r][columnName])
        return coldata
        
    def getColumns(self, colnames, filterby=None, allowempty=True):
        """Get column data for multiple cols, with given filter options,
            filterby: a tuple of key value that allows filtering
            allowempty: boolean if false means rows with empty vals for any
            required fields are not returned
            returns: lists of column data"""
        def evaluate(l):
            for i in l:
                if i == '' or i == None:
                    return False
            return True        
          
        coldata=[]       
        for c in colnames:
            vals = self.getColumnData(columnName=c, filterby=filterby)            
            coldata.append(vals)          
        if allowempty == False:  
            result = [i for i in zip(*coldata) if evaluate(i) == True]               
            coldata = zip(*result)        
        return coldata
        
    def getDict(self, colnames, filterby=None):  
        """Get the model data as a dict for given columns with filter options"""
        data={}
        names, = self.getColumns(['name'], filterby)        
        cols = self.getColumns(colnames, filterby)
        coldata = zip(*cols)   
        for name,cdata in zip(names, coldata):           
            data[name] = dict(zip(colnames,cdata))
        return data
        
    def filterBy(self, filtercol, value, op='contains', userecnames=False,
                     progresscallback=None):
        """Filter recs according to one column"""
        
        from Filtering import Operators        
        funcs = {'contains':Operators.contains,'=':Operators.equals,
                   '>':Operators.greaterthan,'<':Operators.lessthan,
                   'starts with':Operators.startswith,
                   'ends with':Operators.endswith,
                   'has length':Operators.haslength}
        floatops = ['=','>','<']                   
        func = funcs[op]
        data = self.data
        #coltype = self.columntypes[filtercol]
        names=[]
        for rec in self.reclist:
            if data[rec].has_key(filtercol):  
                #try to do float comparisons if required
                if op in floatops:
                    try:
                        item = float(data[rec][filtercol])
                        v=float(value)
                        if func(v, item) == True:                    
                            names.append(rec)
                        continue 
                    except:
                        pass
                if filtercol == 'name' and userecnames == True:
                    item = rec
                else:
                    item = str(data[rec][filtercol])                
                if func(value, item):                    
                    names.append(rec)               
            
        return names        
        
    def filterByExpr(self, expr):
        """Filter recs using simple expression"""
        import shlex
        criteria=[]
        if filterby != None:
            lexer = shlex.shlex(filterby)
            expr = list(lexer)
            print expr
            i=0
            for e in range(0,len(expr),3):
                print e
                try:
                    criteria.append((expr[e],expr[e+2]))
                except:
                    pass

            print criteria

    def getRowCount(self):
         """Returns the number of rows in the table model."""
         return len(self.reclist)

    def getValueAt(self, rowIndex, columnIndex):
         """Returns the cell value at location specified
             by columnIndex and rowIndex."""        
         value = self.getRecordAttributeAtColumn(rowIndex, columnIndex)
         return value

    def setValueAt(self, value, rowIndex, columnIndex):
        """Changed the dictionary when cell is updated by user"""
        name = self.getRecName(rowIndex)
        colname = self.getColumnName(columnIndex)
        coltype = self.columntypes[colname]
      
        if coltype == 'number':
            try:
                if value == '': #need this to allow deletion of values
                    self.data[name][colname] = ''
                else:
                    self.data[name][colname] = float(value)
            except:
                pass
        else:
            self.data[name][colname] = value       
        return

    def setFormulaAt(self, f, rowIndex, columnIndex):
        """Set a formula at cell given"""
        name = self.getRecName(rowIndex)
        colname = self.getColumnName(columnIndex)
        coltype = self.columntypes[colname]
        rec = {}
        rec['formula'] = f
        self.data[name][colname] = rec
        return

    def getColorAt(self, rowIndex, columnIndex, key='bg'):
        """Return color of that record field for the table"""
        name = self.getRecName(rowIndex)
        colname = self.getColumnName(columnIndex)
        if self.colors[key].has_key(name) and self.colors[key][name].has_key(colname):
            return self.colors[key][name][colname]
        else:
            return None

    def setColorAt(self, rowIndex, columnIndex, color, key='bg'):
        """Set color"""
        name = self.getRecName(rowIndex)
        colname = self.getColumnName(columnIndex)
        if not self.colors[key].has_key(name):
            self.colors[key][name] = {}
        self.colors[key][name][colname] = str(color)        
        return

    def resetcolors(self):
        """Remove all color formatting"""
        self.colors={}
        self.colors['fg']={}
        self.colors['bg']={}
        return

    def getRecColNames(self, rowIndex, ColIndex):
        """Returns the rec and col name as a tuple"""
        recname = self.getRecName(rowIndex)
        colname = self.getColumnName(ColIndex)
        return (recname, colname)

    def getRecAtRow(self, recname, colname, offset=1, dim='y'):
        """Get the record name at a specified offset in the current
           table from the record given, by using the current sort order"""
        thisrow = self.getRecordIndex(recname)
        thiscol = self.getColumnIndex(colname)
        #table goto next row
        if dim == 'y':
            nrow = thisrow + offset
            ncol = thiscol
        else:
            nrow = thisrow
            ncol = thiscol + offset

        newrecname, newcolname = self.getRecColNames(nrow, ncol)
        print 'recname, colname', recname, colname
        print 'thisrow, col', thisrow, thiscol
        return newrecname, newcolname

    def appendtoFormula(self, formula, rowIndex, colIndex):
        """Add the input cell to the formula"""
        cellRec = getRecColNames(rowIndex, colIndex)
        formula.append(cellRec)
        return

    def doFormula(self, cellformula):
        """Evaluate the formula for a cell and return the result"""
        value = Formula.doFormula(cellformula, self.data)
        return value

    def copyFormula(self, cellval, row, col, offset=1, dim='y'):
        """Copy a formula down or across, using the provided offset"""
        import re
        frmla = Formula.getFormula(cellval)
        #print 'formula', frmla

        newcells=[]
        cells, ops = Formula.readExpression(frmla)

        for c in cells:
            print c
            if type(c) is not ListType:
                nc = c
            else:
                recname = c[0]
                colname = c[1]
                nc = list(self.getRecAtRow(recname, colname, offset, dim=dim))
            newcells.append(nc)
        newformula = Formula.doExpression(newcells, ops, getvalues=False)
        return newformula

    def merge(self, model, key='name', fields=None):
        """Merge another table model with this one based on a key field, 
           we only add nrecords from the new model where the key is present
           in both models"""
        if fields == None: fields = model.columnNames
        for rec in self.reclist:
            if not self.data[rec].has_key(key):
                continue
            for new in model.reclist:
                if not model.data[new].has_key(key):
                    continue
                if self.data[rec][key] == model.data[new][key]:
                #if new == rec:                                   
                    for f in fields:
                        if not model.data[rec].has_key(f):
                            continue
                        if not f in self.columnNames:
                            self.addColumn(f)                        
                        self.data[rec][f] = model.data[rec][f]       
        return
        
    def addRecord(self, name, **kwargs):
        """Add new rec with all fields in kwargs dict"""
        #name = kwargs[kwargs.keys()[0]]
        self.addRow(name) 
        self.addColumn('name')
        self.data[name]['name'] = name          
        for k in kwargs:            
            if not k in self.columnNames:
                self.addColumn(k)
            self.data[name][k] = str(kwargs[k])
        return

    def save(self, filename=None):
        """Save model to file"""
        if filename == None:
            return
        data = self.getData()
        fd=open(filename,'w')
        import pickle
        pickle.dump(data,fd)
        fd.close()
        return

    def load(self):
        """Load model from pickle file"""        
        return

    def copy(self):
        """Return a copy of this model"""
        M=TableModel()
        data = self.getData()
        M.setupModel(data)
        return M
    
    def exportCSV(self, filename=None, sep=','):
        """Export directly from model"""
        import csv
        if filename==None: return
        #take column labels as field names 
        writer = csv.writer(file(filename, "w"), delimiter=sep)
        colnames = self.columnNames
        collabels = self.columnlabels
        recs = self.getAllCells()
        row=[]
        for c in colnames:
            row.append(collabels[c])
        writer.writerow(row)
        for row in recs.keys():            
            writer.writerow(recs[row])            
        return
        
    def __repr__(self):
        return 'Table Model with %s rows' %len(self.reclist)
