#!/usr/bin/env python
#
# DataPipeline - A data import and fitting tool
# Copyright (C) 2011 Damien Farrell
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
# Email: damien.farrell_at_ucd.ie
# Normal mail:
# Damien Farrell
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
#

import string, time, datetime
import Base, Utilities

class BaseImporter(object):
    """Base Importer class, sub-class this to define methods specific to each kind of
       import format. At minimum we override the doImport method to get specific
       format functionality. This class should not be instantiated directly."""

    name = 'base'
    def __init__(self, cp):
        """Arguments:
            cp - a ConfigParser object that has been loaded in the parent app"""
        Utilities.setAttributesfromConfigParser(self, cp)
        if self.delimeter=='': self.delimeter=' '
        elif self.delimeter=='tab': self.delimeter='\t'
        return

    def guessRowStart(self, lines):
        """If rowstart is not specified in config, it might be non-zero"""
        s = 0
        for line in lines:
            if not self.checkEmptyRow(line):
                s+=1
            else:
                self.rowstart = s
                return

    def checkEmptyRow(self, line):
        if line.startswith('#') or line=='' or line.startswith('\r'):
            return False
        return True

    def getRow(self, lines, row, start=None, grouped=False):
        """Return values in a row"""
        #if we have repeating cols we return multiple lists
        if self.ignorecomments==True and lines[row].startswith('#'):
            return None
        if not self.checkEmptyRow(lines[row]):
            return None
        if start == None:
            start = self.colstart
        vals = string.strip(lines[row]).split(self.delimeter)
        vals = vals[start:]
        if grouped == False:
            return vals
        else:
            if self.colrepeat == 0:
                return [vals]
            elif self.colrepeat > 1:
                return self.groupList(self.colrepeat, vals)
            else:
                return None

    def getColumn(self, lines, col, grouped=False):
        """Return values in a column"""
        vals = []
        for row in range(self.rowstart, self.rowend):
            if self.ignorecomments==True and lines[row].startswith('#'):
                continue
            rowdata = string.strip(lines[row]).split(self.delimeter)
            vals.append(rowdata[col])
        if grouped == False:
            return vals
        else:
            if self.rowrepeat == 0:
                return [vals]
            elif self.rowrepeat > 1:
                return self.groupList(self.rowrepeat, vals)

    def getColumnHeader(self, lines, grouped=False):
        """Column headers are taken from relevant row"""
        if self.colheader == '':
            row = self.rowstart
        else:
            row = self.colheader
        if self.colheaderstart != '':
            start = self.colheaderstart
        else:
            start = self.colstart
        return self.getRow(lines, row, start=start, grouped=grouped)

    def getRowHeader(self, lines, grouped=False):
        """Return values in header row"""
        if self.colheader == '':
            col = self.colstart
        else:
            col = self.rowheader
        return self.getColumn(lines, col, grouped)

    def groupList(self, n, l, padvalue=None):
        """group a list into chunks of n size"""
        return [l[i:i+n] for i in range(0, len(l), n)]

    def getXYValues(self, xd, yd, start=0, xformat='', yformat=''):
        """Return a lists of floats from lists of vals whilst doing
           various checks to remove errors etc."""

        x=[]; y=[]
        for i in range(start,len(xd)):
            if i>=len(yd): break
            xval = self.checkValue(xd[i], xformat)
            yval = self.checkValue(yd[i], yformat)
            if xval==None or yval==None:
                continue
            x.append(xval)
            y.append(yval)
        if len(x)<=1:
            return None, None
        return x,y

    def checkValue(self, val, format=''):
        """Coerce a string to float if possible"""
        #add code to handle commas in thousand separators
        dec = self.decimalsymbol

        if format != '':
            return self.checkTime(val, format)
        if dec == '.':
            try:
                return float(val)
            except:
                return None
        else:
            try:
                return float(val.replace(".","").replace(dec,"."))
            except ValueError, e:
                #print e
                return None

    def checkTime(self, val, timeformat):
        """Coerce to a datetime object and return a value in seconds
           from some reference"""
        try:
            dt = datetime.datetime.strptime(val,timeformat)
            return dt
        except ValueError, e:
            print e
            return None

    def convertTimeValues(self, vals):
        """Convert datetime values to decimal"""
        if vals == None or len(vals)<1:
            return None
        result=[]
        ref = vals[0]
        for v in vals:
            delta = v-ref
            result.append(delta.seconds)
        return result

    def checkUnicode(self, s):
        """Check for unicode string"""
        try:
            s.decode('ascii')
        except UnicodeDecodeError:
            s = unicode(s)
        return s

    def doImport(self, lines):
        """Should be overrriden"""
        return

class DatabyColImporter(BaseImporter):
    """This importer handles data formatted in columns with common x values
       in a specified column, it also handles multiple sets of data grouped
       in evenly spaced distances down each column"""

    name = 'databycolumn'

    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):

        data = {}
        if self.rowend == 0:
            self.rowend=len(lines)
        xdata = self.getRowHeader(lines,grouped=True)
        header = self.getColumnHeader(lines)
        if self.colend == 0:
            self.colend = len(header)
        if xdata == None:
            return

        for col in range(self.colstart+1, self.colend):
            coldata = self.getColumn(lines, col, grouped=True)
            #print xdata, coldata
            if coldata == None: continue
            for xd,yd in zip(xdata,coldata):
                if len(xd)<=1 or len(yd)<=1: continue
                name = yd[0]
                x,y = self.getXYValues(xd[1:],yd[1:])
                data[name] = [x,y]
        return data

class DatabyRowImporter(BaseImporter):
    """This importer handles data formatted in rows with common x values
       along the top (or specified in a specific row), it also handles
       multiple sets of data grouped in evenly spaced distances along the row"""

    name = 'databyrow'

    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):
        
        data = {}        
        self.guessRowStart(lines)
        xdata = self.getColumnHeader(lines,grouped=True)
        if xdata == None:
            return
        if self.rowend == 0:
            self.rowend=len(lines)
            
        for row in range(self.rowstart+1, self.rowend):
            if row>=len(lines):
                break
            rowdata = self.getRow(lines,row,grouped=True)            
            #print xdata, rowdata
            if rowdata == None: continue
            for xd,yd in zip(xdata,rowdata):
                #print xd, yd
                if len(xd)<=1 or len(yd)<=1: continue
                name = yd[0]
                x,y = self.getXYValues(xd[1:],yd[1:])
                data[name] = [x,y]
        return data

class PairedDatabyColImporter(BaseImporter):
    """This importer handles data formatted in rows with paired x-y values,
       there are therefore no common x-values in the header """

    name = 'paireddatabycolumn'

    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):

        data = {}
        if self.rowend == 0:
            self.rowend=len(lines)
        header = self.getColumnHeader(lines)
        if self.colend == 0:
            self.colend = len(header)+1

        i=0
        for col in range(self.colstart, self.colend):
            coldata = self.getColumn(lines, col, grouped=True)
            name = header[i]
            for xyd in coldata:
                x = xyd[1:len(xyd):2]
                y = xyd[2:len(xyd):2]
                data[name] = [x,y]
            i+=1
        return data

class PairedDatabyDoubleColImporter(BaseImporter):

    name = 'paireddatabydoublecolumn'

    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):

        data = {}
        if self.rowend == 0:
            self.rowend=len(lines)
        header = self.getColumnHeader(lines)
        if self.colend == 0:
            self.colend = len(header)

        i=0
        for col in range(self.colstart, self.colend, 2):
            x = self.getColumn(lines, col, grouped=True)[0]
            y = self.getColumn(lines, col+1, grouped=True)[0]
            name = header[i]
            data[name] = [x,y]
            i+=2
        return data

class PairedDatabyRowImporter(BaseImporter):
    """This importer handles data formatted in rows with paired x-y values,
       there are therefore no common x-values in the header """

    name = 'paireddatabyrow'

    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):

        data = {}
        if self.rowend == 0:
            self.rowend=len(lines)

        for row in range(self.rowstart+1, self.rowend):
            if row>=len(lines):
                break
            rowdata = self.getRow(lines,row, grouped=True)
            for xyd in rowdata:
                name = xyd[0]
                x = xyd[1:len(xyd):2]
                y = xyd[2:len(xyd):2]
                data[name] = [x,y]

        return data

class PairedDatabyDoubleRowImporter(BaseImporter):
    """This importer handles data formatted in rows with paired x-y values,
       there are therefore no common x-values in the header """

    name = 'paireddatabydoublerow'

    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):

        data = {}
        if self.rowend == 0:
            self.rowend=len(lines)
        header = self.getRowHeader(lines)
        i=0
        for row in range(self.rowstart, self.rowend, 2):
            if row>=len(lines):
                break
            x = self.getRow(lines,row, grouped=True)[0]
            y = self.getRow(lines,row+1, grouped=True)[0]
            name = header[i]
            data[name] = [x,y]
            i+=2
        return data

class GroupedDatabyRowImporter(BaseImporter):
    """This importer handles data formatted in rows with multiple independent x values in
       each column, each dataset is then repeated in groups every x rows, specified in
       the rowrepeat option. The importer therefore returns dictionary with multiple sets of
       x-y values for each label/dataset"""

    name = 'groupeddatabyrow'

    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):

        data = {}
        #assumes the column header has labels for each set of xy vals
        labels = self.getColumnHeader(lines)
        print labels
        if self.rowend == 0:
            self.rowend=len(lines)
        if self.colend == 0:
            self.colend = len(labels)

        step = self.rowrepeat
        for d in range(1, self.rowrepeat):
            for row in range(self.rowstart, self.rowend, step):
                if row>=len(lines):
                    break
                xdata = self.getRow(lines, row)
                rowdata = self.getRow(lines, row+d)
                name = rowdata[0]
                if not data.has_key(name):
                    data[name] = {}

                for v in zip(labels,xdata,rowdata)[1:]:
                    label = v[0]
                    x = self.checkValue(v[1])
                    y = self.checkValue(v[2])
                    if x==None or y==None:
                        continue
                    if not data[name].has_key(label):
                        data[name][label]=[]
                    l = data[name][label]
                    l.append((x,y))

        #reformat paired vals into x and y lists
        for d in data:
            for lbl in data[d]:
                data[d][lbl] = zip(*data[d][lbl])

        return data

class GroupedDatabyColImporter(BaseImporter):
    """This importer handles data formatted in cols with multiple independent x values in
       each row, each dataset is then repeated in groups every x columns, specified in
       the colrepeat option. The importer therefore returns dictionary with multiple sets of
       x-y values for each label/dataset"""

    name = 'groupeddatabycolumn'

    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):

        data = {}
        #assumes the row header has labels for each set of xy vals

        if self.rowend == 0:
            self.rowend=len(lines)
        labels = self.getRowHeader(lines)
        header  = self.getColumnHeader(lines)
        if self.colend == 0:
            self.colend = len(header)
        grouplen = len(self.getColumnHeader(lines, grouped=True)[0])
        step = self.colrepeat

        for d in range(1,grouplen):
            for col in range(self.colstart, self.colend, step):

                xdata = self.getColumn(lines, col)
                coldata = self.getColumn(lines, col+d)
                name = coldata[0]
                if not data.has_key(name):
                    data[name] = {}
                #print name, xdata,coldata
                for v in zip(labels,xdata,coldata)[1:]:
                    label = v[0]
                    x = self.checkValue(v[1])
                    y = self.checkValue(v[2])
                    if x==None or y==None:
                        continue
                    if not data[name].has_key(label):
                        data[name][label]=[]
                    l = data[name][label]
                    l.append((x,y))

        #reformat paired vals into x and y lists
        #by convention we put secondary labels into nested dicts
        for d in data:
            for lbl in data[d]:
                data[d][lbl] = zip(*data[d][lbl])
        return data
