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

import os, sys, math, string, types
import inspect
from datetime import datetime
import ConfigParser, csv
from Base import BaseImporter

"""Custom Importers should be added here. These will usually sub-classes of
BaseImporter but can also inherit from any of the importers in Base.py
Users should must add an entry to the dictionary below so that the class can
be identified from the format keyword in the config file"""

class KineticsDataImporter(BaseImporter):
    """This is a custom importer to handle the kinetics data supplied as part of the
       case study. This data are kinetic assays measured over time intervals per substrate
       concentration. The data are grouped in rows per time point, each row is represents a
       specific concentration. Each column is one variant.
       The importer returns a nested dictionary with variants as keys"""

    name = 'kineticsdata'   #use as format keyword in conf file

    def __init__(self, cp):
        BaseImporter.__init__(self, cp)
        return

    def doImport(self, lines):
        """Common x values for every substrate concentration"""

        reftime = None
        data = {}
        #assumes the column header has labels for each set of xy vals
        #labels = self.getColumnHeader(lines)
        names = self.colheader.split(',')
        rowlabels = self.rowheader.split(',')
        #print names

        if self.rowend == 0:
            self.rowend=len(lines)-12
        if self.colend == 0:
            self.colend = len(labels)

        #grouplen = (self.rowend - self.rowstart) / self.rowrepeat
        rowstep = self.rowrepeat
        #print self.rowend,self.rowstart,self.rowrepeat
        for col in range(self.colstart, self.colend):
            name = names[col]
            if not data.has_key(name):
                data[name] = {}
            #print name, col
            for d in range(0, rowstep):
                xdata=[]; ydata=[]
                if d < len(rowlabels):
                    label = rowlabels[d]
                else:
                    label = d

                for row in range(self.rowstart, self.rowend, rowstep):
                    xval = self.getRow(lines, row)[0]
                    rowdata = self.getRow(lines, row+d)
                    if len(rowdata) <= col: continue
                    xdata.append(xval)
                    #print d,row+d, len(rowdata), col
                    if d==0: ind=col+2
                    else: ind = col
                    ydata.append(rowdata[ind])

                #print xdata, ydata
                if len(xdata)<=1: continue
                x,y = self.getXYValues(xdata,ydata,xformat='time')
                #if xformat=='time':
                x = self.convertTimeValues(x)
                if x== None or y==None:
                    continue
                #print x,y
                if not data[name].has_key(label):
                    data[name][label]=[x,y]

        return data

    def postProcess(self):
        """Post process raw data"""
        return

