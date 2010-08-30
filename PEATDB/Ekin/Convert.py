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

"""Class methods for converting ekin data"""

import sys,os
import List_Utils
from Tkinter import *
import re

class EkinConvert(object):

    @classmethod
    def ekin2xy(cls, ekindataset, getall=0, geterrors=False):
        """Convert ekin datapoint storage to seprate x-y lists, can be used
           for plotting in gnuplot or Plotter for example"""
        #newer format uses 2 sub-fields for point, old is just a string
        #active field is ignored in final lists
        import copy
        import re
        s=re.compile('[0-9]')
        xlist=[]
        ylist=[]
        activelist=[]
        xerrors=[]
        yerrors=[]
        #to store which points are inactive
        inactivelist=[]
        if not isinstance(ekindataset, dict):
            return None

        # Get the number of data points from the first key
        ks=ekindataset.keys()
        ks.sort()
        first_key=ks[0]
        dps=len(ekindataset[first_key].keys())

        xdatapoints=ekindataset[0].keys()
        ydatapoints=ekindataset[1].keys()
        xdatapoints.sort()
        ydatapoints.sort()
        count=0
        for datapoint in xdatapoints:
            if datapoint=='label' or datapoint=='label_widget':
                continue
            else:
                if ekindataset[0][datapoint].has_key('var'):
                    if ekindataset[0][datapoint]['var'] != '':
                        if not ekindataset[1].has_key(datapoint):
                            continue

                        x = ekindataset[0][datapoint]['var']
                        y = ekindataset[1][datapoint]['var']
                                    
                        if s.search(ekindataset[0][datapoint]['var']) and s.search(ekindataset[1][datapoint]['var']):
                            xlist.append(float(ekindataset[0][datapoint]['var']))
                            ylist.append(float(ekindataset[1][datapoint]['var']))
                            activelist.append(ekindataset[0][datapoint]['active'])
                    else:
                        if getall==1:
                            #we include even empty points
                            xlist.append(ekindataset[0][datapoint]['var'])
                            ylist.append(ekindataset[1][datapoint]['var'])
                            activelist.append(ekindataset[0][datapoint]['active'])
                        else:
                            continue

                    if geterrors == True:
                        if ekindataset[0][datapoint]['var'] != '' and ekindataset[1][datapoint]['var'] != '':
                            if ekindataset[0][datapoint].has_key('error'):
                                xerrors.append(float(ekindataset[0][datapoint]['error']))
                            else:
                                xerrors.append(float(0))
                            if ekindataset[1][datapoint].has_key('error'):
                                yerrors.append(float(ekindataset[1][datapoint]['error']))
                            else:
                                yerrors.append(float(0))

        if len(xerrors) == 0:
            xerrors=None
        if len(yerrors) == 0:
            yerrors=None
        
        if geterrors == True:
            return xlist, ylist, activelist, xerrors, yerrors
        else:
            return xlist, ylist, activelist

    @classmethod
    def xy2ekin(cls, inputlists, errorlists=None, activelist=None, labels=None):
        """Convert xyz lists to ekin dataset, accepts lists tuple"""
       
        if labels == None:
            labels = ('x','y','z')
        if len(inputlists[0]) != len(inputlists[1]):
            return
        ekindataset = {0:{},1:{}}
        for n in range(len(inputlists)):
            ekindataset[n]['label'] = labels[n]
        for i in range(len(inputlists[0])):
            for n in range(len(inputlists)):
                ekindataset[n][i]={}
                ekindataset[n][i]['var']= str(inputlists[n][i])
                if activelist != None:
                    ekindataset[n][i]['active']=activelist[i]
                else:
                    ekindataset[n][i]['active']=1
                if errorlists!=None and errorlists[n] != None:
                    ekindataset[n][i]['error']=str(errorlists[n][i])
        return ekindataset


    @classmethod
    def check_ekin_dataset(cls, ekindatasets):
        """Check if an ekin dataset has an data points in it at all"""
        for k in ekindatasets.keys():
            ekindataset = ekindatasets[k]
            xdatapoints=ekindataset[0].keys()
            count=0
            for datapoint in xdatapoints:
                if datapoint=='label' or datapoint=='label_widget':
                    continue
                else:
                    if not isinstance(ekindataset[0][datapoint], dict):
                        if ekindataset[0][datapoint] != '':
                            return 1
                    elif ekindataset[0][datapoint].has_key('var'):
                        if ekindataset[0][datapoint]['var']!='':
                            return 1
        return 0

    @classmethod
    def set_error(cls, ekindataset, dp, col):
        """Add empty errors for any points with values"""
        if ekindataset[col][dp]['var'] == '':
            return
        if not ekindataset[col][dp].has_key('error'):
            ekindataset[col][dp]['error'] = 0
        return

    @classmethod
    def add_offset(cls, ekindataset, offset=0, column=0):
        """Add an offset value to one of the ekin data columns"""

        data=ekindataset[0].keys()
        for dp in data:
            if dp =='label' or dp=='label_widget':
                continue
            else:
                if ekindataset[column][dp]['var']!='':
                    p=float(ekindataset[column][dp]['var']) + offset
                    ekindataset[column][dp]['var'] = str(p)
                    print p
        return

    @classmethod
    def getAverageY(cls, ekindataset, offset=0, column=0):
        """Get an average of the y values"""
        x,y = cls.ekin2xy(ekindataset)
        return sum(y) / len(y)

