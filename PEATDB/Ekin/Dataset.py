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

from PEATDB.TableModels import TableModel
import numpy as np
import time
import timeit

class EkinDataset(object):
    """Prototype dataset structure for new ekin datasets"""

    names = ['data', 'active', 'errors', 'labels']

    def __init__(self, x=[], y=[], xy=None, active=[],
                 xerrs=[], yerrs=[], 
                 xlabel='x', ylabel='y',                 
                 data=None):
        
        if data == None:
            self.dims = 2
            if xy != None:
                self.x,self.y = xy
            else:    
                self.x=x
                self.y=y
            self.data = [self.x,self.y]                        
            if len(active)==0:
                active = [1 for i in self.x]
            self.active=active
            if xerrs!=None and len(xerrs)==0:
                xerrs=[0 for i in self.x]
            if yerrs!=None and len(yerrs)==0:
                yerrs=[0 for i in self.y]                
            self.errors=[xerrs,yerrs]
            self.labels=[xlabel,ylabel]
            self.fits={}
        else:
            self.__dict__ = data        
        return

    def add(self, dp=(None,None), a=1, e=(None,None)):
        """Add a new datapoint, dp is a tuple"""        
        self.x.append(dp[0])
        self.y.append(dp[1])
        self.active.append(a)
        for d in range(self.dims):
            self.errors[d].append(e[d])
        return

    def remove(self, i=None):
        """Delete a datapoint at index i"""
        if i >  len(self.x):
            print 'index out of range'
            return
        if i==None:
            self.x.pop()
            self.y.pop()
        else:
            del self.x[i]
            del self.y[i]
            del self.active[i]
            for j in range(self.dims):
                del self.errors[j][i]
        return
   
    def removeMultiple(self, lst):
        """Remove points using list of elements"""       
        rng = range(0,self.length())
        self.x=[self.x[i] for i in rng if i not in lst]
        self.y=[self.y[i] for i in rng if i not in lst]
        self.active=[self.active[i] for i in rng if i not in lst]
        for d in range(self.dims):
            self.errors[d] = [self.errors[d][i] for i in rng if i not in lst]            
        return        
            
    def removeBounded(self, bounds):
        """Remove points within selected x-y bounds"""
        if bounds==None or len(bounds)!=4:
            return
        x1,y1,x2,y2 = bounds
        if x1>x2 :
            temp=x1;x1=x2;x2=temp
        if y1>y2:
            temp=y1;y1=y2;y2=temp
        lst=[]
        for i in range(0,self.length()):
            x=self.x[i]; y=self.y[i]
            if (x>x1 and x<x2) and (y>y1 and y<y2):                
                lst.append(i)
        self.removeMultiple(lst)
        return

    def addFit(self, fitdata, name='default'):
        """Add a fit model"""
        if not hasattr(self, 'fits'):
            self.fits={}
        self.fits[name] = fitdata
        return

    def getFit(self):
        """get default fit data"""
        if self.fits.has_key('default'):
            return self.fits['default']
        else:
            return None
    
    def getxy(self):
        """Get x-y lists"""
        xy=zip(self.x,self.y)
        x=[i[0] for i in xy if i[0]!=None and i[1]!=None]
        y=[i[1] for i in xy if i[0]!=None and i[1]!=None]        
        return x,y

    def getxya(self):
        xya=zip(self.x,self.y,self.active)
        x=[i[0] for i in xya if i[0]!=None and i[1]!=None]
        y=[i[1] for i in xya if i[0]!=None and i[1]!=None]
        a=[i[2] for i in xya if i[0]!=None and i[1]!=None]
        return x,y,a
    
    def getAll(self):
	"""Get all data"""
	return self.x,self.y,self.active,self.errors[0],self.errors[1]
    
    def getActive(self):
        """Get only active points"""
        x=[];y=[]
        for i in zip(self.x,self.y,self.active):
            if i[0]!=None and i[1]!=None and i[2] == 1:
                x.append(i[0])
                y.append(i[1])
        return x,y
        
    def setActive(self, i, a=1):
        """Set a point as active"""
        self.active[i] = a
        return

    def setActiveBounded(self, bounds=None, status=1):
        """Set (in)active from bounded points, a tuple x1,y1,x2,y2"""        
        if bounds==None or len(bounds)!=4:
            return
        x1,y1,x2,y2 = bounds
        if x1>x2 :
            temp=x1;x1=x2;x2=temp
        if y1>y2:
            temp=y1;y1=y2;y2=temp
        for i in range(0,self.length()):
            x=self.x[i]; y=self.y[i]
            if (x>x1 and x<x2) and (y>y1 and y<y2):
                self.active[i]= status     
        return

    def adjust(self, column=0, op='+', val=0):
        """Perform some arithmetic on pts"""
        lst=self.data[column]
        for i in range(0,self.length()):            
            lst[i]=eval(str(lst[i]) + op + str(val))      
        return
        
    def setError(self, col, i, e):       
        if e != None:
            self.errors[col][i] = e                
        return

    def getErrors(self):
        xerrs,yerrs=self.errors
        if xerrs != None:
            xerrs=[i[0] for i in zip(xerrs, self.x) if i[0]!=None and i[1]!=None]
        if yerrs != None:    
            yerrs=[i[0] for i in zip(yerrs, self.y) if i[0]!=None and i[1]!=None]
        return xerrs,yerrs
        
    def xval(self, i):
        """Get x value at index i"""
        return self.x[i]

    def yxal(self, i):
        """Get y value at index i"""
        return self.y[i]

    def xerr(self, i):
        """Get x error at index i"""
        return self.errors[0][i]

    def yerr(self, i):
        """Get y error at index i"""
        return self.errors[1][i]

    def xrange(self):
        x = self.x
        return max(x)-min(x)

    def yrange(self):
        y = self.y
        return max(y)-min(y)

    def length(self, getall=False):
        """Get no. of datapoints"""        
        return len(self.x)

    def minX(self):
        """Get average of x datapoints"""
        return min(self.x)

    def maxX(self):
        """Get average of x datapoints"""
        return max(self.x)

    def minY(self):
        """Get average of y datapoints"""
        return min(self.y)

    def maxY(self):
        """Get average of y datapoints"""
        return max(self.y)

    def avgX(self):
        """Get average of x datapoints"""
        return np.mean(self.x)

    def avgY(self):
        """Get average of y datapoints"""
        return np.mean(self.y)

    def getData(self):
        return self.__dict__

    def printAll(self):
        for i in range(self.x):
            print self.x[i], self.y[i]
        return

    def printXY(self):
        """print paired X-Y points"""
        print zip(self.x, self.y)

    def prettyPrint(self):
        """prints the entire dict"""
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        x=pp.pformat(self.__dict__)
        print x
        return

    def __repr__(self):
        return 'dataset with %s points'  %self.length()

    def len(self):
        return len(self.x)

def simpleTest():
    """Do basic tests"""
    m = EkinDataset()
    for x in range(1,11):
        y=x/3.14
        m.add((x,y), a=1, e=(x/10.0,y*0.02))

    print m.xrange(), m.yrange()    
    print m.minX(), m.maxX()
    print m.minY(), m.maxY()
    print m.avgX(), m.avgY()
    print m.xerr(0), m.yerr(0)
    x=m.getData()
    m.setActive(2, 0)
    m.setError(1, 1,0.5)
    m.prettyPrint()   
    return

def main():
    simpleTest() 

if __name__ == '__main__':
    main()
