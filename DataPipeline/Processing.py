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

import os, sys, copy
import math, random, string, types
import numpy as np
import ConfigParser
import Utilities

class Processor(object):
    """Class that defining pre-processing steps applied to data"""

    def __init__(self, conffile=None):
        self.predefined = ['differentiate','detectpeaks','gaussiansmooth','smooth',
                           'fouriernoisefilter','savitzkygolayfilter',
                           'baselinecorrection','removeoutliers','normalise',
                           'omitrange']
        if conffile==None:
            conffile = os.path.join(os.path.expanduser('~'),'functions.conf')
        if not os.path.exists(conffile):
            self.createConfig(conffile)
        self.parseConfig(conffile)
        return

    def createConfig(self, conffile='functions.conf', **kwargs):
        """Create a separate config file to store """

        c = ConfigParser.ConfigParser()        
        defaults = {'smooth': [('window', 'hanning'),('windowlen', 11)],
                    'gaussiansmooth':[('degree',6)],
                    'baselinecorrection':[('segmentsize',0.3)],
                    'detectpeaks':[('lookahead',5),('delta',20)],
                    'savitzkygolayfilter':[('windowsize',11), ('order',3), ('deriv',0)],
                    'normalise':[('value',1)],
                    'removeoutliers':[('percentile',0.95)],
                    'omitrange':[('xmin',''),('xmax',''),('ymin',''),('ymax','')]}
        cp = Utilities.createConfigParserfromDict(defaults, self.predefined ,**kwargs)
        cp.write(open(conffile,'w'))
        print conffile
        return

    def parseConfig(self, conffile):
        self.cp = ConfigParser.ConfigParser()
        try:
            self.cp.read(conffile)
        except Exception,e:
            print 'failed to read config file! check format'
            print 'Error returned:', e
            return
        return

    def doFunctions(self, names, data):
        """Apply one or more functions to the data, need to account
           for nesting of data here
           names: names of functions, should be in class.predefined
           data: dictionary with paired x,y tuples """

        newdata = copy.deepcopy(data)
        for name in names:
            func = getattr(self, name)
            args = self.getArgs(name)
            #print name, args
            for d in data:               
               x,y = newdata[d]
               newdata[d] = func(x,y, **args)               
               #print d, newdata[d][0]
        return newdata

    def getArgs(self, name):
        """Try to get function arguments from conffile"""
        args={}
        cp = self.cp
        if not name in self.predefined:
            return {}
        if name in cp.sections():
            for i in cp.options(name):
                try: 
                    args[i]=eval(cp.get(name, i))
                except:
                    args[i]=cp.get(name, i)
        return args

    def detectpeaks(self, x_axis, y_axis, lookahead=5, delta=20):
        """
        Converted from/based on a MATLAB script at http://billauer.co.il/peakdet.html
        
        Algorithm for detecting local maximas and minmias in a signal.
        Discovers peaks by searching for values which are surrounded by lower
        or larger values for maximas and minimas respectively
        
        keyword arguments:
        y_axis -- A list containg the signal over which to find peaks
        x_axis -- A x-axis whose values correspond to the 'y_axis' list and is used
            in the return to specify the postion of the peaks. If omitted the index
            of the y_axis is used. (default: None)
        lookahead -- (optional) distance to look ahead from a peak candidate to
            determine if it is the actual peak (default: 500) 
            '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
        delta -- (optional) this specifies a minimum difference between a peak and
            the following points, before a peak may be considered a peak. Useful
            to hinder the algorithm from picking up false peaks towards to end of
            the signal. To work well delta should be set to 'delta >= RMSnoise * 5'.
            (default: 0)
                Delta function causes a 20% decrease in speed, when omitted
                Correctly used it can double the speed of the algorithm
        
        return -- two lists [maxtab, mintab] containing the positive and negative
            peaks respectively. Each cell of the lists contains a tupple of:
            (position, peak_value) 
            to get the average peak value do 'np.mean(maxtab, 0)[1]' on the results
        """
        maxtab = []
        mintab = []
        dump = []   #Used to pop the first hit which always if false
           
        length = len(y_axis)
        if x_axis is None:
            x_axis = range(length)
        
        #perform some checks
        if length != len(x_axis):
            raise ValueError, "Input vectors y_axis and x_axis must have same length"
        if lookahead < 1:
            raise ValueError, "Lookahead must be above '1' in value"
        if not (np.isscalar(delta) and delta >= 0):
            raise ValueError, "delta must be a positive number"
        
        #needs to be a numpy array
        y_axis = np.asarray(y_axis)
        
        #maxima and minima candidates are temporarily stored in
        #mx and mn respectively
        mn, mx = np.Inf, -np.Inf
        
        #Only detect peak if there is 'lookahead' amount of points after it
        for index, (x, y) in enumerate(zip(x_axis[:-lookahead], y_axis[:-lookahead])):
            if y > mx:
                mx = y
                mxpos = x
            if y < mn:
                mn = y
                mnpos = x
            
            ####look for max####
            if y < mx-delta and mx != np.Inf:
                #Maxima peak candidate found
                #look ahead in signal to ensure that this is a peak and not jitter
                if y_axis[index:index+lookahead].max() < mx:
                    maxtab.append((mxpos, mx))
                    dump.append(True)
                    #set algorithm to only find minima now
                    mx = np.Inf
                    mn = np.Inf
            
            ####look for min####
            if y > mn+delta and mn != -np.Inf:
                #Minima peak candidate found 
                #look ahead in signal to ensure that this is a peak and not jitter
                if y_axis[index:index+lookahead].min() > mn:
                    mintab.append((mnpos, mn))
                    dump.append(False)
                    #set algorithm to only find maxima now
                    mn = -np.Inf
                    mx = -np.Inf
            
        #Remove the false hit on the first value of the y_axis
        try:
            if dump[0]:
                maxtab.pop(0)               
            else:
                mintab.pop(0)            
            del dump
        except IndexError:
            #no peaks were found, should the function return empty lists?
            return [], []
       
        x,y = zip(*maxtab)              
        return x, y       

    def baselinecorrection(self, x, y, segmentsize=0.3):
        """Detect and remove baseline by dividing list into segments and
           then interpolating the result over all x values"""        
        seg = int(len(x) * segmentsize)
        bx=[];by=[]
        for i in range(0,len(x),seg):
            mean = np.min(y[i:i+seg])
            bx.append(i)
            by.append(mean)
        by = np.interp(x, bx, by)
        y=y-by
        return x, y

    def differentiate(self, x, y):
        """Get numerical derivative of y vales"""
        dy = list(np.diff(y,1))
        dx = list(x[:len(dy)])
        return dx,dy

    def gaussiansmooth(self, x, y, degree=6):
        """Smooth noisy data with gaussian filter"""

        window=degree*2-1
        weight=np.array([1.0]*window)
        weightGauss=[]
        for i in range(window):
            i=i-degree+1
            frac=i/float(window)
            gauss=1/(np.exp((4*(frac))**2))
            weightGauss.append(gauss)
        weight=np.array(weightGauss)*weight
        smoothed=[0.0]*(len(y)-window)
        for i in range(len(smoothed)):
            smoothed[i]=sum(np.array(y[i:i+window])*weight)/sum(weight)
        #cut x vals so they match y - may introduce some error
        sx = self.getMatchingList(smoothed, x)
        #print d,len(sx), len(smoothed)
        return sx, smoothed

    def getMatchingList(self, a, b):
        d=len(b)-len(a)
        if d>0:
            end=start=int(d/2)
            if start%2!=0: end=end+1
            x = b[start:-end]
        return x

    def smooth(self,x,y,windowlen=11,window='hanning'):
        """smooth the data using a window with requested size. 
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
        """

        if len(x) < windowlen:
            raise ValueError, "Input list needs to be bigger than window size."
        if windowlen<3:
            return y
        if not window in ['hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError, 'wrong filter name'
        s=np.r_[y[windowlen-1:0:-1],y,y[-1:-windowlen:-1]]        
        w=eval('np.'+window+'(windowlen)')
        sy = np.convolve(w/w.sum(),s,mode='same')
        sy = list(sy)
        sx = np.r_[x[windowlen-1:0:-1],x,x[-1:-windowlen:-1]]
        print len(sx), len(sy)
        return sx[windowlen:-windowlen+1],sy[windowlen:-windowlen+1]

    def fouriernoisefilter(self, x, y):
        fft=np.fft.fft(y)
        bp=fft[:]
        for i in range(len(bp)):
            if i<=20:bp[i]=0  
        ibp=np.fft.ifft(bp)
        return x, ibp

    def savitzkygolayfilter(self, x, y, windowsize=11, order=3, deriv=0):
        """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
        The Savitzky-Golay filter removes high frequency noise from data.
        It has the advantage of preserving the original shape and
        features of the signal better than other types of filtering
        approaches, such as moving averages techhniques.
        The main idea behind this
        approach is to make for each point a least-square fit with a
        polynomial of high order over a odd-sized window centered at
        the point.
        """

        try:
            windowsize = np.abs(np.int(windowsize))
            order = np.abs(np.int(order))
        except ValueError, msg:
            raise ValueError("window_size and order have to be of type int")
        if windowsize % 2 != 1 or windowsize < 1:
            raise TypeError("window_size size must be a positive odd number")
        if windowsize < order + 2:
            raise TypeError("window_size is too small for the polynomials order")
        order_range = range(order+1)
        half_window = (windowsize -1) // 2
        # precompute coefficients
        b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        m = np.linalg.pinv(b).A[deriv]
        # pad the signal at the extremes with
        # values taken from the signal itself
        y=np.array(y)
        firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
        lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
        y = np.concatenate((firstvals, y, lastvals))
        return x, np.convolve( m, y, mode='valid')

    def fouriertransform(self, x, y):
        """Find discrete fourier transform of y data"""
        return np.fft(x,y)

    def percentile(self, x, y, percent=0.95):
        """Find the percentile of a list of y values.
        
        @parameter percent - a float value from 0.0 to 1.0.
        @return - the percentile of the values
        adapted from http://code.activestate.com/recipes/511478/ """

        key=lambda x:x
        y = sorted(y)
        if not y:
            return None
        k = (len(y)-1) * percent
        f = math.floor(k)
        c = math.ceil(k)
        if f == c:
            return key(y[int(k)])
        d0 = key(y[int(f)]) * (c-k)
        d1 = key(y[int(c)]) * (k-f)        
        return d0+d1

    def removeoutliers(self, x, y, percentile=0.95):
        """Remove outliers based on nth percentile"""
        p = self.percentile(x,y,percentile)
        print p
        rx=[]; ry=[]
        for i,j in zip(x,y):
            print j,p
            if j<p:
                rx.append(i)
                ry.append(j)
        return rx,ry

    def normalise(self, x, y, value=1):
        """Normalise data"""
        sum = reduce(lambda i,j:i+j, y)
        y = [i/(sum*1.0)*value for i in y]
        return x,y

    def omitrange(self, x, y, xmin='',xmax='', ymin='',ymax=''):
        """Omit range of data points by given range"""
        rx=[]; ry=[]
        if xmin=='':
            xmin=min(x)
        if xmax=='':
            xmax=max(x)
        if ymin=='':
            ymin=min(y)
        if ymax=='':
            ymax=max(y)
        for i,j in zip(x,y)[:]:
            print i,j,xmin,xmax,ymin,ymax         
            if i>xmin and i<xmax and j>ymin and j<ymax:                
                rx.append(i)
                ry.append(j)
        return rx,ry

