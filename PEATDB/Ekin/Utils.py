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

"""General Utility functions, mainly for Ekin related work."""

import math

def meanstdv(x):
     from math import sqrt
     n, mean, std = len(x), 0, 0
     for a in x:
         mean = mean + a
     mean = mean / float(n)
     for a in x:
         std = std + (a - mean)**2
     std = sqrt(std / float(n-1))
     return mean, std

def leadingZeros(value, desired_digits):
    """
    Given an integer, returns a string representation, padded with [desired_digits] zeros.
    http://www.djangosnippets.org/snippets/543/
    """
    num_zeros = int(desired_digits) - len(str(value))
    padded_value = []
    while num_zeros >= 1:
        padded_value.append("0")
        num_zeros = num_zeros - 1
    padded_value.append(str(value))
    return "".join(padded_value)

def getDictMax(a):
    """Get key of max element in a dict of ints or floats"""
    b = dict(map(lambda item: (item[1],item[0]),a.items()))
    max_key = b[max(b.keys())]
    return max_key

def getDictSum(a):
    """Get key of max element in a dict of ints or floats"""
    sum=0
    for i in a: sum+=a[i]
    return sum

def get_bins(data, bins):
    """Create bins from 2 or more sets of lists of unknown ranges,
       so they have the same bin scale when graphed together"""
    import numpy, types
    #print data
    if type(data) is types.ListType:
        mnd = min(data[0])
        mxd = max(data[0])
        for d in data:
            #print min(d), max(d)
            if mnd > min(d):
                mnd = min(d)
            if mxd < max(d):
                mxd = max(d)
    elif type(data) is types.DictType:
        fstkey = data.keys()[0]
        mnd = min(data[fstkey])
        mxd = max(data[fstkey])
        for d in data.keys():
            #print min(data[d]), max(data[d])
            if mnd > min(data[d]):
                mnd = min(data[d])
            if mxd < max(data[d]):
                mxd = max(data[d])

    inc = (mxd-mnd)/bins
    try:
        newbins = numpy.arange(mnd, mxd, inc)
    except:
        newbins = None
    #print newbins
    return newbins

def adjust_bars(patches):
    """Adjust overlapping bars in a pylab plot so that smaller ones
        are in front"""

    l = len(patches[0])
    for i in range(0,l):
        rects=[]
        heights=[]
        for p in patches:
            h = patches[p][i].get_height()
            rects.append((h, patches[p][i]))

        rects.sort()
        i=len(rects)
        for r in rects:
            r[1].set_zorder(i)
            i-=1
    return

def removeEmptyKeys(data):
    """Remove keys from a dict that have no data or zero len lists"""
    import types
    for i in data.keys():
        if type(data[i]) is types.ListType:
            if len(data[i])==0:
                del data[i]

    return data

def getlenDict(data, keys=None):
    """Get len of a dict of lists"""
    l = 0
    if keys == None:
        keys = data.keys()
    for i in keys:
        l += len(data[i])
    return l


def doFtest(model1, model2, numparams1, numparams2, numdps, alpha=0.05):
    """Performs an F test between 2 fitted models.
       model1 represents the simpler null hypothesis model
       the null hypothesis will be that the simpler model fits better
       ss1 : sum of squared error 1
       numdps: number of datapoints """

    def fdist(F,df_n,df_d):

        a = df_n/2.
        b = df_d/2.
        x = df_n*F/float(df_n*F+df_d)
        p = 1 - betai(a,b,x,)
        return p

    ss1 = float(model1['error'])
    ss2 = float(model2['error'])
    if ss1==0 or ss2==0:
        return 0,0
    if ss2>=ss1:
        print 'error on new model is worse - no need for f-test..'
        return 0, 0

    df1 = numdps - numparams1
    df2 = numdps - numparams2
    if df2 <= 0:
        return 0, 0
    #if same degrees of freedom then just compare square diffs and return
    if df1 == df2:
        if ss2>=ss1:
            return 0, 0
        else:
            return 1, 0
    
    #get F value
    F = ((ss1-ss2) / (df1-df2)) / (ss2/df2)    
    #find P value   
    P = fdist(F, df1-df2, df2)
    print 'ss1', 'ss2', 'numparams1', 'numparams2', 'df1', 'df2'
    print ss1, ss2, numparams1, numparams2, df1, df2
    print 'F-ratio:', F, 'P:', P
    #if P value is less than alpha than the null hypothesis is rejected,
    #model2 actually does fit the data better not just due to the extra param
    #we return 1
    if P <= alpha:
        return 1, P
    else:
        return 0, P


def betai(a,b,x):

    ## Numerical Recipes (www.nr.com)
    ## incomplete beta function

    if x < 0. or x > 1.:
        stop
    if x == 0. or x == 1.:
        bt = 0.0
    else:
        bt = math.exp(
            gammln(a+b)-gammln(a)-gammln(b)+a*math.log(x)+b*math.log(1.0-x)
            )

    if x < (a+1.)/(a+b+2.):
        f = bt*betacf(a,b,x)/float(a)
    else:
        f = 1.-bt*betacf(b,a,1.-x)/float(b)

    return f


def gammln(xx):

    ## Numerical Recipes (www.nr.com)
    ## gamma function

    coeff = [76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.120858003e-2, -0.536382e-5,]
    x = xx - 1.0
    tmp = x + 5.5
    tmp = tmp - (x+0.5)*math.log(tmp)
    ser = 1.0
    for j in range(len(coeff)):
        x += 1
        ser += coeff[j]/x

    return -tmp + math.log(2.50662827465*ser)


def betacf(a,b,x):

    ## Numerical Recipes (www.nr.com)

    ITMAX = 200
    EPS = 3.0e-7

    bm = az = am = 1.
    qab = a+b
    qap = a+1.
    qam = a-1.
    bz = 1.-qab*x/qap
    for i in range(ITMAX+1):
        em = float(i+1)
        tem = em + em
        d = em*(b-em)*x/((qam+tem)*(a+tem))
        ap = az + d*am
        bp = bz+d*bm
        d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem))
        app = ap+d*az
        bpp = bp+d*bz
        aold = az
        am = ap/bpp
        bm = bp/bpp
        az = app/bpp
        bz = 1.
        if (abs(az-aold)<(EPS*abs(az))):
            break

    if i == ITMAX:
        print 'a or b too big, or ITMAX too small in Betacf.'
        stop

    return az

def prettyPrint(data):

    import pprint
    pp = pprint.PrettyPrinter(indent=4)
    x=pp.pformat(data)
    return x

