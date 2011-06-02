#!/usr/bin/env python

import numpy as np
import math, random
import operator, os, sys, csv
import pickle
import pylab as plt
import scipy.optimize

"""Prototype for newer fit class that allows user created
models to be added dynamically and can do multivariate fitting"""

class testdata(object):
    def line(self, noise=2.0):
        x=np.random.normal(1,10,500)
        y=[i+np.random.normal(0,noise) for i in x]
        return x,y

    def simpleHH(self, noise=.01):
        x=np.arange(1,10,0.2)
        pKa=6;span=5;offset=0.2
        y=[]
        for i in x:
            val = span / (1 + 10**(- i + pKa)) + offset
            val += np.random.normal(0,9*noise)
            y.append(val)
        return x,y

    def complexHH(self, noise=.02):
        x=np.arange(1,10,0.2)
        pKa1=3;span1=5;pKa2=7;span2=5;offset=0.6
        y=[]
        for i in x:            
            val = span1/ (1+10**(pKa1-i)) + span2/ (1+10**(-i+pKa2)) + offset            
            val += np.random.normal(0,9*noise)            
            y.append(val)          
        return x,y
    
class fitter(object):
    def __init__(self, func, params, x, y):
        self.params = params
        self.func = func
        self.x = x; self.y = y
        return

    def lstsq(self, x, y):
        """DIY lsq"""
        p=self.params
        rounds=range(60)
        for r in rounds:           
            r = self.evaluate(y,fit)
        self.fit = fit        
        return fit

    def residuals(self, p, args=None):
        """Evaluate the func residuals given parameters"""
        r=[]
        x=self.x; y=self.y
        fit=[self.func(i,p) for i in x]        
        r = [math.pow(i[0]-i[1],2) for i in zip(fit,y)]
        return r
    
    def evaluate(self, p, args=None):
        """Evaluate func and get sum sq res for given params""" 
        x=self.x; y=self.y
        fit=[self.func(i,p) for i in x]
        r=0
        for i in zip(fit,y):
            r += math.pow(i[0]-i[1],2)
        return r
    
    def minimize(self):
        return

    def fit(self, method='simplex'):
        """Fit by minimizing r-squared using various algorithms"""
        #downhill simplex algorithm
        if method == 'simplex':
            p = scipy.optimize.fmin(self.evaluate, self.params)
        #using scipy version of levenberg-Marquardt algorithm  
        elif method == 'lm':   
            p,ier = scipy.optimize.leastsq(self.residuals, self.params)
        self.params = p
        fit=[self.func(i,p) for i in self.x]
        self.fit = fit
        return fit       

    def plot(self, ax=None):
        x=self.x; y=self.y
        fit = self.fit
        if ax==None:
            fig=plt.figure(figsize=(6,6))
            ax=fig.add_subplot(111)
            self.fig = fig
        ax.plot(x, y,'o',alpha=0.6)
        inc = abs(max(x)-min(x))/30
        fitx = np.arange(min(x)-inc,max(x)+inc,inc)       
        fity = [self.func(i,self.params) for i in fitx]
        ax.plot(fitx, fity,lw=3,alpha=0.7)
        #ax.set_title(self.params)
        ax.text(0.1,0.8,self.params,fontsize=0.8)        
        return ax

    def estimateUncertainty(self,x,y,p,xerr=0.1,yerr=0.1,runs=10):
        """Generic version of monte carlo parameter uncert, returns
           st dev for each parameter over repeated runs"""
        plist=[]
        for r in range(runs):
            mutx=[];muty=[]
            for i in range(len(x)):
                mutx.append(x[i] + random.uniform(-xerr, xerr))
                muty.append(x[i] + random.uniform(-yerr, yerr))
            F=fitter(self.func,p,mutx,muty)
            F.fit()
            plist.append(F.params)
        result = []    
        for i in range(len(p)):           
            result.append(np.std([v[i] for v in plist]))
        return result
    
    
class fitModel(object):
    """Models created dynamically should use this to inherit from"""
    def __init__(self):
        return

    def guessStart(self):
        return
        
def linear(x,p):
    m,b=p    
    y = m * x + b
    return y

def hh1pka(x,p):
    pKa,span,offset=p
    y = span / (1 + 10**(- x + pKa)) + offset
    return y

def hh2pka(x,p):
    pKa1,span1,pKa2,span2,offset=p
    y = span1/ (1+10**(pKa1-x)) + span2/ (1+10**(-x+pKa2)) + offset 
    return y

def sigmoid(x,p):
    t,bottom,top,slope=p
    y = bottom + (top - bottom) / (1 + math.exp((t-x)/slope))
    return y

def depletion(x, p):
    M,D,x0=p
    y=M * (1 - math.exp(-D*(x-x0)))
    return y
    
def michaelismenten(x,p):
    so,vmax,km=p
    y = vmax*(s0/(km+x))
    return y
    
def test():
    T=testdata()
    x,y=T.line()
    #F=fitter(linear,[0.5,1],x,y)
    x,y=T.simpleHH()
    #x,y=T.complexHH()
    F=fitter(hh1pka,[1,1,1],x,y)
    #F=fitter(sigmoi[1,1,1]d,[6,0,1,1],x,y)
    F.fit()
    F.plot()
    F.estimateUncertainty(x,y,[1,1,1])
    
def test10R():
    """pKa fitting from kcats using substr depletion"""
    path = 'fergal_10R'    
    folders = ['fergal_10R/10RWT','fergal_10R/U33W1']
    pkas=[]
    for path in folders:
        fig=plt.figure(figsize=(8,8))
        i=1
        data = []
        ax1=None      
        for f in os.listdir(path):       
            if os.path.splitext(f)[1] != '.csv': continue
            cr = csv.reader(open(os.path.join(path,f),'r'))
            ph=float(f.split(' ')[1])        
            cols = len(cr.next())-1
            print path, f, ph, '%s cols' %cols
            vals = [r for r in cr]
            #may be several replicates
            for c in range(0,cols,2):
                x = [float(r[c]) for r in vals]
                y = [float(r[c+1]) for r in vals]
            
                #fit
                M = max(y)
                F=fitter(depletion,[M,1,1],x,y)
                F.fit()
                D=F.params[1]
                print 'D',D
                if ph==9.0 and D>6: continue
                data.append((ph,D))
                if c==0:
                    ax=fig.add_subplot(4,4,i,sharey=ax1)
                    i+=1
                if ax1==None: ax1=ax
                F.plot(ax)
                ax.set_title(ph)
            
        #fit pKa
        fig.subplots_adjust(wspace=0.4,hspace=0.4)
        x,y=zip(*data)
        F=fitter(hh1pka,[5,2,0],x,y)
        F.fit()
        pkas.append(F.params[0])
        F.plot()
        #res = F.estimateUncertainty(x,y,[5,2,0],xerr=0.1,yerr=0.2,runs=10)        
        pickle.dump(data,open(os.path.basename(path)+'.pickle','w'))
    print pkas    
    return
    
def parametersTest():
    data = pickle.load(open('10RWT.pickle','r'))
    x,y=zip(*data)
    crossValidate(x,y)
    return
    
def crossValidate(x,y, frac=0.2, num=None):
    """Random sub-sampling removal of points to test effects on 
       fit parameters"""
    l=len(x)
    if num==None:
        num = int(l*(1-frac))
    print 'using %s out of %s points..' %(num,l)
    fig=plt.figure(figsize=(8,8))
    c=0
    pkas=[]
    for n in range(20):
        n1 = random.sample(range(l), num)
        x1 = [x[i] for i in range(l) if i in n1]
        y1 = [y[i] for i in range(l) if i in n1]
        F=fitter(hh1pka,[5,2,0],x1,y1)
        F.fit()
        pka = round(F.params[0],3); pkas.append(pka)
        ax=fig.add_subplot(4,5,c)
        F.plot(ax)
        ax.set_title(pka)
        c+=1
    print 'stdev:', np.std(pkas)    
    return
    
def pltconf():
    #plt.rc('font',family='serif')
    plt.rc('font',size=10)
    plt.rc('legend',fontsize=10)
    #plt.rc('text',usetex=True)
    plt.rc('savefig',dpi=300)


if __name__ == '__main__':
    #test()    
    pltconf()
    test10R()
    #parametersTest()
    plt.show()
