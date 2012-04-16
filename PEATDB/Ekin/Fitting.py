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

from LM_Fitter import *
import os, numpy, math, random, types
import Utils, pickle
  
conv=1e-9
grad=1e-9
fitterClasses = {}

def loadModelsFile(modelsfile=None):
    """Update the current set of available models from a dict"""
    if modelsfile == None:
        path = os.path.abspath(os.path.dirname(__file__))
        modelsfile = os.path.join(path,'models.dict')
    if not os.path.exists(modelsfile):
        path = os.getcwd()
        modelsfile = os.path.join(path,'models.dict')
    if not os.path.exists(modelsfile):
        #write a default file if none available
        modelsdict = createModels()
        writeModelsFile(modelsfile, modelsdict)
    currentmodels = pickle.load(open(modelsfile,'r'))
    return currentmodels, modelsfile

def updateModels(data):
    """Reload the current global set of models from a dict"""    
    global currentmodels
    currentmodels = data
    createFitters()
    return
    
def createModels():
    """Create default models file"""
    modelsdict = {
        'Linear': { 'name': 'Linear',
         'description': 'straight line',
         'equation': 'a*x+b',
         'guess': {'a': 'min(y)', 'b': 'max(y)-min(y)/max(x)-min(x)'},         
         'varnames': 'a,b'},
        'Sigmoid': { 'name': 'Sigmoid',
         'description': 'simple sigmoid',
         'equation': 'bottom+(top-bottom)/(1+exp((tm-x)/slope))',
         'guess': {'slope': '1', 'top': 'max(y)', 'tm': '(max(x)-min(x))/2+min(x)', 'bottom': 'min(y)'},
         'varnames': 'tm,bottom,top,slope'},
        'Gaussian': { 'name': 'Gaussian',
         'description': 'gaussian function, a bell-shaped curve',
         'equation': 'a*exp(-(pow((x-b),2)/(pow(2*c,2))))',
         'guess': {'a':'max(y)'},         
         'varnames': 'a,b,c'},
        'Michaelis-Menten': { 'name': 'Michaelis-Menten',
         'description': 'gaussian function, a bell-shaped curve',
         'equation': '(Vmax * x)/(x + Km)',
         'guess': {'Km':'max(x)/100','Vmax':'max(y)'},         
         'varnames': 'Km,Vmax'},         
    }
    return modelsdict
     
def writeModelsFile(modelsfile, modelsdict):  
    """Write a default file to directory"""  
    pickle.dump(modelsdict, open(modelsfile,'w'))
    return

def getCurrentModels():
    """Returns a list of the current models loaded in the module"""
    return sorted(currentmodels.keys())
    
def createClass(equation, varnames,
                   guess=None, name=None, 
                   labels=['x','y'], **kwargs):                       
    """Dynamically create an instance of LM_fitter class using presets
       this function requires we pass varnames and equation at minimum"""
    class tempfitter(LM_Fitter):
        def __init__(self, variables, exp_data, callback, name=None):
            LM_Fitter.__init__(self,variables,exp_data,callback, name=None)
            
            self.name = name
            self.varnames = varnames
            if type(self.varnames) is not types.ListType:
                self.varnames = self.varnames.split(',')
            
            self.labels = labels
            self.equation = equation
            self.guess = guess
            if variables == None:
                self.variables = [1.0 for i in self.varnames]
            else:
                self.variables = variables            
            self.setChangeVars()

        def get_value(self, variables, data_point):
            """Evaluate string equation"""
            try:
                x=data_point[0]
            except:
                x=data_point
            eq=self.equation                
            for i in range(len(self.varnames)):
                #print i, self.varnames[i],variables[i] 
                globals()[self.varnames[i]] = variables[i]
            
            try:   
                value = eval(eq)
            except Exception, e:
                value = 1.0
                print e
            return value

        def guess_start(self):
            """Guess start vals using provided strings in guess dict"""
            if self.guess==None:
                return
            guess = self.guess
            #print guess
            if type(guess) is not types.DictType:
                guess = eval(guess)
            x=[i[0] for i in self.exp_data]
            y=[i[1] for i in self.exp_data]                
            for i in range(len(self.varnames)):
                var = self.varnames[i]                    
                if var in guess:
                    self.variables[i] = eval(guess[var])
                    #print  i, var, guess[var],self.variables[i]
            return
        
    return tempfitter

def createFitters():
    """Create pre-defined fitter classes"""
    for m in currentmodels:
        data = currentmodels[m]
        fitterClasses[m] = createClass(**data)
    #print fitterClasses
    return

def getFitter(model, vrs=None, expdata=None, callback=None):
    """Return the required LM fit object from the model name"""
    if model not in fitterClasses:
        createFitters()
    try:    
        fitclass = fitterClasses[model]
    except:
        print 'model not found, please check name'
        print 'current available models are: %s' %fitterClasses.keys()
        return
    inst = fitclass(variables=vrs,exp_data=expdata,callback=callback)
    return inst

def makeFitter(fitdata, expdata=None, callback=None):
    """Return a fit object created from the provided ekin fit data"""
    
    if fitdata.has_key('model'):
        model = fitdata['model']
        vrs = getFitVars(fitdata)
        if len(vrs)==0:
            return None
    else:
        print 'no model provided to create fitter'
        return None
    return getFitter(model, vrs=vrs, expdata=expdata, callback=callback)

def getFitVars(fitdata):
    """Returns list of fit variables"""
    fitvars = []
    for var in sorted(fitdata.keys()):
        if var != 'model' and var != 'error':
            fitvars.append(fitdata[var])
    return fitvars

def makeFitData(model, vrs=None, error=None):
    """Get preset ekin fit data for a model, can also use provided
       variables and error vals - returns an ekin fitdata item"""
    fitdata={}
    fitdata['model'] = model
    if error == None:
        fitdata['error']=0.0
    else:
        fitdata['error']=error
    i=0
    if vrs == None:
        #get default vals for model
        X = getFitter(model,None,None)
        vrs = X.variables
        for v in vrs:
           fitdata[i]=v
           i+=1
    else:
        for v in vrs:
           fitdata[i]=v
           i+=1
    return fitdata

def getExpData(ek):
    """Get paired active data points for fitter from ekin dataset"""
    #print data
    if ek == None:            
        return None
    x,y = ek.getActive()
    if len(x)==0 or len(y)==0:
        return None
    expdata = zip(x,y)       
    return expdata

def doFit(ekindata=None, fitdata=None, model=None, expdata=None, vrs=None,
                noiter=50, conv=None, grad=None, LM_damper=1.0,
                guess=True, silent=False, callback=None,
                changevars=None, startvalues=None):
    """General class method for doing a fit of provided ekin data and model
       Returns the ekin fit data and the fitter instance used"""

    if conv == None:
        conv = conv
    if grad == None:
        grad = grad

    def updatenow(diff, vrs, fitvals, c, X):
        """Dummy callback for fit update"""
        pass
    if callback == None:
        callback = updatenow
        
    if expdata == None:            
        expdata = getExpData(ekindata)

    if expdata == None or len(expdata)<=1:
        print 'Not enough data to fit'
        return None, None
   
    if fitdata != None and len(fitdata) > 1:
        X = makeFitter(fitdata, expdata=expdata, callback=callback)
    else:
        X = getFitter(model, expdata=expdata, vrs=vrs, callback=callback)
    if X == None:
        print 'No fitter found!'
        return None, None

    if startvalues!=None:
        X.variables = startvalues
    if guess == True:
        X.guess_start()
        if silent == False:
            print 'guessed good start values', X.variables

    # Set the change vars
    X.setChangeVars(changevars)

    # Now do the fitting
    status, variables = X.fit(rounds = noiter, error_crit = conv, gradient_crit = grad,
                              damper_adjustment_factor = 2.0, LM_damper = LM_damper, silent=silent)
    fitdata = makeFitData(model, X.variables, X.getError())

    if silent == False:
        print X.getResult()
    return fitdata, X

def estimateExpUncertainty(ekindata, fitdata, xuncert=0, yuncert=0,
                                 runs=10, doplots=False, callback=None, guess=True,
                                 **kwargs):
    """Chresten's method for estimating the error on each parameter due
       to the effect of the exp uncertainty, which is user provided.
       If the datapoints have x or y error values, we use that instead"""

    '''we iterate over runs and re-fit the data, each time adjusting x-y points
       then accumulate a fit val for each parameter and get mean/stdev.
       Note: errors are added +/- values'''
    
    from PEATDB.Ekin.Dataset import EkinDataset
    grad=1e-6; conv=1e-6
    for k in kwargs:
        if k=='conv':
            conv=kwargs[k]
        if k=='grad':
            grad=kwargs[k]
    if fitdata == None or len(fitdata) == 0:
        return None
    model = fitdata['model']

    F = getFitter(model=model)
    try:
        varnames = F.getVarNames()
    except:
        return None
    fitvars = {}
    for v in varnames:
        fitvars[v] = []     #group lists by variable name
    
    xd,yd,a,xerrs,yerrs = ekindata.getAll()
    estdata={}; estdata['__datatabs_fits__']={}
    for r in range(runs):
        mut_x=[];mut_y=[]
        #if the datapoint has x or y error values, we use that instead
        for i in range(len(xd)):
            if xerrs != None and xerrs[i] > 0:
                mut_x.append(xd[i] + random.uniform(-xerrs[i], xerrs[i]))
            else:
                mut_x.append(xd[i] + random.uniform(-xuncert, xuncert))
            if yerrs != None and yerrs[i] > 0:
                mut_y.append(yd[i] + random.uniform(-yerrs[i], yerrs[i]))
            else:
                mut_y.append(yd[i] + random.uniform(-yuncert,yuncert))
       
        ek = EkinDataset(xy=[mut_x,mut_y],active=a)
        fitres, X = doFit(ek, model=model, fitdata=fitdata, grad=1e-6, conv=1e-6,
                                silent=True, guess=guess, noiter=100)
        vrs = X.getVariables()

        for n in range(len(varnames)):
            fitvars[varnames[n]].append(vrs[n])
        estdata[str(r)] = ek
        estdata['__datatabs_fits__'][str(r)] = fitres

        if callback != None:
            fitstats = {}
            for v in varnames:
                fitstats[v] = numpy.mean(fitvars[v]), numpy.std(fitvars[v])
            callback(r, fitstats)

    fitstats = {}       
    for v in varnames:
        #print fitvars[v]
        err = numpy.std(fitvars[v])/2
        #print v, numpy.mean(fitvars[v]), err
        #we store the final error as +/- half the stdev.
        fitstats[v] = numpy.mean(fitvars[v]), err
    if doplots == True:           
        from PEATDB.Ekin.Web import EkinWeb
        ew = EkinWeb()
        ew.showEkinPlots(ekindata=estdata, outfile='est_exp_err.html',
                            columns=3, normalise=False, showfitvars=True,
                            imgpath='.')

    return fitstats

def findBestModel(ek, models, fitters=None, checkfunc=None,
                    conv=None, grad=None, alpha=0.05, silent=False):
    """Finds the best fit model using f-testing.
       We do an f-test here between each successive model
       we iterate from simplest to last comparing each model to the more complex ones
       when it finds a better model, this is the best one and is compared with later ones
       in turn.
    """
    import copy, operator
    p=None
    modelfits = {}
    modelstocompare = []

    #sort models by number of parameters first
    sortm=[]
    for m in models:
        X = getFitter(model=m)
        n = len(X.variables)
        sortm.append((m,n))
    sortm = sorted(sortm, key=operator.itemgetter(1))
    models = map(operator.itemgetter(0), sortm)

    #fit with all models first
    print 'Estimating best model..'

    for m in models[:]:
        tempfit, X = doFit(ek, model=m, silent=True, noiter=100, conv=conv, grad=grad)

        if tempfit != None:
            modelfits[m] = copy.deepcopy(tempfit)
            modelstocompare.append(m)

    if len(modelstocompare)<=1:
        return None, None

    best=modelstocompare[0]
    numdps = ek.length()
    modelstocompare.remove(best)
    for n in modelstocompare:
        if silent == False: print best,'vs.',n
        error1 = float(modelfits[best]['error'])
        error2 = float(modelfits[n]['error'])
        if silent == False: print 'error1',error1, 'error2',error2

        numparams1 = len(getFitter(model=best).getVariables())
        numparams2 = len(getFitter(model=n).getVariables())
        if silent == False: print 'numparams1 ',numparams1;print 'numparams2 ',numparams2
        result, p = Utils.doFtest(modelfits[best],modelfits[n],numparams1,numparams2,numdps,alpha=alpha)

        if checkfunc != None and result == 1:
            #use some other checking function, define externally
            result = checkfunc(ek.x, modelfits[n])
        if result == 1:
            if silent == False: print 'Better'
            best=n
        else:
            if silent == False: print 'Worse'
            pass
        if result == 0:
            break
    print 'Best model is', modelfits[best]['model'] , 'with error', modelfits[best]['error']
    print '--------------------------------------------------------- '
    fitdata = modelfits[best]
    return fitdata, p


def checkforOutliers(olddata, oldfitdata):
    """
    It may happen that the simpler fit has a lower error than
    the more complex model because of far outliers, so we check this
    by removing any outliers. (We get stdev of datapoints and see if
    each datapoint is within 3 stdevs of fit curve).
    Then fit again and see if error is better.
    see http://www.graphpad.com/help/prism5/prism5help.html?reg_oultlier_2.htm
    """

    print 'Checking for outliers'
    x=[]
    #we don't want to overwrite the datapoints to avoid confusion
    import copy
    newdata = copy.deepcopy(olddata)

    #xd,yd,a = EkinConvert.ekin2xy(olddata)
    xd,yd = oldata.getxy()
    expdata = zip(xd,yd)
    mean,stdev = Utils.meanstdv(yd)
    print 'meanstdv:',mean,stdev
    for i in range(0,newdata.length()):
        y=newdata.y[i]
    if y>=mean+stdev*2 or y<mean-stdev*2:        
        print y,':point is an outlier?'
        newdata.active[i] = 0

    model = oldfitdata['model']
    olderror = oldfitdata['error']
    #newfit = self.fitDataset(datatofit, oldfitdata, model=model)
    newfit, X = doFit(newdata, model=model)
    newerror = float(newfit['error'])
    return newerror, newfit

def getErrDistribution():
    return

#create module level properties
#Load models if available and create a set of fitters
currentmodels, modelsfile = loadModelsFile()
createFitters()
                 
if __name__=='__main__':
    x=range(1,100); y=range(2,101)
    f=getFitter('Linear',vrs=[2,1.5],expdata=[x,y])    
    f.fit()
    print f.getFitDict()

