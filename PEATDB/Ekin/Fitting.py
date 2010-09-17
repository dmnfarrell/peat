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
import numpy, math, random
import Utils

class Fitting(object):
    """Convenience class to handle LM_Fitter objects for ekin usage,
       replaces some of the FITTER class functions.
       Uses:
       -Provide pre-written fitter subclasses and indexing them
       by model name.
       -Allows us to utilise LM_Fitter objects for creating and fitting
       ekin fit data.
       : no tkvars here please!"""

    def __init__(self):
        return

    #add new models here and then create them as a subclass of LM_fitter
    #dict key is a label and value is the class name
    models = {'Linear':'linear',
              'Power':'power',
              'sigmoid':'sigmoid',
              '1st order decay':'one_decay',
              '1 pKa 2 Chemical shifts':'HH1pKa2shifts',
              '2 pKas, 3 Chemical shifts':'HH2pKa3shifts',
              '3 pKas, 4 Chemical shifts':'HH3pKa4shifts',
              'Bell-shaped pH-act profile (2 pKas)':'pHActivity2pkas',
              'Bell-shaped pH-act profile (3 pKas)':'pHActivity3pkas',
              'Michaelis-Menten':'MichaelisMenten',
              'Modified Hill':'modifiedHill',
              'Unfolding':'generalUnfolding',
              'Chemical Denaturation': 'chemicalDenaturation',
              'diffDenaturation': 'diffDenaturation',
              'DSC2state':'DSC2state',
              'DSCindependent':'DSCindependent',
              'DSC2stateIrreversible':'DSC2stateIrreversible',
              'Residual Activity':'residualActivity'}
    conv=1e-9
    grad=1e-9

    @classmethod
    def getFitter(cls, model, vrs=None, expdata=None, callback=None):
        """Return the required LM fit object from the model name"""
        
        #self.__class__
        if model == 'Linear':
            return linear(variables=vrs,exp_data=expdata,callback=callback,name=model)
        if model == 'Power':
            return power(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == 'Sigmoid':
            return sigmoid(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == '1 pKa 2 Chemical shifts':
            return HH1pKa2shifts(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == '2 pKas, 3 Chemical shifts':
            return HH2pKa3shifts(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == '3 pKas, 4 Chemical shifts':
            return HH3pKa4shifts(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == '3 pKas, 2 Chemical shifts':
            return HH3pKa2shifts(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == 'Bell-shaped pH-act profile (2 pKas)':
            return pHActivity2pkas(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == 'Bell-shaped pH-act profile (3 pKas)':
            return pHActivity3pkas(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == 'Michaelis-Menten':
            return MichaelisMenten(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == 'singlepKa':
            """no offset"""
            return singlepKa(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == 'Modified Hill':
            return modifiedHill(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == 'Chemical Denaturation':
            return chemicalDenaturation(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == 'Unfolding':
            return generalUnfolding(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == 'diffDenaturation':
            return diffDenaturation(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == 'elwellschellman':
            return elwellschellman(variables=vrs,exp_data=expdata,callback=callback,name=model)        
        elif model == 'DSC2state':
            return DSC2state(variables=vrs,exp_data=expdata,callback=callback,name=model)           
        elif model == 'DSCindependent':
            return DSCindependent(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == 'DSC2stateIrreversible':
            return DSC2stateIrreversible(variables=vrs,exp_data=expdata,callback=callback,name=model)    
        elif model == 'DSC2stateIrreversibleII':
            return DSC2stateIrreversibleII(variables=vrs,exp_data=expdata,callback=callback,name=model)             
        elif model == 'Residual Activity':
            return residualActivity(variables=vrs,exp_data=expdata,callback=callback,name=model)
        elif model == 'Arrhenius':
            return Arrhenius(variables=vrs,exp_data=expdata,callback=callback,name=model)          
            
        else:
            return None

    @classmethod
    def makeFitter(cls, fitdata, expdata=None, callback=None):
        """Return a fit object created from the provided ekin fit data"""
        if fitdata.has_key('model'):
            model = fitdata['model']
            vrs = cls.getFitVars(fitdata)
            if len(vrs)==0:
                return None
        else:
            print 'no model provided to create fitter'
            return None
        return cls.getFitter(model, vrs=vrs, expdata=expdata, callback=callback)

    @classmethod
    def getFitVars(cls, fitdata):
        """Returns list of fit variables"""
        fitvars = []
        for var in fitdata.keys():
            if var != 'model' and var != 'error':
                fitvars.append(fitdata[var])
        return fitvars

    @classmethod
    def makeFitData(cls, model, vrs=None, error=None):
        """Get preset ekin fit data for a model, can also use provided
           variables and error vals - returns an ekin fitdata item"""
        fitdata={}
        fitdata['model']=model
        if error == None:
            fitdata['error']=0.0
        else:
            fitdata['error']=error
        i=0
        if vrs == None:
            #get default vals for model
            X=cls.getFitter(model,None,None)
            vrs=X.variables
            for v in vrs:
               fitdata[i]=v
               i+=1
        else:
            for v in vrs:
               fitdata[i]=v
               i+=1
        return fitdata

    @classmethod
    def getExpData(cls, ek):
        """Get paired active data points for fitter from ekin dataset"""
        #print data
        if ek == None:            
            return None 
               
        x,y = ek.getActive()
        if len(x)==0 or len(y)==0:
            return None
        expdata = zip(x,y)       
        return expdata

    @classmethod
    def doFit(cls, ekindata=None, fitdata=None, model=None, expdata=None, vrs=None,
                    noiter=50, conv=None, grad=None, LM_damper=1.0,
                    guess=True, silent=False, callback=None,
                    changevars=None, startvalues=None):
        """General class method for doing a fit of provided ekin data and model
           Returns the ekin fit data and the fitter instance used"""

        if conv==None:
            conv=cls.conv
        if grad==None:
            grad=cls.grad

        def updatenow(diff, vrs, fitvals, c, X):
            """Dummy callback for fit update"""
            pass
        if callback == None:
            callback = updatenow
            
        if expdata == None:            
            expdata = cls.getExpData(ekindata)

        if expdata == None or len(expdata)<=1:
            print 'Not enough data to fit'
            return None, None

        if fitdata != None and len(fitdata) > 1:
            X = cls.makeFitter(fitdata, expdata=expdata, callback=callback)
        else:
            X = cls.getFitter(model, expdata=expdata, vrs=vrs, callback=callback)
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
        fitdata = cls.makeFitData(model, X.variables, X.getError())

        if silent == False:
            print X.getResult()
        return fitdata, X
        
    @classmethod
    def estimateExpUncertainty(cls, ekindata, fitdata, xuncert=0, yuncert=0,
                                     runs=10, doplots=False, callback=None, guess=True,
                                     **kwargs):
        """Chresten's method for estimating the error on each parameter due
           to the effect of the exp uncertainty, which is user provided.
           If the datapoints have x or y error values, we use that instead"""

        '''we iterate over runs and re-fit the data, each time adjusting x-y points
           then accumulate a fit val for each parameter and get mean/stdev.
           Note: errors are added +/- values'''

        grad=1e-6; conv=1e-6
        for k in kwargs:
            if k=='conv':
                conv=kwargs[k]
            if k=='grad':
                grad=kwargs[k]
        if fitdata == None or len(fitdata) == 0:
            return None
        model = fitdata['model']

        F = cls.getFitter(model=model)
        try:
            varnames = F.getVarNames()
        except:
            return None
        fitvars = {}
        for v in varnames:
            fitvars[v] = []     #group lists by variable name

        #xd,yd,a,xerrs,yerrs = EkinConvert.ekin2xy(ekindata, geterrors=True)
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

            #edata = EkinConvert.xy2ekin([mut_x, mut_y], activelist=a)
	    ek = EkinDataset(xy=[mut_x,mut_y],active=a)

            fitres, X = cls.doFit(ek, model=model, fitdata=fitdata, grad=1e-6, conv=1e-6,
                                    silent=True, guess=guess, noiter=100)
            vrs = X.getVariables()
            #print vrs, varnames
            for n in range(len(varnames)):
                fitvars[varnames[n]].append(vrs[n])
            estdata[str(r)] = edata
            estdata['__datatabs_fits__'][str(r)] = fitres

            if callback != None:
                fitstats = {}
                for v in varnames:
                    fitstats[v] = numpy.mean(fitvars[v]), numpy.std(fitvars[v])
                callback(r, fitstats)

        fitstats = {}
        #print 'results'
        for v in varnames:
            #print fitvars[v]
            err = numpy.std(fitvars[v])/2
            #print v, numpy.mean(fitvars[v]), err
            #we store the final error as +/- half the stdev.
            fitstats[v] = numpy.mean(fitvars[v]), err
        if doplots == True:
            #from PEATDB.Ekin.Base import EkinProject
            #Es = EkinProject(data=estdata)
            #Es.plotDatasets()            
            from PEATDB.Ekin.Web import EkinWeb
            ew = EkinWeb()
            ew.showEkinPlots(ekindata=estdata, outfile='est_exp_err.html',
                                columns=3, normalise=False, showfitvars=True,
                                imgpath='.')

        return fitstats

    @classmethod
    def findBestModel(cls, ekindata, models, checkfunc=None,
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
            X=cls.getFitter(model=m)
            n=len(X.variables)
            sortm.append((m,n))
        sortm = sorted(sortm, key=operator.itemgetter(1))
        models = map(operator.itemgetter(0), sortm)

        #fit with all models first
        print 'Estimating best model..'

        for m in models[:]:
            tempfit, X = cls.doFit(ekindata, model=m, silent=True, noiter=100, conv=conv, grad=grad)

            if tempfit != None:
                modelfits[m] = copy.deepcopy(tempfit)
                modelstocompare.append(m)

        if len(modelstocompare)<=1:
            return None, None

        best=modelstocompare[0]
        numdps = len(ekindata[0].keys())-1
        modelstocompare.remove(best)
        for n in modelstocompare:
            if silent == False: print best,'vs.',n
            error1 = float(modelfits[best]['error'])
            error2 = float(modelfits[n]['error'])
            if silent == False: print 'error1',error1, 'error2',error2

            numparams1 = len(Fitting.getFitter(model=best).getVariables())
            numparams2 = len(Fitting.getFitter(model=n).getVariables())
            if silent == False: print 'numparams1 ',numparams1;print 'numparams2 ',numparams2
            result, p = Utils.doFtest(modelfits[best],modelfits[n],numparams1,numparams2,numdps,alpha=alpha)

            if checkfunc != None and result == 1:
                #use some other checking function, define externally
                result = checkfunc(ekindata, modelfits[n])
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

    @classmethod
    def checkforOutliers(cls, olddata, oldfitdata):
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
        newfit, X = cls.doFit(newdata, model=model)
        newerror = float(newfit['error'])
        return newerror, newfit

    def getErrDistribution(self):

        return

#
# -- Fitting subclasses --
#

class HH1pKa2shifts(LM_Fitter):
    """1 pKa 2 Chemical shifts"""
    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['pKa','span','offset']
        self.labels = ['pH','Chem. shift']
        if variables == None:
            self.variables = [7.0, 2.0, 100.0]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data] 
        self.variables[0]=min(x)+(max(x)-min(x))/2.0
        self.variables[1]=max(y)-min(y)
        self.variables[2]=min(y)
        return self.variables

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = 'span1/(1+10**(- x + pKa1))+offset'
        return eq

    def get_value(self, variables, data_point):
        """1 pKa, 2 shifts"""
        try:
            x=data_point[0]
        except:
            x=data_point
        pKa1=variables[0]
        span1=variables[1]
        offset=variables[2]
        value = span1/(1+10**(- x + pKa1))+offset
        return value

class HH2pKa3shifts(LM_Fitter):
    """2 pKa 3 Chemical shifts
       Two groups titrating independently, but influencing each others chemical shift"""
    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['pKa1','pKa2','offset','span1','span2']
        self.labels = ['pH','Chem. shift']
        if variables == None:
            self.variables = [5.0,7.0,100,10.0,10.0]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data]
        if self.variables[0]==0:
            self.variables[0]=min(x)+(max(x)-min(x))*0.25
        if self.variables[0]==0:
            self.variables[1]=min(x)+(max(x)-min(x))*0.75
        self.variables[3]=(max(y)-min(y))/2
        self.variables[4]=(max(y)-min(y))/2
        self.variables[2]=min(y)
        return

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = 'span1/(1+10**(pKa1-pH))+span2/(1+10**(-pH+pKa2))+offset'
        return eq

    def get_value(self, variables, data_point):
        """1 pKa, 2 shifts"""
        try:
            x=data_point[0]
        except:
            x=data_point
        pKa1=variables[0]; pKa2=variables[1]
        span1=variables[3]; span2=variables[4]
        offset=variables[2]
        value = span1/(1+10**(pKa1-x))+span2/(1+10**(-x+pKa2))+offset
        return value

class HH3pKa4shifts(LM_Fitter):
    """3 pKa 4 Chemical shifts
       Two groups titrating independently, but influencing each others chemical shift"""
    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['pKa1','pKa2','pKa3','offset','span1','span2','span3']
        self.labels = ['pH','Chem. shift']
        if variables == None:
            self.variables = [5.0,7.0,9.0,100,10.0,10.0,1.0]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data]
        a=max(x)-min(x)
        b=max(y)-min(y)
        self.variables[0]=min(x)+(a*0.25)
        self.variables[1]=min(x)+(a*0.5)
        self.variables[2]=min(x)+(a*0.75)
        self.variables[4]=b/3
        self.variables[5]=b/3
        self.variables[6]=b/3
        self.variables[3]=min(y)
        return

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = 'span1/(1+10**(pKa1-pH))+span2/(1+10**(-pH+pKa2))+span3/(1+10**(-pH+pKa3))+offset'
        return eq

    def get_value(self, variables, data_point):
        """1 pKa, 2 shifts"""
        try:
            x=data_point[0]
        except:
            x=data_point
        pKa1=variables[0]; pKa2=variables[1]; pKa3=variables[2]
        span1=variables[4]; span2=variables[5]; span3=variables[6]
        offset=variables[3]
        value = span1/(1+10**(pKa1-x))+span2/(1+10**(-x+pKa2))+span3/(1+10**(-x+pKa3))+offset
        return value

class HH3pKa2shifts(LM_Fitter):
    """3 pKa 2 Chemical shifts
       Two groups titrating independently, but influencing each others chemical shift"""
    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['pKa1','pKa2','pKa3','offset','span1']
        self.labels = ['pH','Chem. shift']
        if variables == None:
            self.variables = [5.0,7.0,7.0,100,1.0]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data]
        a=max(x)-min(x)
        b=max(y)-min(y)
        self.variables[0]=min(x)+(a*0.25)
        self.variables[1]=min(x)+(a*0.5)
        self.variables[2]=min(x)+(a*0.75)
        self.variables[4]=b/2
        self.variables[3]=min(y)
        return

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = 'span*(pH/pK4 + 1)/(1+pH/pK3+pH/pK4+pH*pH/(pK1*pK4))+offset'
        return eq

    def get_value(self, variables, data_point):
        """1 pKa, 2 shifts"""
        try:
            x=data_point[0]
        except:
            x=data_point
        pKa1=variables[0]; pKa2=variables[1]; pKa3=variables[2]
        span=variables[2]; offset=variables[3]
        value = span*(x/pKa3 + 1)/(1+x/pKa2+x/pKa3+x*x/(pKa1*pKa3))+offset
        return value


class singlepKa(LM_Fitter):
    """1 pKa 2 Chemical shifts, no offset param"""
    def __init__(self, variables, exp_data, callback, name=None):
        LM_Fitter.__init__(self,variables,exp_data,callback,name=None)
        self.names = ['pKa','span']
        self.labels = ['pH','Chem. shift']
        if variables == None:
            self.variables = [7.0, 2.0]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data]
        self.variables[1]=max(y)-min(y)
        self.variables[0]=min(x)+(max(x)-min(x))/2.0
        return self.variables

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = 'span1/(1+10**(- x + pKa1))'
        return eq

    def get_value(self, variables, data_point):
        """1 pKa, 2 shifts"""
        try:
            x=data_point[0]
        except:
            x=data_point
        pKa1=variables[0]
        span1=variables[1]
        value = span1/(1+10**(- x + pKa1))
        return value

class pHActivity2pkas(LM_Fitter):
    """bell shaped ph activity profile, 2pkas"""
    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['pka1', 'pka2', 'scale']
        self.labels = ['ph', 'kcat']
        if self.variables == None:
            self.variables = [4.0, 7.0, 100.0]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data] 
        peakx = x[y.index(max(y))]
        print peakx
        self.variables[0] = peakx - (max(x)-min(x))/4
        self.variables[1] = peakx + (max(x)-min(x))/4
        if self.variables[1] <= peakx:
            self.variables[1] = max(x)+(max(x)-min(x))/2
        return

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = 'scale*(1.0/(math.pow(10,-pH)/math.pow(10,-pKa2)+1+math.pow(10,-pKa1)/math.pow(10,-pH)))'
        return eq

    def get_value(self, variables, data_point):
        """Bell shaped pH activity profile"""
        try:
            x=data_point[0]
        except:
            x=data_point
        pka1=variables[0]; pka2=variables[1]
        scale=variables[2];

        value = scale*(1.0/(math.pow(10,-x)/math.pow(10,-pka1)+1+math.pow(10,-pka2)/math.pow(10,-x)))
        return value

class pHActivity3pkas(LM_Fitter):
    """bell shaped ph activity profile, 3 pkas"""
    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data, callback,name)
        self.name = name
        self.names = ['pKa1','pKa2','scale','pKa3']
        self.labels = ['ph', 'kcat']
        if self.variables == None:
            self.variables = [4.0,7.0,100.0,7.0]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data] 
        peakx = x[y.index(max(y))]
        if self.variables[1] <= peakx:
            self.variables[1] = max(x)+(max(x)-min(x))/2
        return

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = 'scale*(1.0/(math.pow(10,-pH)/math.pow(10,-pKa1)+1+math.pow(10,-pKa2)/math.pow(10,-pKa3)+math.pow(10,-pKa2)/math.pow(10,-pH)))'
        return eq

    def get_value(self, variables, data_point):
        """Bell shaped pH activity profile"""
        try:
            x=data_point[0]
        except:
            x=data_point
        pKa1=variables[0]; pKa2=variables[1]
        scale=variables[2]; pKa3=variables[3]
        value = scale*(1.0/(math.pow(10,-x)/math.pow(10,-pKa1)+1+math.pow(10,-pKa2)/math.pow(10,-pKa3)+math.pow(10,-pKa2)/math.pow(10,-x)))
        return value

class linear(LM_Fitter):
    """Fit to a line.
       Construct a class like this for each function you want to fit"""
    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data, callback, name)
        self.name = name
        self.names = ['a','b']
        self.labels = ['x', 'y']
        if self.variables == None:
            self.variables = [1.0,1.0]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data] 
        if self.variables[0]==0:
            print 'Guesssed 0'
            self.variables[0] = max(y)-min(y)/max(x)-min(x)
        if self.variables[1]==0:
            self.variables[1] = min(y)
        return

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = 'a * x + b'
        return eq

    def get_value(self, variables, data_point):
        """Sigmoid: tm, bottom, top, slope"""
        try:
            x=data_point[0]
        except:
            x=data_point
        a=variables[0]; b=variables[1]
        value = a * x + b
        return value

class MichaelisMenten(LM_Fitter):
    """Michaelis Menten kinetics"""
    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['Km','Vmax']
        self.labels = ['[S]','v']
        if variables == None:
            self.variables = [1.0, 1.0]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data] 
        self.variables[0] = max(x)/2
        self.variables[1] = max(y)
        return

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = '(Vmax * x )/( x + Km)'
        return eq

    def get_value(self, variables, data_point):
        """Michaelis-Menten"""
        try:
            x=data_point[0]
        except:
            x=data_point
        Km=variables[0]
        Vmax=variables[1]
        value = (Vmax * x )/( x + Km)
        return value


class sigmoid(LM_Fitter):
    """This is an example of a fitter that fits to a sigmoid
       Construct a class like this for each function you want to fit"""
    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data, callback,name)
        self.name = name
        self.names = ['tm', 'bottom', 'top', 'slope']
        self.labels = ['mdeg', 'temp']
        if self.variables == None:
            self.variables = [100, 2, 20, 2]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data]
        try:
            self.variables[0]=(max(x)-min(x))/2+min(x)
            self.variables[1]=min(y)
            self.variables[2]=max(y)
        except:
            pass
        return

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq ='b + (t-b)/(1+exp((Tm-x)/slope))'
        return eq

    def get_value(self, variables, data_point):
        """Sigmoid: tm, bottom, top, slope"""
        #value = self.function(variables, data_point[0])
        try:
            x=data_point[0]
        except:
            x=data_point
        tm=variables[0]; bottom=variables[1]
        top=variables[2]; slope=variables[3]
        if slope <0.03: slope = 0.03   #range error
        try:
            value = bottom + (top - bottom) / (1 + math.exp((tm-x)/slope))
        except:
            value = 0.0
        return value

class modifiedHill(LM_Fitter):
    """Modified Hill Equation for multiple titrating pKas"""
    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['pKa','span','offset', 'n']
        self.labels = ['pH','Chem. shift']
        if variables == None:
            self.variables = [7.0, 2.0, 100.0, 1]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data] 
        self.variables[0]=min(x)+(max(x)-min(x))/2.0
        self.variables[1]=max(y)-min(y)
        self.variables[2]=min(y)
        return self.variables

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = 'span1/(1+10**(n*(- x + pKa1)))+offset'
        return eq

    def get_value(self, variables, data_point):
        """Hill"""
        try:
            x=data_point[0]
        except:
            x=data_point
        pKa1=variables[0]
        span1=variables[1]
        offset=variables[2]
        n=variables[3]
        value = span1/(1+10**(n*(- x + pKa1)))+offset
        return value

class generalUnfolding(LM_Fitter):
    """Fit to a denaturation curve. with 6 parameters, including baseline slopes"""

    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['an','bn','ad','bd','m','d50']
        self.labels = ['U','e']
        if variables == None:
            self.variables = [0, 1, 1, 1, 1,1]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """This guess works best for data normalised to 1"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data] 
        #ad and bd are the slopes of the baselines        
        #an and bn are the intrinsic signals at native & denatured states
        an=min(y); bn=max(y)
        ad=1;bd=1;d50=1;m=1
        halfy=an+(bn-an)/2
        cl=100
        for i in sorted(y):
            diff = abs(abs(i) - abs(halfy))
            if diff < cl:
                cl = diff
                y50=i
                                      
        d50=x[y.index(y50)]
        #print an,bn,halfy,y50,d50
        self.variables = [an,bn,ad,bd,m,d50]
        return self.variables

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = '((an + bn * x) + (ad + bd * x) * math.exp(m * (x-d50))) / (1+math.exp(m * (x-d50)))'
        return eq

    def get_value(self, variables, data_point):       
        try:
            x=data_point[0]
        except:
            x=data_point
        an=variables[0]
        bn=variables[1]
        ad=variables[2]
        bd=variables[3]
        m=variables[4]
        t50=variables[5]
        R=8.3144
        term1 = an + bn * x
        term2 = ad + bd * x
        #term3 = math.exp(m * (x-t50)/R)        
        term3 = math.exp(m*(x-t50))
        value = (term1 + term2 * term3) / (1+term3)
        return value
    
class chemicalDenaturation(generalUnfolding):
    """Fit to a denaturation curve. with 6 parameters,
        inclusing baseline slopes, ref Ferguson/Fersht"""

    def __init__(self, variables, exp_data, callback, name):
        generalUnfolding.__init__(self,variables,exp_data,callback,name)
        return

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = '((an + bn * x) + (ad + bd * x) * math.exp(m * (x-d50)/R*T)) / (1+math.exp(m * (x-d50)/R*T))'
        return eq

    def get_value(self, variables, data_point):
        try:
            x=data_point[0]
        except:
            x=data_point
        an=variables[0]
        bn=variables[1]
        ad=variables[2]
        bd=variables[3]
        m=variables[4]
        d50=variables[5]
        R=8.3144
        T=300
        term1 = an + bn * x
        term2 = ad + bd * x       
        term3 = math.exp(m * (x-d50)/R*T)  
        value = (term1 + term2 * term3) / (1+term3)
        return value

class diffDenaturation(LM_Fitter):
    """Fit to a 3 parameter differential denaturation curve,
        ref. John/Weeks Protein Science 2000"""

    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['A', 'Tm','deltaH']
        self.labels = ['T','e']
        if variables == None:
            self.variables = [1, 340, 40]
        else:
            self.variables = variables
        self.setChangeVars()

        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[];y=[]
        for a in self.exp_data:
            x.append(float(a[0]))
            y.append(float(a[1]))
            
        R=8.3144e-3
        Tm = x[y.index(max(y))]
        deltaH = Tm
        A=1      
        self.variables = [A,Tm,deltaH]
        return self.variables

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = 'A* (deltaH / (R*math.pow(x,2)) ) *f * (1-f)'
        return eq

    def get_value(self, variables, data_point):
        """3 parameter differential denaturation """
        try:
            x=data_point[0]
        except:
            x=data_point
        A=variables[0]
        Tm=variables[1]
        deltaH=variables[2]

        R=8.3144e-3
        K = math.exp((deltaH/R) * (1/Tm - 1/x))
        f = K / (K+1)
        #value = A * (deltaH / (R*math.pow(x,2)) ) *f * (1-f)
        value = A * f * (1-f) * math.pow(x,2)
        
        return value

class elwellschellman(LM_Fitter):
    """Fit raw data to all thermodynamic params at same time
       Elwell/Schellman, Biochim Biophys Acta 1977"""

    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['Tm','deltaH','deltacp']
        self.labels = ['T','e']
        if variables == None:
            self.variables = [340, 100, 1]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[];y=[]
        for a in self.exp_data:
            x.append(float(a[0]))
            y.append(float(a[1]))
            
        R=8.3144e-3
        Tm = x[y.index(max(y))]
        deltaH = Tm; deltacp=1       
        self.variables = [Tm,deltaH,deltacp]
        return self.variables

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = 'deltaH * (1- (t/tm)) - deltacp( (tm-t)+ (t*ln(t/tm)) )'       
        return eq

    def get_value(self, variables, data_point):
        """Fit raw data to all thermodynamic params at same time """
        try:
            x=data_point[0]
        except:
            x=data_point
        Tm=variables[0]
        deltaH=variables[1]
        deltacp=variables[2]
        R=8.3144e-3
        value = deltaH * (1- (x/Tm)) - deltacp * ((Tm-x)+ (x * math.log(x/Tm)))
        return value
    
class DSCindependent(diffDenaturation):
    """Fit to a 3 parameter DSC curve"""
    def __init__(self, variables, exp_data, callback, name):
        diffDenaturation.__init__(self,variables,exp_data,callback,name) 
        return
        
    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data]             
        R=8.3144e-3
        Tm = x[y.index(max(y))]
        deltaH = Tm-80
        A=1e3
        self.variables = [A,Tm,deltaH]
        return self.variables
        
    def get_value(self, variables, data_point):
        """3 parameter differential denaturation """
        try:
            x=data_point[0]
        except:
            x=data_point
        A=variables[0]
        Tm=variables[1]
        deltaH=variables[2]
        R=8.3144e-3
        K = math.exp((deltaH/R) * (1/Tm - 1/x))
        f = K / (K+1)
        value = A* (math.pow(deltaH,2) / (R*math.pow(x,2)) ) *f * (1-f)     
        return value
        
class DSC2state(LM_Fitter):
    """Fit to a 4 parameter DSC curve"""

    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['Tm','deltaH']
        self.labels = ['T','e']
        if variables == None:
            self.variables = [340, 40]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data]            
        R=8.3144e-3
        Tm = x[y.index(max(y))]
        deltaH = Tm-80    
        self.variables = [Tm,deltaH]
        return self.variables
        
    def get_value(self, variables, data_point):
        """simple 2 state DSC"""
        try:
            x=data_point[0]
        except:
            x=data_point      
        Tm=variables[0]
        deltaH=variables[1]
        R=8.3144e-3
        K = math.exp((deltaH/R) * (1/Tm - 1/x))
        f = K / (K+1)
        value = deltaH / (R*math.pow(x,2)) * f * (1-f)
        return value
        
class DSCmultipleindependent(LM_Fitter):
    """Fit to multiple independent reactions for DSC curve"""

    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['Tm','deltaH']
        self.labels = ['T','e']
        if variables == None:
            self.variables = [340, 40]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data]            
        R=8.3144e-3
        Tm = x[y.index(max(y))]
        deltaH = Tm-80    
        self.variables = [Tm,deltaH]
        return self.variables
        
    def get_value(self, variables, data_point):
        """DSC"""
        try:
            x=data_point[0]
        except:
            x=data_point      
        Tm=variables[0]
        deltaH=variables[1]
        R=8.3144e-3
        K = math.exp((deltaH/R) * (1/Tm - 1/x))
        f = K / (K+1)
        value = deltaH / ((R*math.pow(x,2)) *f * (1-f))
        return value
        
class DSC2stateIrreversible(LM_Fitter):
    """Fit to 2 state irreversible model for DSC"""

    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['Tm','deltaH','E']
        self.labels = ['T','e']
        if variables == None:
            self.variables = [350, 100, 300]
        else:
            self.variables = variables
        self.setChangeVars()
        return
        
    def get_value(self, variables, data_point):
        """DSC"""
        try:
            x=data_point[0]
        except:
            x=data_point      
        Tm=variables[0]
        deltaH=variables[1]       
        E=variables[2]
        R=8.3144e-3
        uT=(E/R) * (1/Tm - 1/x)       
        value = (deltaH * E /(R*math.pow(x,2))) * math.exp(uT) * (math.exp(-math.exp(uT)))
        return value
        
class DSC2stateIrreversibleII(LM_Fitter):
    """Fit to 2 state irreversible extended model for DSC"""

    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data,callback,name)
        self.name = name
        self.names = ['Tm','deltaH','deltacp','E']
        self.labels = ['T','e']
        if variables == None:
            self.variables = [340, 40,1,1]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data]        
        Tm = x[y.index(max(y))]
        deltaH = Tm-80
        self.variables = [Tm,deltaH,1,1]
        return self.variables
        
    def get_value(self, variables, data_point):
        """DSC"""
        try:
            x=data_point[0]
        except:
            x=data_point      
        Tm=variables[0]
        deltaH=variables[1]
        deltacp=variables[2]
        E=variables[3]
        #R=8.3144e-3
        R=0.00194
        b=2*deltacp/deltaH * (R*math.pow(Tm,2)/E)
        uT=(E/R) * (1/Tm - 1/x)
        z=-(1+b)*math.exp(uT)
        value = deltaH * (E/(R*math.pow(x,2))) * (1+b) * math.exp(uT) * math.exp(z) + (deltacp* (1-math.exp(z)))
        return value
        
class power(LM_Fitter):
    """Fit to a power law"""
    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data, callback, name)
        self.name = name
        self.names = ['a','b']
        self.labels = ['x', 'y']
        if self.variables == None:
            self.variables = [.1,0.0]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq = '(x **a) + b'
        return eq

    def get_value(self, variables, data_point):
        """power"""
        try:
            x=data_point[0]
        except:
            x=data_point
        a=variables[0]; b=variables[1]
        try:
            value = math.pow(x,a) + b
        except:
            value = 0
        return value

class residualActivity(LM_Fitter):
    """For getting Tm  """
    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data, callback,name)
        self.name = name
        self.names = ['t50', 'slope', 'top']
        self.labels = ['temp', 'mdeg']
        if self.variables == None:
            self.variables = [2, .5, 1]
        else:
            self.variables = variables
        self.setChangeVars()
        return

    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data] 
        try:
            self.variables[0]=(max(x)-min(x))/2+min(x)
            self.variables[2]=max(y)
        except:
            pass
        return

    def get_equation(self):
        """Return a text form of the model for printing"""
        eq ='top/((1+(T/T50)^slope))'
        return eq

    def get_value(self, variables, data_point):
        """Res act: t50, slope, top"""
        try:
            x=data_point[0]
        except:
            x=data_point
        t50=variables[0]; slope=variables[1]
        top=variables[2]
        t50=abs(t50)
        if slope >300: slope = 300
        if slope<0.001: slope = 0.001  #prevent a range errors

        value = top/(1+math.pow((x/t50),slope))
        return value

class Arrhenius(LM_Fitter):
    """Fit Ea to derived rate constants vs T using Arrhenius equation"""
    def __init__(self, variables, exp_data, callback, name):
        LM_Fitter.__init__(self,variables,exp_data, callback,name)
        self.name = name
        self.names = ['A', 'Ea']
        self.labels = ['temp', 'k']
        if self.variables == None:
            self.variables = [10,1e-3]
        else:
            self.variables = variables
        self.setChangeVars()
        return
        
    def guess_start(self):
        """Guess start vals for this model"""
        x=[i[0] for i in self.exp_data]
        y=[i[1] for i in self.exp_data]      
        try:
            Ea=self.variables[1]
            self.variables[0] = y[0]/math.exp(-Ea/8.3144e-3*x[0])          
        except:
            pass
        return
        
    def get_equation(self):
        """Return a text form of the model for printing"""
        eq ='A*exp(-Ea/RT)'
        return eq

    def get_value(self, variables, data_point):
        """Arrhenius"""
        try:
            x=data_point[0]
        except:
            x=data_point
        A=variables[0]
        Ea=variables[1]     
        R=8.3144e-3
        value = A * math.exp(-(Ea/R*x))
        return value  
        
