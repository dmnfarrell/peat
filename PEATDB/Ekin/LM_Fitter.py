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

import numpy, math
from math import *

class LM_Fitter:

    def __init__(self, variables, exp_data, callback, change_vars=None, name=None):
        """Initialise the fitter with a list of variables
        The get_difference function must be overwritten in the derived class
        with a function that calculates the difference to be minimized in the fit
        when given two sets of variable dictionaries"""
        
        self.variables = variables        
        if not change_vars:
            self.setChangeVars()
        else:
            self.change_vars=change_vars

        # Store the rest of the data
        self.equation = ''
        self.exp_data=exp_data
        self.callback_function=callback
        self.name = name
        self.residuals = numpy.array(0,float)
        self.stop_fit=None
        self.singular_matrix = 0
        return

    def setChangeVars(self, changevars=None):
        """Change vars is a mask array of 1's and 0's that specifies which
            variables should be fitted if change_vars[varnum] is None,
            then it is not optimised """

        if self.variables==None: return
        self.change_vars=[]
        if changevars == None:
            self.change_vars=[True for v in self.variables]
        else:
            i=0
            for var in self.variables:
                if changevars[i] == 1:
                    self.change_vars.append(True)
                else:
                    self.change_vars.append(False)
                i+=1
        return

    def evaluate(self, x):
        """Evaluate the model at x with the current variables"""
        return self.get_value(self.variables, x)

    def guess_start(self):
        """Guess start vals for this model - optional overrride"""
        return

    def getEquation(self):
        """Return a text form of the model - optional override"""
        return self.equation

    def getName(self):
        """Return model name"""
        return self.name

    def getVariables(self):
        """Get variables"""
        return self.variables

    def getError(self):
        """Return the current error"""
        diff = self.get_difference(self.variables)
        return diff

    def getpercSqDiff(self):
        """Normalised sum of differences for comparison"""
        diff = self.get_difference(self.variables)
        y = tuple(zip(*self.exp_data))[1]
        yrange = max(y) - min(y)
        percsqdiff = diff/yrange
        #print len(y), yrange, diff, percsqdiff
        return percsqdiff

    def getVarNames(self):
        """Get variable names"""
        return self.varnames

    def getFitDict(self):
        """Return the fitted variables formatted in a dict"""
        fitdict={}
        names = self.varnames; vrs = self.variables
        i=0
        for i in range(len(names)):
            fitdict[names[i]] = vrs[i]
            i+=1       
        return fitdict
        
    def getFitLine(self, xdata):
        """Get fit data for Fitter object over range xdata"""
        e=[]
        for i in xdata: e.append(i)
        fity=[]
        #inc = (max(e)-min(e))/50 
        #xx = numpy.arange(min(e)-inc, max(e)+inc, inc).tolist()
        for i in xdata:
            fity.append(self.evaluate(i))            
        return fity
        
    def getResult(self):
        """Return the fit results for printing"""       
        fd = self.getFitDict()
        fd['error'] = self.getError() 
        fd['model'] = self.name
        fd['rmse'] = self.getRMSE()
        fd['rmsr'] = self.getRMSR()
        return fd 
        
    def printResult(self):
        """Return the fit results for printing"""       
        res=self.getResult()
        for v in res.keys():
            result+= v +':'+ "%E" %fd[v] +' '
        return result

    def get_difference(self,function_variables,return_values=None):
        """Calculate r2 (R-squared)"""       
        diff=0.0
        fit_values=[]
        for datapoint_num in range(len(self.exp_data)):
            datapoint=self.exp_data[datapoint_num]
            exp_value=datapoint[-1]
            #fit_value=self.get_value(self.variables,datapoint)
            fit_value=self.get_value(function_variables,datapoint)
            diff=diff+math.pow(exp_value-fit_value,2)
            fit_values.append([datapoint,fit_value])    
       
        if return_values:
            return diff,fit_values
        else:
            return diff

    def getRMSE(self):
        """Return root mean squared error"""
        errs=[]
        for i in range(len(self.exp_data)):
            datapoint=self.exp_data[i]
            exp=datapoint[-1]            
            fit=self.get_value(self.variables,datapoint)
            errs.append(math.pow(exp-fit,2))
        rmse = math.sqrt(numpy.mean(errs))           
        return rmse
        
    def getRMSR(self):
        """Reduced chi squared for fit"""        
        v = len(self.exp_data) - len(self.varnames) - 1
        if v<=0:
            v=len(self.exp_data)
        if v<1:
            v=len(self.variables)
        x = math.sqrt(self.get_difference(self.variables)/v)      
        return x
        
    def get_value(self,function_variables,data_point):
        """To be overridden
        Function should return the value of the function with function_variables at data point
        """
        raise Exception,'You must override this class'

    def callback(self,difference,variables,count,fitter):
        if self.callback_function:

            # Calculate the current fit
            diff,fit_values=self.get_difference(variables,return_values=True)
            #print 'Calling callback function',self.callback_function
            self.callback_function(difference, variables, fit_values, count, fitter)
        return

    def callback_scipy(self,variables):
        if not hasattr(self,'count'):
            self.count=0
        self.count=self.count+1
        if self.callback_function:
            diff,fit_values=self.get_difference(variables,return_values=True)
            self.callback_function(diff,variables,fit_values)
        #print variables
        return

    def simplex_fit(self):
        """Use the Scipy simple algorithm to minimize R^2"""
        try:
            import scipy.optimize
        except:
            return False,[]
        #print self.variables
        #solution=scipy.optimize.fmin(self.get_difference,self.variables,callback=self.callback_scipy)
        #print 'Simple solution',solution
        solution=scipy.optimize.fmin(self.get_difference,self.variables,callback=self.callback_scipy)
        return True,solution

    def fit(self, rounds=500, error_crit=0.00005, gradient_crit = 1e-6, step=1E-8,
                LM_damper = 0.00001, damper_adjustment_factor = 2.0, silent=False):
        """Do 500 rounds of fitting, or until converged"""
        if silent == False:
            print 'Fitting..'
            print 'start vars:', self.variables
        # Set damper
        self.step=step
        self.LM_damper = LM_damper
        
        # Try to set error_crit using exp_data range
        
        if error_crit == None:
            #get y vals
            y=[]
            for i in self.exp_data:
                y.append(i[1])
            yrange = max(y) - min(y)            
            error_crit = yrange/1e6
            print 'setting error_crit:', error_crit

        # Start iterating
        old_diff=self.get_difference(self.variables)
        status = 'MI'
        count=0
        #print 'Step    Diff   LM_damper'
        for x in range(1,rounds):
            self.fit_LM()
            # check for singluar matrix
            if self.singular_matrix == 1:
                status = 'Stopped - singular matrix'
                if silent == False:
                    print 'Stopped (Singular matrix)'
                break

            now_diff=self.get_difference(self.variables)
            count+=1
            self.callback(difference=now_diff,variables=self.variables,count=count,fitter=self)

            # Check convergence
            s = ''
            if now_diff < old_diff:
                s = 'step accepted'
                self.LM_damper = self.LM_damper / damper_adjustment_factor
                old_diff=now_diff
            else:
                s = 'cancelling'
                self.cancel_step()
                self.LM_damper = self.LM_damper * damper_adjustment_factor

            #print '%5d  %6.4f  %6.4e' %(x, now_diff, self.LM_damper), self.variables, s

            # check error criterium
            #print 'Diff: %5.3e, error_crit: %5.3e,| gradient: %5.3e, grad_crit: %5.3e' %(now_diff,error_crit,self.error_gradient,gradient_crit)
            if abs(now_diff)<=error_crit:
                status = 'Converged - error criterium'
                if silent == False:
                    print 'Converged. (Error criterium) Sum of differences: %7.2e' %(now_diff)
                break

            # check gradient criterium
            if self.error_gradient <= gradient_crit:
                status = 'Stopped - gradient criterium'
                if silent == False:
                    print 'Stopped (Gradient criterium) Error gradient: %7.2e' %(self.error_gradient)
                break

            # check for singluar matrix
            if self.singular_matrix == 1:
                status = 'Stopped - singular matrix'
                if silent == False:
                    print 'Stopped (Singular matrix)'
                break

            # Were we told to stop?
            if self.stop_fit==1:
                break
        return status,self.variables


    def fit_LM(self,silent=0):
        """Do Levenberg-Marquardt fitting"""

        J,E =self.get_jacobian()
        self.residuals = E

        JT = numpy.transpose(J)
        JTE = numpy.dot(JT,E)

        JTJ = numpy.dot(JT,J)
        JTJd = JTJ + self.LM_damper*numpy.identity(numpy.shape(JTJ)[0])

        count=0
        while abs(numpy.linalg.det(JTJd)) == 0:
            ## Determinant is zero => matrix is singular, try to fix this by adding the damper again
            #print 'SINGULAR MATRIX - Adding damper again'
            #print numpy.linalg.det(JTJd)
            JTJd = JTJd + self.LM_damper*numpy.identity(numpy.shape(JTJ)[0])
            count=count+1
            if count>100:
                self.singular_matrix = 1
                return

        invJTJd = numpy.linalg.inv(JTJd)
        self.q = -numpy.dot(JTE,invJTJd)        
        for varnum in range(len(self.variables)):
            if self.change_vars[varnum]:
                self.variables[varnum]=self.variables[varnum]+self.q[varnum]
        self.error_gradient = numpy.linalg.norm(JTE)
        return

    def cancel_step(self):
        """Cancel a step"""
        for varnum in range(len(self.variables)):
            if self.change_vars[varnum]:
                self.variables[varnum]=self.variables[varnum]-self.q[varnum]
        return

    def get_jacobian(self,silent=0):
        """Get the Jacobian matrix and errors of the data points"""
        
        # Get the number of data points        
        no_data_points = len(self.exp_data)

        errors = numpy.resize(numpy.array(0,float),[no_data_points])
        jacobian = numpy.resize(numpy.array(0,float),[no_data_points,len(self.variables)])
        
        # calculate the variation of all parameters
        variations=[]
        step =self.step
        for var in range(len(self.variables)):
            if self.change_vars[var]:
                self.variables[var]=self.variables[var]+step
            variations.append(self.variables[:])
            if self.change_vars[var]:
                self.variables[var]=self.variables[var]-step
        
        # construct jacobian        
        data_id=0
        for datapoint_num in range(len(self.exp_data)):
            if datapoint_num%100==0:
                pass
            data_point = self.exp_data[datapoint_num]
            # The value to fit is always the last one
            exp_value=float(data_point[-1])
            errors[datapoint_num] = exp_value-self.get_value(self.variables,data_point)
            
            # Find the derivative for this variable for this data point            
            diff=numpy.resize(numpy.array(0,float),[len(self.variables)])
            count=0
            for variation in variations:
                diff[count]=(self.get_value(self.variables,data_point)-
                             self.get_value(variation,data_point))/step
                count=count+1
            jacobian[data_id]=diff
            data_id=data_id+1

        return jacobian,errors

    def do_statistics(self):
        residualsT = numpy.transpose(self.residuals)
        m = len(self.exp_data)
        n = len(self.variables)
        diagonal = numpy.diag(numpy.resize(numpy.array(1.0,float),m))
        #print 'self.residuals',self.residuals
        covr = numpy.dot(residualsT, self.residuals)
        covr = numpy.dot(covr, diagonal)/(m-n)
        covr = numpy.diag(covr)

        print 'covr', covr
        #covr=resid'*resid/(m-n);  #               %covariance of residuals
        #Vy=1/(1-n/m)*covr;  % Eq. 7-13-22, Bard         %covariance of the data

        # get Jacobian at optimal parameter values
        J,E =self.get_jacobian()

        JT = numpy.transpose(J)
        JTdiagonal = numpy.dot(JT, diagonal)
        JTdiagonalJ = numpy.dot(JTdiagonal, J)

        residualsTdiagonal = numpy.dot(residualsT, diagonal)
        residualsTdiagonalresiduals = numpy.dot(residualsTdiagonal,self.residuals)
        Z = (m-n) *JTdiagonalJ/ (n*residualsTdiagonalresiduals)
        #print 'Z',Z
        #((m-n)*jac'*Qinv*jac)/(n*resid'*Qinv*resid);
        return


class myfitter(LM_Fitter):
    """This is an example of a fitter that fits to a y=ax+b
    Construct a class like this for each function you want to fit"""
    def __init__(self,variables,exp_data):
        LM_Fitter.__init__(self,variables,exp_data,self.callback)
        return

    def get_value(self,function_variables,data_point):
        value=function_variables[0]*data_point[0]+function_variables[1] #a*x+b
        return value

    def callback(self,difference,variables,count):
        #print 'I am called and called and called'
        #print 'variables',variables
        return


class difficult(LM_Fitter):
    def __init__(self,variables,exp_data):
        LM_Fitter.__init__(self,variables,exp_data,self.callback)
        return

    def get_value(self,function_variables,data_point):
        import math
        value = math.sin(function_variables[0] * data_point[0]) + math.cos(function_variables[1] * data_point[0]) ##sin(3*x) + cos(0.5*x)
        return value

    def callback(self,difference,variables,count,fitter):
        #print 'I am called and called and called'
        #print 'variables',variables
        return


if __name__=='__main__':
    #
    # Test of LM_Fitter class
    #

    exp_data=[[0,1],[0.5,1.966],[0.8,1.5965],[1,1.0187],[1.8,-0.15115],[2.6,1.266]]
    X=difficult(variables=[1.8,1.],exp_data=exp_data)

    #exp_data=[[0,1],[0.5,1.966],[0.8,1.5965],[1,1.0187],[1.8,-0.15115],[2.6,1.266]]
    #X=myfitter(variables=[2,1],exp_data=exp_data) ##sin(3*x) + cos(0.5*x)


    status,variables=X.fit(rounds = 10000, error_crit = 0.000001, LM_damper =1.0, damper_adjustment_factor = 1.0, step = 1e-8)
    print
    print 'Done. Status',status
    print 'Fitted variables',variables

    X.do_statistics()


