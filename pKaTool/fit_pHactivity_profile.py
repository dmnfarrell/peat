#!/usr/bin/env python
#
# pKaTool - analysis of systems of titratable groups
# Copyright (C) 2010 Jens Erik Nielsen
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

import LM_Fitter

from numpy import *
from numpy.linalg import * 
inverse=inv
Float=float
class pKaTool_fitter:

    def __init__(self,system_spec,exp_data,CCPS,normalize=True,max_steps=1000,conv_crit=0.000005):
        #
        # Initialize
        #
        self.CCPS=CCPS
        self.normalize=normalize 
        self.exp_pHs=[]
        self.profiles={}
        self.exp_data=exp_data
        for pH,pop in exp_data:
            self.exp_pHs.append(pH)
        #
        # Reformat the variables
        #
        self.groups=system_spec.keys()
        self.groups.sort()
        self.vars=[]
        # intpkas
        for group in self.groups:
            self.vars.append(system_spec[group]['intpka'])
        # Matrix
        for group1 in self.groups:
            for group2 in self.groups:
                if system_spec[group1]['matrix'].has_key(group2):
                    self.vars.append(system_spec[group1]['matrix'][group2])
        #
        # Set up the pKaTool stuff
        #
        self.pHstart=0.0
        self.pHend=20.0
        self.pHstep=1.0
        self.numgroups=len(self.groups)
        self.CCPS=CCPS
        import pKa_calc
        self.PK=pKa_calc.Boltzmann()
        #
        # Set damper
        #
        self.LM_damper = 0.1
        #
        # Start fitting
        #
        old_diff=0.0000001
        #print 'Step    Diff   LM_damper'
        now_diff=self.activity_diff()
        #print '%5d  %6.4f  %6.4f' %(0, now_diff, self.LM_damper)
        for x in range(1,max_steps):
            self.fit_LM_activity()
            now_diff=self.activity_diff()
            #
            # Check convergence
            #
            #print '%5d  %6.4f  %6.4f' %(x, now_diff, self.LM_damper)
            if abs(now_diff-old_diff)<conv_crit:
                #print 'Converged',now_diff
                break
            else:
                old_diff=now_diff
            now_diff=self.activity_diff()
        self.diff=now_diff
        return 

    #
    # ------
    #

    def fit_LM_activity(self):
        """Do Levenberg-Marquardt fitting"""
        J,E =self.get_jacobian_activity()
        JT = transpose(J)
        JTE = dot(JT,E)
        JTJ = dot(JT,J)
        JTJd = JTJ + self.LM_damper*identity(shape(JTJ)[0])
        invJTJd = inv(JTJd)
        q = -dot(JTE,invJTJd)
        #out1 ='bv '
        #out2 ='av '
        for var in range(len(self.vars)):
            #out1 += '%4.6f '%self.vars[var]
            self.vars[var]=max(0.0,self.vars[var]+q[var]) # We don't allow any parameter to drop below 0.0
            #out2 += '%4.6f '%self.vars[var]
        return



    def get_jacobian_activity(self):
        """Get the Jacobian matrix and errors of the data points"""
        #
        # Get the number of data points
        #
        no_data_points = len(self.exp_data)
        errors = resize(array(0,float),[no_data_points])        
        jacobian = resize(array(0,float),[no_data_points,len(self.vars)])
        #
        # Precalculate the variation of all parameters
        #
        now=self.get_act_prof()
        variations=[]
        step = 1e-8
        for var in range(len(self.vars)):
            self.vars[var]=self.vars[var]+step
            variations.append(self.get_act_prof())
            self.vars[var]=self.vars[var]-step
        #
        # construct jacobian
        #
        data_id=0
        #
        # activity data
        #
        for ph,act in self.exp_data:
            x=float(ph)
            y=float(act)
            errors[data_id] = y-now[ph]
            #
            # Find the derivative at this ph (x) for this data point
            #
            diff=resize(array(0,float),[len(self.vars)])
            count=0
            for variation in variations:
                diff[count]=(now[ph]-variation[ph])/step
                count=count+1
            jacobian[data_id]=diff
            data_id=data_id+1
        return jacobian,errors
        
    #
    # -------
    #
        
    def activity_diff(self, just_activity_fit = 0):
        act_prof_and_tc = self.get_act_prof()
        diff = 0.0
        for ph,act in self.exp_data:
            diff = diff + abs(act_prof_and_tc[ph]-act)
        return diff
        
    #
    # -------
    #
    
    def get_act_prof(self):
        """Get the pH-activity profile and titration curves"""
        #
        # titration data
        #
        self.setup_system(self.vars)
        # Calculate the CCPSs
        pKa_values,prot_states=self.PK._calc_pKas(mcsteps=0,phstep=self.pHstep,phstart=self.pHstart,phend=self.pHend,exp_pHs=self.exp_pHs,verbose=1)
        #
        # Primary activity data
        #
        act_prof1={}
        states=0
        for pH in self.exp_pHs:
            act=0.0
            for state in self.PK.all_states[pH].keys():
                if self.PK.all_states[pH][state]['def'] in self.CCPS:
                    act=act+self.PK.all_states[pH][state]['pop']
                    states=states+1
            act_prof1[pH]=act
        if states>0:
            max_act=max(act_prof1.values())
            if max_act==0.0:
                fact=1.0
            else:
                fact=1.0/max_act
            for ph in act_prof1.keys():
                act_prof1[ph] = fact*act_prof1[ph]
        else:
            act_prof1={}
        return act_prof1
        
    #
    # -----
    #
        
    def setup_system(self,variables):
        """Setup the system as per the system_spec dictionary"""
        # 
        # Get the intrinsic pKas
        #
        count=0
        self.PK.intrinsic_pKa={}
        self.PK.matrix={}
        for group in self.groups:
            self.PK.intrinsic_pKa[group]=variables[count]
            count=count+1
        #
        # Matrix
        #
        for group1 in self.groups:
            self.PK.matrix[group1]={}
        for group1 in self.groups:
            for group2 in self.groups:
                if group1<group2:
                    intene=self.PK.acid_base(group1)*self.PK.acid_base(group2)*variables[count]
                    self.PK.matrix[group1][group2]=[intene,0.0,0.0,0.0]
                    self.PK.matrix[group2][group1]=[intene,0.0,0.0,0.0]
                    count=count+1
                elif group1==group2:
                    self.PK.matrix[group1][group2]=[0.0,0.0,0.0,0.0]
        self.PK.prepare_matrix()
        #
        self.PK.groups=self.groups
        return


#
# --------------
#

class profile_fitter:

    def __init__(self,options,args):
        """Instantiate a class that fits an experimental pH-activity profile to a combination of CCPSes"""
        self.options=options
        #
        # read the file with experimental data
        #
        fd=open(options.expdata)
        lines=fd.readlines()
        fd.close()
        self.exp_data=[]
        for line in lines:
            if line.strip()[0]=='#' or line.split()[0]=='pH':
                continue
            sp=line.split()
            self.exp_data.append([float(sp[0]),float(sp[1])])
        #
        # What should we do
        #
        if self.options.analyze:
            self.analyze()
        else:
            self.do_fitting()
        return
            
    def do_fitting(self):
        #
        # Set the group definition
        #
        import string, random
        groups=[]
        for group in range(options.ngroups):
            groups.append('%s:%s' %(string.zfill(group,4),'ASP'))
        #
        CCPS_def=[[0,1,0],[0,1,1]]
        #
        print 'Doing %5d fits with maxstep: %4d with groups: %s and CCPS definition: %s' %(options.nfits,options.max_steps,str(groups),str(CCPS_def))
        solutions=[]
        for fitno in range(options.nfits):
            self.set_random_state(groups)
            score,variables,profile=self.fit_system(CCPS_def)
            print 'Fit# %4d, score: %6.3f' %(fitno,score)
            solutions.append([score,variables,profile,CCPS_def,groups])
        #
        # Save the solutions
        #
        fd=open(options.outfile,'w')
        import pickle
        pickle.dump(solutions,fd)
        fd.close()
        return
        
    #
    # -----
    #
    
    def align_parms(self,solutions):
        """Sort the solutions so parameters are shifted around to give optimal overlap"""
        solutions.sort()
        ref_sol=solutions[0][1]
        aligned_sols=[solutions[0]]
        for score,variables,profile,CCPS_def,groups in solutions:
            best_vars=variables[:]
            diff=self.get_similarity(variables,ref_sol)
            for permutation in range(len(groups)-1):
                nv=len(variables)*[0.0]
                nv[0]=variables[1]
                nv[1]=variables[2]
                nv[2]=variables[0]
                nv[3]=variables[5]
                nv[4]=variables[3]
                nv[5]=variables[4]
                # I should write something general here
                variables=nv[:]
                if self.get_similarity(variables,ref_sol)<diff:
                    #print 'I found a better alignment',self.get_similarity(variables,ref_sol),diff
                    diff=self.get_similarity(variables,ref_sol)
                    best_vars=variables[:]
            aligned_sols.append([score,best_vars,profile,CCPS_def,groups])
        return aligned_sols
        
    #
    # ------
    #
        
    def get_similarity(self,solution1,solution2):
        """Get the similarity of two solutions"""
        import math, numpy
        s1=numpy.array(solution1)
        s2=numpy.array(solution2)
        I=numpy.inner(s1,s2)
        costerm=I/(numpy.linalg.norm(s1)*numpy.linalg.norm(s2))
        #
        # Maybe better to get the length of the difference vector
        #
        return numpy.linalg.norm(s1-s2)
    
    #
    # -------
    #
        
    def analyze(self):
        #
        # Load the pickle file and do the analysis
        #
        fd=open(self.options.outfile)
        import pickle
        solutions=pickle.load(fd)
        fd.close()
        #
        # Sort the parameters - make sure that redundant solutions appear the same
        #
        solutions=self.align_parms(solutions)
        #
        # Print the solutions
        #
        solutions.sort()
        #
        # Keep only solutions that are less than 10% worse than the best solution
        #
        best_score=solutions[0][0]
        ok_sols=[]
        count=1
        for sol in solutions:
            if abs((sol[0]-best_score)/best_score)<=self.options.cluster_cutoff:
                print count,sol[0],sol[1]
                count=count+1
                ok_sols.append(sol)
        #
        # Plot the ranges of the variables
        #
        xs=[]
        ys=[]
        for score,variables,profile,CCPS_def,groups in ok_sols:
            count=1
            for var in variables:
                xs.append(count)
                ys.append(var)
                count=count+1
        import pylab
        pylab.plot(xs,ys,'ro')
        pylab.xlabel('Variable number')
        pylab.ylabel('Value')
        pylab.show()
        #
        # Plot all the fits
        #
        import pylab
        count=0
        for score,variables,profile,CCPS_def,groups in ok_sols:
            pHs=sorted(profile.keys())
            xs=[]
            ys=[]
            for pH in pHs:
                xs.append(pH)
                ys.append(profile[pH])
            if count==0:
                pylab.plot(xs,ys,'r-',label='Fitted solutions')
            else:
                pylab.plot(xs,ys,'r-')
            count=count+1
        #
        # Plot the exp data
        #
        xs=[]
        ys=[]
        for pH,act in self.exp_data:
            xs.append(pH)
            ys.append(act)
        pylab.plot(xs,ys,'bo',label='Experimental data')
        pylab.legend()
        pylab.xlabel('pH')
        pylab.ylabel('Normalized population')
        pylab.show()

        return
        
    #
    # ------
    #
            
    def set_random_state(self,groups):
        #
        # Construct random starting values for the system
        #
        import pKa_calc
        X=pKa_calc.Boltzmann()
        self.system_spec={}
        for group in groups:
            import random
            self.system_spec[group]={'intpka':random.uniform(0,8),'matrix':{}}
        for group1 in groups:
            for group2 in groups:
                if group1<group2:
                    self.system_spec[group1]['matrix'][group2]=random.uniform(0,10.0)
        return
        
    #
    # -----
    #
        
    def fit_system(self,CCPS_def):
        #
        # Fit
        #
        self.FIT=pKaTool_fitter(self.system_spec,self.exp_data,CCPS_def,max_steps=self.options.max_steps,conv_crit=self.options.conv_crit)
        prof=self.FIT.get_act_prof()
        return self.FIT.diff,self.FIT.vars,prof
    
    #
    # -----
    #
        
    def plot_current_fit(self):
        """Plot the current fit"""
        import pylab
        xs=sorted(prof.keys())
        ys=[]
        for x in xs:
            ys.append(prof[x])
        pylab.plot(xs,ys,'r-',label='fitted')
        #
        xs=[]
        ys=[]
        for x,y in self.exp_data:
            xs.append(x)
            ys.append(y)
        pylab.plot(xs,ys,'bo',label='experimental')
        pylab.legend()
        pylab.show()
        print
        print 'Done.'
        print 'Fitted variables',self.FIT.variables
        return
        
    #
    # -------
    #
    

def main():
    print
    print 'Fit enzymatic pH-activity profiles to pKaTool CCPS populations'
    print 'Jens Erik Nielsen, 2010'
    print
    import sys, os
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <file>',version='%prog 1.0')
    parser.add_option('-e',"--expdata",type='string',dest='expdata',action='store',
                      help='File with experimental pH-activity profile. Space separation. Default: %default',default='pH_act.txt')
    parser.add_option('-o',"--outfile",type='string',dest='outfile',action='store',
                      help='output file for solutions. Default: %default',default='solutions.pickle')
    parser.add_option('-n','--nfits',dest='nfits',action='store',type='int',
                        help='Number of random starting positions to evaluate. Default: %default',default=10)
    parser.add_option('-m','--maxsteps',dest='max_steps',action='store',type='int',
                        help='Max number of LM steps in any one fitting run. Default: %default',default=100)
    parser.add_option('-g','--ngroups',dest='ngroups',action='store',type='int',
                        help='Number of titratable groups. Default: %default',default=3)
    parser.add_option('-c','--convcrit',dest='conv_crit',type='float',action='store',
                        help='LM gradient convergence criterion. Default: %default',default=0.000005)
    parser.add_option('-d','--cluster',dest='cluster_cutoff',type='float',action='store',
                        help='Fractional deviation of solution scores from best score that is accepted when filtering solutions. Default: %default',default=0.1)
    parser.add_option('-a','--analyze',dest='analyze',action='store_true',default=False,
                        help='Do only analysis (i.e. skip fitting and use old output file). Default: %default')
    
    (options, args) = parser.parse_args()
    #
    # Call main
    #
    X=profile_fitter(options,args)
    return


if __name__=="__main__":
    main()
    #import cProfile
    #cProfile.run('main()','profiler.out')
    #import pstats
    #p = pstats.Stats('profiler.out')
    #p.sort_stats('cumulative').print_stats(10)


