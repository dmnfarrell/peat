#!/usr/bin/env python

"""Fit stability profile data to extract energies"""

def average(l):
    sum=0.0
    for i in l:
        sum=sum+i
    #print 'Sum: %f, dividing by %d' %(sum,len(l))
    avg=sum/float(len(l))
    #
    # Variance
    #
    N=len(l)
    sum2=0.0
    for i in l:
        sum2=sum2+math.pow((i-avg),2)
    if N>1:
        variance=sum2/float(N-1)
    else:
        variance=0
        print 'Warning - N was 1'
    std_dev=math.sqrt(variance)
    return avg,variance,std_dev

import numpy as np
import os, sys, csv, math
sys.path.append('../')
from myfitter import fitter
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.rc('font',size=10)
plt.rc('font',family='monospace')
plt.rc('axes',linewidth=.5)

R=8.314 #(J/mol K)
T=298.15 # K

import Double_mutant

class fit_stability(Double_mutant.double_mutant):

    def __init__(self,PDB,csvfile,exp_pKa_file,options):
        """Load the PDB file and load the experimental data"""
        self.options=options
        import Protool
        self.PI=Protool.structureIO()
        self.PI.readpdb(PDB)
        self.TGs=self.PI.get_titratable_groups()
        #
        # Read the experimental data
        #
        stab_data=self.read_stability(csvfile)
        #
        # Read the pKa values
        #
        pKa_values=self.read_pKa_values(exp_pKa_file)
        #
        # Should we test the predictions?
        #
        if options.test:
            #self.test_charge()
            self.test_pHstab(pKa_values,stab_data)
        #
        # Do the fitting
        # 
        if options.singlefits:
            pKa_values=self.singlefits(stab_data,pKa_values)
            self.fit_pKa_values(stab_data,pKa_values)
        elif options.doublemuts:
            self.doublemuts(stab_data,pKa_values)
        else:
            #
            # Just fit the full pH-stability profile
            #
            self.fit_pKa_values(stab_data,pKa_values)
        return
        
    #
    # ------
    #
    
    def singlefits(self,stab_data,pKa_values):
        """Perform fits of single pKa values from mutant data"""
        mutants=stab_data.keys()
        mutants.remove('wt')
        mutants.remove('ph')
        results={}
        pkas=[]
        for mutant in mutants:
            results[mutant]=self.fit_mutant_pKa(mutant,pKa_values=pKa_values,stabdata=stab_data)
            for res in results[mutant].keys():
                if not res in pkas:
                    pkas.append(res)
        #
        # Print a nice table with the values
        #
        fit_pKas={}
        pkas.sort()
        txt='%15s' %('Mutant')
        for pka in pkas:
            txt=txt+' %13s' %pka
        print txt
        for mutant in sorted(results.keys()):
            txt='%15s' %mutant
            for pka in pkas:
                if results[mutant].has_key(pka):
                    if not fit_pKas.has_key(pka):
                        fit_pKas[pka]=[]
                    fit_pKas[pka].append(results[mutant][pka][0])
                    avgpka=results[mutant][pka][0]
                    SDpka=results[mutant][pka][2]
                    txt=txt+' %5.1f (%5.1f)' %(avgpka,SDpka)
                else:
                    txt=txt+' %13s' %('-')
            print txt
        #
        # Get the average fitted pKa values
        #
        txt='%15s' %('Average pKa')
        avg_pkas={}
        for pka in pkas:
            avg_pKa=sum(fit_pKas[pka])/float(len(fit_pKas[pka]))
            avg_pkas[pka]=avg_pKa
            txt=txt+' %13.1f' %avg_pKa
        print txt
        #
        # Insert these pKa values into the pKa value matrix and re-fit the complete pH-activity profile
        #
        for pka in avg_pkas.keys():
            if pka=='score':
                continue
            sp=pka.split()
            group=sp[0]
            state='U'
            pKa_values[group][state]=avg_pkas[pka]
            print 'Inserting fitted pKa value for %s state %s' %(group,state)
        return pKa_values
    
    #
    # ------
    #    
        
    def test_charge(self):
        import scipy
        vals=[[],[],[]]
        for pH in scipy.arange(0.1,12,0.1):
            acid=self.get_charge(':0001:ASP',pH,4.0)
            base=self.get_charge(':0002:LYS',pH,8.0)
            vals[0].append(pH)
            vals[1].append(acid)
            vals[2].append(base)
        import pylab
        pylab.plot(vals[0],vals[1],'r-')
        pylab.plot(vals[0],vals[2],'b-')
        pylab.show()
        return
            
    
    #
    # ------
    #
    
    def test_pHstab(self,pKa_values,exp_stab):
        variables=[]
        for group in sorted(self.TGs):
            for state in ['F','U']:
                pka=pKa_values[group][state]
                if pka=='fit':
                    import pKaTool.pKadata
                    group_type=group.split(':')[-1]
                    pka=pKaTool.pKadata.modelpKas[group_type]
                    pKa_values[group][state]=pka
        #
        self.pKa_values=pKa_values.copy()
        stabprof,group_contribs=self.get_pH_stab_profile('wt',variables)
        xs=[]
        ys=[]
        for pH in sorted(stabprof.keys()):
            xs.append(pH)
            ys.append(stabprof[pH])
        import pylab
        pylab.plot(xs,ys,'r.-')
        #
        # Plot the wild type pH-activity profile
        #
        xs=[]
        ys=[]
        for pH,ddG in zip(exp_stab['ph'],exp_stab['wt']):
            xs.append(pH)
            ys.append(ddG)
        pylab.plot(xs,ys,'bo-',label='wt exp')
        pylab.legend()
        pylab.xlabel('pH')
        pylab.ylabel('ddGfold (kJ/mol)')
        pylab.show()
        return
    #
    # ------
    #
    
    def read_stability(self,csvfile):    
        """Load the experimental data"""
        f=open(csvfile,'r')
        cr = csv.reader(f)
        tmp_names = cr.next()
        names=[]
        for name in tmp_names:
            names.append(name.lower().strip())
        data={}
        #
        for n in names:
            data[n.lower()]=[]
        for row in cr:
            for i in range(len(row)):
                if row[i]=='':
                    break
                if options.ddGunfold and names[i]!='ph':
                    x=-float(row[i])
                else:
                    x=float(row[i])
                if options.kcal and names[i]!='ph':
                    x=x*4.19
                if names[i]=='ph':
                    if abs(x-9.0)<0.1:
                        break
                    x='%5.1f' %x
                n=names[i]
                data[n].append(x)
        return data
        

    #
    # ------
    #
    
    def read_pKa_values(self,filename):
        """Read the known and unknown pKa values"""
        pKa_values={}
        fd=open(filename,'rU')
        import csv
        cr=csv.reader(fd)
        names=cr.next()
        pKa_values={}
        for row in cr:
            grp=[]
            for col in range(len(row)):
                grp.append(row[col])
            vals=[]
            for val in row[1:]:
                try:
                    val=float(val)
                except:
                    val=val
                vals.append(val)
            pKa_values[grp[0]]={'F':vals[0],'U':vals[1]}
        #
        # Fill in the pKa values we assume are unchanged
        #
        import pKaTool.pKadata
        print 'I read the following pKa values'
        print '%-15s %7s %7s  %8s' %('Group','pKa(F)','pKa(U)','pKa diff')
        fitcount=0
        totpkas=0
        for group in sorted(pKa_values.keys()):
            #
            # Check that we know this group
            #
            if not group in self.TGs:
                raise Exception('The PDB file does not hold this titratable group: %s' %group)
            #
            txt='%-15s' %group
            vals=[]
            for state in ['F','U']:
                value=str(pKa_values[group][state])
                if value.lower()=='model' or value.lower=='model pka' or value=='m':
                    group_type=group.split(':')[-1]
                    model_pKa=pKaTool.pKadata.modelpKas[group_type]
                    pKa_values[group][state]=model_pKa
                    txt=txt+'%7.2f' %(float(model_pKa))
                    vals.append(model_pKa)
                    totpkas=totpkas+1
                elif value.lower()=='fit':
                    txt=txt+'%7s' %('fit')
                else:
                    txt=txt+'%7.2f' %float(value)
                    vals.append(float(value))
                    totpkas=totpkas+1
            #
            # ----
            #
            if len(vals)==2:
                txt=txt+'  %7.2f' %(vals[0]-vals[1])
            elif len(vals)==1:
                txt=txt+'  %8s' %('from fit')
                fitcount=fitcount+1
            else:
                txt=txt+'  %8s' %('fit')
            print txt
        print 'Total %3d known pKa values and %d values to fit' %(totpkas,fitcount)
        return pKa_values


    #
    # -----
    #

    def getdQ(self,ph, pkas):
        """Get fractional charge for specific set of pka vals at a ph
           uses eq 7 from linse paper"""   
        total=0 
        q=0 
        for pka,Q in zip(pkas, hewlchrg):
            try:
                q += Q/(1+math.pow(10,Q*(ph-pka)))
            except:
                return 1
        total = q*2.3*R*T
        return total
    
    #
    # -----
    #
    
    def get_charge(self,group,pH,pKa):
        """Get the charge of this group at this pH"""
        import pKaTool.pKadata
        ab=pKaTool.pKadata.acidbase[group.split(':')[-1]]
        if ab==-1:
            charge=-1.0/(1.0+math.pow(10,pKa-pH))
        else:
            charge=1.0/(1.0+math.pow(10.0,pH-pKa))
        return charge
        
    #
    # ------
    #
    
    def get_pH_stab_profile(self,protein,function_variables):
        """Get the pH_stability profile of a specific protein (mutant)"""
        import scipy
        pH_values=scipy.arange(self.options.pH_start,16.0,options.pH_step)
        calc_stab={}
        group_contrib={}
        calc_stab={}
        stab=0.0
        pH_count=1
        for pH in pH_values[1:]:
            # Get the pH step for the integration
            pka_count=0
            group_contrib[pH]={}
            thispH=0.0
            for group in sorted(self.TGs):
                #
                # See if this group was mutated
                #
                use_group=True
                if protein!='wt':
                    muts=protein.split('+')
                    for mutation in muts:
                        number=int(mutation[1:-1])
                        group_number=group.split(':')[1]
                        if number==int(group_number):
                            use_group=False
                if use_group:
                    #
                    # Get the dQ for this group at this pH
                    #
                    pka_folded=self.pKa_values[group]['F']
                    #
                    # Change in pKa value of E35 when D52 is removed
                    #
                    pka_changes={'d52a':[[':0035:GLU',5.0]],
                                'd52n':[[':0035:GLU',5.0]],
                                'd66a':[[':0035:GLU',5.8],[':0052:ASP',3.1]],
                                'd48a':[[':0035:GLU',5.8],[':0052:ASP',3.3]]}
                    for mutation in protein.split('+'):
                        if pka_changes.has_key(mutation):
                            for group_change,newpka in pka_changes[mutation]:
                                if group==group_change:
                                    pka_folded=newpka
                                    break
                    if pka_folded=='fit':
                        pka_folded=function_variables[pka_count]
                        pka_count=pka_count+1
                    Q_folded=self.get_charge(group,pH,pka_folded)
                    # Unfolded state
                    pka_unfolded=self.pKa_values[group]['U']
                    if pka_unfolded=='fit':
                        pka_unfolded=function_variables[pka_count]
                        pka_count=pka_count+1
                    Q_unfolded=self.get_charge(group,pH,pka_unfolded)
                    #
                    Qdiff=R*T*math.log(10)*(Q_folded-Q_unfolded)/1000.0
                    thispH=thispH+Qdiff
                    group_contrib[pH][group]=Qdiff
            #
            stab=stab+thispH*options.pH_step
            calc_stab['%5.1f' %pH]=stab 
        return calc_stab,group_contrib
        
    #
    # -----
    #
    
    def get_difference_all(self,function_variables):
        #
        # Calculate the total difference value between predicted and observed pH-stability profiles
        #
        # First calculate the current predicted stability profile
        #
        calc_stab={}
        stab_prof,group_contrib=self.get_pH_stab_profile('wt',function_variables)
        calc_stab['wt']=stab_prof
        #
        # Now we have the calculated pH-stability profile. Compare it to the experimental one
        #
        plot_xs=[]
        plot_ys=[]
        #
        import math
        diffs=[]
        for protein in ['wt']:
            for pH,ddG_exp in zip(self.stab_data['ph'],self.stab_data[protein]):
                diff=ddG_exp-calc_stab[protein][pH]
                #print pH, ddG_exp-calc_stab[protein][pH],diff
                diff=math.pow(diff,2)
                diffs.append(diff)
            self.sumdiffs=sum(diffs)
        return sum(diffs)
    
    def print_score(self,vars):
        #print self.sumdiffs
        return
        
    #
    # -----
    #
        
    def get_difference_mutant(self,function_variables):
        """Calculate the difference between two stability profiles"""
        calc_stab={}
        stab_prof,group_contrib=self.get_pH_stab_profile('wt',function_variables)
        calc_stab['wt']=stab_prof
        # Mutant
        stab_prof,group_contrib=self.get_pH_stab_profile(self.mutant,function_variables)
        calc_stab[self.mutant]=stab_prof
        #
        # Subtract the two calculated profiles
        #
        calc_diff={}
        for pH in calc_stab['wt'].keys():
            diff=calc_stab[self.mutant][pH]-calc_stab['wt'][pH]
            calc_diff[pH]=diff
        #

        exp_diff={}
        for pH,ddG_wt,ddG_mut in zip(self.stab_data['ph'],self.stab_data['wt'],self.stab_data[self.mutant]):
            exp_diff[pH]=ddG_mut-ddG_wt
        #
        # Make sure that the difference starts at 0 at pH 1.5
        #
        diffs=[]
        pHs=sorted(exp_diff.keys())
        offset=exp_diff[pHs[0]]
        for pH in pHs:
            diffs.append(math.pow(exp_diff[pH]-offset-calc_diff[pH],2))
        RMSD=math.sqrt(sum(diffs)/float(len(diffs)))
        return RMSD
        
    #
    # ------
    #
    
    def fit_mutant_pKa(self,mutant,pKa_values=None,stabdata=None):
        """Use the mutant data to determine the pKa values of the unfolded form"""
        self.mutant=mutant
        self.stab_data=stabdata
        import copy
        self.pKa_values=copy.deepcopy(pKa_values)
        #
        fit_variables=[]
        fit_names=[]
        for group in sorted(self.TGs):
            for state in ['F','U']:
                pka=self.pKa_values[group][state]
                if pka=='fit':
                    import pKaTool.pKadata
                    group_type=group.split(':')[-1]
                    pka=pKaTool.pKadata.modelpKas[group_type]
                    if state=='U':
                        pka=self.pKa_values[group]['F']
                    #
                    # Only fit this variable if it's removed in the mutant
                    #
                    if self.group_removed(self.mutant,group):
                        fit_variables.append(pka)
                        fit_names.append('%s %s' %(group,state))
                    else:
                        self.pKa_values[group][state]=pka 
        #
        # Print info and fit
        #
        print '---------------------------------'
        print 'Fitting %d pKa values for mutant %10s' %(len(fit_variables),mutant)
        print 'Groups being fit and their starting values'
        stats={}
        for group,value in zip(fit_names,fit_variables):
            print '%15s %5.1f' %(group,value)
            group=group.split()[0]
            stats[group]=[]
        stats['score']=[]
        print 'Doing %3d fits to get reliable errors' %self.options.numfits
        for numfit in range(self.options.numfits):
            import scipy.optimize
            solution,score,numits,func_calls,warnflag=scipy.optimize.fmin(self.get_difference_mutant,fit_variables,full_output=True)
            print 'Solution for fit #%3d:' %(numfit)
            
            print '%15s %5s %5s' %('Group','start','fitted')
            for group,value,value2 in zip(fit_names,fit_variables,solution):
                print '%15s %5.1f %5.1f' %(group,value,value2)
                #result[group]=value2
                group=group.split()[0]
                stats[group].append(value2)
            stats['score'].append(score)
            #
            # Randomize the starting values
            #
            import random
            new_start=[]
            for val in fit_variables:
                new_start.append(random.uniform(0,10.0))
            fit_variables=new_start[:]
                
        #
        # Calculate average and std. deviation
        #
        result={}
        for group in stats.keys():
            result[group]=average(stats[group])
        #
        # Show the fit of the difference in pH-stability profiles
        #
        print result.keys()
        fit_count=0
        for group in sorted(self.TGs):
            for state in ['F','U']:
                pka=self.pKa_values[group][state]
                if pka=='fit' and self.group_removed(self.mutant,group):
                    pka=result[group][0] #solution[fit_count]
                    self.pKa_values[group][state]=pka 
                    fit_count=fit_count+1
        wt_stabprof,group_contrib=self.get_pH_stab_profile('wt',[])
        mut_stabprof,group_contrib=self.get_pH_stab_profile(self.mutant,[])   
        #
        # Subtract the two calculated profiles
        #
        xs=[]
        ys=[]
        for pH in sorted(wt_stabprof.keys()):
            diff=mut_stabprof[pH]-wt_stabprof[pH]
            xs.append(pH)
            ys.append(diff)
        import pylab
        pylab.clf()
        pylab.plot(xs,ys,'ro-',label='Calculated')
        #
        xs=[]
        ys=[]
        offset=self.stab_data[self.mutant][0]-self.stab_data['wt'][0]
        for pH,ddG_wt,ddG_mut in zip(self.stab_data['ph'],self.stab_data['wt'],self.stab_data[self.mutant]):
            exp_mut_diff=ddG_mut-ddG_wt-offset
            xs.append(float(pH))
            ys.append(float(exp_mut_diff))
        pylab.plot(xs,ys,'bo-',label='Experimental')
        pylab.xlabel('pH')
        pylab.ylabel('ddG (kJ/mol)')
        pylab.title(self.mutant.upper())
        pylab.savefig('fits/%s.png' %self.mutant.upper(),dpi=300)
        return result 
    #
    # ------
    #
    
    def group_removed(self,mutation,group):
        """Find out if a group was removed in the mutation"""
        group_removed=False
        muts=mutation.split('+')
        for mutation in muts:
            number=int(mutation[1:-1])
            group_number=group.split(':')[1]
            if number==int(group_number):
                group_removed=True
                break
        return group_removed
    
    #
    # -------
    #
    
    def fit_pKa_values(self,stabdata,pKa_values):
        """Fit the pH-activity profiles"""
        self.stab_data=stabdata
        self.pKa_values=pKa_values.copy()
        #
        fit_variables=[]
        fit_names=[]
        for group in sorted(self.TGs):
            for state in ['F','U']:
                pka=pKa_values[group][state]
                if pka=='fit':
                    import pKaTool.pKadata
                    group_type=group.split(':')[-1]
                    pka=pKaTool.pKadata.modelpKas[group_type]
                    if state=='U':
                        pka=self.pKa_values[group]['F']
                    fit_variables.append(pka)
                    fit_names.append('%s %s' %(group,state))
        #
        # Print info and optimize
        #
        results={}
        print 'Groups being fit and their starting values'
        for group,value in zip(fit_names,fit_variables):
            print '%15s %5.1f' %(group,value)
            group=group.split()[0]
            results[group]=[]
        results['score']=[]
        #
        # Do the fits
        #
        for numfit in range(self.options.numfits):
            import scipy.optimize
            solution,score,numits,func_calls,warnflag=scipy.optimize.fmin(self.get_difference_all,fit_variables,full_output=True)
            print 'Solution for fit #%3d:' %(numfit)
            print 'Groups and fitted variables'
            for group,value,value2 in zip(fit_names,fit_variables,solution):
                print '%15s %5.1f %5.1f' %(group,value,value2)
                group=group.split()[0]
                results[group].append(value2)
            results['score'].append(score)
            #
            # Randomize the starting values
            #
            print fit_variables
            import random
            new_start=[]
            for var in fit_variables:
                new_start.append(random.uniform(0,14))
            fit_variables=new_start[:]
        #
        # Calculate average and std. deviation
        #
        result={}
        for group in results.keys():
            results[group]=average(results[group])
            print '%13s avg: %7.1f SD: %7.1f' %(group,results[group][0],results[group][2])
        #
        # Show the currently fitted pH-stability profile
        #
        fit_count=0
        for group in sorted(self.TGs):
            for state in ['F','U']:
                pka=pKa_values[group][state]
                if pka=='fit':
                    pka=results[group][0] #solution[fit_count]
                    pKa_values[group][state]=pka 
                    fit_count=fit_count+1
        self.pKa_values=pKa_values.copy()
        stabprof,group_contrib=self.get_pH_stab_profile('wt',[])
        #
        xs=[]
        ys=[]
        for pH in sorted(stabprof.keys()):
            xs.append(pH)
            ys.append(stabprof[pH])
        import pylab
        pylab.clf()
        pylab.plot(xs,ys,'r.-')
        #
        # Plot the wild type pH-activity profile
        #
        xs=[]
        ys=[]
        for pH,ddG in zip(stabdata['ph'],stabdata['wt']):
            xs.append(pH)
            ys.append(ddG)
        pylab.plot(xs,ys,'bo-',label='wt exp')
        pylab.legend()
        pylab.xlabel('pH')
        pylab.ylabel('ddGfold (kJ/mol)')
        pylab.show()    
        return     
        
#
# ------
#

if __name__=='__main__':
    print
    print 'Fitting multiple pH-stability profiles'
    print
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]',version='%prog 1.0')
    parser.add_option('-p','--pdb',dest='pdbfile',action='store',type='string',default='2lzt.pka.pdb',
                      help='The PDB file. Default: %default')
    parser.add_option('-x','--expdata',dest='expdata',action='store',type='string',default='ddG_ref.csv',help='The CSV file containing the stability data')
    parser.add_option('-k','--pkadata',dest='pkadata',action='store',type='string',default='pKa_values.csv',help='The CSV file containing the pKa data')
    parser.add_option('-d','--pHstep',dest='pH_step',action='store',type='float',default=0.1,help='pH integration step. Default: %default')
    parser.add_option('--pHstart',dest='pH_start',action='store',type='float',default=1.4,help='pH start (lowest value). Default: %default')
    

    parser.add_option('-u','--unfold',dest='ddGunfold',action='store_true',default=False,help='Stability values are ddG(unfold) values. Default: %default')
    parser.add_option('--kcal',dest='kcal',action='store_true',default=False,help='Stability values are in kcal/mol. Default: %default')
    #
    # Execution options
    #
    parser.add_option('--test',dest='test',action='store_true',default=False,help='Test charge calculation routines. Default: %default')
    parser.add_option('--singlefits',dest='singlefits',action='store_true',default=False,help='Fit unfolded pKa values from single mutants. Default: %default')
    parser.add_option('--doublemuts',dest='doublemuts',action='store_true',default=False,help='Do double mutant cycle analysis. Default: %default') 
    parser.add_option('--numfits',dest='numfits',action='store',type='int',default=100,help='Number of refits to perform to get standard deviation. Default: %default')    
    
    (options, args) = parser.parse_args()
    # 
    # Fit with the current parameters
    #
    X=fit_stability(options.pdbfile,options.expdata,options.pkadata,options)
