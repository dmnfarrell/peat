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

class sensitivity:

    def __init__(self,pIO,selected_groups,pKa_MCsteps=20000):
        """
        # Setup the system
        """
        self.wtpIO=pIO
        self.pHstart=0.0
        self.pHend=20.0
        self.pHstep=0.1
        self.pKa_MCsteps=pKa_MCsteps
        self.selected_groups=selected_groups
        #
        import pKa_calc
        self.PK=pKa_calc.Monte_Carlo_CPP()
        self.set_wt_energies()
        #
        # Get the pKa values
        #
        print 'Calculating wild type pKa values'
        pKa_values,prot_states=self.PK.calc_pKas(mcsteps=self.pKa_MCsteps,
                                            phstep=self.pHstep,
                                            phstart=self.pHstart,
                                            phend=self.pHend,
                                            exp_pHs=[],
                                            verbose=1)
        self.wt_pKas=pKa_values.copy()
        self.wt_prot_states=prot_states.copy()
        #
        # Get the HH fit for each group
        #
        import pKarun.pKa_general
        self.wt_HHfits={}
        for group in self.wt_prot_states.keys():
            ANA=pKarun.pKa_general.pKanalyse()
            tc=self.wt_prot_states[group]
            del tc['pKa']
            solution,sq=ANA.fit_to_henderson(tc)
            self.wt_HHfits[group]=[solution,sq]
            #print group,abs(solution[1]),abs(solution[0]),sq
        return 
        
    #
    # ----
    #
    
    def set_wt_energies(self):
        """Setup the MC_CPP system for the wild type"""
        # 
        # Get the intrinsic pKas
        #
        intpkas={}
        backgr={}
        desolv={}
        for group in self.wtpIO.pka.keys():
            if group in self.selected_groups:
                intpkas[group]=self.wtpIO.pka[group]['intpka']
                backgr[group]=self.wtpIO.pka[group]['backgr']
                desolv[group]=self.wtpIO.pka[group]['desolv']
        self.PK.groups=intpkas.keys()
        self.PK.groups.sort()
        self.PK.desolv=desolv.copy()
        self.PK.backgr=backgr.copy()
        #
        matrix=self.wtpIO.read_matrix()
        prune_matrix={}
        for group1 in self.selected_groups:
            prune_matrix[group1]={}
            for group2 in self.selected_groups:
                prune_matrix[group1][group2]=matrix[group1][group2]
        self.PK.matrix=prune_matrix.copy()
        return
    
    #
    # ------
    #
    
    def random_sensitivity_analysis(self,monitor_groups,MCsteps):
        """Randomly vary the setup"""
        #
        # Start looping for all intrinsic pKa values
        #
        num_groups=len(self.selected_groups)
        self.do_states=MCsteps
        self.stab_test_on=1
        #
        self.curv_diffs=[]
        self.counter=0
        import time
        self.starttime=time.time()
        self.intpka_var=0.5
        self.intene_var=20
        diffs=[]
        while self.counter<self.do_states:
            for group in self.selected_groups:
                #
                # Change intrinsic pKa
                #
                # Vary as fraction of desolv+backgr (=intpka-modelpk)
                #
                intpka_var=self.intpka_var/100.0
                #print self.wt_pKas
                this_var=(self.wtpIO.pka[group]['intpka']-self.wtpIO.pka[group]['modelpK'])*intpka_var
                #
                import random
                change=random.uniform(-this_var,this_var)
                #
                # Now change the interaction energies for this group
                #
                intene_var=self.intene_var/100.0
                groups2=self.PK.matrix[group].keys()
                groups2.sort()
                for igroup in groups2:
                    change=random.uniform( 1.00-intene_var,1.00+intene_var)
                    ene=self.PK.matrix[group][igroup][0]-self.PK.matrix[group][igroup][1]-self.PK.matrix[group][igroup][2]+self.PK.matrix[group][igroup][3]
                    self.PK.matrix[group][igroup]=[change*ene,0.0,0.0,0.]
            diffs.append(self.measure_change(monitor_groups))
            self.set_wt_energies() # Restore the wt state
        #
        # normal mode
        #
        self.stab_test_on=1
        return diffs
        
    #
    # -----
    #

    def measure_change(self,monitor_groups,reference=None,absolut=False):
        """Measure the difference in the titration curve as compared to the wt situration"""
        #
        # Get the reference
        #
        if not reference:
            reference=self.wt_prot_states
        #
        new_pKa_values,new_prot_states=self.PK.calc_pKas(mcsteps=self.pKa_MCsteps,
                                            phstep=self.pHstep,
                                            phstart=self.pHstart,
                                            phend=self.pHend,
                                            exp_pHs=[],
                                            verbose=1)
        #
        this_diff=0.0
        self.counter=self.counter+1
        import time
        timenow=time.time()
        time_passed=timenow-self.starttime
        fraction_completed=float(self.counter)/float(self.do_states)
        time_to_finish=(1.0-fraction_completed)/(fraction_completed/time_passed)
        txt='%d seconds' %(time_to_finish)
        if time_to_finish/60.0/60.0/24.0>365:
            txt='%5.2f years' %(time_to_finish/60.0/60.0/24.0/365.0)  
        elif time_to_finish/60.0/60.0>24:
            txt='%5.2f days' %(time_to_finish/60.0/60.0/24.0)    
        elif time_to_finish/60.0>60:
            txt='%5.2f hours' %(time_to_finish/60.0/60.0)
        elif time_to_finish>60:
            txt='%5.2f minutes' %(time_to_finish/60.0)
        print '%.2f %% done (%5.2e of %5.2e). Finish in %s' %(fraction_completed*100.0,self.counter,self.do_states,txt)
        import sys
        sys.stdout.flush()
        #
        # Measure the change
        #
        for group_test in monitor_groups:
            pH_values=new_prot_states[group_test].keys()
            if 'pKa' in pH_values:
                pH_values.remove('pKa')
            for pH in pH_values:
                new_crg=new_prot_states[group_test][pH]
                old_crg=reference[group_test][pH]
                if abs(new_crg)>1.0 or abs(old_crg)>1.0:
                    print 'Charge outside range',new_crg,old_crg
                    stop
                this_diff=this_diff+abs(new_crg-old_crg)*self.pHstep
        print 'Difference is',this_diff, self.pHstep
        if abs(this_diff)>10.0:
            print this_diff
            print 'Something is wrong'
            stop
        return this_diff

#
# -----
#

def average(l):
    """Calculate the average, variance and standard deviation"""
    import math
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

#
# ----
#

def perform_ID(PD_pKa,Nuc_pKa):
    """Can we identify the proton donor"""
    if not PD_pKa is None and not Nuc_pKa is None:
        if PD_pKa>5.0 and PD_pKa-Nuc_pKa>=1.5:
            ID_result='Correct'
        elif  Nuc_pKa>5.0 and Nuc_pKa-PD_pKa>=1.5:
            ID_result='Incorrect'
        else:
            ID_result='Inconclusive'
    else:
        ID_result='Missing pKa values!'
        raw_input('Did you fix it?')
        #raise Exception('Fix this one')
    return ID_result

#
# -----
#

def print_result(calc,results,PD,Nuc,resultdir=None):
    """Print the result"""
    groups=results[calc].keys()
    groups.sort()
    print 'Calculation',calc
    if results[calc]=={}:
        print 'Calculation not finisihed'
        return None,None
    print '%10s\t%5s\t%22s\t%7s\t%22s\t%5s' %('Group','pKa','HH fit','subpKa','subHH fit','sens(MC)')
    SD={}
    csvline='%s' %calc
    for group in [PD,Nuc]:
        if not group:
            continue
        SD[group]={}
        
        this=results[calc][group]
        sens_single=''
        sens_MC=''
        if this.has_key('sensitivity'):
            #sens_single='%5.2f' %(this['sensitivity']['single'][1])
            sens_MC='%5.2f' %(this['sensitivity']['random'][1])
        pka=this['pKa']
        if pka:
            pka='%5.2f' %pka
        subpka='-'
        if this.has_key('subpKa'):
            subpka=this['subpKa']
            if subpka<-1000.0:
                subpka='ND '
            elif subpka:
                subpka='%5.2f' %subpka

        #
        subHHfit='-'
        if this.has_key('subHHfit'):
            HHfit=this['subHHfit']
            if type(HHfit) is type([]):
                HHfit=HHfit[0][0]
            if type(HHfit) is type(' x'):
                subHHfit='-'
            else:
                subHHfit='%5.2f' %(abs(HHfit))
        HHfit='-'
        if this.has_key('HHfit'):
            HHfit='%5.2f' %(abs(this['HHfit'][0][0]))
        txt= '%10s\t%5s\t%22s\t%7s\t%22s\t%8s' %(group,get_pKa(group,resultdir,calc),HHfit,subpka,subHHfit,sens_MC)
        print txt,
        csvline=csvline+',%10s,%5s,%22s,%7s,%22s,%8s' %(group,get_pKa(group,resultdir,calc),HHfit,subpka,subHHfit,sens_MC)
        if this.has_key('stats'):
            stat=this['stats']
            txt= '%5.2f (%5.2f)' %(stat['random'][0],stat['random'][2])
            print txt,
            csvline=csvline+',%5.2f, (%5.2f)' %(stat['random'][0],stat['random'][2])
            SD[group]['random']=stat['random'][2]    
        print
    #
    # Can we identify the PD?
    #
    import string
    print 'PD-NUC',PD,Nuc
    for value in ['pKa']: #['pKa','subpKa']:
        text='%s correct PD identification: ' %value
        ID_result=None
        print text,
        #
        # Get pKa values and perform ID
        #   
        PD_pKa=get_pKa(PD,resultdir,calc)
        Nuc_pKa=get_pKa(Nuc,resultdir,calc)
        ID_result=perform_ID(PD_pKa,Nuc_pKa)
        #print 'pKa value',PD_pKa,Nuc_pKa
        #
        # Qualify the proton donor ID
        #
        print ID_result+',',
        if ID_result!='Missing pKa values!':
            success=0
            for test in ['random']:
                for PD_ch in [-SD[PD][test],SD[PD][test]]:
                    for Nuc_ch in [-SD[Nuc][test],SD[Nuc][test]]:
                        #print perform_ID(PD_pKa+PD_ch,Nuc_pKa+Nuc_ch),ID_result
                        if perform_ID(PD_pKa+PD_ch,Nuc_pKa+Nuc_ch)==ID_result:
                            success=success+1
            qual=['Very low','Low','Medium','High','Solid']
            txt= 'Confidence: %s (%d/4)' %(qual[success],success)
            print txt,
            csvline=csvline+',%s,Confidence: %s (%d/4)' %(ID_result,qual[success],success)
        elif PD and Nuc:
            stop
            print
        else:
            print
        print
    print '---------------------'
    #
    # Convert this into a score. 1.0 is correct/solid, -1.0 is incorrect/solid
    #
    if ID_result=='Correct':
        score=success*1.0/4
    elif ID_result=='Incorrect':
        score=success*-0.25
    elif ID_result=='Missing pKa values!':
        score=None
    else:
        score=0.0
    return score,csvline

#
# -----
#
def get_intpKa(PD,resultdir,calc):
    return get_pKa(PD,resultdir,calc,intpka=True)

def get_pKa(PD,resultdir,calc,intpka=False):
    """Get the pKa value even if the group doesn't titrate"""
    import os
    pdbfile=os.path.join(resultdir,calc)
    if not PD:
        return None
    if PD.lower()=='mutated':
        return None
    print 'PD',PD
    import os
    if len(PD.split(':'))!=3:
        import Protool
        P=Protool.structureIO()
        P.readpdb(pdbfile)
        PD='%s:%s' %(PD,P.resname(PD))
        print 'Found resname',P.resname(PD)
        import pKarun.pKa_utility_functions as util
        if not util.is_titratable(P.resname(PD)):
            print 'NOT titratable',P.resname(PD)
            return None
    #
    # Read the pKa values
    #
    import pKaTool.pKaIO
    X=pKaTool.pKaIO.pKaIO(pdbfile)
    if not X.assess_status():
        return None
    tcs=X.readtitcurv()
    pkas=X.readpka()
    #
    # Do we want the real pka value or just the intrinsic pKa?
    #
    if intpka:      
        return pkas[PD]['intpka']
    #
    # Da real pKa value
    #
    pKa_value=pkas[PD]['pKa']
    if pKa_value:
        PD_pKa=pKa_value
    else:
        tc=tcs[PD]
        pos=True
        neu=True
        neg=True
        del tc['pKa']
        pHs=tc.keys()
        pHs.sort()
        #print PD
        for pH in pHs:
            charge=tc[pH]
            #print pH,charge
            if charge>-0.40:
                neg=False
            if charge<-0.4 or charge>0.4:
                neu=False
            if charge<0.4:
                pos=False
        if pos:
            PD_pKa=20.0
        elif neu: 
            if pkas[PD]['acidbase']==-1:
                PD_pKa=20.0
            else:
                PD_pKa=0.0
        elif neg:
            PD_pKa=0.0
        else:
            PD_pKa=None
            raise Exception('This really happens')
    return PD_pKa


#
# ----
#
def find_interacting_groups(matrix,selected_groups,int_cutoff=1.0):
        """Select groups that interact strongly with the selected groups"""
        append_groups={}
        #
        # Loop over all selected groups and check for interaction energies
        #
        interaction_energy_cutoff=int_cutoff
        for group_sel in selected_groups:
            #
            # Loop over all groups
            #
            for group in matrix[group_sel].keys():
                record=matrix[group_sel][group]
                intene=record[0]-record[1]-record[2]+record[3]
                if abs(intene)>=interaction_energy_cutoff:
                    append_groups[group]=1
        #
        # Found all, now add them
        #
        for group in append_groups.keys():
            if not group in selected_groups:
                selected_groups.append(group)
        return

#
# -----
#
    
def sensitivity_analysis(pdbfile,center_residues,cut_off,monitor_residues,MCstates):
    """Perform Sensitivity analsysis using pKaTool"""
    data={}
    #
    # Start pkaTool
    #
    import pKaTool.pKaIO
    X=pKaTool.pKaIO.pKaIO(pdbfile)
    if not X.assess_status():
        return None,None,None
    X.readtitcurv()
    pkavals=X.readpka()
    matrix=X.read_matrix()
    #
    # Get HH fits of from full pKa calculation
    #
    
    #
    # Find strongly coupled groups
    #
    cut_off_adjustment=0.0
    num_groups=2000
    while num_groups>10:
        selected_groups=[]
        #
        # Add the residues at the center
        #
        for residue in center_residues:
            selected_groups.append(residue)
        #
        # Select them
        #
        int_cutoff=cut_off+cut_off_adjustment
        find_interacting_groups(matrix,selected_groups,int_cutoff)
        #
        # Remove Tyrs with high pKa values (i.e. pKa of -999) no check for
        # completely negatively charged TYR is done..
        #
        remove=[]
        for group in selected_groups:
            data[group]={'pKa':pkavals[group]['pKa']}
            if group.find('TYR')!=-1 and (pkavals[group]['pKa'] is None or pkavals[group]['pKa']>13):
                remove.append(group)
        for rem in remove:
            selected_groups.remove(rem)
            del data[rem]
        #
        # How many groups?
        #
        num_groups=len(selected_groups)
        print 'Number of selected groups is %d with %5.3f kT cut_off' %(num_groups,cut_off+cut_off_adjustment)
        cut_off_adjustment=cut_off_adjustment+0.1
    #
    # Get the HHfits from the full pKa calculation
    #
    import pKarun.pKa_general
    titcurvs=X.read_titration_curve()
    for group in selected_groups:
        ANA=pKarun.pKa_general.pKanalyse()
        tc=titcurvs[group]
        del tc['pKa']
        solution,sq=ANA.fit_to_henderson(tc)
        data[group]['HHfit']=[solution,sq]
    #
    # Instantiate sensitivity class
    #
    SENS=sensitivity(X,selected_groups)
    #
    # Get the pKa values in the reduced system
    #
    groups=SENS.wt_pKas.keys()
    groups.sort()
    for group in groups:
        data[group]['subpKa']=float(SENS.wt_pKas[group])
        data[group]['subHHfit']=SENS.wt_HHfits[group]
    #
    # Start the sensitivity test
    #
    results={}
    diffset={}
    tests={'random':SENS.random_sensitivity_analysis}
    for do_test in tests.keys():
        #
        # Do the test
        #
        sens_diffs=tests[do_test](monitor_groups=monitor_residues,MCsteps=MCstates)
        #
        # Get the differences
        #
        diffs=[]
        rdiffs=[]
        sum=0.0
        for diff in sens_diffs: #0 is the wt state
            diffs.append(abs(diff))
            rdiffs.append(diff)
            #print 'doff',diff
            sum=sum+abs(diff)
        avg_change=sum/float(len(diffs))
        max_change=diffs[0]
        print 'TEST type',do_test
        #print 'Total number of possible permutatinos: %6.2e' %X.pKasys.tot_perms
        print 'Max. change: %5.2f' %max_change
        print 'Average change in pKa value; %5.2f for %d states' %(avg_change,len(diffs))
        print 'Number of groups in sub_system',len(groups)
        print '--------------------'
        import sys
        sys.stdout.flush()
        results[do_test]=[max_change,avg_change]
        #
        diffset[do_test]=rdiffs
    return results,data,diffset

#
# -----
#

def main(options,args):
    """Do automatic sensitivity analysis with pKaTool"""
    #
    # Parse the arguments
    #
    file=args[0]
    center_residues=options.center_residues
    cut_off=options.cut_off
    monitor_residues=options.monitor_residues
    #
    # single file
    #
    if not options.many_files:
        return sensitivity_analysis(file,center_residues,cut_off,monitor_residues,options.MCstates)
    #
    # Multiple files
    #
    if options.filestructure=='dirs':
        fd=open(file)
        lines=fd.readlines()
        fd.close()
        #
        # Array for storing data
        #
        count=0
        for line in lines:
            results={}
            print 'Reading: %s' %(line.strip())
            #
            # skip comments
            #
            if line[0]=='#':
                continue
            #
            # Do some work for this line
            #
            import os
            data=line.split()
            ID=data[0]
            dirname=data[0]
            filename=data[1]
            filename=os.path.join(os.getcwd(),dirname,filename)
            resultfile=filename+'.sensresult_%d' %options.MCstates
            if os.path.isfile(resultfile):
                continue
            import pickle
            #
            # get the residue name
            #
            import Protool
            X=Protool.structureIO()
            X.readpdb(filename)
            X.RemoveALT()
            #
            # Get the Nucleophile and Proton donor residues
            #
            PD=data[2]
            Nuc=data[3]
            if len(data)>4:
                status=data[4]
            else:
                status='Yes'
            #
            c_res=[]
            if status=='Yes':
                if PD.lower()=='mutated' or PD.lower()==':mutated' or X.resname(PD) not in ['ASP','GLU','HIS','LYS','ARG']:
                    PD=None
                else:
                    c_res.append('%s:%s' %(PD,X.resname(PD)))
                    PD='%s:%s' %(PD,X.resname(PD))
                if Nuc.lower()=='mutated' or Nuc.lower()==':mutated' or X.resname(Nuc) not in ['ASP','GLU','HIS','LYS','ARG']:
                    Nuc=None
                else:
                    c_res.append('%s:%s' %(Nuc,X.resname(Nuc)))
                    Nuc='%s:%s' %(Nuc,X.resname(Nuc))
            else:
                print 'Skipping file',line
                continue
            #
            # Do the sens analysis
            #
            calc=data[1]
            results={}
            results[calc]={}
            for res in c_res:
                print 'Doing sensitivity analysis for:',filename,c_res,[res]
                import sys
                sys.stdout.flush()
                sens_results,data,diffs=sensitivity_analysis(filename,c_res,cut_off,[res],options.MCstates)
                if diffs is None:
                    continue
                #
                # Use a temp array for the large data chunks
                #
                temp_results={}
                for k in data.keys():
                    if not results[calc].has_key(k):
                        results[calc][k]=data[k]
                temp_results['sensitivity']=sens_results
                #
                # Calculate standard deviation on the pKa values
                #
                results[calc][res]['stats']={}
                pka=results[calc][res]['pKa']
                if not pka:
                    pka=0.0
                for diffset in diffs.keys():
                    ndiffs=[]
                    for diff in diffs[diffset]:
                        ndiffs.append(pka+diff)
                    avg,variance,sdev=average(ndiffs)
                    results[calc][res]['stats'][diffset]=[avg,variance,sdev]
            #
            # Save the results
            #
            fd=open(resultfile,'w')
            data=[calc,results,PD,Nuc]
            pickle.dump(data,fd)
            fd.close()
            #
            # Print the results
            #
            import os
            res_dir=os.path.join(os.getcwd(),dirname)
            print_result(calc,results,PD,Nuc,resultdir=res_dir)
            count=count+1
            print 'Finished %d of %d' %(count,len(lines[1:]))
            data=None
            results=None
    return

if __name__=='__main__':
    print
    print 'Perform calculated pKa value sensitivity analysis'
    print 'Jens Erik Nielsen, 2009'
    print
    import sys, os
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <file>',version='%prog 1.0')
    parser.add_option('-r',"--residues",type='string',dest='center_residues',action='append',
                      help='Residues to include in sub-system')
    parser.add_option('-s',"--MCstates",type='int',dest='MCstates',action='store',
                      help='Number of Monte Carlo steps for MC sensitivity analysis. Default: %default',default=1000)
    parser.add_option('-m',"--monitor",type='string',dest='monitor_residues',action='append',
                      help='Residues that will be monitored in the sensitivity analysis')
    parser.add_option('-c',"--cut_off",type='float',dest='cut_off',default=1.0,action='store',
                    help='Cut-off for including stronlgy interacting groups. Default: cut_off=%default',metavar='cut_off')
    parser.add_option('-n','--many_files',dest='many_files',action='store_true',default=False,
                    help='Loop over many files. Default: %default')
    parser.add_option('-f',"--filestructure",type='string',dest='filestructure',default='dirs',action='store',
                      help='Are the PDB files in individual dirs (dirs) or in the same dir (files)? Default: %default')
    (options, args) = parser.parse_args()
    if len(args)!=1:
        parser.error('You must specify a PDB file')
    #
    # Call main
    #
    main(options,args)
