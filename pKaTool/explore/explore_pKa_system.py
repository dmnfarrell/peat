#!/usr/bin/env python
# %Id$
"""Routines for exploring the electrostatic sub-space of a particular system"""
#
# Copyright (C) Jens Erik Nielsen 2005 University College Dublin
#
import sys
import numpy as np

class explore_system:
    
    def __init__(self,options):
        """Set up the arrays for a full exploration of the behaviour of a system of titratable groups"""
        self.options=options
        self.number_of_groups=options.numgroups
        self.pka_resolution=options.step 
        #
        # -----
        #
        import math
        self.intpka_range=np.arange(7.0,13.0,self.pka_resolution)
        self.intene_range=np.arange(0,15.0,self.pka_resolution*math.log(10.0)) # in pKa units
        self.acid_base_range=['Acid','Base']
        #
        perms=[]
        #
        # Number of free parameters in the system are (N*N-1)/2 interaction energies
        # N intrinsic pKa values
        # and N acid-base specifications
        #
        self.num_levels=self.number_of_groups*(self.number_of_groups-1)/2+self.number_of_groups+self.number_of_groups
        self.levels=range(self.num_levels)
        print 'number of parameters',self.num_levels
        #
        # Add all the permuters 
        #
        counters=[]
        self.values=[]
        self.max_count=[]
        for l in self.levels:
            counters.append(0)
        #
        # Intrinsic pKa values and Acid-Base
        #
        for l in range(self.number_of_groups):
            self.max_count.append(len(self.acid_base_range))
            self.values.append(self.acid_base_range)
            self.max_count.append(len(self.intpka_range))
            self.values.append(self.intpka_range)
            
        #
        # Interaction energies
        #
        for l in range(self.num_levels-2*self.number_of_groups):
            self.max_count.append(len(self.intene_range))
            self.values.append(self.intene_range)
        #
        # How many permutations in total?
        #
        self.tot_perms=1
        for l in self.levels:
            self.tot_perms=self.tot_perms*self.max_count[l]
        print 'Total number of permutations: %.2e' %self.tot_perms
        #
        # Set some pKa calculation parameters
        #
        self.pHstep=0.1
        self.pHstart=0.1
        self.pHend=20.0
        return

    #
    # -----
    #

    def construct_permutations(self):
        """Construct the specification for all possible permutations specified in __init__"""
        #
        # Make pickled files for each configuration
        #
        import large_dataset_manager, os
        self.options.dir=os.path.join(os.getcwd(),options.dir)
        if not os.path.isdir(self.options.dir):
            os.mkdir(options.dir)
        self.Data=large_dataset_manager.data_manager(self.options.dir)
        print 'Constructing system permutations'
        print '000.00 perc done',
        perm=0
        last_perc=0
        #
        # Init the counters
        #
        import time
        start=time.time()
        #
        # Start the loop
        #
        while perm<self.tot_perms:
            #
            # Get the specification
            #
            spec=self.get_specs(perm)
            #print spec
            self.Data.put(perm,{'configuration':spec,'numgroups':self.number_of_groups})
            #
            # Get the time to completion
            #
            now=time.time()
            time_passed=now-start
            perc=float(perm)/float(self.tot_perms)*100.0+0.00001 # To avoid zerodivision error
            total_time=time_passed/(perc/100.0)
            time_to_completion=total_time-time_passed
            text='%5.2f perc done. Complete in %.2f minutes (= %.2f hours or %.2f days)          ' %(perc,time_to_completion/60.0,time_to_completion/60.0/60.0,
                                                                                           time_to_completion/60.0/60.0/24.0)
            text=(len(text)+4)*'\b'+text
            print text,
            #
            # Calculate titration curves for this permutation
            #
            new_data=self.do_one_pKa_calc(perm)
            self.Data.put(perm,new_data)
            perm=perm+1
        return


    #
    # ---------
    #

    def get_specs(self,perm_number):
        """Given a permutation number get the configuration"""
        #
        base_number=[]
        base=1
        for level in range(self.num_levels-1,-1,-1):
            #print level,'base',base
            base_number.append(base)
            base=base*self.max_count[level]
        base_number.reverse()
        #
        # Find the configuration
        #
        conf=[]
        for base in base_number:
            frac=float(perm_number)/float(base)
            if frac>=1.0:
                num=int(frac)
                conf.append(num)
                perm_number=perm_number-num*base
            else:
                conf.append(0)
        #
        # Fill in the real values
        #
        configuration=[]
        count=0
        for index in conf:
            value=self.values[count][index]
            import types
            if type(value) is types.IntType:
                value=float(value)
            configuration.append(value)
            count=count+1
        return configuration

    #
    # ------
    #

    def calculate_system(self):
        """Calculate titration curves for all permutations of the system"""
        print 'Solving system'
        count=0
        import time
        start=time.time()
        for perm in self.Data.index.keys():
            count=count+1
            perc=float(count)/float(self.tot_perms)*100.0
            #
            #
            now=time.time()
            time_passed=now-start
            total_time=time_passed/(perc/100.0)
            time_to_completion=total_time-time_passed
            text='%5.2f perc done. Complete in %.2f minutes' %(perc,time_to_completion/60.0)
            text=(len(text)+4)*'\b'+text
            print text,
            new_data=self.do_one_pKa_calc(perm)
            self.Data.put(perm,new_data)
        return
    

    #
    # -----
    #

    def do_one_pKa_calc(self,permutation_name):
        #
        # Set up the pKa calculation
        #
        old_data=self.Data(permutation_name)
        if old_data.has_key('pKas'):
            return old_data
        #
        # Haven't done this one yet
        #
        configuration=old_data['configuration']
        num_groups=old_data['numgroups']
        import pKa
        #
        # Fill instance with data
        #
        import pKaTool.pKa_calc, math
        X=pKaTool.pKa_calc.Boltzmann()
        matrix_dummy=self.setup_system(configuration,num_groups,X)
        #
        # Set the pKa value variables
        # 
        X.groups=X.intrinsic_pKa.keys()
        X.groups.sort()
        #
        # Print the setup for verification
        #
        #print X.groups
        #print X.intrinsic_pKa
        #print X.intene
        pKa_values=X._calc_pKas(0,self.pHstep,self.pHstart,self.pHend,1)
        #print pKa_values
        #print X.prot_states
        #print '=========================\n'
        old_data['pKas']=pKa_values
        old_data['titcurvs']=X.prot_states
        old_data['All_states']=X.all_states
        #
        return old_data

    #
    # -----
    #

    def setup_system(self,configuration,num_groups,X):
        #
        # Create the description of the system
        #
        self.names={}
        self.intpka={}
        self.type={}
        count=0
        for position in range(0,2*num_groups,2):
            #
            # Set the acid-base flag
            #
            a_b=configuration[position]
            intpka=configuration[position+1]
            #
            count=count+1
            number=str(count).zfill(4)
            name=':%s:' %number
            if a_b=='Acid':
                name=name+'ASP'
                self.type[count]=-1
            elif a_b=='Base':
                name=name+'ARG'
                self.type[count]=1
            else:
                raise 'Neither acid nor base. Bug!'
            #
            # set the name and the id
            #
            self.names[count]=name
            self.intpka[count]=intpka
        #
        # Store the interaction energies
        #
        intenes=configuration[2*num_groups:]
        #
        # Add all data
        #
        matrix={}
        X.intene={}
        X.intrinsic_pKa={}
        int_ene_count=0
        for group in self.names.keys():
            #
            # Set everything
            #
            name1=self.names[group]
            # Int pKa
            intpka=self.intpka[group]
            X.intrinsic_pKa[name1]=intpka
            type=self.type[group]
            #
            # Set the interaction energies
            #
            if not X.intene.has_key(name1):
                X.intene[name1]={}
            for group2 in self.names.keys():
                type2=self.type[group2]
                name2=self.names[group2]
                #
                # We do the symmetry immediately
                #
                if not X.intene.has_key(name2):
                    X.intene[name2]={}
                #
                # Now do the assignments
                #
                if not X.intene[name1].has_key(name2):
                    if name1==name2:
                        X.intene[name1][name2]=0.0
                        X.intene[name2][name1]=0.0
                    elif type==type2:
                        X.intene[name1][name2]=intenes[int_ene_count]
                        X.intene[name2][name1]=intenes[int_ene_count]
                        int_ene_count=int_ene_count+1
                    else:
                        X.intene[name1][name2]=-intenes[int_ene_count]
                        X.intene[name2][name1]=-intenes[int_ene_count]
                        int_ene_count=int_ene_count+1
        return X.intene



if __name__=="__main__":
    print
    print 'Explore the parameter space for a system of titratable groups'
    print
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]',version='%prog 1.0')
    parser.add_option('-n','--groupnum',dest='numgroups',type='int',action='store',default=3,help='Number of titratable groups. Default: %default')
    parser.add_option('-s','--step',dest='step',action='store',type='float',default=1.0,help='Stepsize for each parameter in pKa units. Default: %default')
    parser.add_option('-d','--dir',dest='dir',action='store',type='string',default='resultdir',help='Directory where results will be stored. Default: %default')
    (options, args) = parser.parse_args()
    #
    # Call the class
    #
    X=explore_system(options)
    X.construct_permutations()
    #X.calculate_system()
    
        
    
