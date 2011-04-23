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

import sys, os
import pKarun.pKarun_main as pKarun
from pKarun.pKa_utility_functions import *


import pKaIO

class recordstate:
    """Class for keeping track of states in f.ex. MC sampling"""

    def __init__(self):
        self.states=[]
        return

    def recordstate(self,state):
        self.states.append((state.copy()))
        return

    
#
# Dummy function for accessing pKarun routines
#
class pKacalc(pKarun.pKarun):

    def dummy_function(self):
        return

#
# ------
#



class Monte_Carlo(pKaIO.pKaIO):
    #
    # Monte Carlo algorithm and base class for the other pKa calc routines
    #

    def acid_base(self,group):
        #
        # Return 1 for a base, return -1 for an acid
        #
        import string
        return self.acidbase[string.split(group,':')[-1]]

    #
    # -----------------------
    #

    def get_modelpKa(self,group):
        #
        # Return the model pKa value for the group
        #
        import string 
        group=string.split(group,':')[-1]
        return self.modelpKas[group]
        
    #
    # -------------------------
    #

    def prepare_matrix(self):
        """
        # Prepare the matrix
        #
        # Precompute full term
        """
        self.intene={}
        for key in self.matrix.keys():
            self.intene[key]={}
            for key2 in self.matrix[key].keys():
                #print key,key2,self.matrix[key][key2]
                self.intene[key][key2]=  self.matrix[key][key2][0] \
                                        -self.matrix[key][key2][1] \
                                        -self.matrix[key][key2][2] \
                                        +self.matrix[key][key2][3]

        #
        # Make sure that the matrix is symmetrical
        #
        residues=self.intene.keys()
        residues.sort()
        for key in residues:
            for key2 in residues:
                if not self.intene[key].has_key(key2):
                    print 'prepare matrix failed'
                    print 'Matrix[%s] is missing key; %s' %(key,key2)
                    print key
                    print self.intene[key].keys()
                    raise Exception('Matrix[%s] is missing key; %s' %(key,key2))
                E12=self.intene[key][key2]
                E21=self.intene[key2][key]
                new_ene=min(E12,E21)
                self.intene[key][key2]=new_ene
                self.intene[key2][key]=new_ene
        return

    #
    # ----------------------
    #

    def calc_intpKas(self):
        #
        # Set a few constants
  
        #
        # Check that all dictionaries have been filled
        #
        self.groups=self.desolv.keys()
        if self.groups!=self.backgr.keys():
            print
            print 'Inconsistent desolv and backgr'
            print
            raise Exception('Error in Python Monte Carlo routine')
        #
        # Calculate the intrinsic pKa values
        #
        import math
        self.ln10=math.log(10)
        self.intrinsic_pKa={}
        self.groups.sort()
        for group in self.groups:
            resname=get_resname(group)
            self.intrinsic_pKa[group]=self.get_modelpKa(group)+ \
                                       -float(self.acid_base(group))*self.desolv[group]/self.ln10 + \
                                       float(self.acid_base(group))*self.backgr[group]/self.ln10
        return

    #
    # -----------------
    #
    
    def calc_pKas(self,mcsteps=2000,phstep=0.1,phstart=2.0,phend=14.0,verbose=1,complete_pka=None,exp_pHs=[]):
        """Calculate pKa values for the system"""
        #
        # Init
        #
        # Set a few constants
        #
        # KBol and Temp are equal to 1.0 everywhere in this derived class!
        #
        # Check that all dictionaries have been filled
        #
        self.groups=self.desolv.keys()
        self.groups.sort()
        b_groups=self.backgr.keys()
        b_groups.sort()
        m_groups=self.matrix.keys()
        m_groups.sort()
        #
        if self.groups!=b_groups or self.groups!=m_groups:
            print
            print 'Inconsistent desolv, backgr and matrix dictionaries'
            print
            ndes=len(self.desolv.keys())
            nback=len(self.backgr.keys())
            nmat=len(self.matrix.keys())
            print 'Groups in desolv: %3d, groups in backgr: %3d, groups in matrix: %3d \n' %(ndes,nback,nmat)
            groups=self.backgr.keys()+self.desolv.keys()+self.matrix.keys()
            g_dict={}
            for group in groups:
                g_dict[group]=1
            groups=g_dict.keys()
            groups.sort()
            for group in groups:
                has=['-','-','-']
                if self.desolv.has_key(group):
                    has[0]='+'
                if self.backgr.has_key(group):
                    has[1]='+'
                if self.matrix.has_key(group):
                    has[2]='+'
                print '%14s Desolv: %s, Backgr: %s, Matrix: %s' %(group,has[0],has[1],has[2])
                print 'Done with ',group
            print 'Totally done'
            import sys
            sys.stdout.flush()
            raise 'Error in Python Monte Carlo routine'
        #
        # Prepare the matrix
        #
        self.prepare_matrix()
        #
        # Calculate the intrinsic pKa values
        #
        self.calc_intpKas()
        #
        # Calculate pKa values
        #
        return self._calc_pKas(mcsteps,phstep,phstart,phend,verbose,complete_pka,exp_pHs=exp_pHs)

    #
    # ------
    #
        
    def _calc_pKas(self,mcsteps=200000,phstep=0.1,phstart=1.0,phend=20.0,verbose=1,complete_pka=None,exp_pHs=[]):
        """Calculate pKa values from intrinsic pKa values and interaction energies. No checking done"""
        #
        # KBol and Temp are equal to 1.0 everywhere
        #
        import math
        self.ln10=math.log(10)
        #
        # Start calculating protonation states
        #
        self.mcsteps=mcsteps
        pHs=range(int(phstart*100.0),int(phend*100.0),int(phstep*100.0))
        pHvalues=[]
        for pH in pHs:
            pHvalues.append(float(pH)/100.0)

        pHvalues.extend(exp_pHs)

        self.pHvalues=pHvalues
        self.prot_states_tmp={}
        #
        # Calculate protonation states at each pH value
        #
        self.all_states={}
        import copy
        for pH in pHvalues:
            tmp,states=self.calc_fractional_charge(pH)
            if tmp=={}:
                print pH
                print 'No result'
                raise Exception('I dont believe it')
            self.prot_states_tmp[pH]=copy.deepcopy(tmp)
            self.all_states[pH]=copy.deepcopy(states)
        if verbose>1:
            print
        #
        # Determine pKas
        #
        pkavalues=self.determine_pKa_values()
        #
        # Reformat the titration data
        #
        self.prot_states={}
        for group in self.groups:
            self.prot_states[group]={}
            for ph in pHvalues:
                self.prot_states[group][ph]=self.prot_states_tmp[ph][group]
            self.prot_states[group]['pKa']=pkavalues[group]
        self.prot_states_tmp=None
        #
        # ----------
        #
        self.pka={}
        for group in pkavalues.keys():
            self.pka[group]={'pKa':pkavalues[group]}
        if complete_pka:
            self.complete_pka()

        return pkavalues,self.prot_states

    #
    # --------------------------
    #

    def determine_pKa_values(self):
        """Determine pKa values as half-points of titration from titration data"""
        pkavalues={}
        pHvalues=self.pHvalues
        pHvalues.sort()
        for group in self.groups:
            pka=-99.9
            last_crg=self.prot_states_tmp[pHvalues[0]][group]
            phstep=float(pHvalues[1]-pHvalues[0])
            for ph in pHvalues:
                try:
                    crg=self.prot_states_tmp[ph][group]
                except:
                    grps=self.prot_states_tmp[ph].keys()
                    grps.sort()
                    print grps
                    print group
                    raise 'same error'
                #
                # ----
                #
                if crg<last_crg:
                    if self.acid_base(group)==1:
                        if crg<=0.5 and last_crg>0.5:
                            pka=(last_crg-0.5)/(last_crg-crg)*phstep+(ph-phstep)
                            break
                    else:
                        if crg<=-0.5 and last_crg>-0.5:
                            pka=(last_crg-(-0.5))/(last_crg-crg)*phstep+(ph-phstep)
                            break
                last_crg=crg
            if pka<-90.0:
                if self.acid_base(group)==1:
                    if last_crg>0.5:
                        pka=99.9
                    else:
                        pka=-99.9
                else:
                    if last_crg>-0.5:
                        pka=99.9
                    else:
                        pka=-99.9
            pkavalues[group]=pka
        return pkavalues
    
    #
    # --------------------------
    #

    def calc_fractional_charge(self,pH):
        #
        # Calculate the fractional charge for all residues
        # at this pH
        #
        # Get all the groups
        #
        # Define the Monte Carlo parameters
        #
        eqsteps=self.mcsteps/10
        #
        # Initialise the random number generator
        #
        import random, math
        rand=random.Random(198984)
        #
        # Initialise helper MC class
        #
        X=recordstate()
        #
        # Construct the starting state
        # State is a dictionary. For each group the value
        # is either 1 (charged) or 0 (neutral)
        #
        state={}
        old_cha={}
        for group in self.groups:
            state[group]=rand.randint(0,1)
        curE=self.get_energy(pH,state)
        #
        # Start the MC steps
        #
        for step in range(self.mcsteps):
            #
            # Construct the new state
            #
            change_group=rand.choice(self.groups)
            new_state=state.copy()
            new_state[change_group]=abs(new_state[change_group]-1)
            #
            # Calculate the new energy
            #
            newE=self.get_energy(pH,new_state)
            if newE<=curE:
                state=new_state.copy()
                curE=newE
            else:
                deltaE=newE-curE
                if deltaE<50.0:
                    if rand.random()<=math.exp(-deltaE):
                        state=new_state.copy()
                        curE=newE
                    else:
                        pass
                else:
                    pass
            if step>eqsteps:
                X.recordstate(state)
        #
        # Find the fractional degree of protonation
        #
        sumstate={}
        for group in self.groups:
            sum=0
            for state in X.states:
                sum=sum+state[group]
            sumstate[group]=float(sum)/float(len(X.states))
            if isacid(group):
                sumstate[group]=-sumstate[group]
        #
        # Done
        #
        return sumstate,{}
    
    #
    # --------------------
    #
    
    def get_energy(self,pH,state):
        #
        # Get the energy for this state
        #
        energy=0.0
        for group in self.groups:
            #
            # Add the effect of the non-titratable environment
            #
            if state[group]==1:
                energy=energy+float(self.acid_base(group))*self.ln10* \
                        (pH-self.intrinsic_pKa[group])
                #
                # Add the effect of all other titratable groups in the system
                #
                for group2 in self.groups:
                    if state[group2]==1 and group!=group2:
                        energy=energy+self.intene[group][group2]/2.0
                #
                # If we have a non-system groups to take into account, then we include
                # that here
                #
                if hasattr(self,'non_system_groups'):
                    energy=energy+self.non_system_groups[group][round(pH,1)]
        return energy

    #
    # -----
    #

    def complete_pka(self):
        """Complete the self.pka dictionary. Insert delec,ddesolv,dbackgr,dpka and intpka values"""
        for group in self.pka.keys():
            self.pka[group]['intpka']=self.intrinsic_pKa[group]
            self.pka[group]['modelpK']=self.get_modelpKa(group)
            self.pka[group]['desolv']=-float(self.acid_base(group))*self.desolv[group]/self.ln10 
            self.pka[group]['backgr']=float(self.acid_base(group))*self.backgr[group]/self.ln10
            self.pka[group]['delec']=self.pka[group]['pKa']-self.pka[group]['intpka']
        return

#
# -------------------------------
#



class Monte_Carlo_CPP(Monte_Carlo):
    """
    # C++ implementation of the Monte Carlo alg
    """

    def test(self):
        """Test if we can import the C++ module"""
        import pMC
        return

    #
    # --------
    # 
    
    def make_matrix_linear(self):
        """Change the matrix to linear form - this makes it easy to pass it to the C++ code"""
        linear=[]
        residues=self.intene.keys()
        residues.sort()
        for group1 in residues:
            for group2 in residues:
                linear.append(self.intene[group1][group2])
        return linear

    #
    # -----------------
    #

    def calc_pKas(self,mcsteps=200000,phstep=0.1,phstart=0.0,phend=14,verbose=1,complete_pka=None,exp_pHs=[],monitor_states=None):
        """
        # Calculate pKa values
        """
        # Init
        #
        # Set constants
        # KBol and Temp are equal to 1.0 everywhere in this derived class!
        #
        # Check that all dictionaries have been filled
        #
        self.groups=self.desolv.keys()
        self.groups.sort()
        b_groups=self.backgr.keys()
        b_groups.sort()
        m_groups=self.matrix.keys()
        m_groups.sort()
        if self.groups!=b_groups or self.groups!=m_groups:
            print
            print 'Inconsistent desolv, backgr and matrix dictionaries'
            print
            ndes=len(self.desolv.keys())
            nback=len(self.backgr.keys())
            nmat=len(self.matrix.keys())
            print 'Groups in desolv: %3d, groups in backgr: %3d, groups in matrix: %3d \n' %(ndes,nback,nmat)
            groups=self.backgr.keys()+self.desolv.keys()+self.matrix.keys()
            g_dict={}
            for group in groups:
                g_dict[group]=1
            groups=g_dict.keys()
            groups.sort()
            for group in groups:
                has=['-','-','-']
                if self.desolv.has_key(group):
                    has[0]='+'
                if self.backgr.has_key(group):
                    has[1]='+'
                if self.matrix.has_key(group):
                    has[2]='+'
                print '%14s Desolv: %s, Backgr: %s, Matrix: %s' %(group,has[0],has[1],has[2])
            print 'Totall done here'
            import sys
            sys.stdout.flush()
            raise Exception('Error in C++ Monte Carlo module')
        #
        # Prepare the matrix
        #
        self.prepare_matrix()

        #
        # Calculate the intrinsic pKa values
        #
        self.calc_intpKas()

        return self._calc_pKas(mcsteps,phstep,phstart,phend,verbose,complete_pka,monitor_states=monitor_states)

    def allok(self,list):
        for value in list:
            if not value and value!=0.0:
                return None

        return 1

    
    #
    # ----
    #
    
    def _calc_pKas(self,mcsteps=200000,phstep=0.1,phstart=1.0,phend=20.0,verbose=1,complete_pka=None,exp_pHs=[],
        monitor_groups=None,monitor_states=None):
        """Do the pKa calculation with the CPP module"""
        #
        # Do specific CPP setup
        #
        import time
        starttime=time.time()
        residues=self.intrinsic_pKa.keys()
        residues.sort()
        intpkas=[]
        acidbase=[]
        for residue in residues:
            intpkas.append(self.intrinsic_pKa[residue])
            acidbase.append(float(self.acid_base(residue)))
            #print residue,self.intrinsic_pKa[residue],self.acid_base(residue)
            # Works ok with mutations since self.acid_base works on the last residue
        linear_matrix=self.make_matrix_linear()
        #raw_input('ok?')
        #
        # Import the C++ module and run pKa calc from there
        #
        import time, math
        import pMC
        reload(pMC)
        use_MC=1
        FAST=pMC.MC(intpkas,linear_matrix,acidbase,use_MC)
        loadtime=time.time()
        FAST.set_MCsteps(int(mcsteps))
        #
        # If we monitor a CCPS then define it
        #
        if monitor_states:
            FAST.set_monitor_states(monitor_states)
        #
        #
        # Calculate pKa values
        #
        pKa_values=FAST.calc_pKas(phstart,phend,phstep)
        #
        # Get the CCPS population back
        #
        if monitor_states:
            self.pH=[]
            self.CCPS_pop=[]
            CCPS_pop=FAST._CCPS_population
            for count in range(0,len(CCPS_pop),2):
                #print CCPS_pop[count],CCPS_pop[count+1]
                self.pH.append(CCPS_pop[count])
                self.CCPS_pop.append(CCPS_pop[count+1])
        #
        # Put the pKa values back in a dictionary
        #
        pkavalues={}
        count=0
        for residue in residues:
            pkavalues[residue]=pKa_values[count]
            if pKa_values[count]<-100.0:
                pkavalues[residue]=-9999.9
            count=count+1
        #
        # Now get the titration curves
        #
        # pMC returns a specification of how many values it returns
        #
        pH_start=pKa_values[count]
        pH_step=pKa_values[count+1]
        pH_num=pKa_values[count+2]
        count=count+3
        #
        #
        #
        self.pHvalues=[]
        self.prot_states={}
        for residue in residues:
            self.prot_states[residue]={}
            while 1:
                pH=float(pKa_values[count])
                count=count+1
                charge=float(pKa_values[count])
                count=count+1
                #
                # If we find the terminator value
                #
                if pKa_values[count]==999.0:
                    if pKa_values[count+1]==-999.0:
                        count=count+2
                        break
                self.prot_states[residue][pH]=charge
            #
            #
            if self.pHvalues==[]:
                self.pHvalues=self.prot_states[residue].keys()
                self.pHvalues.sort()
            #
            # prot_states also contain the pKa value
            #
            self.prot_states[residue]['pKa']=pkavalues[residue]
          
        #
        # Construct the proper pKa array
        #
        self.pka={}
        for group in pkavalues.keys():
            self.pka[group]={'pKa':pkavalues[group]}
        if complete_pka:
            self.complete_pka()
        #
        # Construct the pHvalues array
        #
        
        #
        # All done
        #
        loadtid=(loadtime-starttime)/60.0
        tid=time.time()-starttime
        tid=tid/60.0
        return pkavalues,self.prot_states


    #
    # ---------------------
    #

    def calc_fractional_charge(self,pH):
        #
        # Not defined
        #
        raise 'Function not defined for C++ module'


class Boltzmann_CPP(Monte_Carlo_CPP):
    """ C++ implementation of the Boltzmann algorithm"""

    def _calc_pKas(self,mcsteps=200000,phstep=0.1,phstart=1.0,phend=20.0,verbose=1,complete_pka=None,exp_pHs=[]):
        """Do the pKa calculation with the CPP module for Boltzmann stats"""
        #
        # Do specific CPP setup
        #
        import time
        starttime=time.time()
        residues=self.intrinsic_pKa.keys()
        residues.sort()
        intpkas=[]
        acidbase=[]
        for residue in residues:
            intpkas.append(self.intrinsic_pKa[residue])
            acidbase.append(float(self.acid_base(residue)))
            #print residue,self.intrinsic_pKa[residue],self.acid_base(residue)
            # Works ok with mutations since self.acid_base works on the last residue
        linear_matrix=self.make_matrix_linear()
        #raw_input('ok?')
        #
        # Import the C++ module and run pKa calc from there
        #
        import time, math
        import pMC
        reload(pMC)
        use_MC=0 # This is Boltzmann
        FAST=pMC.MC(intpkas,linear_matrix,acidbase,use_MC)
        loadtime=time.time()
        FAST.set_MCsteps(int(mcsteps))
        pKa_values=FAST.calc_pKas(phstart,phend,phstep)
        #
        # Put the pKa values back in a dictionary
        #
        pkavalues={}
        count=0
        for residue in residues:
            pkavalues[residue]=pKa_values[count]
            if pKa_values[count]<-100.0:
                pkavalues[residue]=-9999.9
            count=count+1
        #
        # Now get the titration curves
        #
        # pMC returns a specification of how many values it returns
        #
        pH_start=pKa_values[count]
        pH_step=pKa_values[count+1]
        pH_num=pKa_values[count+2]
        count=count+3
        #
        #
        #
        self.pHvalues=[]
        self.prot_states={}
        for residue in residues:
            self.prot_states[residue]={}
            while 1:
                pH=float(pKa_values[count])
                count=count+1
                charge=float(pKa_values[count])
                count=count+1
                #
                # If we find the terminator value
                #
                if pKa_values[count]==999.0:
                    if pKa_values[count+1]==-999.0:
                        count=count+2
                        break
                self.prot_states[residue][pH]=charge
            #
            #
            if self.pHvalues==[]:
                self.pHvalues=self.prot_states[residue].keys()
                self.pHvalues.sort()
            #
            # prot_states also contain the pKa value
            #
            self.prot_states[residue]['pKa']=pkavalues[residue]
          
        #
        # Construct the proper pKa array
        #
        self.pka={}
        for group in pkavalues.keys():
            self.pka[group]={'pKa':pkavalues[group]}
        if complete_pka:
            self.complete_pka()
        #
        # Construct the pHvalues array
        #
        
        #
        # All done
        #
        loadtid=(loadtime-starttime)/60.0
        tid=time.time()-starttime
        tid=tid/60.0
        return pkavalues,self.prot_states


#
# ---------------------------------
#

class Monte_Carlo_Mult_CPP(Monte_Carlo):
    """Derived class for calling the mult version of the C++ code"""

    def test(self):
        """Test if we can import the C++ module"""
        import pMC_mult
        return

    #
    # --------
    # 
    
    def make_matrix_linear(self):
        """Change the matrix to linear form - this makes it easy to pass it to the C++ code
        In the mult form self.intene is called self.intene_mult, and each self.intene_mult[group1][group2] entry
        consists of a list that defines the interaction energies for each of the possible states of
        group1 and group2 with each other"""
        linear=[]
        residues=self.intene_mult.keys()
        residues.sort()
        for group1 in residues:
            for group2 in residues:
                for grp1_state in self.intene_mult[group1][group2]:
                    for intene in grp1_state:
                        linear.append(intene)
        return linear

    #
    # -----------------
    #

    def prepare_matrix(self):
        #
        # Prepare the matrix
        #
        # Precompute full term
        #
        self.intene={}
        for key in self.matrix.keys():
            self.intene[key]={}
            for key2 in self.matrix[key].keys():
                new_vals=[]
                for state1 in self.matrix[key][key2]:
                    row=[]
                    for state2 in state2:
                        row.append(state2[0]\
                                   -state2[1] \
                                   -state2[2] \
                                   +state2[3])
                new_vals.append(row)
                self.matrix[key][key2]=new_vals[:]
        self.intene=self.matrix.copy()
        #
        # Make sure that the matrix is symmetrical
        #
        #residues=self.intene.keys()
        #residues.sort()
        #for key in residues:
        #    for key2 in residues:
        #        E12=self.intene[key][key2]
        #        E21=self.intene[key2][key]
        #        ave_ene=(E12+E21)/2.0
        #        self.intene[key][key2]=ave_ene
        #        self.intene[key2][key]=ave_ene
        return

    #
    # ---------
    #

    def calc_pKas(self,mcsteps=200000,phstep=0.1,phstart=0.0,phend=14,verbose=1,complete_pka=None,exp_pHs=[]):
        """Calculate pKa values from the defined arrays - this is a class specific for the mult CPP code
        THIS CLASS STILL HAS TO BE MODIFIED FOR THE MULT CODE!!!!"""
        #
        # Calculate pKa values
        #
        # Init
        #
        # Set constants
        # KBol and Temp are equal to 1.0 everywhere in this derived class!
        #
        # Check that all dictionaries have been filled
        #
        self.groups=self.desolv.keys()
        self.groups.sort()
        b_groups=self.backgr.keys()
        b_groups.sort()
        m_groups=self.matrix.keys()
        m_groups.sort()
        if self.groups!=b_groups or self.groups!=m_groups:
            print
            print 'Inconsistent desolv, backgr and matrix dictionaries'
            print
            ndes=len(self.desolv.keys())
            nback=len(self.backgr.keys())
            nmat=len(self.matrix.keys())
            print 'Groups in desolv: %3d, groups in backgr: %3d, groups in matrix: %3d \n' %(ndes,nback,nmat)
            groups=self.backgr.keys()+self.desolv.keys()+self.matrix.keys()
            g_dict={}
            for group in groups:
                g_dict[group]=1
            groups=g_dict.keys()
            groups.sort()
            for group in groups:
                has=['-','-','-']
                if self.desolv.has_key(group):
                    has[0]='+'
                if self.backgr.has_key(group):
                    has[1]='+'
                if self.matrix.has_key(group):
                    has[2]='+'
                print '%14s Desolv: %s, Backgr: %s, Matrix: %s' %(group,has[0],has[1],has[2])
            print 'Totally done here'
            import sys
            sys.stdout.flush()
            raise 'Error in C++ Monte Carlo module'
        #
        # Prepare the matrix
        #
        self.prepare_matrix()

        #
        # Calculate the intrinsic pKa values
        #
        self.calc_intpKas()

        return self._calc_pKas(mcsteps,phstep,phstart,phend,verbose,complete_pka)

    #
    # ----
    #
    
    def _calc_pKas(self,mcsteps=200000,phstep=0.1,phstart=1.0,phend=20.0,verbose=1,complete_pka=None,exp_pHs=[]):
        """Do the pKa calculation with the CPP module"""
        #
        # Do specific CPP_Mult setup
        #
        import time
        starttime=time.time()
        #
        # Get the number of groups
        #
        residues=self.intrinsic_pKa.keys()
        residues.sort()
        #
        # num_states holds the number of states for each group
        # charged_state identified whether a state is charged or neutral
        #
        intpkas=[]
        acidbase=[]
        num_states=[]
        linear_charged_state=[]
        for residue in residues:
            num_states.append(len(self.intrinsic_pKa[residue]))
            #this_state=[]
            state_count=0
            for state in self.intrinsic_pKa[residue]:
                intpkas.append(float(state))
                acidbase.append(int(self.acid_base[residue]))
                linear_charged_state.append(self.charged_state[residue][state_count])
                state_count=state_count+1
        #
        # Linearise the matrix
        #
        linear_matrix=self.make_matrix_linear()
        #
        # Import the C++ module and run pKa calc from there
        #
        import time, math
        import pMC_mult
        FAST=pMC_mult.MC(intpkas,linear_matrix,acidbase,num_states,linear_charged_state)
        loadtime=time.time()
        #print 'CPP calculations starting at',time.strftime('%d/%m-%Y %H:%M',time.localtime(loadtime))
        FAST.set_MCsteps(int(mcsteps))
        #print 'phstart: %f, phend: %f, phstep: %f' %(phstart,phend,phstep)
        pKa_values=FAST.calc_pKas(phstart,phend,phstep)
        #
        # Put the pKa values back in a dictionary
        #
        pkavalues={}
        count=0
        for residue in residues:
            pkavalues[residue]=pKa_values[count]
            if pKa_values[count]<-100.0:
                pkavalues[residue]=-9999.9
            count=count+1
        #
        # Now get the titration curves
        #
        # pMC returns a specification of how many values it returns
        #
        pH_start=pKa_values[count]
        pH_step=pKa_values[count+1]
        pH_num=pKa_values[count+2]
        count=count+3
        #
        #
        #
        self.pHvalues=[]
        self.prot_states={}
        for residue in residues:
            self.prot_states[residue]={}
            while 1:
                pH=float(pKa_values[count])
                count=count+1
                charge=float(pKa_values[count])
                count=count+1
                #
                # If we find the terminator value
                #
                if pKa_values[count]==999.0:
                    if pKa_values[count+1]==-999.0:
                        count=count+2
                        break
                self.prot_states[residue][pH]=charge
            #
            #
            if self.pHvalues==[]:
                self.pHvalues=self.prot_states[residue].keys()
                self.pHvalues.sort()
            #
            # prot_states also contain the pKa value
            #
            self.prot_states[residue]['pKa']=pkavalues[residue]
          
        #
        # Construct the proper pKa array
        #
        self.pka={}
        for group in pkavalues.keys():
            self.pka[group]={'pKa':pkavalues[group]}
        if complete_pka:
            self.complete_pka()
        #
        # Construct the pHvalues array
        #
        
        #
        # All done
        #
        loadtid=(loadtime-starttime)/60.0
        tid=time.time()-starttime
        tid=tid/60.0
        #print 'CPP calculations took %.4f minutes' %tid
        #print 'setup time %.4f minutes' %loadtid
        #import sys
        #sys.stdout.flush()
        return pkavalues

    
#
# ---------------------------------
#

class Boltzmann(Monte_Carlo):
    #
    # Everything the same as the Monte Carlo routine, except the way the
    # fractional charges are calculated
    #

    def calc_fractional_charge(self,pH):
        #
        # Calculate the fractional charge for all residues
        # at this pH
        #
        # Get all the groups
        #
        groups=self.groups[:]
        groups.sort()
        #
        # Define the states
        #
        states={}
        import math
        number_of_states=int(math.pow(2,len(groups)))
        #
        # Create all states
        #
        for x in range(number_of_states):
            states[x]={'def':self.get_binary_number(x,len(groups)),'E':None}
        #
        # Calculate the energy for all states 
        #
        for state_number in states.keys():
            state={}
            count=0
            for group in groups:
                state[group]=states[state_number]['def'][count]
                count=count+1
            states[state_number]['E']=self.get_energy(pH,state)
        #
        # Find the fractional degree of protonation
        #
        # denominator
        nomi=0.0
        for state_number in states.keys():
            Energy=states[state_number]['E']
            nomi=nomi+math.exp(-Energy)
        #
        # For each group get the fractional degree of protonation
        #
        sumstate={}

        count=0
        for group in groups:
            numerator=0.0
            #
            # Sum energy for all states where this group is charged
            #
            for state_number in states.keys():
                if states[state_number]['def'][count]==1:
                    Energy=states[state_number]['E']
                    numerator=numerator+math.exp(-Energy)
            #
            # Calc fractional charge
            #
            sumstate[group]=numerator/nomi
            if isacid(group):
                sumstate[group]=-sumstate[group]
            #
            # Update count
            #
            count=count+1
        #
        # Calculate the population of each state
        #
        for state_number in states.keys():
            E_state=states[state_number]['E']
            numerator=math.exp(-E_state)
            states[state_number]['pop']=numerator/nomi
            #if states[state_number]['pop']<0.0000001:
            #    states[state_number]['pop']=0.0
            
        #
        # Done
        #
        return sumstate,states

    #
    # -----------------------
    #

    def get_binary_number(self,number,digits):
        #
        # Decompose the number
        #
        binary=[]
        import math
        for position in range(digits-1,-1,-1):
            if math.pow(2,position)>number:
                binary.append(0)
            else:
                binary.append(1)
                number=number-math.pow(2,position)
        #
        # Done
        #
        return binary

#
# ---------------------------------
#

class Tanford_Roxby(Monte_Carlo):

    def calc_fractional_charge(self,pH):
        #
        # Do a quick Tanford-Roxby run
        #
        TRmatrix=self.intene.copy()
        self.groups=TRmatrix.keys()
        pkas={}
        for group in TRmatrix.keys():
            if self.intrinsic_pKa.has_key(group):
                pkas[group]=self.intrinsic_pKa[group]
            else:
                 pkas[group]=self.get_modelpKa(group)
        #
        # Start looping
        #
        temppkas=pkas.copy()
        intpkas=pkas.copy()
        converged=None
        import math
        count=0
        charges_recorded={}
        while not converged and count<50:
            #
            # Calculate charges
            #
            count=count+1
            charge={}
            for group in pkas.keys():
                acibas=1.0
                if isacid(group):
                    acibas=-1.0
                charge[group]=acibas/(1.0+math.pow(10.0,acibas*(pH-temppkas[group])))
            #
            # Calculate energies
            #
            converged=1
            for group in pkas.keys():
                acibas=1.0
                if isacid(group):
                    acibas=-1.0
                #
                energy=0.0
                #
                # Add intpka term
                #
                energy=energy+acibas*math.log(10.0)*(pH-intpkas[group])
                #
                # site-site interaction energies
                #
                for group2 in pkas.keys():
                    import os
                    if not TRmatrix.has_key(group):
                        print 'TRmatrix does not have %s' %group
                        os._exit(0)
                    if not TRmatrix[group].has_key(group2):
                        print 'TRmatrix[%s] does not have %s' %(group,group2)
                        print TRmatrix[group]
                        os._exit(0)
                    if not charge.has_key(group2):
                        print 'Charge does not have %s' %group2
                        os._exit(0)
                    if not group==group2:
                        intene=TRmatrix[group][group2]
                        energy=energy+intene*abs(charge[group2])
  
                #
                # Will we converge this loop?
                ##
                charges_recorded=charge.copy()
                if abs(temppkas[group]-(pH-energy/(acibas*math.log(10.0))))>0.01:
                    converged=None
                #
                # Calculate new temppka
                #
                temppkas[group]=pH-energy/(acibas*math.log(10.0))
        #
        # Tadaa!
        #
        return charges_recorded,{}

    #
    # ------------------------
    #

    def get_energy(self,pH,state):
        #
        # Not defined
        #
        raise 'Function not defined in Tanford-Roxby'
    
