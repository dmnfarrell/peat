#!/usr/bin/env python
#
# pKa - various programs and scripts for pKa value analysis, calculation and redesign
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


import pKarun
import pKaTool

from output import *

#
# --------------------
#

def fact(number):
    #
    # Return fact(number)
    #
    if number>1:
        return number*fact(number-1)
    else:
        return number

#
# ----------------------
#

class pKa_calc:

    """This class is inherited by the main Design_pKa class"""

    def calculate_pKa_change(self,mutations):
        """ Calculate the dpKa value(s) for the residue(s) in self.desired_pKa_changes
            and calculate the score for this set of mutations"""
        #
        # Which method are we going to use?
        #
        dpKas=None
        if not self.params['tabulated']:
            #
            # Explicit (each combination of mutations is evaluated separately)
            #
            dpKas=self.calculate_dpKas_explicit(mutations,return_protstates=self.params['stability_mode'])
            if self.params['stability_mode']:
                self.prot_states=dpKas[1]
                dpKas=dpKas[0]
        else:
            #
            # Tabulated (dpKa values from combinations of mutations is calculated
            # by adding the dpKa values for each individual mutation)
            #
            dpKas=self.calculate_dpKas_tabulated(mutations)
        #
        # If dpKas is None then we used a mutation that cannot be modelled
        #
        if not dpKas:
            return None,None
        #
        # Remove pKa changes that are too small
        #
        for group in dpKas.keys():
            #
            # Is the shift significant?
            #
            if dpKas[group]:
                if abs(dpKas[group])<self.params['sens_limit'] and self.params['calc_dpka']!=1:
                    dpKas[group]=0.0
        #
        # Score the final result
        # Score is calculated as sum of (desired_change-achieved_change)*weight
        #
        score=0.0
        import math
        desired_residues=self.desired_pKa_changes.keys()
        desired_residues.sort()
        for dpKa_residue in desired_residues:
            #
            # Score each delta pKa
            #
            weight=self.desired_pKa_changes[dpKa_residue]['weight']
            #
            # Did we get a pKa shift?
            #
            if dpKas.has_key(dpKa_residue):
                if dpKas[dpKa_residue]!=None:
                    # Yes
                    score=score+abs(self.desired_pKa_changes[dpKa_residue]['dpKa']-dpKas[dpKa_residue])*weight
                else:
                    # No - this means an extreme shift
                    # We give an unfavourable score and warn the user
                    #
                    score=score+200
                    print 'A set of mutations gave an undefined pKa value for %s\n The set of mutations was %s' %(residue,str(mutations))
            else:
                print 'No dpKa for %s with these mutations: %s' %(dpKa_residue,str(mutations))                
        #
        # Add term penalising ineffective mutations
        #
        for mutation in mutations:
            if mutation:
                score=score+self.params['mutation_penalty']
        return score,dpKas

    #
    # ---------------------
    #

    def find_reporter_groups(self,mutations,check_groups=False):
        """Find all the reporter groups that are valid with this set of mutations

        If check_groups is True then check that the pKa values of the groups are
        inside pH_start and pH_stop

        """
        if not check_groups:
            reporter_groups=self.desired_pKa_changes.keys()
        else:
            #
            # If so requested get dpKa values for all wild type residues
            #
            reporter_groups=[]
            groups=self.wt_pKas.keys()
            groups.sort()
            for group in groups:
                if self.params['use_titration_curves']:
                    reporter_groups.append(group)
                else:
                    if self.wt_pKas[group]>self.params['pHstart'] and self.wt_pKas[group]<self.params['pHstop']:
                        reporter_groups.append(group)
        #
        # Make sure that we remove the groups that are included in the mutation
        #
        delete=[]
        for mutation in mutations:
            if mutation:
                import pKD_tools
                orgres=pKD_tools.get_oldres_from_mut(mutation)
                if orgres in reporter_groups:
                    delete.append(orgres)
        #
        rep_tmp=[]
        for group in reporter_groups:
            if group in delete:
                pass
            else:
                rep_tmp.append(group)
        reporter_groups=rep_tmp
        return reporter_groups

    #
    # ------
    #
            
    def calculate_dpKas_explicit(self,mutations,find_reporter_groups=None,return_protstates=False):
        #
        # Calculate the change in pKa value given the mutations in mutations
        #
        reporter_groups=self.find_reporter_groups(mutations,check_groups=find_reporter_groups)
        #
        # ========================================================================
        # Done finding reporter groups
        #
        # Calculate the pKa values for the mutant
        #
        if self.pKa_calc_alg=='MC' or self.pKa_calc_alg=='Tanford-Roxby':
            #
            # Calculate the pKas explicit for the mutations
            #
            # get the wild type pKa values
            #
            wt_pKas=[]
            wt_pKa_names={}
            for group in reporter_groups:
                if self.wt_pKas.has_key(group):
                    wt_pKas.append(self.wt_pKas[group])
                    wt_pKa_names[group]=group
                else:
                    wt_pKas.append(None)
            # ==============================================================
            #
            # Check if we already know this solution
            #
            real_muts=[]
            for mut in mutations:
                if mut:
                    real_muts.append(mut)
            real_muts.sort()
            found_all=1
            sol_key=str(real_muts)
            #
            # See if we have a solution
            #
            import os
            solfile=os.path.join(self.sol_dir,sol_key)
            if os.path.isfile(solfile):
                fd=open(solfile)
                import pickle
                solution=pickle.load(solfile)
                #print 'Found solution in sol_dict'
                for group in reporter_groups:
                    if not solution.has_key(group):
                        print 'Did not find all reporter groups. First miss is',group
                        found_all=None
                        residues=solution.keys()
                        residues.sort()
                        print residues
                        print
                        break
            else:
                found_all=None
            #
            # We have to calculate for this mutation
            #
            if not found_all:
                #
                # See if we have enough data to do this
                #
                for test_mut in mutations:
                    if test_mut:
                        if not self.data.has_key(test_mut):
                            #
                            # This means we have to get interaction energies and a rotamer score for this mutation
                            #
                            import Design_pKa_help
                            rotamer_score=self.get_rotamer_score(test_mut)
                            self.data[test_mut]={'Rotamer quality':rotamer_score}
                        elif self.data[test_mut]==None:
                            import Design_pKa_help
                            rotamer_score=self.get_rotamer_score(test_mut)
                            self.data[test_mut]={'Rotamer quality':rotamer_score}
                            print 'Calculated rotamer quality for %s: %5.3f' %(test_mut,rotamer_score)
                #
                # Call the calculation routine
                #
                self.O._print('Calling C++ routine with these mutations: %s' %(str(mutations)),1)
                mut_pKas,prot_states=self.pKaCALC.calc_MC_mut_pKas(mutations=mutations,
                                                                   data=self.data,
                                                                   atomdata=self.atomdata,
                                                                   wt_pKas=wt_pKas,
                                                                   require_pKa=self.desired_pKa_changes.keys(),
                                                                   for_tabulated=find_reporter_groups)
                if mut_pKas=='SS-bridge':
                    return {'Failed':'Possible SS-bridge'}
                #
                # Store this solution
                #
                if mut_pKas:
                    sol_file=os.path.join(self.sol_dir,sol_key)
                    if os.path.isfile(sol_file):
                        import pickle
                        fd=open(sol_file)
                        solution=pickle.load(sol_file)
                        fd.close()
                    else:
                        solution={}
                        
                    if solution.has_key(sol_key):
                        for res in mut_pKas.keys():
                            solution[sol_key][res]=mut_pKas[res]
                        solution[sol_key]['titcurves']=prot_states
                    else:
                        solution[sol_key]=mut_pKas.copy()
                        solution[sol_key]['titcurves']=prot_states
            else:
                #
                # Use the old results
                #
                mut_pKas=solution[sol_key].copy()
                del mut_pKas['titcurves']
                prot_states=solution[sol_key]['titcurves'].copy()
            #
            # ..
            #
        elif self.pKa_calc_alg=='phidiff':
            #
            # Here for the simple routine
            #
            mut_pKas=self.pKas_philn10(mutations)
            prot_states={}
        else:
            print 'Unknown pKa calculation algorithm'
            print self.pKa_calc_alg
            raise pKD_error('bug')
        #
        # If we get no pKa values and we're not using titration curves then return
        #
        if not mut_pKas and not self.params['use_titration_curves']:
            #
            # Sometimes we get a None back
            #
            return {}
        elif not prot_states and self.params['use_titration_curves']:
            return {}
        #
        # Calculate the delta pKa values
        #
        import string
        dpKas={}
        for group in reporter_groups:
            wt_group=wt_pKa_names[group]
            if not self.wt_pKas.has_key(wt_group):
                print 'The wild type does not have this reporter group: %s' %group
                continue
            if not self.wt_pKas[wt_group] and not self.params['use_titration_curves']:
                #
                # If we don't have a wt pKa value then there's no
                # point in continuing..
                #
                print group,self.wt_pKas[group]
                raise pKD_error('Wild type pKa value not defined')
            #
            #
            if not self.params['use_titration_curves']:
                if not mut_pKas.has_key(group):
                    #
                    # No dpKa value for this group
                    #
                    print 'I could not get a dpKa value for %s' %group
                    mut_pKas[group]=None
            #
            # Is it a real pKa value?
            #
            if self.params['use_titration_curves']:
                #
                # Subtract the two titration curves from each other in the integration range
                #
                if not self.wt_prot_states or not prot_states:
                    print 'No pKa calculation result at all...'
                    raise Exception()
                diff=0.0
                pH_values=self.wt_prot_states[wt_group].keys()
                pH_values.sort()
                for pH in pH_values:

                    if pH>=self.params['pHstart_integrate'] and pH<=self.params['pHstop_integrate']:
                        diff=diff+(prot_states[group][pH]-self.wt_prot_states[wt_group][pH])*self.params['pHstep']
                        #print 'Integrating',pH,group,diff,self.params['pHstep'],prot_states[group][pH],self.wt_prot_states[wt_group][pH]
                dpKas[group]=diff
                #print '--------------'
            else:
                #
                # We're using pKa 1/2 values
                #
                if not mut_pKas[group]:
                    self.muts_excluded.append(mutations)
                    mut_pKas[group]=None
                    print 'Mutant pKa value not defined for target %s. Setting dpKa to None' %group
                    dpKas[group]=None
                else:
                    dpKas[group]=mut_pKas[group]-self.wt_pKas[wt_group]
        #
        # Return the dpKa values for all the reporter groups
        #
        print 'I am here and return_protstates is',return_protstates
        if return_protstates:
            return dpKas,prot_states
        return dpKas

    #
    # -----------------------
    #

    def calculate_dpKas_tabulated(self,mutations):
        #
        # Calculate the score (abs. cumulated sum of pKa differences) for a particular
        # set of mutations.
        #
        # This is the subroutine that's used for phi/ln(10) and MC
        #
        dpKas={}
        reporter_groups=self.find_reporter_groups(mutations)
        #
        # Loop over all reporter_groups
        #
        for target_residue in reporter_groups:
            dpKas[target_residue]=0.0
        #
        for mutation in mutations:
            if not mutation:
                continue
            #
            # Check if this mutation is over the sensitivity limit for any target residue
            #
            temp_target={}
            ok=None
            for target_residue in reporter_groups:
                #
                # Get the dpKa value for the mutation
                #
                dpKa=self.tab_data.get_dpKa(mutation,target_residue,self)
                try:
                    if abs(dpKa)<self.params['sens_limit']:
                        dpKa=0.0
                except:
                    print 'dpKa for %s was %s. This mutation can not be used!!' %(str(mutation),str(dpKa))
                    dpKa=None
                    #if len(reporter_groups)>1:
                    #    raise Exception()
                    #else:
                    return None
                if abs(dpKa)>0.0:
                    ok=1
                temp_target[target_residue]=dpKa
            #
            # Did we have any non-zero dpKa values? (this is to make sure that we
            # only add dpKas if they are higher than that sensitivity limit)
            #
            if ok:
                for target in temp_target.keys():
                    #
                    # Add it to the total dpKa if it is not None (a dpKa[target] can
                    # be set to None if one of the mutations has mutated it)
                    #
                    if not dpKas[target] is None:
                        dpKas[target]=dpKas[target]+temp_target[target]
        #
        # Return the total dpKas
        #
        return dpKas

    #
    # -------------------
    #

    def pKas_philn10(self,mutations):
        #
        # Get the dpKas for all titratable groups and add it to their current pKa
        #
        pKas={}
        reporter_groups=self.desired_pKa_changes.keys()
        for group in reporter_groups:
            dpKa_tot=0.0
            for mutation in mutations:
                dpKa=self.dpKa_philn10(mutation,group)
                if dpKa:
                    dpKa_tot=dpKa_tot+dpKa
            pKas[group]=dpKa_tot+self.wt_pKas[group]
        return pKas

    #
    # ---------------------
    #       

    def dpKa_philn10(self,mutation,target_residue,verbose=None):
        """
        # For a given mutation, calculate the change in the pKa value if the target using phi/ln(10)
        """
        # If we have no mutation, then we have no change
        #
        if not mutation:
            return 0.0
        #
        # Check that target_residue and mutation is not the same residue
        #
        import pKD_tools
        if pKD_tools.get_resnum_from_mut(mutation)==pKD_tools.get_resnum_from_res(target_residue):
            return None
        #
        # Get the interaction energy 
        #
        int_ene=None
        if self.data.has_key(mutation):
            if self.data[mutation].has_key('neutral mutation'):
                #
                # We have a neutral mutation e.g. ARG:0112:GLN
                # In this case we get the interaction energy from the matrix file
                #
                oldres=pKD_tools.get_oldres_from_mut(mutation)
                newres=pKD_tools.get_newres_from_mut(mutation)
                #
                # Check that we have a neutral mutation
                #
                if pKarun.is_titratable(oldres) and not pKarun.is_titratable(newres):
                    int_ene=self.matrix[oldres][target_residue]
                    int_ene=int_ene[0]-int_ene[1]-int_ene[2]+int_ene[3]
                else:
                    print 'Could not find interaction energy'
                    print 'oldres',oldres
                    print 'target',target_residue
                    int_ene=0.0
                    print 'This should not happen'
                    raise pKD_error('Something is wrong')
            elif self.data[mutation].has_key(target_residue):
                int_ene=self.data[mutation][target_residue]
            else:
                print 'Could not find residue: %s ' %mutation
                print 'When trying to get the interaction energy for mutation %s with target residue %s' %(str(mutation),str(target_residue))
                keys=self.data[mutation].keys()
                keys.sort()
                print keys
                raise pKD_error('Could not find interaction energy')
        else:
            print 'Could not find mutation: %s' %mutation
            keys=self.data.keys()
            keys.sort()
            print keys
            raise Exception()
        if not int_ene:
            raise pKD_error('Could not get interaction energy')
        #
        # Calculate the dpKa 
        #
        dpKa=0.0
        import math
        #
        # The sign of the interaction energy tells us if we have
        # an upward or a downward shift
        #
        if pKarun.isacid(target_residue):
            # Acid
            dpKa=int_ene/math.log(10)
        else:
            # Base
            dpKa=-int_ene/math.log(10)
        #
        # Was the wt residue (that of the mutation) a charged group?
        #
        import string
        number=string.split(mutation,':')[1]
        new_residue=':'+number+':'+string.split(mutation,':')[-1]
        old_residue=':'+number+':'+string.split(mutation,':')[0]
        if pKarun.is_charged(old_residue):
                #
                # Calculate effect of the mutation
                #
                #if pKa.istitratable(new_residue):
                if pKarun.is_charged(new_residue):
                    if pKarun.isacid(old_residue)==pKarun.isacid(mutation):
                        dpKa=0.0
                    else:
                        #
                        # If it was the opposite sign, then we get twice the effect
                        #
                        dpKa=2.0*dpKa
                else:
                    #
                    # Not titratable, - a neutral mutation
                    # Opposite effect of old residue
                    #
                    if pKarun.isacid(target_residue):
                        # Acid
                        dpKa=-int_ene/math.log(10)
                    else:
                        # Base
                        dpKa=int_ene/math.log(10)
        #
        # ok, all done. Return dpKa
        #
        return dpKa

    # -------------------
    #
    #
    # Routines for handling the calculation of dpKa values from an
    # already specified set of mutations
    #

    def calc_dpKa(self):
        """
        # Calculate the dpKa for a particular set of mutations
        """
        # Initialise
        #
        self.initialise('calcdpka')
        #
        # Calculate starting state (wt pKa vals)
        #
        none_muts=[None]
        #
        # Set the mutations
        #
        mutations=self.params['mutations']
        if not mutations:
            self.O._print('You have to specify mutations with the -mutations flag when calculating dpKas',0)
            print 'I got these',mutations
            return

        #
        # Do we have multiple sets of multiple mutations?
        #
        results={}
        import types 
        self.O._print('Residue   Old pKa   New pKa   Difference   Energy     Mutation(s)',1)
        self.O._print('-------------------------------------------------------------------',1)
        if type(mutations[0]) is types.ListType:
            #
            # Multiple sets of mutations
            #
            results['target']=self.desired_pKa_changes.keys()
            for muts in mutations:
                do_it=1
                for mut in muts:
                    for key in self.desired_pKa_changes.keys():
                        if key[:-3]==mut[:-3]:
                            do_it=None
                if do_it:
                    import string
                    score,dpKas=self.calc_dpKa_sub(muts)
                    muts.sort()
                    results[str(muts)]=dpKas
            #
            # Set the outputfile extension
            #
            rtype='phidiff'
            if self.params['MC']:
                rtype='MC'
            #
            # How did we combine mutations?
            #
            comb_type='Explicit'
            if self.params['tabulated']:
                comb_type='Tabulated'
            #
            # Get the target residue
            #
            target_res=self.desired_pKa_changes.keys()
            target_res.sort()
            import string
            target_res=string.join(target_res,'+')
            #
            # Write the file
            #
            fd=open('performance_%s_%s_%s.pickle' %(rtype,comb_type,target_res),'w')
            import cPickle
            cPickle.dump(results,fd)
            fd.close()
        else:
            #
            # Just a single set of mutations
            #
            # Make sure that we don't mutate the residue that we are measuring!
            #
            do_it=1
            for mut in mutations:
                for key in self.desired_pKa_changes.keys():
                    if key[:-3]==mut[:-3]:
                        do_it=None      
            if do_it:
                score,results[str(mutations)]=self.calc_dpKa_sub(mutations)
        return results

    #
    # ----------------------
    #

    def calc_dpKa_sub(self,mutations):
        #
        # Calculate the new pKas
        #
        import string, sys
        score,dpKas=self.calculate_pKa_change(mutations)
        if self.params['stability_mode']:
            import pylab
            for titcurv,color in [[self.wt_prot_states,'b'],[self.prot_states,'r']]:
                import pKarun.pKa_general
                X=pKarun.pKa_general.pKanalyse()
                stab,sums=X.analyse_stability(titcurv,write_file=False)
                curve={}
                for residue in stab.keys():
                    for pH in stab[residue].keys():
                        if not curve.has_key(pH):
                            curve[pH]=0.0
                        curve[pH]=curve[pH]+stab[residue][pH]
                print
                
                xs=curve.keys()
                xs.sort()
                ys=[]
                for x in xs:
                    ys.append(curve[x])
                pylab.plot(xs,ys,'%s-' %color)
            pylab.show()
            stop
            print self.wt_prot_states.keys()
            print 'In stability mode'
            print self.prot_states.keys()
            print
            stop
        #
        # Print summary
        #
        sys.stdout.flush()
        residues=self.desired_pKa_changes.keys()
        residues.sort()
        for dpKa_name in residues:
            wt_pka=self.wt_pKas[dpKa_name]
            if wt_pka<-8000:
                wt_pka=0.0
            #
            # Get the name of this residue for this set of mutations
            #
            if not dpKas:
                print '%s: dpKa calculation cannot be performed. Probably incorrect modelling of the mutation' %(dpKa_name)
                print dpKas
                continue
            if not dpKas.has_key(dpKa_name):
                dpKa_name=None
            #
            # Print the info
            #
            if dpKas is None or dpKa_name is None:
                self.O._print('%12s  %6.3f       %s    %s      %6.1f   %s' %(dpKa_name,wt_pka,
                                                                            'ERROR ',
                                                                            'ERROR ',
                                                                            -999.9,string.join(mutations,'+')),1) 
            else:
                self.O._print('%12s  %6.3f       %6.3f    %6.3f      %6.1f   %s' %(dpKa_name,wt_pka,
                                                                                  wt_pka+dpKas[dpKa_name],
                                                                                  dpKas[dpKa_name],
                                                                                  score,string.join(mutations,'+')),1)
        return score,dpKas
    

#
# ---------------------
#

class Design_pKa(pKa_calc):

    def __init__(self,params):
        """
        Constructor for the Design pKa class.
        This is the place where we define the mutations we can use

        """
        self.params=params
        self.parse_parameters()
        #
        # Define the pdb files
        #
        import os
        self.pdbfile=self.params['pdb']
        #
        # Init mutfile_names
        #
        self.mutfile_names={}
        #
        #
        # Make sure we have an absolute path
        #
        import os
        self.pdbfile=os.path.join(os.getcwd(),self.pdbfile)
        if not self.pdbfile:
            self.pdbfile_missing()
        if not os.path.isfile(self.pdbfile):
            self.pdbfile_missing(self.pdbfile)
        #
        # If we have a ligand then add that to the pdb file
        #
        if self.params['ligands']!=[]:
            import Protool
            X=Protool.structureIO()
            X.readpdb(self.pdbfile)
            orgpdbfile=self.pdbfile
            #
            # Now add the ligand(s) and save the combined file as a different pdbfile
            #
            L=Protool.ligand(X)
            ligand_name=''
            for ligandfile in self.params['ligands']:
                print 'Adding ligandfile: %s' %ligandfile
                L.readmol2(ligandfile,tag='STRIP')
                ligand_name=ligand_name+ligandfile
            self.pdbfile='%s_%s.pdb' %(self.pdbfile[:-4],ligand_name)
            X.writepdb(self.pdbfile)
            #
            # Copy the pKa calc files - they don't change
            #
            files=['.PKA.DAT','.MATRIX.DAT','.DESOLV.DAT','.BACKGR.DAT','.TITCURV.DAT']
            for filename in files:
                real_file=orgpdbfile+filename
                new_file=self.pdbfile+filename
                if os.path.isfile(real_file):
                    import shutil
                    shutil.copy(real_file,new_file)
            
        # -----------------------------------------------
        #
        # Load the PDB file into Protool
        #
        import Protool
        self.PI=Protool.structureIO()
        self.PI.readpdb(self.pdbfile)
        #
        # Set the topdir arg - this if for storing datafiles
        #
        self.topdir=os.getcwd()
        #
        # Create a special parameter dictionary for pKa.pKarun
        #
        self.pKarun_params={}
        include=['indi','dbcrit','ion','exdi']#,'allow_unknown_atoms','unknown_crg','unknown_rad']
        import copy
        for key in include:
            self.pKarun_params[key]=copy.copy(self.params[key])
        #
        # Array definitions
        #
        self.wt_pKas=None
        self.muts_excluded=[]
        #
        # Instantiate pKa_info
        #
        import Design_pKa_help
        self.pKa_info=Design_pKa_help.pKa_info(self.pdbfile,parent=self,PBEsolver=self.params['PBEsolver'])
        #
        # Stupid
        #
        methods={'TR':'Tanford-Roxby',
                 'MC':'MC',
                 'BM':'Boltzmann',
                 'phidiff':'phidiff'}
        self.pKa_calc_alg=methods[self.params['dpKa_method']]
        print 'pKa calculation algorithm is "%s"' %self.pKa_calc_alg
        #if self.params['MC']:
        #    self.pKa_calc_alg='MC'
        #elif self.params['TR']:
        #    self.pKa_calc_alg='Tanford-Roxby'
        #elif self.params['BM']:
        #    self.pKa_calc_alg='Boltzmann'
        #else:
        #    self.pKa_calc_alg='phidiff'
        #
        # Instantiate the tabulated data
        #
        recalc_intpka=None
        if self.params['recalc_intpka']:
            recalc_intpka=self.params['recalc_intpka_dist']
        self.tab_data=Design_pKa_help.pKa_tabulated(self.pdbfile,
                                                    self.pKa_calc_alg,
                                                    self.params['use_titration_curves'],
                                                    recalc_intpka)
        #
        # Set the pKa solver
        #
        import pKa_MC
        self.pKaCALC=pKa_MC.pKa_calculation_class(self.pdbfile,self.pKa_info,self.params,self)
        self.pKa_calc_alg=None
        if self.params['dpKa_method']=='MC':
            self.pKaCALC.set_MC_CPP()
            self.pKa_calc_alg='MC'
        elif self.params['dpKa_method']=='TR':
            self.pKaCALC.set_Tanford_Roxby()
            self.pKa_calc_alg='Tanford-Roxby'
        elif self.params['dpKa_method']=='BM':
            self.pKaCALC.set_Boltzmann()
            self.pKa_calc_alg='Boltzmann'
        elif self.params['dpKa_method']=='phidiff':
            self.pKa_calc_alg='phidiff'
        else:
            raise Exception('Incorrect pKa calculation method: %s' %self.params['dpKa_method'])
        #
        # Get the intrinsic pKa values
        #
        self.intrinsic_pKa=self.pKaCALC.MC.intrinsic_pKa.copy()
        #
        # Open the PKA.DAT file to get basic information
        #
        import pKaTool.pKaIO
        IO=pKaTool.pKaIO.pKaIO(self.pdbfile)
        self.pkavals=IO.readpka()
        #
        # I also need the matrix file to estimate the effect of neutral
        # mutations using phi/ln(10)
        #
        self.matrix=IO.read_matrix()
        #
        # ================================================
        #
        # Get the interaction energies
        #
        self.get_interaction_energies()
        #
        # Score the accessibilty the mutations
        #
        self.ok_mutations=self.score_accessibility_of_parent_residues()
        #
        # Set up a directory for storing the solutions
        #
        self.sol_dir=os.path.join(os.getcwd(),self.pdbfile+'.'+self.pKa_calc_alg+'.sol_dict')
        print 'Directory for solutions',self.sol_dir,self.pKa_calc_alg
        if self.params['use_titration_curves']:
            self.sol_dir=self.sol_dir+'_titration_curves'
        if self.params['recalc_intpka']:
            self.sol_dir=self.sol_dir+'_reintpka_%d' %self.params['recalc_intpka_dist']
        #
        # Create the dir if it doesn't exist
        #
        if not os.path.isdir(self.sol_dir):
            if os.path.isfile(self.sol_dir):
                os.unlink(self.sol_dir)
            os.mkdir(self.sol_dir)
        #
        # Done with init
        #
        return

    def parse_parameters(self):
        #
        # Output level
        #
        self.O=verbose(self.params['verbose'])
        #
        # Reformat some of the parameters
        #
        import string
        #if self.params['mutate_res']:
        #    self.getreslist('mutate_res')
        if self.params['mutations']:
            self.getreslist('mutations')
        return
    #
    # -----------------
    #
    def get_wild_type_pKas(self):
        #
        # Get the pKas for the wild type
        #
        import os, string
        tit_curve_text='dpKas'
        if self.params['use_titration_curves']:
            tit_curve_text='titration_curves'
        filename=os.path.join(self.topdir,'WT_pKas_%s_%s_%s' %(self.pKa_calc_alg,
                                                               tit_curve_text,
                                                               string.join(self.desired_pKa_changes.keys(),'_')))
        #
        # Open the list of filenames
        #
        catalog=os.path.join(self.topdir,'filenames.pickle')
        if os.path.isfile(catalog):
            fd=open(catalog)
            import pickle
            try:
                dict=pickle.load(fd)
            except EOFError:
                dict={}
                fd.close()
            fd.close()
        else:
            dict={}
        #
        # Look for the real filename
        #
        if dict.has_key(filename):
            real_filename=dict[filename]
        else:
            import tempfile
            real_filename=tempfile.mkstemp(suffix='pKD_wtpKas',dir=self.topdir)[1]
            dict[filename]=real_filename
            import pickle
            fd=open(catalog,'w')
            pickle.dump(dict,fd)
            fd.close()
        #
        # Now open the file
        #
        readpkas=False
        if os.path.isfile(real_filename):
            #
            # We have a file with the results
            #
            self.O._print('',1)
            self.O._print('Loading wild type pKa values........',1)
            import sys
            import cPickle
            fd=open(real_filename)
            try:
                wtpkas,wt_prot_states=cPickle.load(fd)
                fd.close()
                readpkas=True
            except EOFError:
                fd.close()
                os.unlink(real_filename)
        #
        # If the file doesn't exist or we couldn't read it
        #
        if not readpkas:
            #
            # Dang - we have to calculate them
            #
            self.O._print('',1)
            self.O._print('Calculating wild type pKa values........',1)
            wtpkas,wt_prot_states=self.pKaCALC.calc_wt_pkas(self.desired_pKa_changes.keys())
            fd=open(real_filename,'w')
            import cPickle
            cPickle.dump([wtpkas,wt_prot_states],fd)
            fd.close()
        #
        # Return the values
        #
        self.O._print('done',1)
        return wtpkas,wt_prot_states

    #
    # ------------------
    #

    def initialise(self,name):
        """
         Initialise the pKa Design class
        """
        self.name=name
        #
        # Parse the design argument
        #
        self.parse_pKa_arguments(self.params['pKas'])
        #
        # Get the wild type pKa values
        #
        if self.params['dpKa_method']=='MC' or self.params['dpKa_method']=='TR':
            #
            # For the MC + TR methods we have to calculate wt pKa values
            #
            self.pKaCALC.set_reporter_groups(self.desired_pKa_changes.keys())
            #
            # Calculate the wild type pKa values
            #
            self.wt_pKas,self.wt_prot_states=self.get_wild_type_pKas()
            #
            # Print the pKa values that we calculated for the wild type
            #
            wt_residues=self.wt_pKas.keys()
            wt_residues.sort()
            for residue in wt_residues:
                if self.wt_pKas[residue]:
                    self.O._print('%15s %5.3f' %(residue,self.wt_pKas[residue]),1)
                else:
                    self.O._print('%15s Undefined' %(residue),1)
        else:
            self.wt_pKas={}
            for residue in self.pkavals.keys():
                self.wt_pKas[residue]=self.pkavals[residue]['pKa']
        #
        # -----------------------------------------------------------------------
        #
        # Which shifts do we want? (parse arguments)
        #
        import string
        self.parse_pKa_arguments(self.params['pKas'])
        #
        # Filter the mutations
        #
        ok_mutations=self.ok_mutations.copy()
        self.filter_mutations(ok_mutations)
        #
        # Done
        #
        return
    #
    # ------------
    #

    def getreslist(self,parameter):
        #
        # Get a list of residues from the input parameter
        #
        import string, os
        tmp_split=string.split(self.params[parameter],',')
        split=[]
        for entry in tmp_split:
            split=split+entry.split('+')
        if split[0][0]==':' or not os.path.isfile(self.params[parameter]):
            self.params[parameter]=split
        else:
            #
            # If we do not have a colon at the start, then it is probably a filename (or an error)
            #
            # The format of the file is either:
            # <1-lt restyp><resnum> e.g R10
            # or
            # <resnum>
            # or, required for -mutations:
            # <3lt restyp_old>:<resnum>:<3lt restyp_new>[,/+]<3ltrestyp_old>:<resnum2>:<3lt restyp_new> ... etc.
            # 
            fd=open(self.params[parameter])
            self.params[parameter]=[]
            line=string.strip(fd.readline())
            while line:
                if len(line)>2:
                    number=None
                    if line[0]=='#':
                        pass
                    elif line[0] in string.letters and not line[1] in string.letters:
                        number=line[1:] 
                    elif line[0]==':' or line[1]==':':
                        #
                        # Is the line of the form O:NNNN:ORG:NEW or :NNNN:ORG:NEW ?
                        #
                        tmpsplit=string.split(line,',')
                        split=[]
                        for entry in tmpsplit:
                            split=split+entry.split('+')
                        self.params[parameter].append(split)
                    else:
                        number=line
                    if number:
                        self.params[parameter].append(string.zfill(number,4))
                line=string.strip(fd.readline())
            fd.close()
        return

    #
    # -------------------
    #

    def pdbfile_missing(self,pdbfile):
        raise Exception('PDB file not found: %s' %pdbfile)

    #
    # ---------------------
    #

    def pkas_missing(self):
        import os
        print
        print 'You did not give any pKa values to re-design.'
        print 'The PDB file contains the following titratable groups:'
        groups= self.pkavals.keys()
        groups.sort()
        print groups
        print 'Exiting...'
        print
        raise Exception()
    
    #
    # ----------------------
    #

    def get_interaction_energies(self):
        """
        Calculate the interaction energies for all the mutations that change the net charge on the protein
        This function runs sugelm in WHAT IF
        """
        import os, cPickle, sys
        stdout=sys.stdout
        import tempfile
        fd_tmp=open(tempfile.mktemp(),'w')
        sys.stdout=fd_tmp
        fd_tmp.close()
        sys.stdout=stdout
        resultfile=os.path.join(self.topdir,self.pdbfile+'.sugelm_data')
        if not os.path.isfile(resultfile) and self.params['generate_mutations']:
            self.O._print('--------------------------------------',1)
            self.O._print('Running SUGELM',1)
            self.O._print('--------------------------------------',1)
            if self.params['PBEsolver']=='DelPhi':
                #
                # Create the directory
                #
                sugdir=os.path.join(self.topdir,'sugelm')
                if not os.path.isdir(sugdir):
                    import os
                    os.mkdir(sugdir)
                #
                # Generate the sugelm data in the new way
                #
                import pKD_dict
                self.data=pKD_dict.pKD_dict()
                self.atomdata=pKD_dict.pKD_dict()
                import Protool
                P=Protool.structureIO()
                P.parse_terms=None # Make sure we don't invent new chain IDs
                P.readpdb(self.pdbfile)
                #
                residues=P.residues.keys()
                residues.sort()
                charge_mutations=['ASP','GLU','ARG','LYS','HIS','TYR','CYS']
                total=len(residues)*len(charge_mutations)
                count=0
                Failed=[]
                for residue in residues:
                    if not P.isaa(residue):
                        count=count+1
                        continue
                    resname=P.resname(residue)
                    for new_resname in charge_mutations:
                        print '-----------------------------------------'
                        print '%5d of %5d mutations. %5.2f%% done' %(count,total,float(count)/total*100.0)
                        print '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
                        print
                        import sys
                        sys.stdout.flush()
                        count=count+1
                        if resname!=new_resname:
                            import Design_pKa_help
                            mutation='%s:%s:%s' %(residue,resname,new_resname)
                            if self.mutfile_names.has_key(mutation):
                                newfilename=self.mutfile_names[mutation]
                                bump_score=0.0
                            else:
                                newfilename,bump_score=Design_pKa_help.make_mutation(pdbfile=self.pdbfile,mutation=mutation,topdir=self.topdir)
                            #
                            # Define the filename for this PBE files
                            #
                            matrix_file=os.path.join(sugdir,mutation+'.matrix')
                            atomdata_file=os.path.join(sugdir,mutation+'.atomdata')
                            #
                            # Is the bump score low enough?
                            #
                            if bump_score<=self.params['mutation_quality'] and not bump_score is False and not bump_score is None:
                                #
                                # Do we have the matrix saved?
                                #
                                if os.path.isfile(matrix_file) and os.path.isfile(atomdata_file):
                                    import pickle
                                    fd=open(matrix_file)
                                    matrix=pickle.load(fd)
                                    fd.close()
                                    #
                                    fd=open(atomdata_file)
                                    atomdata=pickle.load(fd)
                                    fd.close()
                                    if atomdata=='PBE Failed' or atomdata=='Modelling Failed':
                                        Failed.append([mutation,atomdata])
                                else:
                                    #
                                    # Re-calculate
                                    #
                                    print 'Calculating PBE charge-charge matrix for %s' %mutation
                                    self.close_stdout()
                                    try:
                                        matrix,atomdata=Design_pKa_help.run_WI_pKa_calc(pdbfile=newfilename,
                                                                                        mutated_res='%s:%s' %(residue,new_resname),
                                                                                        pKarun_params=self.pKarun_params,
                                                                                        do_matrix=True)
                                    except:
                                        self.open_stdout()
                                        print '!!!!\n\nCould not calculate PBE charge-charge matrix for mutation: %s. Skipping this mutation\n\n!!!!!' %mutation
                                        fd2=open(atomdata_file,'w')
                                        fd=open(matrix_file,'w')
                                        import pickle
                                        pickle.dump('PBE Failed',fd)
                                        pickle.dump('PBE Failed',fd2)
                                        fd.close()
                                        fd2.close()
                                        Failed.append([mutation,'PBE Failed'])
                                        continue
        
                                    #
                                    self.open_stdout()
                                    matrix['Rotamer quality']=bump_score
                                    fd=open(matrix_file,'w')
                                    import pickle
                                    pickle.dump(matrix,fd)
                                    fd.close()
                                    #
                                    fd=open(atomdata_file,'w')
                                    import pickle
                                    pickle.dump(atomdata,fd)
                                    fd.close()
                                #
                                # Put the result in the self.data array
                                #
                                if type(matrix)==type({}):
                                    self.data[mutation]=matrix.copy()
                                    self.atomdata[mutation]=atomdata.copy()
                            else:
                                print 'Mutation cannot be modelled - not using %s' %mutation
                                Failed.append([mutation,'Modelling Failed'])
                                fd=open(matrix_file,'w')
                                import pickle
                                pickle.dump('Modelling Failed',fd)
                                fd.close()
                                #
                                fd=open(atomdata_file,'w')
                                import pickle
                                pickle.dump('Modelling Failed',fd)
                                fd.close()
                        else:
                            print 'No mutation - skipping mutation to wt residue'
                print 'Failed mutations'
                print Failed
            elif self.params['PBEsolver']=='APBS':
                import pdb2pka_interface
                P2P=pdb2pka_interface.P2pI(os.getcwd(),self.pdbfile,self.pKarun_params)
                self.data,self.atomdata=P2P.pdb2pka_sugelm()
            else:
                raise Exception,'Unknown PBEsolver:'+self.params['PBEsolver']
            #
            # Read the original PDB file and insert the org residue
            # so we get the mutation in the new format: ChainID:residue number:Orgres:Newres
            # e.g. A:0010:ASP:ASN
            #
            import Protool
            P=Protool.structureIO()
            P.parse_terms=None # Make sure we don't invent new chain IDs
            P.readpdb(self.pdbfile)
            #
            # Save the results
            #
            import copy, pickle
            data={}
            for key in self.data.keys():
                data[key]=self.data[key]
            atomdata={}
            for key in self.atomdata.keys():
                atomdata[key]=self.atomdata[key]
            savedict={'data':data,'atomdata':atomdata}
        else:
            import pKD_dict
            self.data=pKD_dict.pKD_dict()
            self.atomdata=pKD_dict.pKD_dict()
        return



    #
    # -----------------------
    #

    def map_on_Bfactors(self,pdbfile,array):
        #
        # Map the scaled differences onto the PDB file
        #
        import Protool, pKD_tools
        X=Protool.structureIO()
        X.readpdb(pdbfile)
        xresidues=X.residues.keys()
        xresidues.sort()
        for residue in xresidues:
            found=None
            resname=None
            for res2 in array.keys():
                resnumber=pKD_tools.get_resid_from_mut(res2)
                if residue==resnumber:
                    found=1
                    resname=res2
                    break
            
            if not found:
                print '%s was not found' %residue
                for atom in X.residues[residue]:
                    X.atoms[atom]['B-FACTOR']=0.0
            else:
                for atom in X.residues[residue]:
                    X.atoms[atom]['B-FACTOR']=array[resname]*40.0
        X.writepdb('%s.mapped.pdb' %pdbfile)
        return

    #
    # ----------------------
    #

    def parse_pKa_arguments(self,pkas):
        #
        # Parse the -pKas arguments
        #
        import string
        if not pkas:
            self.pkas_missing()
        pkas=string.split(pkas,',')
        self.desired_pKas={}
        self.desired_pKa_changes={}
        for pkaset in pkas:
            #
            # We can have three kinds of designs: =, <, >
            #
            if string.find(pkaset,'=')!=-1:
                pkaset_split=string.split(pkaset,'=')
            elif string.find(pkaset,'<'):
                pkaset_split=string.split(pkaset,'<')
            elif string.find(pkaset,'>'):
                pkaset_split=string.split(pkaset,'>')
            else:
                #
                # No valid qualifier
                #
                print
                print 'Valid qualifiers for pKa changes are: =, < and >'
                print
                raise Exception()
            #
            # Get the residue, its number and the weight
            #
            print pkaset_split
            residue=pkaset_split[0]
            dpKa_weight=pkaset_split[1].split('%')
            self.desired_pKas[residue]={'pKa':dpKa_weight[0],'weight':1.0}
            if len(dpKa_weight)==2:
                self.desired_pKas[residue]['weight']=float(dpKa_weight[1])
        #
        # Test the validity of the resiudes
        #
        self.O._print('',1)
        self.O._print('                Target pKa values',1)
        self.O._print('%12s   %5s     %5s   %5s    %6s'   %('residue','cur. pKa','new pKa','change','weight'),1)
        self.O._print('-------------------------------------------------------',1)
        error=None
        desired_residues=self.desired_pKas.keys()
        desired_residues.sort()
        for residue in desired_residues:
            if self.pkavals.has_key(residue):
                if self.pkavals[residue]['pKa'] or self.params['use_titration_curves']:
                    #
                    # Find the wild type pKa value
                    #
                    pka_value=self.pkavals[residue]['pKa']
                    if not pka_value:
                        pka_value=0.0
                    pka_value=float(pka_value)
                    #
                    # Calculate the desired pKa - in the case of undefined pKa values, this is the dpKa
                    #
                    if self.desired_pKas[residue]['pKa'][0]=='+':
                        self.desired_pKas[residue]['pKa']=pka_value+ \
                                                    float(self.desired_pKas[residue]['pKa'][1:])
                    elif self.desired_pKas[residue]['pKa'][0]=='m':
                        self.desired_pKas[residue]['pKa']=pka_value- \
                                           float(self.desired_pKas[residue]['pKa'][1:])
                    else:
                        self.desired_pKas[residue]['pKa']=float(pka_value)
                    #
                    # Print
                    #
                    self.desired_pKa_changes[residue]={'dpKa':self.desired_pKas[residue]['pKa']-pka_value}
                    self.desired_pKa_changes[residue]['weight']=self.desired_pKas[residue]['weight']
                    self.O._print('%8s    %5.1f     %5.1f     %5.1f      %5.1f'  %(residue,
                                                                                   pka_value,
                                                                                   self.desired_pKas[residue]['pKa'],
                                                                                   self.desired_pKa_changes[residue]['dpKa'],
                                                                                   self.desired_pKa_changes[residue]['weight']),1)
                else:
                    error=1
                    self.O._print('%s does not have a defined pKa values and cannot be redesigned' %residue,1)
            else:
                error=1
                self.O._print('%8s   %5s   %5s  %5s %10s' %(residue,'---','--','--','Residue not found'),1)
        #
        self.O._print('',1)
        self.O._print('',1)
        if error:
            print 'Please correct errors in the desired pKa list'
            print
            print 'Known residues are:'
            residues=self.pkavals.keys()
            residues.sort()
            print residues
            raise Exception()
        return

    #
    # ---------------------
    #



    def filter_mutations(self,ok_mutations):
        """Filter the mutations"""
        #
        # Add neutral mutations
        #
        if self.params['generate_mutations']:
            self.neutral_muts={}
            self.neutral_muts['ASP']='ASN'
            self.neutral_muts['GLU']='GLN'
            self.neutral_muts['HIS']='PHE'
            self.neutral_muts['LYS']='MET'
            self.neutral_muts['ARG']='MET'
            import string
            residues=self.pkavals.keys()
            residues.sort()
            import pKD_tools
            for residue in residues:
                restype=pKD_tools.get_restype_from_titgroup(residue)
                if self.neutral_muts.has_key(restype):
                    neutral_mut='%s:%s' %(residue,self.neutral_muts[restype])
                    self.data[neutral_mut]={'neutral mutation':1}
        #
        # Store the total number of mutations
        #
        deleted={'unhelpful':[],'distance':[],'accessibility':[],'target':[],'total':len(self.data.keys())}
        #
        # Filter the mutations
        #
        #
        # 1.   Do not consider mutations that involve the residues in self.desired_pKa_changes
        # 1.1  Do not allow Asp -> Glu, Lys -> Arg or similar mutations
        # 2.   If mutate_res is true, then keep only residues in that list
        # 3.   Only keep mutations that are in ok_mutations 
        # 4.   See if min_target_dist is set, and exlude mutations accordingly
        #
        self.O._print('Filtering mutations.....',1)
        #
        # Instantiate the pKa_dist class
        #
        import Design_pKa_help
        self.DIST=Design_pKa_help.pKa_dist(self.pdbfile,parent=self,save_file=self.params['save_temp_files'])
        #
        # Prepare for loop
        #
        import string
        delete=[]
        mutations=self.data.keys()
        mutations.sort()
        for mutation in mutations:
            #
            # Make sure we don't mutate the residue we want to design
            #
            import pKD_tools
            new_res=pKD_tools.get_newrestyp_from_mut(mutation)
            mut_res_number=pKD_tools.get_resnum_from_mut(mutation)
            old_res=pKD_tools.get_oldrestyp_from_mut(mutation)
            for residue in self.desired_pKa_changes.keys():
                des_res_number=pKD_tools.get_resid_from_res(residue)
                if des_res_number==mut_res_number:
                    delete.append(mutation)
                    deleted['target'].append(mutation)
                    break
            #
            # Disallow obviously unhelpful mutations
            #
            unhelpful={}
            unhelpful['ASP']='GLU'
            unhelpful['GLU']='ASP'
            unhelpful['LYS']='ARG'
            unhelpful['ARG']='LYS'
            if unhelpful.has_key(old_res):
                if unhelpful[old_res]==new_res:
                    delete.append(mutation)
            #
            # If -mutate_res is set, then we can only mutate those residues
            #
            #if self.params['mutate_res']:
            #    keep=None
            #    for keep_res in self.params['mutate_res']:
            #        keep_res_num=pKD_tools.get_resid_from_res(keep_res)
            #        if keep_res_num==mut_res_number:
            #            keep=1
            #            break
            #    if not keep:
            #        deleted['unhelpful'].append(mutation)
            #        delete.append(mutation)
            #
            # if -accessibility is set, then we can only mutate if the parent residue is in ok_mutations
            #
            if self.params['accessibility']:
                keep=None
                for res in ok_mutations.keys():
                    resnumber=pKD_tools.get_resnum_from_res(res)
                    if resnumber==mut_res_number:
                        keep=1
                    if new_res=='ASN' or new_res=='GLN' or new_res=='PHE':
                        keep=1
                if not keep:
                    deleted['accessibility'].append(mutation)
                    delete.append(mutation)
            #
            # if -distactivesite is set then include that criteria
            #
            if self.params['min_target_dist']>0.0:
                #
                # Make sure that all mutations are at least
                # self.params['min_target_dist'] from all target residues
                #
                dists_ok=1
                for target_res in self.desired_pKa_changes.keys():
                    min_dist=self.DIST.get_min_dist(target_res,mutation)
                    #
                    # If we cannot get the minimum distance then we cannot model the mutation for some reason...
                    #
                    if not min_dist:
                        dists_ok=None
                        break
                    if min_dist<self.params['min_target_dist']:
                        dists_ok=None
                        break
                #
                # Keep or trash?
                #
                if dists_ok:
                    pass
                else:
                    deleted['distance'].append(mutation)
                    delete.append(mutation)
            #
            # We don't want to mutate the residue that we are designing
            #
            import pKD_tools
            for target_res in self.desired_pKa_changes.keys():
                if pKD_tools.get_resid_from_mut(mutation)==pKD_tools.get_resid_from_res(target_res):
                    delete.append(mutation)
        # -----------------------------------------------------
        #
        # Delete the mutations we don't want
        #
        self.mutations={}
        for mut in self.data.keys():
            self.mutations[mut]='ok to use'
        #
        for purge_res in delete:
            if self.mutations.has_key(purge_res):
                del self.mutations[purge_res]
        self.O._print('done filtering',1)
        ok_mutations=self.mutations.keys()
        ok_mutations.sort()
        #
        # If the user specified -list_mutations then do that and exit here
        #
        if self.params['list_mutations']:
            print
            print 'Listing mutations that would have been considered for this design problem'
            print 'Number of mutations selected: %4d' %len(ok_mutations)
            print 'Number of mutations before selection: %4d' %(deleted['total'])
            print 'Cutoff for accessibility: %3d' %(self.params['acclevel'])
            print 'Reasons for exclusion'
            for reason in deleted.keys():
                if reason=='total':
                    continue
                print reason,len(deleted[reason]),deleted[reason]
                print
            print 'Listing all mutations that can be used'
            print ok_mutations
            unique={}
            for mutation in ok_mutations:
                sp=mutation.split(':')
                orgres=sp[0]+':'+sp[1]
                unique[orgres]=1
            print 'Number of wt residues mutated',len(unique.keys())
            #
            # Load the pdbfile
            #
            pdbfilename=self.params['pdb']
            import Protool
            X=Protool.structureIO()
            X.readpdb(pdbfilename)
            perc=float(len(unique.keys()))/float(len(X.residues.keys()))*100.0
            print 'Percent of residues mutated: %5.1f%%' %perc
            raise Exception()
        #
        # Return
        #
        return
    

    #
    # ---------------------
    #

    def score_accessibility_of_parent_residues(self):
        #
        # Get the accessibilities of all residues and
        # use only the ones with a high relative accessibility
        #
        # Notice that this routine only provides a list of ok_mutations. It does
        # not remove any mutations...
        #
        self.O._print('Calculating accessibilities with WHAT IF',2)
        import pKarun.WI_tools, os
        pdbfile=os.path.join(os.getcwd(),self.pdbfile)
        acc_file=pdbfile+'.accs'
        if os.path.isfile(acc_file):
            import cPickle
            fd=open(acc_file)
            accs=cPickle.load(fd)
            fd.close()
        else:
            done=None
            count=0
            accs={}
            accs=pKarun.WI_tools.relative_accessibility(pdbfile)
            fd=open(acc_file,'w')
            import cPickle
            cPickle.dump(accs,fd)
            fd.close()
        #
        # Continue the analysis
        #
        residues=accs.keys()
        residues.sort()
        ok_to_mutate={}
        self.O._print ('done',2)
        for residue in residues:
            if accs[residue]['sum']['rel']>=float(self.params['acclevel']):
                ok_to_mutate[residue]=1
        #
        # Map on B-factors
        #
        if self.params['list_mutations']:
            print ok_to_mutate, len(ok_to_mutate.keys())
            #self.map_on_Bfactors(pdbfile,ok_to_mutate)
            #stop
        self.O._print('Returning from WI acc function',4)
        return ok_to_mutate


    #
    # ---------------------
    #

    def design(self,name):
        #
        # Design everything
        #
        self.initialise(name)
        #
        # Find best solution
        #
        best_solutions,dpKa_dict=self.find_solution(self.params['max_mutations'])
        #
        # Print the summary
        #
        if len(best_solutions)==0:
            print 'No solutions'                
            return {},{}
        #
        # Final ranking
        #
        self.O._print('',1)
        self.O._print('Ranking of best solutions: (lower score = better solution, Emin=0.0)',1)
        printed={}
        count=0
        import string
        filename=''
        if filename=='' and self.params['save_solutions']:
            import tempfile, os
            design_target=string.replace(string.replace(self.params['pKas'],',','|',),':','_')
            fd,filename=tempfile.mkstemp('_pKaDesign',design_target+string.join(self.name,'_'),os.getcwd())
            self.O._print('Saving results in %s' %filename,1)
            os.close(fd)
        #
        # Print results, write the file - and construct an array to hold the results
        #
        results={}
        if filename!='':
            fd=open(filename,'w')
        #
        # Loop over all the solutions
        #
        for sol in best_solutions:
            val=sol[0]
            if abs(val-sol[0])<0.000001:
                if not printed.has_key(str(sol[1])):
                    count=count+1
                    #
                    # Add entry to array
                    #
                    results[count]={'score':val,'mutations':sol[1]}
                    #
                    # Print the result
                    #
                    text='Score: %5.2f. Mutations: ' %val
                    for mutation in sol[1]:
                        if mutation:
                            text=text+mutation+','
                    text=text[:-1]
                    self.O._print(text,1)
                    #
                    # Write the file
                    #
                    printed[str(sol[1])]=1
                    count2=0
                    if filename!='':
                        for mut in sol[1]:
                            if mut:
                                if count2==0:
                                    count2=1
                                    fd.write(mut)
                                else:
                                    fd.write(','+mut)
                            fd.write('\n')
                        for mut in sol[1]:
                            if mut:
                                fd.write(mut+'\n')
            if count==10:
                break
        if filename!='':
            fd.close()
        return results,dpKa_dict

    #
    # ---------------------------------
    #


    def find_solution(self,number_of_mutations):
        """Find the best solution for the design problem"""
        #
        # See if we have a dictionary of previously calculated solutions
        #
        import os
        #
        # Find the best solution
        #
        number_of_mutations=int(number_of_mutations)
        best_solutions=[]
        #
        # Find best solution(s) in order of increasing mutations
        #
        import random
        self.rand=random.Random(7)
        #
        # We test for a maximum of number_of_mutations mutations
        # and we allow no (None) mutations
        #
        self.O._print('',1)
        self.O._print('Searching for the best way to design the pKa value(s) using a maxmimum of %d mutations' %number_of_mutations,1)
        self.O._print('All mutations must be at least %4.1f A from the target residue' %self.params['min_target_dist'],1) 
        self.O._print('',1)
        self.O._print('Statistics',1)
        self.O._print('----------',1)
        #
        # Calculate statistics
        #
        self.O._print('Maximum number of allowed mutations in solution: %d' %(number_of_mutations),1)
        num_mutations=len(self.data.keys())
        self.O._print('Number of different mutations that can be tried: %d' %(num_mutations),1)
        #
        # If no mutations then we can't do anything...
        #
        if num_mutations==0:
            self.O._print('',1)
            self.O._print('No mutations available for this position with the current parameters',1)
            self.O._print('',1)
            return [],{}
        #
        # Calculate the number of permutations
        #
        numbers={}
        import string
        for mut in self.data.keys():
            resnum=string.join(string.split(mut,':')[:-1],':')
            if not numbers.has_key(resnum):
                numbers[resnum]=1
            else:
                numbers[resnum]=numbers[resnum]+1
        self.O._print('Number of different positions: %d' %(len(numbers.keys())),1)
        perms=1
        for x in range(number_of_mutations):
            vals=numbers.values()
            vals.sort()
            sum=0
            for val in vals:
                sum=sum+val
            perms=perms*sum
            #
            # Find out how many muts to remove
            #
            for number in numbers.keys():
                if numbers[number]==vals[-1]:
                    del numbers[number]
                    break
        #
        # Remove redundant solutions
        #
        perms=perms/fact(number_of_mutations)
        self.O._print('Min. number of permutations: %3.2e for %d mutations' %(perms,number_of_mutations),1)
        self.O._print('------------',1)
        #
        # Get a starting state
        #
        self.O._print('Phase I: Getting starting states',1)
        current_mutations,other_solutions=self.get_starting_state()
        self.O._print('Starting state:'+str(current_mutations),1)
        #
        # If we don't have any real mutations in the starting state then
        # we wont get a solution...
        #
        real=None
        for mut in current_mutations:
            if mut:
                real=1
                break
        if not real:
            #
            # No solution
            #
            print 'No single mutation gave an improvement'
            print 'Cannot find a solution'
            print
            return [],{}
        #
        # In Monte Carlo we never, ever use the tabulated data
        # since we already know what the best solution is from
        # the sorting of the scores
        #
        if self.params['tabulated']:
            self.params['tabulated']=None
            #
            # Here we should also filter mutations that do not give a primary or secondary effect
            # if we did tabulated optimisation
            #
        #
        # Get starting energy
        #
        import sys
        if current_mutations!=['wt']:
            oldE,dpKas=self.calculate_pKa_change(current_mutations)
            self.O._print('Energy of starting state: %5.3f' %oldE,3)
        else:
            oldE=None
            print 'The best starting state was the wild type!!!'
            print 'Trying to recover..... (I probably wont)...'
        #
        # Start Monte Carlo loop
        #
        self.O._print('------------',1) 
        import copy
        step=self.params['MCsteps']
        if step>100:
            print 'Step is',self.params
            raise Exception('Too many MC steps')
        self.O._print('Doing %5d MC steps' %step,1)
        #
        # Construct the graveyard from the other solutions
        #
        dpKa_dict={}
        self.O._print('Phase II: Calculating energies for all starting states',3)
        graveyard=[]
        count=len(other_solutions)
        for solution in other_solutions:
            if self.O.level<3:
                streng='\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bSol : %8d' %count
                print streng,
                sys.stdout.flush()
            count=count-1
            if solution==['wt']:
                continue
            self.O._print ('-------------------------------------------------',4)
            XE,dpKasX=self.calculate_pKa_change(solution)
            graveyard=[[XE,solution[:]]]+graveyard
            dpKa_dict[str(solution)]=dpKasX.copy()
        if graveyard==[] and oldE==None:
            return [],{}
        if oldE==None:
            print graveyard
            oldE=graveyard[0][0]
            grv_yrd_sol=str(graveyard[0][1])
            dpKas=dpKa_dict[grv_yrd_sol]
        #
        # Record the starting state
        #
        self.O._print('Phase III: Entering Monte Carlo loop with %3d steps' %step,3)
        dpKa_dict[str(current_mutations)]=dpKas.copy()
        graveyard=[[oldE,current_mutations[:]]]+graveyard
        import random, sys
        while step>0:
            if self.O.level>3:
                streng='\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bStep: %8d' %step
                print streng,
                sys.stdout.flush()
            #
            # Update step counter
            #
            step=step-1
            #
            # Change a random mutation
            #
            new_state=self.change_random_mutation(current_mutations)
            #
            # Clean the mutations - make sure that they are all real
            # mutations, and not just reversions of unsuccessful muts
            #
            new_state=self.clean_state(new_state)
            #
            # Calculate dpKas
            #
            newE,dpKas=self.calculate_pKa_change(new_state)
            #
            # Accept better scores, and worse ones with a probability
            # 
            accept=None
            diff=newE-oldE
            if diff<0.0:
                accept=1
                self.O._print('Accepting diff %5.3f, new state: \n%s\n' %(diff,str(new_state)),4)
            elif diff<1.0:
                # We accept a worse solution based on the randomisation
                tal=0.5*random.random()
                if tal<diff:
                    accept=1
                else:
                    pass
            else:
                pass
            #
            # Accept this state?
            #
            if accept:
                current_mutations=new_state[:]
                oldE=newE 
                #
                # Record the state
                #
                current_mutations.sort()
                graveyard=[[oldE,current_mutations[:]]]+graveyard
                dpKa_dict[str(current_mutations)]=dpKas.copy()
        self.O._print('Finished Monte Carlo loop with %d iterations' %(self.params['MCsteps']),1)
        self.O._print('------------',1)
        self.O._print('Sorting solutions',3)
        #
        # Done with MC steps
        #
        current_mutations.sort()
        graveyard=[[oldE,current_mutations[:]]]+graveyard
        dpKa_dict[str(current_mutations)]=dpKas.copy()
        #
        # Get Energy of the wild type
        #
        wt_score=self.get_wild_type_score()
        self.O._print('Wild type score: %5.3f' %wt_score,1)
        #
        # Sort the graveyard
        #
        sorted_graveyard=[]
        values=[]
        for solution in graveyard:
            values.append(solution[0])
        values.sort()
        recorded={}
        for value in values:
            #
            # Get rid of solutions that are very similar to the wild type
            #
            if abs(value-wt_score)<0.05:
                continue
            #
            # Get rid of worse solutions
            #
            if value>wt_score:
                continue
            #
            # Find the solution
            #
            same_Es=[]
            for solution in graveyard:
                E=solution[0]
                mutations=solution[1]
                if abs(value-E)<0.1:
                    mut_str=str(mutations)
                    if not recorded.has_key(mut_str):
                        recorded[mut_str]=1
                        same_Es.append(solution)
            #
            # Get rid of all the 'None's and sort by number of mutations
            #
            real_sols=[]
            for solution in same_Es:
                clean_sol,num_muts=self.clean_solution(solution)
                real_sols.append([clean_sol,num_muts])
            if real_sols!=[]:
                real_sols.sort(cmp=lambda x, y: cmp(x[1],y[1]))
                #
                # Add the solutions to the graveyard
                #
                for sol in real_sols:
                    sorted_graveyard.append(sol[0])
        #
        # 
        #
        # done
        #
        self.O._print('Done with design job',1)
        return sorted_graveyard,dpKa_dict

    #
    # ----
    #

    def clean_solution(self,solution):
        """Get rid of all non-real mutations"""
        c=[]
        count=0
        for s in solution[1]:
            c.append(s)
            if s:
                count=count+1
        return [solution[0],c],count

    #
    #--------------------------
    #
    
    def get_starting_state(self):
        #
        # Get the starting state for the MC calculations
        #
        # Calculate the net effect of each mutation, - either by phi/ln(10) or
        # by using tabulated values from the MC runs
        #
        scores={}
        possibilities=[]
        #
        # Get dpKa values for each individual mutation
        #
        if not self.params['tabulated']:
            self.params['tabulated']='tempset'
        #
        #
        tmp_dpKas={}
        for mutation in self.mutations.keys():
            if self.data[mutation] is None:
                #
                # No data for this mutation, so even if it is a crg->neutral mutation we should
                # get a rotamer score
                #
                self.data[mutation]={'Rotamer quality':self.get_rotamer_score(mutation)}
            if not self.data[mutation].has_key('Rotamer quality'):
                self.data[mutation]['Rotamer quality']=self.get_rotamer_score(mutation)
            #
            # Now we can check if the rotamer score is good enough
            #
            if not self.data[mutation].has_key('Rotamer quality'):
                print '\n\nNo rotamer quality in sugelm data.\nYou should update this run\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n'
                raise Exception('No rotamer quality score for %s' %mutation)
            if self.data[mutation]['Rotamer quality']<=self.params['mutation_quality']:
                #
                # Include the mutation only if the rotamer quality (the debump score from Protool)
                # is low enough
                #
                score,dpKas=self.calculate_pKa_change([mutation])
                scores[mutation]=score
            import copy
            tmp_dpKas[mutation]=copy.copy(dpKas)
        #
        # Dump all data to tempfiles
        #
        self.pKa_info.update_data()
        #
        # Make sure we go back to normal
        #
        if self.params['tabulated']=='tempset':
            self.params['tabulated']=None
        # --------------------------------------------------------------------------------
        # Remove all mutations that don't give dpKa values for one of the desired residues
        #
        deleted=[]
        for mutation in self.mutations.keys():
            if not tmp_dpKas[mutation]:
                deleted.append(mutation)
                del self.mutations[mutation]
                del scores[mutation]
        if deleted!=[]:
            print
            print 'WARNING!'
            print 'I did not use the following mutations because I could not calculate dpKa values for them'
            print 'This could mean several things:'
            print '1. The mutation is not feasible at this position. Lower the cut-off value that you use for rotamer_quality'
            print '2. The mutation induced a pKa shift so large that the pKa value of your target residue is outside the pH range considered'
            print '   To fix this, use the -use_titration_curves flag'
            print '3. This program is broken and you should submit an error report'
            print
            print 'Mutations excluded:'
            x=0
            for mut in deleted:
                print mut+',',
                x=x+1
                if x>5:
                    print
                    x=0
            print
            print
        # ---------------------------------------------------------------------------------
        #
        # To get the starting state we simply pick the x mutations with the best scores
        #
        # We never pick mutations that are worse than the wild type
        #
        wt_score=self.get_wild_type_score()
        #
        # We try to predict 10 different staring states - each time exluding the n best mutations
        # we do this to get a wider range of starting states to look at
        #
        for round_num in range(10):
            #
            # Initialise the starting state
            #
            these_scores=scores.copy()
            starting_state=[]
            last_energy=wt_score
            while len(starting_state)<self.params['max_mutations']:
                #
                # Get all the scores
                #
                vals=these_scores.values()           
                #
                # If we don't have enough positions then we have to fill up with Nones
                #
                # Sort the scores
                #
                vals.sort()
                vals=vals[round_num:] # Exlude one or more mutations
                if len(vals)==0:
                    starting_state.append(None)
                    continue
                best_score=vals[0]
                #
                # If the best score is very close to the score for the wild type then we cannot
                # get a better solution (anything we do will make it worse)
                #
                if best_score-wt_score>-0.05:
                    starting_state.append(None)
                    continue
                #
                # ok, we can actually get a better score - find the mutation
                #
                found=None
                for mut in these_scores.keys():
                    if these_scores[mut]==best_score:
                        found=mut
                        break
                if not found:
                    raise 'panic - mutation not found'
                #
                # Make sure that we did not already mutate this position
                #
                illegal=None
                import pKD_tools
                for exist_mut in starting_state:
                    if pKD_tools.get_resid_from_mut(exist_mut)==pKD_tools.get_resid_from_mut(found):
                        illegal=1
                        break
                if not illegal:
                    #
                    # Check if we get closer to the preferred energy
                    #
                    new_energy,dpKas=self.calculate_pKa_change(starting_state+[mut])
                    #
                    # We store the state to get more solutions to look at
                    #
                    tmp=starting_state[:]+[mut]
                    tmp.sort()
                    while len(tmp)<self.params['max_mutations']:
                        tmp.append(None)
                    if not tmp in possibilities:
                        possibilities.append([tmp,new_energy])
                    #
                    # Is this a better solution?
                    #
                    if new_energy<last_energy:
                        last_energy=new_energy
                        #
                        # Add the mutation
                        #
                        starting_state.append(mut)
                    #
                    # Make sure we don't add the mutation again
                    #
                    del these_scores[mut]

                else:
                    del these_scores[mut]
        # --------------------------------------------------------
        #
        # Sort all the possibilities
        #
        possibilities.sort(cmp=lambda x, y: cmp(x[1],y[1]))
        refm=[]
        for posi in possibilities:
            refm.append(posi[0])
        #
        # Check that we can actually get a solution
        #
        if len(refm)>0:
            starting_state=refm[0][:]
            starting_state.sort()
        else:
            starting_state=[]
            refm=[]
        #
        # If we design for more than one group, then we need a better starting state
        #
        if len(self.desired_pKa_changes.keys())>1:
            #
            # Make sure we have a starting state
            #
            if len(starting_state)==0:
                for x in range(self.params['max_mutations']):
                    starting_state.append(None)
            return self.get_multdesign_starting_state(starting_state)
        #
        # Finally return everything
        #
        return starting_state,refm[1:6]

    #
    # ---------------------
    #

    def get_multdesign_starting_state(self,starting_state):
        """Get a number of solutions for the design of more than one group"""
        #
        # In the MC - we always use the tabulated version
        #
        if not self.params['tabulated']:
            self.params['tabulated']='tempset'
        #
        # Identify all the mutations that potentially could contribute to the design
        # Those are the mutations that change any target pKa value by at least 0.1
        #
        contributing_mutations=[]
        mutations=self.mutations.keys()
        mutations.sort()
        for mutation in mutations:
            score,dpKas=self.calculate_pKa_change([mutation])
            for val in dpKas.values():
                if abs(val)>0.1:
                    contributing_mutations.append(mutation)
        contributing_mutations.sort()
        print 'Choosing multiple design starting state from',contributing_mutations
        #
        # Do a Monte Carlo run to find the best starting states
        #
        step=25000
        dpKa_dict={}
        current_mutations=starting_state[:]
        current_mutations.sort()
        oldE,dpKas=self.calculate_pKa_change(current_mutations)
        dpKa_dict[str(current_mutations)]=dpKas.copy()
        graveyard=[[oldE,current_mutations[:]]]
        solutions=[str(current_mutations)]
        import random, sys
        print 'Steps: %8d' %step
        while step>0:
            #
            # Update step counter
            #
            step=step-1
            if self.O.level>2:
                streng='\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bStep: %8d' %step
                print streng,
                sys.stdout.flush()
            #
            # Change a random mutation
            #
            new_state=self.change_random_mutation(current_mutations,contributing_mutations)
            #
            # Clean the mutations - make sure that they are all real
            # mutations, and not just reversions of unsuccessful muts
            #
            new_state=self.clean_state(new_state)
            #
            # Calculate dpKas
            #
            newE,dpKas=self.calculate_pKa_change(new_state)
            #
            # Accept better scores, and worse ones with a probability
            # 
            accept=None
            diff=newE-oldE
            if diff<0.0:
                accept=1
            elif diff<1.0:
                # We accept a worse solution based on the randomisation
                tal=0.5*random.random()
                if tal<diff:
                    accept=1
                else:
                    pass
            else:
                pass
            #
            # Accept this state?
            #
            if accept:
                current_mutations=new_state[:]
                if newE<oldE:
                    #
                    # Record the state
                    #
                    current_mutations.sort()
                    if not str(current_mutations) in solutions:
                        graveyard=[[oldE,current_mutations[:]]]+graveyard
                        solutions.append(str(current_mutations))
                #
                # Update the energy
                #
                oldE=newE 
        #
        # Make sure we go back to normal tabulated state
        #
        if self.params['tabulated']=='tempset':
            self.params['tabulated']=None
        #
        # Sort the solutions
        #
        current_mutations.sort()
        if not str(current_mutations) in solutions:
            graveyard=[[oldE,current_mutations[:]]]+graveyard
        #dpKa_dict[str(current_mutations)]=dpKas.copy()
        #
        # Get Energy of the wild type
        #
        wt_score=self.get_wild_type_score()
        graveyard=[[wt_score,['wt']]]+graveyard
        #
        # Sort
        #
        graveyard.sort(cmp=lambda x, y: cmp(x[0],y[0]))
        #
        # Get rid of the scores
        #
        posib=[]
        for score,sol in graveyard:
            if sol=='wt':
                break
            posib.append(sol)
        #
        #
        starting_state=[]
        if len(posib)>0:
            starting_state=posib[0]
        other_sols=posib[min(1,len(posib)):min(20,len(posib))]
        return starting_state,other_sols

    #
    #-----------------------
    #

    def save_summary(self,current_mutations):
        #
        # Saves a summary of the current mutations
        #
        current_mutations.sort()
        #
        # Print the overall summary
        #
        score,dpKas=self.calculate_pKa_change(current_mutations)
        #
        # Print the breakdown of the pKa changes
        #
        des_res=self.desired_pKa_changes.keys()
        des_res.sort()
        text=[]
        text.append('Score :%5.3f' %score)
        text.append('wt pKa \t mut pKa \t dpKa \t desired dpKa \t Error')
        for residue in des_res:
            text2='%6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f' %(self.wt_pKas[residue],
                                                                self.wt_pKas[residue]+dpKas[residue],
                                                                dpKas[residue],
                                                                self.desired_pKa_changes[residue],
                                                                dpKas[residue]-
                                                                self.desired_pKa_changes[residue])
            text.append(text2)
        filename=string.join(current_mutations,'_')
        fd=open(filename,'w')
        import cPickle
        cPickle.dump(text,fd)
        fd.close()
        #
        # Done
        #
        return
       

    #
    # -------------------
    #

    def get_wild_type_score(self):
        """Get the score of the wild type"""
        score=0.0
        import math
        for residue in self.desired_pKa_changes.keys():
            score=score+abs(self.desired_pKa_changes[residue]['dpKa']-0.0)*self.desired_pKa_changes[residue]['weight']
        return score


    #
    # ----------------
    #

    def change_random_mutation(self,current_mutations,mutations_considered=None):
        #
        # Perform a random change
        #
        number_of_mutations=len(current_mutations)
        test_mutations=[]
        done=None
        #
        # Generate the list of mutations
        #
        list_of_mutations=self.mutations.keys()
        list_of_mutations.sort()
        while not done:
            #
            # Change a random mutation
            #
            if number_of_mutations==1:
                position=0
            else:
                position=self.rand.randint(0,number_of_mutations-1)
            #
            # We change this mutation to None or to a different mutation
            # in 20% of cases we go to None
            #
            new_mutation=''
            change_to_None=self.rand.randint(1,10)
            thismut=current_mutations[position]
            if change_to_None>8 and thismut: # Make sure that we are not going from None to None
                newmutation=None
            else:
                #
                # Pick a mutation to change to
                #
                found=None
                while not found:
                    #print
                    #print 'Need to include rotamer scores in Monte Carlo steps'
                    #print
                    newmut_pos=self.rand.randint(0,len(list_of_mutations)-1)
                    newmutation=list_of_mutations[newmut_pos]
                    if mutations_considered:
                        if newmutation in mutations_considered:
                            found=1
                        else:
                            # Not using this mutation
                            pass
                    else:
                        found=1
            #
            # Get the new list
            #
            test_mutations=current_mutations[:]
            test_mutations[position]=newmutation
            #
            # See if configuration is realistic
            # (i.e. that the same residue isn't mutated twice)
            #
            import pKD_tools
            error=None
            for s1 in range(len(test_mutations)-1):
                for s2 in range(s1+1,len(test_mutations)):
                    if test_mutations[s1] and test_mutations[s2]:
                        if pKD_tools.get_resid_from_mut(test_mutations[s1])==pKD_tools.get_resid_from_mut(test_mutations[s2]):
                            error=1
            #
            # Are we done?
            #
            if not error:
                done=1
        #
        # Yep, done
        #
        return test_mutations

    #
    # --------------------------
    #

    def clean_state(self,state):
        #
        # Remove any reversions (i.e. D7D)
        #
        revertants=[]
        new_state=[]
        for mut in state:
            if not self.pkavals.has_key(mut):
                new_state.append(mut)
            else:
                new_state.append(None)
        return new_state

    #
    # ------------------------------
    #

    def get_rotamer_score(self,mutation):
        """Get the rotamer score for the mutation"""
        import Design_pKa_help
        if self.mutfile_names.has_key(mutation):
            mutant_pdbfile=self.mutfile_names[mutation]
            rotamer_score=0.0
        else:
            mutant_pdbfile,rotamer_score=Design_pKa_help.make_mutation(self.pdbfile,mutation)
        if not mutant_pdbfile or rotamer_score is None:
            print '\nModelling failed'
            raise Exception()
        return rotamer_score

    #
    # -----
    #

    def close_stdout(self):
        """Closes std out and dumps output to file"""
        import tempfile, sys
        self.newstdout=tempfile.TemporaryFile()
        self.old_stdout=sys.stdout
        sys.stdout=self.newstdout
        return

    def open_stdout(self):
        """Reopens stdout"""
        import sys
        sys.stdout=self.old_stdout
        self.newstdout.close()
        return
        

#
# -------------------------------
#

class pKD_error(Exception):

    def __init__(self,value,txt=''):
        self.value=value
        self.txt=txt

    def __str__(self):
        return repr(self.value)+repr(self.txt)

#
# ----
#

def print_setup(dict):
    #
    # Print the setup
    #
    keys=dict.keys()
    keys.sort()
    name=[]
    O=verbose(dict['verbose'])
    #
    # Subset selection and calculation of pKa values for higher order mutants
    #
    O._print('**** Selection of titratable groups ****',1)
    O._print('All titratable groups are included at the same level',1)
    O._print('',1)
    O._print('**** Combination of dpKa values ****',1)
    if dict['dpKa_method']=='tabulated':
        O._print('dpKa values are evaluated explicitly for both single, double and higher order mutants',1)
        O._print('',1)
        name.append('Excplicit')
    else:
        O._print('All residues are included in calculation. Higher order dpKa values are found by linear combination',1)
        O._print('',1)
        name.append('Tabulated')
    #
    # Intrinsic pKa recalculation
    #
    O._print('**** Intrinsic pKa recalculation *****',1)
    O._print('Intrinsic pKa values are calculated for mutated residues (if the new residue is titratable)',1)
    if dict['recalc_intpka']:
        O._print('Recalculating intrinsic pKa for target residues if mutation is closer than %5.3fA' %dict['recalc_intpka_dist'],1)
    O._print(' ',1)
    #
    # pKa calculation method
    #
    pKa_calc_alg=None
    O._print('**** pKa calculation method *****',1)
    if dict['dpKa_method']=='TR':
        O._print('Using Tanford-Roxby method for evaluating pKa values',1)
        O._print('pKa values are evaluated using the Tanford-Roxby iteration scheme. pKa values are evaluated for all residues',1)
        O._print('Desolvation and background interaction energies of mutated residues are taken into account',1)
        O._print('',1)
        name.append('Tanford-Roxby')
    elif dict['dpKa_method']=='phidiff':
        O._print('Using dPhi/ln10',1)
        O._print('Delta pKa values are calculated using Phi/ln(10)',1)
        O._print('Desolvation energies and background interaction energies are *not* taken into account',1)
        O._print('',1)
        name.append('phidiff')
    elif dict['dpKa_method']=='MC':
        O._print('Using Monte Carlo sampling for calculating pKa values. pKa values are evaluated for all residues',1)
        O._print('Desolvation and background interaction energies of mutated residues are taken into account',1)
        #O._print('Parameters:',1)
        O._print
        name.append('MC')
    else:
        O._print('No pKa calculation method set - this is a bug',1)
        raise pKD_error('This is a bug')
    return name
    
#
G_Design_instance = None

def run_opt(defaults,new_instance=None,mutant_PDB_files={}):
    """
    Do post-processing of defaults and branch into correct routine
    """
    #
    # Did we get a new design instance?
    #
    global G_Design_instance
    if new_instance:
        G_Design_instance=None
    #
    # Reformat parameters first though...
    #
    params={}
    for key in defaults.keys():
        params[key]=defaults[key][0]
    #
    # Print setup, and ask if ok
    #
    name=print_setup(params)
    if params['calc_dpka'] and params['list_mutations']:
        print
        print 'You cannot specify both -calc_dpka and -list_mutations'
        print
        raise Exception()
    #
    # Start the design process
    #
    if not G_Design_instance:
        print 'Instantiating Design_pKa anew'
        G_Design_instance=Design_pKa(params)
    else:
        # If we have an instance then just copy parameters
        G_Design_instance.params=params.copy()
        G_Design_instance.parse_parameters()
    #
    # Make sure we use pre-modelled PDB files if available
    #
    G_Design_instance.mutfile_names=mutant_PDB_files
    try:
        if params['calc_dpka']:
            G_Design_instance.name=name
            
            results=G_Design_instance.calc_dpKa()
            dpKa_dict=None
        else:
            results,dpKa_dict=G_Design_instance.design(name)
    finally:
        G_Design_instance.pKa_info.update_data()
    #
    # Done
    #
    return results,dpKa_dict

#
# -----------------------------------------
#

def get_defaults():
    #
    # Types for parameter passing:
    # text: text. res_num_list: List of residu numbers + associated number for each
    # number: number
    # res_list: List of residue numbers
    # T/F: True/False: If you give the argument then set to true, otherwise false
    #
    defaults={}
    options,args=parse_options([])
    dicte=eval(str(options))
    #
    for option in dicte.keys():
        defaults[option]=[dicte[option],'junk']
    import string
    defaults['mutations'][0]=string.join(defaults['mutations'][0],',')
    defaults['pKas'][0]=string.join(defaults['pKas'][0],',')
    defaults['sens_limit']=[2.0*defaults['pHstep'][0],'junk']
    defaults['main']=[False,'junk']
    if __name__=="__main__":
        defaults['main'][0]=True
    return defaults.copy()

#
# -----------------------------------------
#

def delete_old_files(pdbfile,delete=None,keep_sugelm=None):
    #
    # Should we delete old files
    #
    if not pdbfile:
        return
    if delete=='N':
        #
        # Not deleting anything
        #
        return
    import os, string
    oldfiles=['.accs','.MC.tabdata','.MC.tabdata.titcurves','.sugelm_data','sugelm''.intpka_data','.Tanford-Roxby.tabdata','.distances','.pdbs']
    for filename in oldfiles:
        rfile=os.path.join(os.getcwd(),pdbfile+filename)
        
        if os.path.isfile(rfile) or os.path.isdir(rfile):
            if not delete:
                while delete!='Y' and delete!='N':
                    print 'Files from a previous run exist in this dir.'
                    delete=string.upper(raw_input('Do you want to delete them? (Y/N) [N] '))
                    if delete=='':
                        delete='N'
            if delete=='Y':
                #
                # We make an exception for the sugelm data
                #
                if filename=='.sugelm_data':
                    if not keep_sugelm:
                        print 'The sugelm data can often be resused if the PDB file is the same.'
                        keep_sugelm=raw_input(' Do you want to keep the SUGELM data? (Y/N) [Y] ')
                    keep_sugelm=string.upper(keep_sugelm)
                    #
                    # Default
                    #
                    if keep_sugelm=='':
                        keep_sugelm='Y'
                    if keep_sugelm=='Y':
                        # Dont delete the sugelm data
                        pass
                    else:
                        # Delete teh sugelm data
                        pKarun.pKa_ostools.rmr(rfile)
                        # Delete sugelm datadir
                        if os.path.isdir('sugelm'):
                            pKarun.pKa_ostools.rmr('sugelm')
                        # Delete models
                        pdbdir=os.path.join(os.getcwd(),pdbfile+'.pdbs')
                        if os.path.isdir(pdbdir):
                            pKarun.pKa_ostools.rmr(pdbdir)
                else:
                    import pKarun.pKa_ostools
                    print 'Deleting: %s' %rfile
                    pKarun.pKa_ostools.rmr(rfile)

            elif delete=='N' or delete=='n':
                #
                # Not deleting
                #
                pass
    #
    # Delete all the sol_dicts and tabdata
    #
    if delete:
        if delete.lower()=='y':
            thisdir=os.getcwd()
            for filename in os.listdir(thisdir):
                rname=os.path.join(thisdir,filename)
                if rname.find('sol_dict')!=-1 or rname.find('tabdata')!=-1:
                    pKarun.pKa_ostools.rmr(rname)
        #
        # If we deleted old files, then also delete the Wild type pKa calcs
        #
        if delete=='Y':
            files=os.listdir(os.getcwd())
            for file in files:
                if file[:7]=='WT_pKas':
                    import pKarun.pKa_ostools
                    pKarun.pKa_ostools.rmr(file)
    return

#
# --------------------------------
#


def parse_options(inputargs):
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]',version='%prog 1.0')
    #
    # Selection of design targets
    #
    parser.add_option('-p','--pdb',dest='pdb',action='store',type='string',default='2lzt.pka.pdb',
                      help='The PDB file. Default: %default')
    parser.add_option('--pKas',dest='pKas',action='append',type='string',
                      help='Target specification of form A:0035:GLU=+2.0 (use m for -). You can define multiple pKa values to redefine. Default=%default',default=[])

    #
    # Options controlling the mode
    #
    parser.add_option('-l','--calc_dpKa',dest='calc_dpka',action='store_true',default=False,
                      help='Calculate dpKa values for mutations specified. Do not design. Default: %default')
    parser.add_option('-b',"--donot_generate_mutations",dest='generate_mutations',action='store_false',
                      help='Do not generate all charge-changing mutations at onset')
    parser.add_option('-g',"--generate_mutations",dest='generate_mutations',action='store_true',
                      help='Generate all charge-changing mutations at onset. Default: %default',default=True)
    parser.add_option('--stability_mode',action='store_true',dest='stability_mode',
                      help='Calculate and design changes in pH-stability profile. Default: %default.\nDesign not yet implemented.',default=False)
    #
    # Options controlling the selection of mutations
    #
    parser.add_option('-d',"--min_target_dist",type='int',dest='min_target_dist',action='store',
                      help='Minimum distance between a mutated atom and the target residue',default=10)
    parser.add_option('-n',"--max_mutations",type='int',dest='max_mutations',action='store',
                      help='Maximum number of mutations to employ',default=6)
    parser.add_option('-m','--mutations',type='string',dest='mutations',action='append',
                      help="Mutations to use in design process of dpKa calculation",default=[])
    parser.add_option('--use_access',dest='accessibility',action='store_true',
                      help='Use accessibility when selecting mutants. Default= %default',default=True)
    parser.add_option('--acc_level',dest='acclevel',type='float',action='store',
                      help='Only mutate residues that have a surface accessibility higher than this value. Deafult=%default',default=30)
    parser.add_option('--MCsteps',dest='MCsteps',type='int',action='store',
                        help='Number of MC steps to perform when combining mutations. Default: %default',default=20)
    #
    # Options controlling recalculation of the intrinsic pKa
    #
    parser.add_option('-i','--recalc_intpka',dest='recalc_intpka',action='store_true',
                      help='Recalculate intrinsic pKa for target residues. Default= %default',default=True)
    parser.add_option('-j','--no_recalc_intpka',dest='recalc_intpka',action='store_false',
                      help='Do not recalculate intrinsic pKa value for target')
    parser.add_option('-e','--recalc_intpka_dist',dest='recalc_intpka_dist',action='store',type='float',
                      help='Mutations closer than this distance to the target group will force a recalculation of the intrinsic pKa of the target. Default: %default A',
                      default=10.0)
    parser.add_option('--mutation_quality',dest='mutation_quality',type='float',default=0.0,
                      help='Maximum bump value accepted for mutations. Default: %default')
    #
    # Options for calculating dpKa values
    #
    parser.add_option('--dpKa_method',dest='dpKa_method',type='string',action='store',
                      help='Method to use for calculating dpKa values.(MC,BM,TR,phidiff) Default= %default',default='MC')
    parser.add_option('--use_titration_curves',dest='use_titration_curves',action='store_true',
                      help='Use integration of titration curve differences to calculate dpKa values. Default= %default',default=True)
    parser.add_option('--pHstart_integrate',dest='pHstart_integrate',type='float',action='store',
                      help='Lowest pH value to use for titration curve difference integration. Default= %default',default=0.1)

    parser.add_option('--pHstop_integrate',dest='pHstop_integrate',type='float',action='store',
                      help='Highest pH value to use for titration curve difference integration. Default= %default',default=12.0)
    parser.add_option('--pKMCsteps',dest='pKMCsteps',type='int',default=200000,
                      help='Number of Monte Carlo steps used when calculating titration curves. Default: %default')
    #'pHstart':[0.1,'number'],'pHstop':[12.0,'number'],'pHstep':[0.05,'number']
    parser.add_option('--pHstart',dest='pHstart',type='float',default=0.1,
                      help='Starting pH value when calculating titration curves. Default: %default')
    parser.add_option('--pHstop',dest='pHstop',type='float',default=12.0,
                      help='End pH value when calculating titration curves. Default: %default')
    parser.add_option('--pHstep',dest='pHstep',type='float',default=0.05,
                      help='pH value step when calculating titration curves. Default: %default')
    #
    # Combining dpKa values
    #
    parser.add_option('--tabulate_dpKas',dest='tabulated',action='store_true',
                       help='Tabulate dpKa values for individual mutations and add these to arrive at final dpKa (when searching). Default=%default',default=False)
    #
    # PBE options
    #
    parser.add_option('--PBEsolver',dest='PBEsolver',type='string',action='store',
                      help='PBE solver to use (DelPhi/APBS/pKaTool[still to be implemented]). Default= %default',default='APBS')
    
    parser.add_option('--ionic_strength',dest='ion',type='float',action='store',
                      help='ionic strength to use in the calculations. Default= %default',default=0.144)
    
    parser.add_option('--indi',dest='indi',type='float',action='store',
                      help='Protein dielectric constant Default= %default',default=8)

    parser.add_option('--exdi',dest='exdi',type='float',action='store',
                      help='solution dielectric constant Default= %default',default=80)
    
    #
    # Other pKa calc options
    #
    parser.add_option('--dbcrit',dest='dbcrit',type='int',action='store',
                      help='Use different dielectric constant if modified B-factor is above this value. Default= %default',default=1000)
    #
    # Ligands
    #
    parser.add_option('--ligands',dest='ligands',type='string',action='append',
                      help='Mol2 file for ligand. Can be specified multiple times.',default=[])
    #
    # Options for calculating the final score and finding the best solution
    #
    parser.add_option('--mutation_penalty',dest='mutation_penalty',type='float',default=0.0,
                      help='Add this penalty multiplied by the number of mutations to a solution. This discourages solutions with many mutations. Default: %default') 
    #
    # General options
    #
    parser.add_option('-v',"--verbose",type='int',dest='verbose',action='store',
                      help='verbosity level. Default: %default',default=1)
    parser.add_option('-t',"--temp_files",dest='save_temp_files',action='store_true',
                      help='Save temporary files. Default: %default',default=True)
    parser.add_option('-y','--no_temp_files',dest='save_temp_files',action='store_false',
                      help='Do not save temporary files.')
    parser.add_option("--save_solutions",dest='save_solutions',action='store_true',
                      help='Save solutions. Default: %default',default=True)
    parser.add_option("--list_mutations",dest='list_mutations',action='store_true',
                      help='List all mutations then exit',default=False)

    (options, args) = parser.parse_args(inputargs)
    return options,args

#
# ----
#

def main(options,args):
    #
    # Start the preprocessing
    #
    defaults={}
    #print options
    options_list=dir(options)
    for option in options_list:
        if option[0]!='_':
            defaults[option]=[getattr(options,option),'junk']
    import string
    defaults['mutations'][0]=string.join(defaults['mutations'][0],',')
    defaults['pKas'][0]=string.join(defaults['pKas'][0],',')
    defaults['sens_limit']=[2.0*defaults['pHstep'][0],'junk']
    #
    # Delete old files if we are the main program
    #
    if __name__=='__main__':
        delete_old_files(defaults['pdb'][0])
    #
    # Start the program
    #
    results,dpKa_dict=run_opt(defaults)
    if not results:
        return
    if defaults['calc_dpka'][0]:
        return
    #
    # Show the final dpKa values for each solution
    #
    print 'Calculating delta pKa values for solutions'
    defaults['calc_dpka']='T'
    sols=results.keys()
    sols.sort()
    for sol in sols:
        mutations=results[sol]['mutations']
        clean=[]
        for mut in mutations:
            if mut:
                clean.append(mut)
        mutations=string.join(clean,'+')
        defaults['mutations']=[mutations]
        defaults['verbose']=[0]
        dpKas,dpKa_dict=run_opt(defaults)
        text='#%2d ' %sol
        this_sol=dpKas.keys()[0]
        targets=dpKas[this_sol].keys()
        for target in targets:
            text=text+'dpKa %9s:%5.2f ' %(target,dpKas[this_sol][target])
        print text
    #
    # Done
    #
    print 'Normal program end.'
    return


#
# ------------------------
#

if __name__=="__main__":
    print
    print 'pKD - redesign of protein pKa values'
    print 'Copyright Jens Erik Nielsen, 2005-2010 All rights reserved'
    print 'Email: Jens.Nielsen@ucd.ie http://enzyme.ucd.ie'
    print
    import sys, os
    options,args=parse_options(sys.argv[1:])
    main(options,args)
    print
    print
    print '---------------------------------------------------------------'
    print 'Please cite:'
    print 'Re-designing protein pKa values'
    print 'Tynan-Connolly BM, Nielsen JE'
    print 'Protein science 2007 Feb;16(2):239-49. Epub 2006 Dec 22.'
    print ' '
    print 'pKD: Re-Designing protein pKa values'
    print 'Tynan-Connolly BM, Nielsen JE'
    print 'Nucleic Acids Research 2006 Jul 1;34(Web Server issue):W48-51.'
