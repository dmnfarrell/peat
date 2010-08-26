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

import pKaTool.pKa_calc
import pKarun
from pKarun.pKa_utility_functions import *
from output import *

class pKa_calculation_class:
    """
    # Manages many small pKa calculations with the pKa_MC code
    """
    
    def __init__(self,pdbfile,pKa_info=None,params=None,parent=None):
        #
        # Init
        #
        # Set parameters
        #
        self.verbose=params['verbose']
        self.O=verbose(self.verbose)
        self.MCsteps=params['pKMCsteps']
        self.pHstart=params['pHstart']
        self.pHstop=params['pHstop']
        self.pHstep=params['pHstep']
        #
        # Instance of pKa_info class from Design_pKa_help.py
        #
        self.pKa_info=pKa_info 
        #
        # Default pKa calc method is Monte_Carlo with C++
        #
        self.pdbfile=pdbfile
        self.MC=pKaTool.pKa_calc.Monte_Carlo_CPP(self.pdbfile)
        #
        # Initialise all
        #
        self.initall() # Calls routine that reads all the wt files from disk (should be made to read from memory)
        self.parent=parent

        return

    #
    # Choose your pKa calculation method
    #
    import pKaTool.pKa_calc
    def set_MC_CPP(self):
        self.MC=pKaTool.pKa_calc.Monte_Carlo_CPP(self.pdbfile)
        self.initall()
        return

    def set_MC_Python(self):
        self.MC=pKaTool.pKa_calc.Monte_Carlo(self.pdbfile)
        self.initall()
        return

    def set_Tanford_Roxby(self):
        self.MC=pKaTool.pKa_calc.Tanford_Roxby(self.pdbfile)
        self.initall()
        return

    def set_Boltzmann(self):
        self.MC=pKaTool.pKa_calc.Boltzmann(self.pdbfile)
        self.initall()
        return

    #
    # -------------------------
    #

    def initall(self):
        #
        # Set variables
        #
        self.groups=[]
        self.essential_groups={}
        #
        # Load data and calculate intrinsic pKa values
        #
        self.loadall()
        self.MC.calc_intpKas()
        return

    #
    # -------------------------
    #

    def loadall(self):
        #
        # Load all datafiles
        #
        self.MC.read_desolv()
        self.MC.read_backgr()
        self.MC.read_matrix()
        self.groups=self.MC.desolv.keys()
        self.groups.sort()
        self.reporter_groups=[]
        return

    #
    # -------------------------
    #

    def set_reporter_groups(self,residues):
        self.reporter_groups=residues[:]
        return

    #
    # ---------------
    #

    def no_network(self):
        self.deleted_groups={}
        return

    #
    # ----
    #

    def remove_interactions_with(self,groups=[]):
        """Remove all interactions with the groups"""
        for group1 in self.MC.matrix.keys():
            for group2 in self.MC.matrix[group1].keys():
                if group1 in groups or group2 in groups:
                    self.MC.matrix[group1][group2]=[0.0,0.0,0.0,0.0]
        return

    #
    # -------------
    #

    def find_network(self,residues):
        #
        # Find the essential network for designing our pKa values
        #
        # First we check that we have all of them in the wt
        #
        for group in residues:
            if not self.MC.matrix.has_key(group):
                print
                print 'I do not know this group: %s' %group
                print
                raise Exception()
            self.essential_groups[group]=1
        #
        # Find two layers like in Nielsen & McCammon, Pro.Sci. 2004
        #
        # Here we set the selection cut-off
        #
        for int_cutoff in [1,2]:
            # Selection cutoff is 0.1 and 0.2
            cutoff=float(int_cutoff)/10.0
            addgroups={}
            for group in self.essential_groups.keys():
                for group2 in self.MC.matrix[group].keys():
                    intene=self.MC.matrix[group][group2][0]
                    #if abs(intene)>=cutoff:
                    addgroups[group2]=1
            #
            # Add the groups to essential groups
            #
            for add_group in addgroups.keys():
                self.essential_groups[add_group]=1
        #
        # Construct a new matrix with only the essential residues
        # Update the backgr and desolv matrices too
        #
        delete_groups=[]
        for group1 in self.MC.matrix.keys():
            if not self.essential_groups.has_key(group1):
                delete_groups.append(group1)
        #
        # Keep a log of the groups, so that we can add them back again
        # if needed 
        #
        self.deleted_groups={}
        for group in delete_groups:
            self.deleted_groups[group]=self.MC.matrix[group].copy()
            del self.MC.matrix[group]
            del self.MC.backgr[group]
            del self.MC.desolv[group]
        #
        # Delete the groups in the second layer
        #
        for group1 in self.MC.matrix.keys():
            delete_groups=[]
            for group2 in self.MC.matrix.keys():
                if not self.essential_groups.has_key(group2):
                    delete_groups.append(group2)
            for group2 in delete_groups:
                del self.MC.matrix[group][group2]
        #
        # Done
        #
        return self.MC.matrix.keys()

    #
    # -------------
    #

    def calc_wt_pkas(self,dummy=None):
        #
        # Calculate the pKa values for the wt
        #
        print 'pHrange: step: %.2f start: %.2f, end: %.2f. MCsteps: %d ' %(self.pHstep,self.pHstart,self.pHstop,self.MCsteps)
        wtpkas,populations=self.MC.calc_pKas(self.MCsteps,self.pHstep,self.pHstart,self.pHstop,self.verbose)
        return wtpkas,self.MC.prot_states.copy()

    #
    # --------------
    #

    def calc_MC_mut_pKas(self,mutations=[],data={},atomdata={},wt_pKas=[],reporter_residues=[],require_pKa=[],for_tabulated=None):
        """
        Calculate pKa values for a mutant
        """
        #
        # Copy mutfile_names from parent
        #
        if hasattr(self.parent,'mutfile_names'):
            self.mutfile_names=self.parent.mutfile_names.copy()
        #
        # Make backups of files 
        #
        self.backgr=self.MC.backgr.copy()
        self.desolv=self.MC.desolv.copy()
        import copy
        self.matrix=copy.deepcopy(self.MC.matrix)
        self.org_matrix=copy.deepcopy(self.MC.matrix)
        mkeys=self.matrix.keys()
        mkeys.sort()
        #
        # Fill in the new interaction energies
        #
        self.MC.matrix,delete_groups=self.add_mutations_to_matrix(mutations=mutations,
                                                                  target_array=self.MC.matrix,
                                                                  source_array=data,
                                                                  source_atom_array=atomdata)
        if self.MC.matrix=='SS-bridge':
            #
            # Before we return we have to restore the arrays again..
            #
            import copy
            self.MC.backgr=self.backgr.copy()
            self.MC.desolv=self.desolv.copy()
            self.MC.matrix=copy.deepcopy(self.matrix)
            return 'SS-bridge','SS-bridge'
        #
        # Delete groups from matrix and backgr
        #
        for group in delete_groups:
            if self.MC.backgr.has_key(group):
                del self.MC.backgr[group]
                del self.MC.desolv[group]
       
        #
        # Calculate the desolvation and background interaction energies for the mutations
        # The pKa_info class keeps track of these...
        # Fill in the new values
        #
        self.O._print('Getting intrinsic pKa data',4)
        real_mutations=0
        real_mutation_list=[]
        #
        # Array for holding dbackgr and ddesolv values
        #
        d_intpka_vals={}
        #
        # Loop over all the mutations
        #
        for mutation in mutations:
            #
            # Get the new residue
            #
            if not mutation or mutation=='None':
                continue
            import pKD_tools
            new_residue=pKD_tools.get_newres_from_mut(mutation)
            if new_residue:
                real_mutations=real_mutations+1
                real_mutation_list.append(mutation)
                if is_normally_titratable(new_residue):
                    desolv,backgr=self.pKa_info.get_desolv_and_backgr_mut(mutation)
                    if desolv is None or backgr is None:
                        #
                        # Could not get background or desolvation information
                        #
                        # Remember to restore arrays...
                        #
                        self.MC.backgr=self.backgr.copy()
                        self.MC.desolv=self.desolv.copy()
                        import copy
                        self.MC.matrix=copy.deepcopy(self.matrix)
                        self.O._print('Intpka recalculation failed for %s' %mutation,3)
                        return None,None
                    self.O._print('%s, background: %.2f, desolvation: %.2f. net intr dpKa: %.2f' %(mutation,backgr,desolv,(backgr+desolv)/2.302),3)
                    import pKD_tools
                    new_residue=pKD_tools.get_newres_from_mut(mutation)
                    self.MC.backgr[new_residue]=backgr
                    self.MC.desolv[new_residue]=desolv
        #
        # done with getting intpkas for mutated residues
        #
        # Check if we have everything
        #
        need_intpka=[]
        for group in self.MC.matrix.keys():
            if not self.MC.backgr.has_key(group) or not self.MC.desolv.has_key(group):
                need_intpka.append(group)
        if len(need_intpka)>0:
            #
            # ok, something is missing, and hat's not allowed
            #
            raise Exception,'could not find backgr and desolv for one of these groups: %s' %(str(need_intpka))
        #
        # Now recalculate intrinsic pKa values
        #
        for mutation in mutations:
            if not mutation or mutation=='None':
                continue
            #
            #
            # See if the mutation is so close to a reporter group that
            # we have to recalculate the intrinsic pKa value for the reporter group
            #
            if self.parent.params['recalc_intpka']:
                import Design_pKa_help
                DIST=Design_pKa_help.pKa_dist(self.pdbfile,parent=self,save_file=self.parent.params['save_temp_files'])
                for group in require_pKa:
                    #
                    # If we mutate the residue then skip
                    #
                    import pKD_tools
                    if pKD_tools.get_resid_from_mut(mutation)==pKD_tools.get_resid_from_titgroup(group):
                        continue
                    #
                    # Get the distance
                    #
                    dist=DIST.get_min_dist(group,mutation)
                    if dist<=self.parent.params['recalc_intpka_dist']:
                        print '\n\n***Recalculating intrinsic pKa value for %s due to mutation %s, distance: %5.2f A' %(group,mutation,dist)
                        desolv,backgr=self.pKa_info.get_desolv_and_backgr_mut(mutation,group,measuring_mutated_residue=None)
                        #
                        # Insert the energies
                        #
                        if self.MC.backgr.has_key(group):
                            print 'Mutation: %15s denergies for group %s: desolv: %5.2f, backgr: %5.2f' %(mutation,group,desolv-self.MC.desolv[group],backgr-self.MC.backgr[group])
                            if not d_intpka_vals.has_key(mutation):
                                d_intpka_vals[mutation]={}
                            if not d_intpka_vals[mutation].has_key(group):
                                d_intpka_vals[mutation][group]={}
                            d_intpka_vals[mutation][group]['backgr']=backgr-self.MC.backgr[group]
                            d_intpka_vals[mutation][group]['desolv']=desolv-self.MC.desolv[group]
                        else:
                            print 'Arrrg'
                            groups= self.MC.backgr.keys()
                            groups.sort()
                            print 'I know of these groups',groups
                            raise Exception('Could not change intpka for group: %s' %group)
                    else:
                        print 'Mutation: %s is too far away (%5.2f A) from group %s to recalculate intrinsic pKa' %(mutation,dist,group)
                        if not d_intpka_vals.has_key(mutation):
                            d_intpka_vals[mutation]={}
                        if not d_intpka_vals[mutation].has_key(group):
                            d_intpka_vals[mutation][group]={}
                        d_intpka_vals[mutation][group]['backgr']=0.0
                        d_intpka_vals[mutation][group]['desolv']=0.0
        #
        # If we recalculate intrinsic pKa values, then add up all the dbackgr and ddesolv values
        #
        if self.parent.params['recalc_intpka']:
            for group in require_pKa:
                #
                # If we mutated the residue then skip
                #
                skip_group=False
                for mutation in mutations:
                    if mutation:
                        import pKD_tools
                        if pKD_tools.get_resid_from_mut(mutation)==pKD_tools.get_resid_from_titgroup(group):
                            skip_group=True
                            continue
                if skip_group:
                    continue
                #
                # Otherwise adjust background and desolv energies
                #
                ddesolv=0.0
                dbackgr=0.0
                for mutation in mutations:
                    if not mutation or mutation=='None':
                        continue
                    print 'mutation: %s, group: %s, ddesolv: %5.3f, dbackgr: %5.3f' %(mutation,group,d_intpka_vals[mutation][group]['desolv'],d_intpka_vals[mutation][group]['backgr'])
                    ddesolv=ddesolv+d_intpka_vals[mutation][group]['desolv']
                    dbackgr=dbackgr+d_intpka_vals[mutation][group]['backgr']
                self.MC.desolv[group]=self.MC.desolv[group]+ddesolv
                self.MC.backgr[group]=self.MC.backgr[group]+dbackgr
        #
        # Now calculate the new pKa value(s)
        #
        # Loop until we get a pKa value
        #
        got_all_reporter_pKa_values=None
        step=2.0
        while not got_all_reporter_pKa_values:
            #
            # We try to guess the approximate value of the new pKa value and save time by limiting
            # the pH range we're calculating charges in.
            #
            # We do this only if we calculate for more than one mutation, otherwise
            # we want to calculate for the full range so we can store the results in tabdata
            #
            # We also do not do this when calculating dpKa values using titration curve integration
            #
            pHstep=self.pHstep
            pHstop=self.pHstop
            pHstart=self.pHstart
            if len(wt_pKas)==1 and not for_tabulated and not self.parent.params['use_titration_curves']:
                wt_pKa=float(int(10.0*wt_pKas[0]))/10.0
                pHstart=wt_pKa-1.0*step
                pHstop=wt_pKa+1.0*step
                #
                # Make sure we don't go outside 0.1 - 12.0
                #
                pHstart=max(0.1,pHstart)
                pHstop=min(12.0,pHstop)
            else:
                wt_pKa=-1.0
            #
            # If we are not to combine the calculated pKa value with anything, then we can save some time
            #
            if real_mutations>1 and not self.parent.params['use_titration_curves']:
                pHstep=2.0*self.pHstep
            #
            # Do the pKa calc
            #
            self.O._print('Calculating pKa values....',text_level=4,tab=1)
            self.O._print('Calculating for %s' %(str(real_mutation_list)),text_level=4,tab=1)
            self.mut_pKas,populations=self.MC.calc_pKas(mcsteps=self.MCsteps,
                                            phstep=pHstep,
                                            phstart=pHstart,
                                            phend=pHstop,
                                            verbose=self.verbose,
                                            complete_pka=None)
            prot_states=self.MC.prot_states.copy()
            self.O._print('Done',text_level=4,tab=1)
            #
            # Check that we got all the pKa values we need
            # We skip this step if we use titration curves.
            #
            got_all_reporter_pKa_values=1
            if not self.parent.params['use_titration_curves']:
                for group in require_pKa:
                    if not self.mut_pKas[group] or self.mut_pKas[group]<-90.0:
                        got_all_reporter_pKa_values=None
            #
            # If we are still missing pKa values then increase the size of the pH range
            #
            step=step+1.0
            #
            # No more than five steps away before we give up
            #
            if step>5.0:
                got_all_reporter_pKa_values=1
        #
        # Done with pKa calcs
        #
        # Before we return we have to restore the arrays again..
        #
        import copy
        self.MC.backgr=self.backgr.copy()
        self.MC.desolv=self.desolv.copy()

        self.MC.matrix=copy.deepcopy(self.matrix)
        #
        # Done
        #
        return self.mut_pKas,prot_states

    #
    # ==================
    #

    def get_source_array_entry(self,group,source_array,wt_groups):
        """Loop over all the entries in source array, and see if one of them produces this group"""
        #
        # Check if this is a wild type group
        #
        if group in wt_groups:
            return None
        #
        # Not a wild type group. Do we have a mutation
        #
        muts=source_array.keys()
        muts.sort()
        for mutation in muts:
            if len(source_array[mutation].keys())>1: # Do we have more than the Rotamer quality entry
                import pKD_tools
                new_residue=pKD_tools.get_newres_from_mut(mutation)
                if new_residue==group:
                    return mutation
        return 'Missing'

    #
    # ------
    #
    def calc_matrix(self,group):
        """Calculate the charge-charge interaction matrix for a mutation"""
        #
        # Find the mutation
        #
        import pKD_tools
        resid=pKD_tools.get_resid_from_res(group)
        new_resname=pKD_tools.get_restype_from_titgroup(group)
        pdb_resname=self.parent.PI.resname(resid)
        if pdb_resname==new_resname or new_resname in ['CTERM','NTERM',None]:
            return None,None,None
        mutation='%s:%s:%s' %(resid,pdb_resname,new_resname)
        #
        # Calculate the matrix
        #
        import Design_pKa_help
        if self.parent.mutfile_names.has_key(mutation):
            newfilename=self.parent.mutfile_names[mutation]
            bump_score=0.0
            print 'Using PDB file:',newfilename
        else:
            newfilename,bump_score=Design_pKa_help.make_mutation(pdbfile=self.pdbfile,mutation=mutation,topdir=self.parent.topdir)
        if bump_score<=self.parent.params['mutation_quality'] and not bump_score is False and not bump_score is None:
            #
            # Re-calculate
            #
            print 'Calculating matrix for %s' %mutation
            self.parent.close_stdout()
            matrix,atomdata=Design_pKa_help.run_WI_pKa_calc(pdbfile=newfilename,
                                                            mutated_res='%s:%s' %(resid,new_resname),
                                                            pKarun_params=self.parent.pKarun_params,
                                                            do_matrix=True)
            self.parent.open_stdout()
            return matrix,atomdata,mutation
        else:
            return None,None,None

    #
    # ----
    #
            
    def add_mutations_to_matrix(self,mutations,target_array,source_array,source_atom_array):
        """
        # Add the mutations in mutations to the matrix "target_array".
        # Find new interactions in source_array
        """
        import string
        #
        # Store the wild type groups
        #
        wt_groups=target_array.keys()
        #
        # Make sure that we print something sensible
        #
        pr_mutations=[]
        for mut in mutations:
            if mut:
                pr_mutations.append(mut)
        if pr_mutations==[]:
            pr_mutations=['No mutations']
        self.O._print('I am calculating for these mutations: '+string.join(pr_mutations,'+')+'...',4,0)
        #
        # Now find the best solution
        #
        delete_groups=[]
        self.O._print('Phase I: Delete original residues if present in matrix',4)
        for mutation in mutations:
            if mutation is None:
                continue
            #
            # Get the old titratable group
            #
            import pKD_tools
            resid=pKD_tools.get_resid_from_mut(mutation)
            new_residue=pKD_tools.get_newres_from_mut(mutation)
            for group in target_array.keys():
                resid_titgroup=pKD_tools.get_resid_from_res(group)
                titgroup_type=pKD_tools.get_titgroup_type_from_titgroup(group)
                restype=pKD_tools.get_restype_from_titgroup(group)
                #
                # Delete if the residue was changed 
                #
                if resid==resid_titgroup: #restype is None for CTERM and NTERM
                    if restype!=new_residue and restype:
                        delete_groups.append(group)
                #
                # Look in all partner groups
                #
                for group2 in target_array[group].keys():
                    resid_titgroup=pKD_tools.get_resid_from_res(group2)
                    titgroup_type=pKD_tools.get_titgroup_type_from_titgroup(group2)
                    restype=pKD_tools.get_restype_from_titgroup(group2)
                    if resid==resid_titgroup:
                        if restype!=new_residue and restype:
                           delete_groups.append(group2)
        #
        # Delete, delete, delete
        #
        for group in delete_groups:
            if target_array.has_key(group):
                del target_array[group]
        for group in target_array.keys():
            for group2 in delete_groups:
                if target_array[group].has_key(group2):
                    del target_array[group][group2]
        #
        # Fill in the new values
        #
        self.O._print('Phase II: Add placeholders for new residues',4)
        for mutation in mutations:
            if mutation:
                import pKD_tools
                newresidue=pKD_tools.get_newres_from_mut(mutation)
                if newresidue:
                    #
                    # Only include charged residues
                    #
                    if is_normally_titratable(newresidue):
                        target_array[newresidue]={}
                        self.O._print('Added to matrix:'+newresidue,4,tab=1)
                else:
                    raise Exception,'something is seriously wrong'
        #
        # Fill in the interactions
        #
        self.O._print('Phase III: Add all missing interactions...',4)
        null=[0.0,0.0,0.0,0.0]
        could_not_fill=[]
        import copy
        target_groups=target_array.keys()
        target_groups.sort()
        #
        # Start finding the interactions
        #
        import types
        self.O._print('Finding missing interactions for:',4)
        for group in target_groups:
            self.O._print(group,4,tab=1)
            for group2 in target_groups:
                if not target_array[group].has_key(group2):
                    source_name1=group
                    source_name2=group2
                    #
                    # See if we can find an interaction
                    #
                    filled=None
                    #
                    # If we don't have the mutation in source array then check if we should make it
                    #
                    mutation_producing_group1=self.get_source_array_entry(group,source_array,wt_groups)
                    mutation_producing_group2=self.get_source_array_entry(group2,source_array,wt_groups)
                    for this_group,mut in [[group,mutation_producing_group1],[group2,mutation_producing_group2]]:
                        if mut=='Missing':
                            print 'Calculating matrix for', this_group
                            matrix,atom_array,mutation=self.calc_matrix(this_group)
                            if matrix=='SS-bridge':
                                return matrix,matrix
                            if matrix and atom_array:
                                rotq=source_array[mutation]['Rotamer quality']
                                source_array[mutation]=matrix
                                source_array[mutation]['Rotamer quality']=rotq
                                source_atom_array[mutation]=atom_array
                    mutation_producing_group1=self.get_source_array_entry(group,source_array,wt_groups)
                    mutation_producing_group2=self.get_source_array_entry(group2,source_array,wt_groups)
                    if not mutation_producing_group1 and not mutation_producing_group2:
                        print 'This interaction should not be missing'
                        print 'G1/G2',group,group2
                        raise Exception('Fatal error in code')
                    #
                    # Copy interaction to target[group][group2]
                    #
                    if source_array.has_key(mutation_producing_group1):
                        if source_array[mutation_producing_group1].has_key(source_name2):
                            source=source_array[mutation_producing_group1][source_name2]
                            if type(source) is types.FloatType:
                                source=[source,0.0,0.0,0.0]
                                #print 'Converted source to list1'
                            target_array[group][group2]=source
                            filled=1
                        elif source_array[mutation_producing_group1].has_key(group2):
                            source=source_array[mutation_producing_group1][group2]
                            if type(source) is types.FloatType:
                                source=[source,0.0,0.0,0.0]
                                #print 'Converted source to list2'
                            target_array[group][group2]=source
                            filled=1
                    #
                    # Copy interaction to target[group][group2]
                    #

                    if source_array.has_key(mutation_producing_group2) and not filled:
                        if source_array[mutation_producing_group2].has_key(source_name1):
                            source=source_array[mutation_producing_group2][source_name1]
                            if type(source) is types.FloatType:
                                source=[source,0.0,0.0,0.0]
                                #print 'Converted source to list3'
                            target_array[group][group2]=source
                            filled=1
                        elif source_array[mutation_producing_group2].has_key(group):
                            source=source_array[mutation_producing_group2][group]
                            import types
                            if type(source) is types.FloatType:
                                source=[source,0.0,0.0,0.0]
                                #print 'Converted source to list4'
                            target_array[group][group2]=source
                            filled=1
                    #
                    # Maybe we can fill this interaction from the original matrix?
                    #
                    if mutation_producing_group1 is None and mutation_producing_group2 is None:
                        if self.org_matrix.has_key(group):
                            if self.org_matrix[group].has_key(source_name2):
                                #
                                # Org matrix has all values
                                #
                                target_array[group][group2]=copy.copy(
                                    self.org_matrix[group][source_name2])
                                filled=1
                    #
                    # Log interactions that could not be filled
                    #
                    if not filled and group!=group2:
                        could_not_fill.append([group,group2])
                    elif not filled and group==group2:
                        #
                        # set diagonal elements to null
                        #
                        target_array[group][group2]=null
                        target_array[group2][group]=null
        #
        # ------------------------------------------------------------------------------
        #
        # Fill unresolved interactions
        #
        self.O._print('\nPhase IV: Filling unresolved interactions',4)
        for unresolved in could_not_fill:
            #
            # If we could not find an interaction energy in the residue-residue interaction data
            # then look for it in the atom_based array
            #
            group=unresolved[0]
            group2=unresolved[1]
            import string
            filled=None
            mutation_producing_group1=self.get_source_array_entry(group,source_atom_array,wt_groups)
            if source_atom_array.has_key(mutation_producing_group1) and not filled:
                #
                # Trying to find interaction [group][group22]
                #
                if group!=group2:
                    self.O._print('trying to fill from atomdata1',4)
                    import pKD_tools
                    resid2=pKD_tools.get_resid_from_res(group2)
                    type1=pKD_tools.get_titgroup_type_from_titgroup(group)
                    type2=pKD_tools.get_titgroup_type_from_titgroup(group2)
                    if source_atom_array[mutation_producing_group1].has_key(resid2):
                        energy=self.get_correct_energy(source_atom_array
                                                       [mutation_producing_group1][resid2],type1,type2)
                        target_array[group][group2]=[energy,0.0,0.0,0.0]
                        self.O._print('Added energy %5.3f for %s and %s' %(energy,group,group2),4)
                        filled=1
            #
            # Array 2
            #
            mutation_producing_group2=self.get_source_array_entry(group2,source_atom_array,wt_groups)
            if source_atom_array.has_key(mutation_producing_group2) and not filled:
                #
                # Trying to find interaction [group2][group]
                #
                if group!=group2:
                    self.O._print('trying to fill from atomdata2',4)
                    import pKD_tools
                    resid1=pKD_tools.get_resid_from_res(group)
                    type1=pKD_tools.get_titgroup_type_from_titgroup(group)
                    type2=pKD_tools.get_titgroup_type_from_titgroup(group2)
                    if source_atom_array[mutation_producing_group2].has_key(resid1):
                        energy=self.get_correct_energy(source_atom_array
                                                       [mutation_producing_group2][resid1],type1,type2)
                        target_array[group2][group]=[energy,0.0,0.0,0.0]
                        self.O._print('Added energy %5.3f for %s and %s' %(energy,group2,group),4)
                        filled=1
            #
            # Is it the same group? - if so then zero energy of course...
            #
            if group==group2:
                target_array[group][group]=[0.0,0.0,0.0,0.0]
                filled=1
            #
            # If we could not fill an interaction then it's because we have
            # never calculated a phimap for the mutation
            #
            if not filled:
                newg1=pKD_tools.get_newrestyp_from_mut(mutation_producing_group1)
                newg2=pKD_tools.get_newrestyp_from_mut(mutation_producing_group2)
                if newg1=='CYS' or newg2=='CYS':
                    return 'SS-bridge','SS-bridge'
                print 'Could not fill: %s-%s' %(group,group2)
                print
                print 'I cannot predict dpKa values for this mutation because modelling of one'
                print 'or more mutations produced bumps in the structure.'
                print 'If you want to allow bumps when modelling mutations, then adjust the value of '
                print 'the mutation_quality parameter.'
                print 'The current maximum bump for a mutation is %5.3f' %(self.parent.params['mutation_quality'])
                print
                print 'In a future version we will include a structure relaxation protocol that can be used to model'
                print '"difficult" mutations'
                print
                print 'Bump values for mutations in this run:'
                import Design_pKa_help
                for mutation in mutations:
                    pdbfile,bump_value=Design_pKa_help.make_mutation(pdbfile=self.parent.pdbfile,mutation=mutation,topdir=self.parent.topdir)
                    print '%10s, bump: %5.3f' %(mutation,bump_value)
                print
                raise Exception('Cannot model mutant %s' %mutation)
        #
        # Do another pass with delete groups to make sure we didn't fill in an interaction we shouldn't have
        #
        self.O._print('Phase V: Making sure that all unwanted groups have been deleted',4)
        for group in delete_groups:
            if target_array.has_key(group):
                del target_array[group]
        for group in target_array.keys():
             for group2 in delete_groups:
                if target_array[group].has_key(group2):
                    del target_array[group][group2]
        #
        # Return results
        #
        self.O._print('Copying return arrays',4)
        return target_array.copy(), delete_groups

    #
    # -------------------
    #
        
    def get_correct_energy(self,dict,type1,type2):
        #
        # Get the interaction energy with another residue
        # We sum over all atoms, but it should be minor error
        #
        count=0.0
        pot_sum=0.0
        #
        # dict is source_atom_array and holds potentials - not energies
        #
        for atom in dict.keys():
            pot_sum=pot_sum+dict[atom]
            count=count+1.0
        #
        # We have the total potential - get the average potential of each atom
        # and multiply that will the charge to get the full energy
        # This is not exactly equal to the full energy but a good approximation
        #
        avg_pot=pot_sum/count
        #
        # Get the correct sign
        #
        ene=avg_pot*float(acibas(type1)*acibas(type2))
        return ene
        
