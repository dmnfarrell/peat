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

class dpKa:

    def __init__(self,pdbfile,user_params={},reporter_groups=None):
        """Read the pdbfile and the wild type pKa values
        Define the reporter groups"""
        # Topdir
        import os
        self.topdir=os.getcwd()
        #
        # Set the PDB file name
        #
        self.pdbfile=os.path.join(self.topdir,pdbfile)
        import Protool
        self.PI=Protool.structureIO()
        self.PI.readpdb(self.pdbfile)
        #
        # Get the defaults from Design_pKa
        #
        import Design_pKa
        defaults=Design_pKa.get_defaults()
        self.params={}
        for key in defaults.keys():
            self.params[key]=defaults[key][0]
        self.params['pKMCsteps']=20000
        self.params['recalc_intpka']=True
        self.params['recalc_intpka_dist']=20.0
        self.params['save_temp_files']=False
        #
        # Change the values that the user specified
        #
        for u_par in user_params.keys():
            self.params[u_par]=user_params[u_par]
        #
        # Set the pKa calc parms as a subset
        #
        self.pKarun_params={}
        include=['indi','dbcrit','ion','allow_unknown_atoms','unknown_crg','unknown_rad']
        import copy
        for key in include:
            self.pKarun_params[key]=copy.copy(self.params[key])
        #
        # Set the output handler
        #
        self.O=Design_pKa.verbose(self.params['verbose'])
        #
        # Instantiate the pKa_info class and instruct it to save nothing
        #
        import Design_pKa_help
        print 'pKa_info, PDB file',self.pdbfile
        self.pKa_info=Design_pKa_help.pKa_info(pdbfile=self.pdbfile,parent=self,PBEsolver=self.params['PBEsolver'],save_file=False)
        #
        # Set the calculation routine to CPP
        #
        import pKa_MC
        self.pKaCALC=pKa_MC.pKa_calculation_class(pdbfile=self.pdbfile,pKa_info=self.pKa_info,params=self.params,parent=self)
        self.pKaCALC.set_MC_CPP()
        #
        # Re-calculate the wild type pKa values
        #
        print 'Calculating wild type pKa values'
        self.wt_pKas,self.wt_prot_states=self.pKaCALC.calc_wt_pkas()
        #
        # Set the reporter groups
        #
        if reporter_groups:
            self.reporter_groups=reporter_groups
        else:
            self.reporter_groups=self.wt_pKas.keys()
        return
    
    def calc_dpKa(self,mutations,mutant_PDB_files={}):
        """Calculate the delta pKa for this set of mutations
        mutations is a list of the form [':0052:ASP:ALA+:0087:ASP:ALA']
        mutant_PDB_files is a dictionary with the path of the mutant PDB files.
        If a mutant PDB file is not found then the routines models it:
        {':0052:ASP:ALA+:0087:ASP:ALA':<path to PDB file>}"""
        #
        # save the dictionary with mutant PDB files
        #
        self.mutfile_names=mutant_PDB_files
        #
        # Remove the mutated groups from reporter groups
        #
        reporter_groups=[]
        for group in self.reporter_groups:
            skip=False
            for mutation in mutations:
                if group in mutation:
                    skip=True
                    break
            if not skip:
                reporter_groups.append(group)
        #
        # Set the mutant PDB file array
        #
        self.pKaCALC.mutfile_names=mutant_PDB_files
        #
        # Call the calculation routine
        #
        print 'Calling calc_MC_mut_pKas'
        mut_pKas,prot_states=self.pKaCALC.calc_MC_mut_pKas(mutations=mutations,
                                                           data={},
                                                           atomdata={},
                                                           wt_pKas=self.wt_pKas,
                                                           require_pKa=reporter_groups)
        if not mut_pKas and not self.params['use_titration_curves']:
            #
            # Sometimes we get a None back
            #
            return {}
        #
        wt_pKa_names={}
        for group in reporter_groups:
            if self.wt_pKas.has_key(group):
                #wt_pKas.append(self.wt_pKas[group])
                wt_pKa_names[group]=group
            #else:
                #wt_pKas.append(None)
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
                diff=0.0
                for pH in self.wt_prot_states[wt_group].keys():
                    if pH>=self.params['pHstart_integrate'] and pH<=self.params['pHstop_integrate']:
                        diff=diff+(prot_states[group][pH]-self.wt_prot_states[wt_group][pH])*self.params['pHstep']
                dpKas[group]=diff
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
        return dpKas

    def close_stdout(self):
        return
        
    def open_stdout(self):
        return
                    
if __name__=='__main__':
    import sys
    X=dpKa(sys.argv[1])
    mutant_PDB_files={':0087:ASP:ALA':'/home/people/nielsen/lib/development/PEAT/trunk/pKa/D87A.pdb',
                      ':0087:ASP:HIS':'/home/people/nielsen/lib/development/PEAT/trunk/pKa/D87H.pdb'}

    for mutation in [[':0087:ASP:HIS']]:
        dpKas=X.calc_dpKa(mutation,mutant_PDB_files=mutant_PDB_files)
        print
        print '---------------------'
        print 'Mutation: %s' %mutation
        print 'dpKas'
        groups=dpKas.keys()
        groups.sort()
        for group in groups:
            print group,dpKas[group]
    
