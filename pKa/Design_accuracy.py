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

class dictIO:

    #
    # Class for managing the results dictionary
    # This class allows us to run several copies of Design_accuracy at once.
    #

    def __init__(self,filename):
        self.filename=filename
        self.lockfile=filename+'.lock'
        import lock
        self.LOCK=lock.lockfile(self.lockfile)
        return

    def load(self):
        #
        # Lock
        #
        count=0
        while self.LOCK.lock()==-1:
            count=count+1
            if count==1000:
                print 'waiting for lock when loading'
        dict=self._load()
        #
        # Unlock
        #
        self.LOCK.unlock()
        return dict

    def _load(self):
        #
        # Load the dict
        #
        fd=open(self.filename)
        import cPickle
        dict=cPickle.load(fd)
        fd.close()
        #
        # Return data
        #
        return dict

    def save(self,dict):
        #
        # Lock
        #
        while self.LOCK.lock()==-1:
            pass
        #
        # Save
        #
        self._save(dict)
        #
        # Unlock
        #
        self.LOCK.unlock()
        return
        
    def _save(self,dict):      
        #
        # Save the dict
        #
        fd=open(self.filename,'w')
        import cPickle
        cPickle.dump(dict,fd)
        fd.close()
        #
        # Return 
        #
        return 

    def update(self,newdict):
        while self.LOCK.lock()==-1:
            pass
        #
        # Locked
        #
        saveddict=self._load()
        #
        # Merge the two dictionaries
        #
        saveddict.update(newdict)
        #
        # Save the new dict
        #
        self._save(saveddict)
        #
        # Unlock
        #
        self.LOCK.unlock()
        return saveddict

#
# ----------------------------
#

def get_this_defaults(pdbfile,target_residue,target_dpKa):
    #
    # Set the parameters that are the same for all permutations
    #
    import Design_pKa

    defaults=Design_pKa.get_defaults()
    # PDB file
    defaults['pdb'][0]=pdbfile
    #
    # pKa calculation parameters
    #
    defaults['pHstart'][0]=0.0
    defaults['pHstop'][0]=14.0
    defaults['pHstep'][0]=0.05
    defaults['MCsteps'][0]=1
    defaults['pKMCsteps'][0]=50000
    #
    # Design settings
    #
    # Target
    #
    target=target_residue+'='+target_dpKa
    defaults['pKas'][0]=target
    #
    # Be not-so-noisy
    #

    defaults['silent'][0]=1

    #
    # Minimum distance between target and mutation
    #
    defaults['min_target_dist'][0]=5.0
    defaults['max_mutations'][0]=20.0
    return defaults

#
# ---------------------------
#

def get_solutions(pdbfile,target_residue,target_dpKa,dpKas):
    #
    # Main routine.
    # Get the solutions with MC_tab and then rescore them with a subset of
    # the different methods
    #
    # Get solutions
    #
    import Design_pKa
    defaults=Design_pKa.get_defaults()
    defaults=get_this_defaults(pdbfile,target_residue,target_dpKa)
    # Method
    defaults['MC'][0]=1
    defaults['tabulated'][0]=1
    #
    # Check if we already did this..
    #
    missing=None
    for x in range(1,11):
        if not dpKas.has_key(x):
            missing=1
            print 'Design solutions not calculated'
            break
    #
    # Do it?
    #
    if missing:
        solutions=Design_pKa.run_opt(defaults)
        solution_keys=solutions.keys()
    else:
        solution_keys=dpKas.keys()
    solution_keys.sort()
    #
    # Calculate pKa values using the other methods
    #
    # Define all the other methods that we will use
    #
    # First element is method, second is value of 'tabulated'
    #
    methods=[['TR',None],[None,None],['MC',None],['TR',1],['MC',1]]
    for solution in solution_keys[:10]:
        #
        # Make sure the dictionary is ready
        #
        if not dpKas.has_key(solution):
            dpKas[solution]={}
            dpKas[solution]['mutations']=solutions[solution]['mutations']
            realmut=None
            for mutation in solutions[solution]['mutations']:
                if mutation:
                    realmut=1
            if not realmut:
                continue
        #
        # Now calc dpKas with different methods
        #
        for method,tabulated in methods:
            #
            # Print info
            #
            print
            print 'Calculating dpKa for %s with %s. Tabulated: %s' %(str(dpKas[solution]['mutations']),method,tabulated)
            #
            # Get defaults
            #
            defaults=get_this_defaults(pdbfile,target_residue,target_dpKa)
            if method:
                # Choose method
                defaults[method][0]=1
            #
            # Set tabulated
            #
            defaults['tabulated'][0]=tabulated
            #
            # Set the general parameters
            #
            mutations=dpKas[solution]['mutations']
            import string
            #
            # Remove all None elements from mutations
            #
            filtered_muts=[]
            for mut in mutations:
                if mut:
                    filtered_muts.append(mut)
            #
            # Set the mutations variable in defaults
            #
            defaults['mutations'][0]=string.join(filtered_muts,',')
            # Calculate delta pkas
            defaults['calc_dpka'][0]=1
            #
            # Set a more meaningful method name
            #
            method_name=method
            if not method:
                method_name='phi/ln10'
            #
            # Is this tabulated mode?
            #
            if tabulated:
                method_name=method_name+'_tab'
            #
            # Get the dpKa(s)
            #
            if not dpKas[solution].has_key(method_name):
                tmp=Design_pKa.run_opt(defaults)
                #
                # We should have only one key
                #
                keys=tmp.keys()
                if len(keys)!=1:
                    print keys
                    raise 'Too many keys when calculating dpkas for one solution'
                dpKas[solution][method_name]=tmp[keys[0]]
            
    #
    # Done
    #
    return dpKas

#
# -------------------------
#


def main():
    #
    # Get the PDB file
    #
    import sys,os
    pdbfile=sys.argv[1]
    suffix=sys.argv[2]
    if os.environ['USER']=='nielsen':
        filename=os.path.join('/enzyme/nielsen/work/pKa_design/accuracy','accuracy_'+os.path.split(pdbfile)[1])
    elif os.environ['USER']=='btconnolly':
        filename=os.path.join('/enzyme/btconnolly/pKa_design',suffix+'accuracy_'+os.path.split(pdbfile)[1])
    #filename=os.path.join(os.getcwd(),'accuracy_'+os.path.split(pdbfile)[1])
    print 'Setting dictionary filename to',filename
    #
    # Make sure that we delete all info from previous runs
    #
    wtpKafile=pdbfile+'_wt_pKavals'
    if os.path.isfile(wtpKafile):
        os.unlink(wtpKafile)
    #
    # Do we have a completed pKa calc for it?
    #
    import pKa
    X=pKa.pKaIO(pdbfile)
    X.usenewnames()
    if not X.assess_status():
        import os
        print 'You have to run a pKa calculation first'
        raise Exception()
    #
    # OK, got the pKa calc. Read results
    #
    wt_pKas=X.readpka()
    groups=wt_pKas.keys()
    groups.sort()
    #
    # Design target: for all groups with pKa value between 2 and 10: +2 and -2
    #
    DB=dictIO(filename)
    import os
    if os.path.isfile(filename):
        #
        # Load old file
        #
        accuracy=DB.load()
    else:
        #
        # No, no old restuls
        #
        accuracy={}
        #
        # Add the PDB file
        #
        fd=open(pdbfile)
        lines=fd.readlines()
        fd.close()
        accuracy['pdbfile']=lines
        #
        # Add the full wt pKa values
        #
        accuracy['wt_full']=wt_pKas
        DB.save(accuracy)
    # ---------------------------------------
    #
    # Delete old files
    #
    import Design_pKa
    Design_pKa.delete_old_files(pdbfile,'N','N')
    #
    # Start looping
    #
    groups.sort()
    for group in groups:
        #if not group==':0129:LEU:CTERM':
        #    continue
        if not accuracy.has_key(group):
            accuracy[group]={}
        #
        # Is somebody working on this group?
        #
        accuracy=DB.update(accuracy)
        if accuracy[group].has_key('locked'):
            if accuracy[group]['locked']==1:
                continue
        #
        # Lock group
        #
        accuracy[group]['locked']=1
        #
        # Save the dictionary
        #
        accuracy=DB.update(accuracy)
        #
        # OK, now we can calculate in peace
        #
        pKaval=wt_pKas[group]['pKa']
        if pKaval>2.0 and pKaval<10.0:
            #
            # First design pKa +2.0
            #
            design='+2.0'
            if not accuracy[group].has_key(design):
                accuracy[group][design]={}
            print 'Designing and evalulating: %s' %(group+design)
            dpKas=get_solutions(pdbfile,group,design,accuracy[group][design])
            accuracy[group][design]=dpKas.copy()
            #
            # Then do pKa -2.0
            #
            design='m2.0'
            if not accuracy[group].has_key(design):
                accuracy[group][design]={}
            print 'Designing and evalulating: %s' %(group+design)
            dpKas=get_solutions(pdbfile,group,design,accuracy[group][design])
            accuracy[group][design]=dpKas.copy()
        else:
            accuracy[group]['pKa out of range']=1
        #
        # Unlock group and merge results
        #
        accuracy[group]['locked']=None
        accuracy=DB.update(accuracy)
    #
    # All done
    #
    return





if __name__=='__main__':
    main()
