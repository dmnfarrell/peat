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

            
from Design_accuracy import dictIO
import Design_pKa

#
# ----------------------------
#

def get_solutions_dist_nummuts(pdbfile=None,
                               target_residue=None,
                               target_dpKa=2.0,
                               num_muts=6,
                               min_target_dist=5.0,
                               dpKas={},
                               method='MC',
                               pKMCsteps=0,
                               X=None):
    #
    # Check if we already did this..
    #
    #
    # Do we have solutions for this problem?
    #
    added_data=None
    missing=None
    if dpKas=={}:
        missing=1
    #
    # Do it?
    #
    if missing:
        defaults=local_defaults(pdbfile,target_residue,target_dpKa,float(min_target_dist),pKMCsteps)
        defaults['max_mutations'][0]=int(num_muts)
        print 'Target: %s, dpKa: %s, Number of mutations: %3d, min_dist: %5.3f' %(
            target_residue,target_dpKa,defaults['max_mutations'][0],defaults['min_target_dist'][0])
        #
        # Set the method default
        #
        allmethods=['TR','MC']
        for m in allmethods:
            defaults[m][0]=None
        if not method=='phidiff':
            defaults[method][0]=1
        #
        # Set the parameters
        #
        if not X:
            params={}
            for key in defaults.keys():
                params[key]=defaults[key][0]
            X=Design_pKa.Design_pKa(params)
        #
        # Get the solutions and the corresponding dpKas
        #
        solutions,dpKa_dict=Design_pKa.run_opt(defaults,X)
        solution_keys=solutions.keys()
        #

        #
        # If there are no solutions then store that info
        #
        if len(solution_keys)==0:
            added_data=1
            dpKas[1]={'no solutions':1}
        #
        # Loop over all the solutions and get the dpKas        
        #
        for solution in solution_keys:
            #
            # Make sure the dictionary is ready
            #
            if not dpKas.has_key(solution):
                added_data=1
                dpKas[solution]={}
            #
            # Check that we have a real mutation (i.e. a non- None) mutation
            #
            dpKas[solution]['mutations']=solutions[solution]['mutations']
            realmut=None
            for mutation in solutions[solution]['mutations']:
                if mutation:
                    realmut=1
            if not realmut:
                continue
            #
            # Get the delta pKa value
            #
            key=str(dpKas[solution]['mutations'])
            if dpKa_dict.has_key(key):
                dpKas[solution][method]=dpKa_dict[key]
    #
    # Did we calculate dpKa values?
    #
    for solution in dpKas.keys():
        if not dpKas[solution].has_key('mutations'):
            #
            # If it's not a real solution then don't check it
            #
            continue
        #
        # Now calculate the dpKas with the method we want, it it's needed
        #
        if method=='phidiff':
            methods=[[None,None]]
        else:
            methods=[[method,None]]
        #
        # Set the method name
        #
        for method,tabulated in methods:
            method_name=method
            if not method:
                method_name='phidiff'
            #
            # Is this tabulated mode?
            #
            if tabulated:
                method_name=method_name+'_tab'
            #
            # Did we calculate this dpKa?
            #

            if not dpKas[solution].has_key(method_name) and not dpKas[solution].has_key(method_name+'_tab'):
                #
                # Get defaults
                #
                import Design_accuracy
                defaults=Design_accuracy.get_this_defaults(pdbfile,target_residue,target_dpKa)
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
                realmut=None
                filtered_muts=[]
                for mut in mutations:
                    if mut:
                        filtered_muts.append(mut)
                        realmut=1
                if not realmut:
                    continue
                #
                # Set the mutations variable in defaults
                #
                defaults['mutations'][0]=string.join(filtered_muts,',')
                # Calculate delta pkas
                defaults['calc_dpka'][0]=1
                #
                # Get the dpKa(s)
                #
                print 'Still calculating here'
                tmp=Design_pKa.run_opt(defaults,X)
                #
                # We should have only one key
                #
                keys=tmp.keys()
                if len(keys)!=1:
                    print keys
                    raise 'Too many keys when calculating dpkas for one solution'
                added_data=1
                dpKas[solution][method_name]=tmp[keys[0]]
    #
    # Done
    #
    return dpKas,added_data,X

#
# -------------------------
#
def local_defaults(pdbfile,target_residue,target_dpKa,min_dist,pKMCsteps=None):
    #
    # Set the parameters that are the same for all permutations
    #
    defaults=Design_pKa.get_defaults()
    # PDB file
    defaults['pdb'][0]=pdbfile
    #
    # pKa calculation parameters
    #
    defaults['pHstart'][0]=0.1
    defaults['pHstop'][0]=12.0
    defaults['pHstep'][0]=0.05
    defaults['pKMCsteps'][0]=200000
    if pKMCsteps:
        defaults['pKMCsteps'][0]=pKMCsteps
    #
    # Design settings
    #
    # Target
    #
    target=target_residue+'='+target_dpKa
    defaults['pKas'][0]=target
    #
    # Method
    #
    defaults['MC'][0]=1
    defaults['tabulated'][0]=1
    defaults['MCsteps'][0]=0
    #
    # Be not-so-noisy
    #
    defaults['silent'][0]=2
    #
    # Minimum distance between target and mutation
    #
    defaults['min_target_dist'][0]=min_dist
    #
    # Do not save the solutions
    #
    defaults['save_solutions'][0]=None
    return defaults

#
# ----
#


def main(pdbfile=None):
    #
    # Calculate the matrix of pKa shifts [number of mutations:min distance from active site]
    #
    # Get the PDB file
    #
    print
    print 'Construct dpKa matrix'
    print
    print 'Usage: Design_dist_nummuts <pdbfile> <database file> <method>'
    print
    #
    # Are we running in parallel?
    #
    import os

    mpi_rank=None
    try:
        import mpi
        mpi_rank=mpi.rank
        mpi_size=mpi.size
        print 'MPI rank',mpi_rank
        print 'MPI size',mpi_size
        if size==1:
            mpi_rank=None
    except:
        #
        # No MPI
        #
        pass
    #
    # start the run
    #
    import sys,os
    if len(sys.argv)<4:
        raise 'Incorrect usage'
    if not pdbfile:
        pdbfile=sys.argv[1]
    if len(sys.argv)>2:
        dir=sys.argv[2]
    else:
        raise 'You have to provide a dir for the output'
    os.chdir(dir)
    
    #
    # Get the method
    #
    method=sys.argv[3]
    #
    # Are we doing this cheap style?
    #
    mcstep_factor=1.00
    if sys.argv[-1]=='cheap':
        mcstep_factor=0.01
    #
    # Set the file name
    #
    filename=os.path.join(dir,'distance_nummuts__'+os.path.split(pdbfile)[1])+'_'+method
    #
    # Keep the cheap results separated from the real results
    #
    if sys.argv[-1]=='cheap':
        filename=filename+'cheap'
    #
    # Print the filename of the dictionary
    #
    print 'Setting dictionary filename to',filename
    DB=dictIO(filename)
    print 'I am running in %s' %os.getcwd()
    #
    # Do we have a completed pKa calc for it?
    #
    import pKaTool.pKaIO
    X=pKaTool.pKaIO.pKaIO(pdbfile)
    if not X.assess_status():
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
    import os
    if os.path.isfile(filename):
        #
        # Load old file
        #
        results=DB.load()
    else:
        #
        # No, no old restuls
        #
        results={}
        #
        # Add the PDB file
        #
        fd=open(pdbfile)
        lines=fd.readlines()
        fd.close()
        results['pdbfile']=lines
        #
        # Add the full wt pKa values
        #
        results['wt_full']=wt_pKas
        DB.save(results)
    # ---------------------------------------
    #
    # Delete old files
    #
    Design_pKa.delete_old_files(pdbfile,'N','Y')
    #
    # Start looping
    #
    dist_range=range(1,26)
    num_muts_range=range(1,21)
    count=0
    design_low=2.0
    design_high=10.0
    #
    # Count the number of designable groups
    #
    design_groups=[]
    for group in groups:
        pKaval=wt_pKas[group]['pKa']
        if pKaval>design_low and pKaval<design_high:
            design_groups.append(group)
    tot_count=float(len(design_groups)*2*len(dist_range)*len(num_muts_range))/100.0
    #
    # Set the number of MC steps
    #
    #
    #if len(groups)<60:
    pKMCsteps=200000
    #else:
    #    #print 'We have %d groups, so increasing pKMCsteps' %(len(groups))
    #    pKMCsteps=200000
    pKMCsteps=int(mcstep_factor*float(pKMCsteps))
    #
    # Start of loop
    #
    groups.sort()
    #
    # Are we running MPI?
    #
    if mpi_rank!=None:
        fraction=1.0/float(mpi_size)
        print 'I will be doing %.2f %% of the job' %(100.0*fraction)
        num_groups=len(groups)
        print 'Specifically I will be doing groups',groups[int(fraction*mpi_rank*num_groups):int(fraction*(mpi_rank+1)*num_groups)]
        print fraction*mpi_rank*num_groups,'to',fraction*(mpi_rank+1)*num_groups
        groups=groups[int(fraction*mpi_rank*num_groups):int(fraction*(mpi_rank+1)*num_groups)]
        #
        # Sleep to desynchronize processes
        #
        import time
        time.sleep(mpi_rank)
        print mpi_rank,'done sleeping!'
        import sys
        sys.stdout.flush()
    #
    # Do the calculation
    #
    for group in groups:
        #
        # Loop over all groups
        #
        if not results.has_key(group):
            results[group]={}
        #
        # Is somebody working on this group?
        #
        #results=DB.update(results)
        if results[group].has_key('locked'):
            if results[group]['locked']==1:
                print '%s is locked' %group
                continue
        #
        # Lock group
        #
        results[group]['locked']=1

        #
        # Force reinitialisation of instance
        #
        X=None
        #
        # Loop over number of mutations and distance cut-off
        #
        pKaval=wt_pKas[group]['pKa']
        if pKaval>design_low and pKaval<design_high:
            for design in ['+20.0','m20.0']:
                if not results[group].has_key(design):
                    results[group][design]={}
                for min_target_dist in dist_range:
                    for num_muts in num_muts_range:
                        if not results[group][design].has_key(num_muts):
                            results[group][design][num_muts]={}
                        #
                        # Print what we are doing
                        #
                        print 'Checking: %20s design: %s, dist: %5.1f, nummuts: %3d %%done %5.2f' %(group,design,
                                                                                  min_target_dist,
                                                                                  num_muts,float(count)/tot_count)
                        count=count+1
                        #
                        # ..
                        #
                        if not results[group][design][num_muts].has_key(min_target_dist):
                            results[group][design][num_muts][min_target_dist]={}
                            #
                            # Get the solutions
                            #
                            dpKas,data_added,X=get_solutions_dist_nummuts(pdbfile,
                                                                          group,
                                                                          design,num_muts,
                                                                          min_target_dist,
                                                                          results[group][design][num_muts][min_target_dist],
                                                                          method,pKMCsteps,X)
                            results[group][design][num_muts][min_target_dist]=dpKas.copy()
                            pass
                    results=DB.update(results)
        else:
            results[group]['pKa out of range']=1
        #
        # Unlock group and merge results
        #
        del results[group]['locked']
    results=DB.update(results)
    #
    # All done
    #
    print
    print 'All done - normal exit'
    print
    return





if __name__=='__main__':
    #try:
    import profile
    #profile.run('main()','profile_out')
    main()
    #except:
    #    print
    #    print 'Usage: Design_dist_nummuts.py <pdbfile> [<result dir>]'
    #    print
