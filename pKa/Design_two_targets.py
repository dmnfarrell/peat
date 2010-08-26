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

def set_parameters(pdbfile=None,design_statement=None,min_dist=None,num_muts=None):
    #
    # Set the parameters that are the same for all permutations
    #
    defaults=Design_pKa.get_defaults()
    #
    # PDB file
    #
    defaults['pdb'][0]=pdbfile
    #
    # pKa calculation parameters
    #
    defaults['pHstart'][0]=0.1
    defaults['pHstop'][0]=12.0
    defaults['pHstep'][0]=0.05
    defaults['pKMCsteps'][0]=200000
    #
    # Design settings
    #
    # Target
    #
    defaults['pKas'][0]=design_statement
    #
    # Method
    #
    defaults['MC'][0]=1
    defaults['tabulated'][0]=1
    defaults['MCsteps'][0]=0
    defaults['TR'][0]=None
    defaults['MC'][0]=1
    #
    # Be not-so-noisy
    #
    defaults['silent'][0]=2
    #
    # Minimum distance between target and mutation
    #
    defaults['min_target_dist'][0]=min_dist
    #
    # Max number of mutations
    #
    defaults['max_mutations'][0]=int(num_muts)
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
    print 'Construct dpKa matrix for a single design criterium'
    print
    print 'Usage: Design_two_targets <pdbfile> <database file> <design statement>'
    print
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
    method='MC'
    #
    # Which design are we doing?
    #
    design_statement=sys.argv[3]
    #
    # Set the file name
    #
    filename=os.path.join(dir,'designtwo__'+os.path.split(pdbfile)[1])+'_'+method+design_statement
    #
    # Print the filename of the dictionary
    #
    print 'Setting dictionary filename to',filename
    DB=dictIO(filename)
    print 'I am running in %s' %os.getcwd()
    #
    # Get wild type pKa values
    #
    import pKaTool.pKaIO
    X=pKaTool.pKaIO.pKaIO(pdbfile)
    if not X.assess_status():
        import os
        print 'You have to run a pKa calculation first'
        raise Exception()
    #
    # OK, got the pKa calc. Read results
    #
    wt_pKas=X.readpka()
    #
    # See if we have an old database file we should work on
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
    import Design_pKa
    Design_pKa.delete_old_files(pdbfile,'N','N')
    #
    # Start looping
    #
    dist_range=range(1,26)
    num_muts_range=range(1,21)
    count=0
    #
    # Set the number of MC steps
    #
    pKMCsteps=200000
    #
    # Start of loop
    #

    #
    # Loop over number of mutations and distance cut-off
    #
    X=None
    if not results.has_key(design_statement):
        results[design_statement]={}
    for min_target_dist in dist_range:
        for num_muts in num_muts_range:
            #
            # Make sure the dictionary entries are there
            #
            if not results[design_statement].has_key(num_muts):
                results[design_statement][num_muts]={}
            if not results[design_statement][num_muts].has_key(min_target_dist):
                results[design_statement][num_muts][min_target_dist]={}
            #
            # Have we done this one yet?
            #
            x='Checking dist: %5.2f #muts: %2d....' %(float(min_target_dist),num_muts)
            print x,
            if results[design_statement][num_muts][min_target_dist]=={}:
                print 'not done. Designing solutions:'
                #
                # Get the parameters
                #
                defaults=set_parameters(pdbfile=pdbfile,
                                        design_statement=design_statement,
                                        min_dist=min_target_dist,
                                        num_muts=num_muts)
                params={}
                for key in defaults.keys():
                    params[key]=defaults[key][0]
                #
                # Call the design routine
                #
                if not X:
                    X=Design_pKa.Design_pKa(params)
                solutions,dpKa_dict=Design_pKa.run_opt(defaults,X)
                res=dpKa_dict.keys()
                res.sort()
                #for r in res:
                #    print r,dpKa_dict[r]
                #
                # Store the solutions
                #
                print 'Found these solutions'
                print dpKa_dict.keys()
                if dpKa_dict.keys()==[]:
                    results[design_statement][num_muts][min_target_dist]={method:{'None':'No solutions'}}
                else:
                    results[design_statement][num_muts][min_target_dist]={method:dpKa_dict.copy()}
                results=DB.update(results)
            else:
                print 'done.'
    #
    # All done
    #
    print
    print 'All done - normal exit'
    print
    return


if __name__=="__main__":
    main()
