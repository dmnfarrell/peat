#!/usr/bin/env python
#
# pKarun - scripts for running pKa calculations
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

import string, os, pKa_general
from pKarun_base import *

#=======================================================

def makenewargv(argv,X,infile):
    #
    # Construct the new Invocation file
    #
    newarg=[]
    subset_found=None
    skipnext=None
    for arg in argv:
        if skipnext:
            skipnext=None
        elif arg=='-groups':
            skipnext=1
            subset_found=1
        elif arg=='-subset':
            subset_found=1
            pass
        else:
            newarg.append(arg)
    if subset_found or X.params['groups']:
        newarg.append('-subset')
        if X.params['groups']:
            newarg.append('-groups')
            newarg.append(X.params['groups'])
    argv=newarg
    print 'New Invocation file written.'
    fd=open(infile,'w')
    # print argv
    fd.write(string.join(argv))
    fd.close()
    #print 'Returning argv:', argv
    return argv

#
# ------------------------------------------------
#

def do_main(action=None,new_cutoff=None,selection_rounds=2):
    from sys import argv
    infile='Invocation'
    if os.path.isfile(infile) and (len(argv)==1 or action):
        fd=open(infile)
        line=fd.readline()
        fd.close()
        argv=string.split(line)
        print 'Using previous command:',string.join(argv)
    else:
        print 'Invoked with command:', string.join(argv)
        fd=open(infile,'w')
        fd.write(string.join(argv))
        fd.close()
    
    if len(argv)<=1:
	print 'Usage:'
	print 'pKarun.py <pdbfile> [-noopt] [-opth <x>] [-super] [-dessoup] [-pair] [-indi <xx>] [-pairene <xx>] '
        print '---------------Hydrogen bond optimisation-------------'
        print '-opth 0: No hydrogen bond optimisation.'
        print '-opth 1: Hydrogen bond optimisation level 1 (only initial optimisation of the H-bond network).'
        print '-opth 2: Hydrogen bond optimisation level 2 (dynamic H-bond optimisation for every protonation state).)'
	print '-dessoup: Force WHAT IF to skip optimisation of the H-bond network when calculating desolvation energies.'
	print '-pair: Calculate also neutral-charged, charged-neutral, and neutral-neutral interaction energies for strongly interacting pairs of titratable groups. (default)'
        print '-pairene <xx>: The energy cut-off (in kT) for strongly interacting pairs of titratable groups (default: 1kT)'
        print '-bumps: Include a simple bump score in the Hbond energy. (disabled by default).'
        print '-----------------Focussing--------------------------'
        print '-yang: Use Yang et al. focussing mode,. (disabled by default)'
        print '-----------------Poisson-Boltzmann solver-----------------------'
        print '-indi <xx>: Protein dielectric constant (default: 8)'
        print '-exdi <xx>: Solvent dielectric constant (default: 80)'
        print '-dbcrit <xx>: Use the alternative dielectric constant for all energy calculations if the B-factor of the residue is above xx (default 30).'
        print '-dbdiel <xx>: Alternative dielectric constant for desolvation and backgr calcs. (default: 15)'
        print '-mcrit <xx>: Use the alternative dielectric constant for the charge-charge energies if the B-factor of the residue is above xx (default 1000).'
        print '-mdiel <xx>: Alternative dielectric constant for the charge-charge energy calculations. (default: 60)'
        print '------------Other options-----------------------------'
        print '-subset : Calculate only for groups a subset of the titratable groups'
        print '-groups <group 1> <group 2> ... : Specify groups to calculate for '
        print '-group_cutoff: The energy requirement for finding groups that interact strongly with the selected subset of groups'
        print '-listgroups: Lists the groups that can be specified with "-groups" '
        #
        # Done - raise and error
        #
	error='Too few arguments'
        print argv
	raise Exception('Too few arguments')
    pdbfile=argv[1]
    x=pKamisc()
    params=x.parse_parameterlist(argv)
    X=pKarun(os.getcwd(),pdbfile,params)
    #
    # Get the titratable groups
    #
    Y=pKa_general.pKatools()
    list=Y.get_titratable_residues(pdbfile)
    okgroups={}
    for x,y in list:
        okgroups[y]=1
    actgroups={}
    #
    # Find out what the sucker wants to do
    #
    while 1:
        if not action:
            print
            print 'Select one of the following:'
            print '1. Run everything that is missing in parallel.'
            print '2. Run all that is missing sequentially'
            print '3. Clear all data files in this directory.'
            print '4. Calculate only desolvation energies'
            print '5. Calculate only background interaction energies'
            print '6. Calculate only the charge-charge interaction matrix'
            print '7. Calculate pKa values from desolvation, background and charge-charge energies'
            print '8. List titratable groups in the PDB file'
            print '9. Select titratable groups'
            print '10. Print the groups currently selected'
            print '11. Find groups that interact strongly with the selected groups'
            print '12. Change the cut-off value for selecting strongly interacting groups'
            print '13. Change the cut-off value for detailed calculation of pairs of titratable groups'
            print '14. List setup'
            print '15. Run the calculation for estimating the autonomy of clusters of titratable groups'
            print '16. Calculate the effect of individual groups on subset pKa values [in network]'
            print '17. Calculate the effect of individual groups on subset pKa values [isolated]'
            print '18. Exit'
            number=-1
            if params['auto']!=0:
                number=int(params['auto'])
                print 'Auto-running option %d' %number
            #
            # If no auto-run, then ask the user for input
            #
            while number<1 or number>18:
                try:
                    number=string.atoi(raw_input('-----> '))
                except ValueError:
                    print 'You have to give a number'
                    number=-1
            print
            print
        else:
            #
            # -------------------------------------------------------------
            #
            # This is for doing the autonomy calculations in an easy way
            #
            number=action
            #
            # Change the cut-off value for selecting strongly interacting groups
            #
            change_next=None
            for x in range(len(argv)):
                if argv[x]=='-group_cutoff':
                    change_next=1
                elif change_next:
                    #print 'Present value: %5.2f' %(string.atof(argv[x]))
                    val=new_cutoff
                    argv[x]=val
                    break
            if not change_next:
                x=new_cutoff
                argv.append('-group_cutoff')
                argv.append(str(x))
            argv=makenewargv(argv,X,infile)
            x=pKamisc()
            params=x.parse_parameterlist(argv)
            params['selection_rounds']=selection_rounds
            X=pKarun(os.getcwd(),pdbfile,params)
        #
        # And now we do whatever we were told to do
        #
	if number==1:
	    X.Runall()
	elif number==2:
	    X.runseq()
	elif number==3:
	    X.clean()
	elif number==4:
	    X.desolv(1)
	elif number==5:
	    X.backgr(1)
	elif number==6:
	    X.matrix(1)
	elif number==7:
	    X.solvepka()
        elif number==8:
            print '%s contains the following titratable groups:' %pdbfile
            print 'Residue\tGroup Number'
            for residue,groupnum in list:
                print '%10s\t%4d' %(residue,groupnum)
            print
            print 'Use "select titratable groups" to select a subset.'
        elif number==9:
            done=None
            while not done:
                group=string.atoi(raw_input('Enter group number [end selection process with 0]: '))
                if group==0:
                    break
                if actgroups.has_key(group):
                    print 'You already selected this group'
                else:
                    if okgroups.has_key(group):
                        actgroups[group]=1
                    else:
                        print 'Invalid selection [group does not exist]'
            X.params['subset']=1
            X.params['groups']=string.replace(str(actgroups.keys())[1:-1],' ','')
            #
            # Construct the new Invocation file
            #
            argv=makenewargv(argv,X,infile)
        elif number==10:
            #
            # List the selected groups
            #
            if not X.params['groups']:
                print 'You have not selected any groups'
            else:
                print 'You have selected the following groups:'
                print '%10s\t%s' %('Residue','Group Number')
                groups_tmp=string.split(X.params['groups'],',')
                groups=[]
                for group in groups_tmp:
                    groups.append(string.atoi(group))
                for group in groups:
                    found=None
                    for residue,groupnum in list:
                        if group==groupnum:
                            print '%10s\t%4d' %(residue,groupnum)
                            found=1
                            break
                    if not found:
                        raise "Invalid group: ",group
            print

        elif number==11:
            #
            # Find strongly interacting groups
            #
            group_list=X.runactsite(list)
            sortlist=[]
            for group in group_list:
                included=None
                split=string.split(group)
                for residue,groupnum in list:
                    if residue==group:
                        sortlist.append(string.atoi(str(groupnum)))
                        included=1
                if not included:
                    print 'Error in including: %s' %(str(group))
            #
            # Sort them
            #
            sortlist.sort()
            X.params['groups']=''   
            for group in sortlist:
                X.params['groups']=X.params['groups']+str(group)+','
            #
            # Insert commas
            #
            X.params['groups']=string.replace(X.params['groups'][:-1],' ','') 
            if X.params['groups']:
                argv=makenewargv(argv,X,infile)
            else:
                print 'Did not find any groups. Something is wrong'
                raise 'error'
        elif number==12:
            #
            # Change the cut-off value for selecting strongly interacting groups
            #
            change_next=None
            for x in range(len(argv)):
                if argv[x]=='-group_cutoff':
                    change_next=1
                elif change_next:
                    print 'Present value: %5.2f' %(string.atof(argv[x]))
                    val=raw_input('Enter new value: ')
                    argv[x]=val
                    break
            if not change_next:
                x=raw_input('Enter value: ')
                argv.append('-group_cutoff')
                argv.append(str(x))
            argv=makenewargv(argv,X,infile)
            x=pKamisc()
            params=x.parse_parameterlist(argv)
            X=pKarun(os.getcwd(),pdbfile,params)
        elif number==13:
            #
            # Change the cut-off value for strongly interacting pairs
            #
            change_next=None
            for x in range(len(argv)):
                if argv[x]=='-pairene':
                    change_next=1
                elif change_next:
                    print 'Present value: %5.2f' %(string.atof(argv[x]))
                    val=raw_input('Enter new value: ')
                    argv[x]=val
                    break
            if not change_next:
                x=raw_input('Enter value: ')
                argv.append('-pairene')
                argv.append(str(x))
            #print argv
            argv=makenewargv(argv,X,infile)
            x=pKamisc()
            params=x.parse_parameterlist(argv)
            X=pKarun(os.getcwd(),pdbfile,params)
        elif number==14:
            #
            # List the setup
            #
            x=pKamisc()
            params=x.parse_parameterlist(argv)
            X=pKarun(os.getcwd(),pdbfile,params)
        elif number==15:
            #
            # Do the autonomy calculation
            #
            X.autonomy()
            if action:
                return
        elif number==16:
            #
            # Autonomy #2 GRPINF
            #
            X.grpinf()
            if action:
                return
        elif number==17:
            #
            # Autonomy #2 GRPINF
            #
            X.ex_system(pdbfile)
        elif number==18:
            #
            # Exit
            #
            break

        ## exit loop if automated choice (e.g. when executing from cluster)
        if params['auto'] != 0:
            params['auto'] = 18
    
if __name__ == "__main__":
    do_main()


