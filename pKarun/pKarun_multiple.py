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

import pKarun_base
import sys
sys.path.append('/home/people/nielsen/lib/development/python')

def corall(pdbfile):
    """Do Corall on the PDBfile"""
    print 'Performing CORALL on %s' %pdbfile
    import pKarun.WI_tools
    log,files=pKarun.WI_tools.corall(pdbfile)
    #for line in log:
    #    print line,
    filename=files.keys()[0]
    fd=open(pdbfile,'w')
    for line in files[filename]:
        fd.write(line)
    fd.close()
    print 'CORALL done'
    print

#
# ----
#

def main(options,args):
    """Start the calculations"""
    import os, shutil
    top=os.getcwd()
    #
    # Read the filelist
    #
    if args[0]!='all':
        fd=open(args[0])
        files=fd.readlines()
        fd.close()
    else:
        files=os.listdir(top)
    #
    # -------------------------------------
    #
    # Loop over all files and do the task
    #
    for filename in sorted(files):
        #
        # clean filename
        #
        if filename[0]=='#':
            continue
        if filename[-4:]!='.pdb':
            continue
        import string
        filename=string.strip(filename)
        filename=filename.split()[0]
        if filename[0]=='#':
            continue
        #
        # If we have a file then create a dir
        #
        if options.filestructure=='files':
            if filename[-4:]=='.pdb':
                dirname=filename[:-4]
            else:
                dirname=filename+'_dir'
            if not os.path.isdir(dirname):
                os.mkdir(dirname)
            #
            # copy the pdb file to the dir
            #
            shutil.copy(filename,dirname)
        else:
            dirname=filename
        #
        # Make the dirname absolute
        #
        dirname=os.path.join(top,dirname)
        #
        # Change dir to the calc dir
        #
        os.chdir(dirname)
        #
        # Delete TOPOLOGY.H and DELRAD.DAT + DELCRG.DAT if they exist
        #
        import shutil
        copyfiles=['DELRAD.DAT','DELCRG.DAT','TOPOLOGY.H']
        for copyfile in copyfiles:
            if os.path.lexists(os.path.join(dirname,copyfile)):
                os.unlink(os.path.join(dirname,copyfile))
        #
        # Find the pdb file
        #
        pdbfile=False
        searchnames=[filename,filename+'.pdb']
        for sname in searchnames:
            rname=os.path.join(dirname,sname)
            if os.path.isfile(rname):
                pdbfile=rname
                break
        if not pdbfile:
            raise Exception('Could not find PDB file in %s' %os.getcwd())
        #
        # EM + MD?
        #
        if options.EM:
            #corall(pdbfile)
            class options:

                def __init__(otherself):
                    otherself.type='pKa'
                    otherself.pdbfile=pdbfile
                    otherself.clean=False
                    return
            Goptions=options()
            #
            import GromacsEM
            pdblines=GromacsEM.EMone(Goptions)
            import Protool
            PI=Protool.structureIO()
            PI.parsepdb(pdblines)
            PI.Remove_All_NonAminoAcids() # Make sure all waters are removed
            pdbfile=pdbfile+'test'
            PI.writepdb(pdbfile)
            corall(pdbfile) # Do final corall
            stop
            
        #
        # Should we do a corall?
        #
        if options.corall:
            corall(pdbfile)
        #
        # Copy the DEL* files and TOPOLOGY.H
        #       
        for copyfile in copyfiles:
            shutil.copy(os.path.join(top,copyfile),os.path.join(dirname,copyfile))
        #
        # Instantiate pKarun
        #
        PM=pKarun_base.pKamisc()
        params=PM.parse_parameterlist(['-dbcrit 1000'],skip2first=False)
        print params
        X=pKarun_base.pKarun(os.getcwd(),pdbfile,params)
        #
        # Carry out the tasks
        #
        if options.tasks=='titration':
            print 'Runing solvepka in ',os.getcwd()
            X.solvepka()
        elif options.tasks=='desolv':
            X.desolv()
        elif options.tasks=='backgr':
            X.backgr()
        elif options.tasks=='matrix':
            X.matrix()
        elif options.tasks=='all':
            print 'Running all'
            print options.tasks
            X.all()
        #
        # Change dir back
        #
        os.chdir(top)
        print 'Back in ',os.getcwd()
        import sys
        sys.stdout.flush()
        


if __name__=='__main__':
    print
    print 'Calculate pKa values for multiple structures'
    print 'Jens Erik Nielsen, 2009'
    print
    import sys, os
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <filelist>',version='%prog 1.0')
    parser.add_option('-t',"--tasks", dest='tasks',type='string',action='store',
                        help="Specify which tasks to perform: all, desolv, backgr, matrix, titration. Default; %default",
                        default='all')
    parser.add_option('-f',"--filestructure",type='string',dest='filestructure',default='dirs',action='store',
                      help='Are the PDB files in individual dirs (dirs) or in the same dir (files)? Default: %default')
    parser.add_option('-c','--corall',dest='corall',action='store_true',help='Perform WHAT IF corall call before calculating pKa values. Default=%default',default=False)
    parser.add_option('--EM',dest='EM',action='store',help='Energy mimimize each snapshot with Gromacs before doing pKa calc',default='None')
    (options, args) = parser.parse_args()
    if len(args)!=1:
        parser.error('You must specify a file holding the PDB files/directories')
    #
    # Call main
    #
    main(options,args)
