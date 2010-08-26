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

gunzip='/usr/sbin/gunzip'

def main():
    #
    # Run pKa calculations in all subdirs given as argument
    #
    import sys
    print
    print 'Run pKa calculations in all subdirs. '
    print
    print 'Looking for pdb files with the name: %s' %sys.argv[1]
    print 'Running pKa calculations for: '
    import os, sys, string
    initial_pdbfile='pdbfile.pdb'
    dirs=sys.argv[1:]
    #
    # if the last argument is '-name', then the name of the dir is part of the pdbfile name
    #
    top=os.getcwd()
    for dir in dirs:
        if os.path.isdir(dir):
            print dir
    print
    print '----------------------'
    for dir in dirs:
        if os.path.isdir(dir):
            initial_pdbfile=dir+'.pdb'
            pkafile=os.path.join(dir,'clean.pka.pdb')
            if not os.path.isfile(pkafile):
                import WI_tools
                logfile,files=WI_tools.delwat_dellig_corall(os.path.join(os.getcwd(),dir,initial_pdbfile),None,1)
                print os.path.join(os.getcwd(),dir,initial_pdbfile)
                for line in logfile:
                    print line,
                fs=files.keys()
                print fs
                if len(fs)!=1:
                    done=None
                    print dir
                    import string
                    pdbfile=string.lower(dir[-4:])+'.pdb.gz'
                    pdbdir='/net/scratch67/jnielsen/pdb'
                    pdbfile=os.path.join(pdbdir,pdbfile)
                    if os.path.isfile(pdbfile):
                        print 'Could not find file'
                        #a=raw_input('Get file from pdb dir? ')
                        #a='y'
                        a='n'
                        if a=='y':
                            destination=os.path.join(os.getcwd(),dir,os.path.split(pdbfile)[1])
                            print destination
                            os.system('cp '+pdbfile+' '+destination)
                            os.system(gunzip+' '+destination)
                            destination=destination[:-3]
                            print destination
                            os.link(destination,os.path.join(os.getcwd(),dir,initial_pdbfile))
                            done=1
                    if not done:
                        print initial_pdbfile
                        raise 'Huh? - could not clean PDB file'
                else:
                    fd=open(pkafile,'w')
                    fd2=open(pkafile+'.backup','w')
                    for line in files[fs[0]]:
                        fd2.write(line)
                        fd.write(line)
                    fd.close()
                    fd2.close()
                    print 'Cleaned PDB file written'
            #
            # Prepare the pKa calculation (write the Invocation file)
            #
            import splitpdb
            splitpdb.prep_pKa(dir,os.path.split(pkafile)[1])
            os.chdir(dir)
            import pKa
            #
            # Read the parameters
            #
            infile='Invocation'
            fd=open(infile)
            line=fd.readline()
            fd.close()
            argv=string.split(line)
            #
            # Do we need to run the pKa calculations for this dir?
            #
            pdbfile=argv[1]
            if not os.path.isfile(pdbfile+'.PKA.DAT'):
                print 'Running pKa calculations for: %s ' %dir
                print
                print 'Using command:',string.join(argv)
                #
                # Parse parameters
                #
                x=pKa.pKamisc()
                params=x.parse_parameterlist(argv)

                X=pKa.pKarun(os.getcwd(),pdbfile,params)
                X.runseq()
                print '============================='
                print
            os.chdir(top)
    return

main()
        
    
