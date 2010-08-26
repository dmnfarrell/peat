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

def main(options,args):
    import sys
    filename=args[0]
    dirname=options.outdir
    MCsteps=options.MCsteps
    #
    fd=open(filename)
    lines=fd.readlines()
    fd.close()
    for line in lines[1:]:
        if line[0]=='#':
            continue
        import string
        print string.strip(line)
        sp=line.split()
        import os
        pdbfilename=os.path.join(os.getcwd(),sp[0],sp[1])
        import pKaTool.get_CCPS_population
        PD=sp[2]
        NU=sp[3]
        if PD.lower().find('mutated')==-1 and NU.lower().find('mutated')==-1:
            #
            # find the residue number of the HIS
            #
            import make_pKaSens_table as mpt
            HIS,other_side=mpt.get_HIS(pdbfilename,PD,NU)
            #
            # Get the CCPS pop
            #
            CCPS=['%s:GLU=0' %PD,'%s:GLU=-1' %NU]
            if HIS:
                CCPS.append('%s:HIS=1' %HIS)
            picklefile=os.path.join(os.getcwd(),dirname,sp[0]+'.pickle')
            if os.path.isfile(picklefile):
                continue
            #
            # Do the calc
            #
            data=pKaTool.get_CCPS_population.get_CCPS_population(pdbfilename,ccps=CCPS,MCsteps=MCsteps)
            #
            # Save the pickle file
            #
            import pickle
            fd=open(picklefile,'w')
            pickle.dump(data,fd)
            fd.close()
            print '---------------------------'
            print

if __name__=='__main__':
    print
    print 'Get pH-activity profiles from a set of pKa calculations'
    print 'Jens Erik Nielsen, 2009'
    print
    import sys, os
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <file>',version='%prog 1.0')
    parser.add_option('-m',"--MCsteps",type='int',dest='MCsteps',action='store',
                      help='Number of Monte Carlo steps to perform when calculating titration curves',default=200000)
    parser.add_option('-o',"--outputdir",type='string',dest='outdir',action='store',
                      help='Dir where result files will be stored')
    (options, args) = parser.parse_args()
    if len(args)!=1:
        parser.error('You must specify a file with pKa run information')
    #
    # Call main
    #
    main(options,args)

