#!/usr/bin/env python
#
# This file is part of the pKaTool package
# Copyright (C) Jens Erik Nielsen 2008
# All rights reserved
#
#comment
def Nonefix(value):
    if value:
        return '%6.3f' %(value)
    else:
        return 'NA'

def main(options,args):
    """Calculate the change in chemical shift due to a full charge on each titratable group"""
    import get_dEF
    method=options.method 
    X=get_dEF.map_struct(args[0])
    X.build_Hs()
    residues=X.PI.residues.keys()
    residues.sort()
    #
    # Get the titratable groups
    #
    titgroups=X.PI.get_titratable_groups()
    titgroups.sort()
    import pKa.pKD_tools
    if options.group_by_titgroup:
        for titgroup in titgroups:
            titgroup_type=pKa.pKD_tools.get_titgroup_type_from_titgroup(titgroup)
            charge=X.PI.titgroups[titgroup_type]['charge']
            print 'TITRATABLE GROUP',titgroup
            print 'Residue  CS Nitrogen    CS Hydrogen'
            
            for residue in residues:
                dCS_N=X.get_dCS(residue+':N',titgroup,charge=charge,method=method)
                dCS_H=X.get_dCS(residue+':H',titgroup,charge=charge,method=method)
                print '%8s,    %s,    %s' %(residue,Nonefix(dCS_N),Nonefix(dCS_H))
    else:
        #
        # Group by atom
        #
        for residue in residues:
            for atom in [':N',':H']:
                changes=[]
                for titgroup in titgroups:
                    titgroup_type=pKa.pKD_tools.get_titgroup_type_from_titgroup(titgroup)
                    charge=X.PI.titgroups[titgroup_type]['charge']
                    dCS=X.get_dCS(residue+atom,titgroup,charge=charge,method=method)
                    changes.append([titgroup,dCS])
                #
                if options.sort:
                    def cmpfunc(x,y):
                        if x[1] is None:
                            return 1
                        if y[1] is None:
                            return -1
                        return cmp(abs(y[1]),abs(x[1]))
                    changes.sort(cmpfunc)
                print 'Residue: %s, Atom: %s' %(residue,atom[1])
                for titgroup,dCS in changes[:options.listnumber]:
                    if dCS:
                        if abs(dCS)<options.limit:
                            continue
                        print titgroup,dCS
                print
    return

if __name__=='__main__':
    print
    print 'Calculate charge-induced Chemical shift changes in protein backbone atoms'
    print 'Jens Erik Nielsen, 2008'
    print
    import sys, os
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <pdbfile>',version='%prog 1.0')
    parser.add_option('-m',"--method", dest='method',
                        help="Choose method for calculating electric field: APBS [PBE solver], PBE [Chresten's PBE solver], Coulomb [Coulomb's law]. Default: %default",
                        default='APBS')
    parser.add_option('-n',"--number",type='int',dest='listnumber',default=100000,action='store',
                    help='Print X values for each grouping. Default: %default',metavar='X')
    parser.add_option('-l',"--limit",type='float',dest='limit',default=0.0,action='store',
                    help='Do not list chemical shift differences where abs(d chem shift)< LIMIT. Default: LIMIT=%default',metavar='LIMIT')
    parser.add_option('-g','--group',dest='group_by_titgroup',action='store_true',default=True,
                    help='Group results by titratable group. Default: %default')
    parser.add_option('-a','--atom',dest='group_by_titgroup',action='store_false',
                        help='Group results by atom. Default: False')
    parser.add_option('-s','--sort',dest='sort',action='store_true',default=False,
                        help='Sort chemical shift differences. Default: %default')                    
    (options, args) = parser.parse_args()
    if len(args)!=1:
        parser.error('You must specify a PDB file')
    if options.sort and options.group_by_titgroup:
        parser.error('Sorting not (yet) implemented when grouping by titratable group')
    #
    # Call main
    #
    main(options,args)
