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

"""This script calculates delta pKa values for one or more residues for a list of mutations"""



#
# -------------------------
#
def local_defaults(pdbfile,target_residues,recalc_intpka):
    """
    # Set the parameters that are the same for all mutations
    """
    import pKa.Design_pKa as Design_pKa
    defaults=Design_pKa.get_defaults()
    # PDB file
    defaults['pdb'][0]=pdbfile
    #
    # pKa calculation parameters
    #
    #defaults['pHstart'][0]=0.1
    #defaults['pHstop'][0]=12.0
    defaults['pHstep'][0]=0.01
    defaults['pKMCsteps'][0]=200000
    #
    # Design settings
    #
    # Target
    #
    #target_residues=target_residues.split(',')
    target_text=''
    for target in target_residues:
        target_text=target_text+target+'=0.0,' # Dummy pKa value
    defaults['pKas'][0]=target_text[:-1]
    #
    # Method
    #
    defaults['dpKa_method']=['MC','junk']
    defaults['tabulated'][0]=0
    defaults['MCsteps'][0]=0
    defaults['stability_mode']=[False,'junk']
    defaults['PBEsolver']=['DelPhi','junk']
    #
    # Be not-so-noisy
    #
    defaults['verbose'][0]=5
    #
    # Minimum distance between target and mutation
    #
    #defaults['min_target_dist'][0]=min_dist
    #
    # Do not save the solutions
    #
    defaults['save_solutions'][0]=None
    #
    #
    defaults['recalc_intpka'][0]=options.recalc_intpka
    defaults['recalc_intpka_dist'][0]=options.recalc_intpka_dist
    defaults['use_titration_curves'][0]=1
    defaults['calc_dpka'][0]=1
    defaults['generate_mutations'][0]=False
    defaults['mutation_quality'][0]=0.5
    return defaults

#
# ----
#

def main(options,args):
    """Load the PDB file and the list of mutations"""
    import sys
    pdbfile=args[0]
    target_residues=options.target_groups
    #
    # Load the PDB file
    #
    import Protool
    P=Protool.structureIO()
    P.readpdb(pdbfile)
    residues=P.residues.keys()
    residues.sort()
    if not options.allmutations:
        #
        # Load the mutations
        #
        mutfile=args[1]
        fd=open(mutfile)
        mutlines=fd.readlines()
        fd.close()
    else:
       mutlines=[]
       aas=P.trueaminoacids.keys()
       aas.sort()
       count=1
       for residue in residues:
           for aa in aas:
               if aa==P.resname(residue):
                   continue
               mutlines.append('clone%d,%s:%s:%s' %(count,residue,P.resname(residue),aa))
               count=count+1
       print 'Created %d mutant proteins each containing 1 mutation' %len(mutlines)
    #
    # Make the resultdir
    #
    import os
    resultdir=os.path.join(os.getcwd(),'pKa_mutscan_results')
    if not os.path.isdir(resultdir):
        os.mkdir(resultdir)
    #
    # which target residues
    #
    if target_residues==[] or target_residues==['ALL']:
        target_residues=P.get_titratable_groups()
        import string
        target_residues=string.join(target_residues,',')
    results={}
    import pickle, os
    for mline in mutlines:
        import string
        if mline[0]=='#' or mline[:2]=='//':
            continue
        #
        line=string.strip(mline)
        sp_line=line.split(',')
        variant_name=sp_line[0]
        mutation=sp_line[1]
        print 'Variant: %s, mutations: %s' %(variant_name,mutation)
        if mutation.find('insert')!=-1:
            print 'Skipping insertions'
            continue
        #
        # Define result filename
        #
        resultfile=os.path.join(resultdir,'mutscan_%s.result' %variant_name)
        if os.path.isfile(resultfile):
            fd=open(resultfile)
            results[variant_name]=pickle.load(fd)
            fd.close()
            #print 'Already did',mutation
        else:
            recalc_intpka=1
            defaults=local_defaults(pdbfile,target_residues,recalc_intpka)
            #
            # Set the mutations
            #
            import string
            defaults['mutations'][0]=string.strip(mutation)
            print 'Calculating for',mutation
            import pKa.Design_pKa as Design_pKa
            #
            # Set other parameters
            #
            defaults['ion'][0]=options.ion
            #
            # Calculate the dpKas
            #
            #try:
            solutions=Design_pKa.run_opt(defaults)
            #except Exception,inst:
            #    if str(inst).find('Cannot model mutant')!=-1:
            #        solutions='Cannot model mutant'
            #        raise Exception('Cannot model mutant')
            #    elif str(inst).find('We cannot model insertions')!=-1:
            #        solutions='Skipping insertions'
            #    else:
            #        print inst
            #        raise Exception(str(inst))
            print
            print
            print 'Results are ',solutions
            results[variant_name]=solutions
            #
            # Save this result
            #
            print 'Saving',results[variant_name],'in',resultfile
            import os
            if len(os.path.split(resultfile)[1])>80:
                continue
            
            fd=open(resultfile,'w')
            pickle.dump(results[variant_name],fd)
            print '*********************'
            fd.close()
    #
    # Save all
    #
    name=os.path.join(os.getcwd(),'%s.mutscan.pickle' %pdbfile)
    fd=open(name,'w')
    import pickle
    pickle.dump(results,fd)
    fd.close()

if __name__=="__main__":
    print
    print 'Calculate changes in pKa values for a list of mutations'
    print '(c) Copyright Jens Erik Nielsen, 2008-2010, All rights reserved'
    print
    import sys, os
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <pdbfile> <mutation_file>',version='%prog 1.0')
    parser.add_option('-t',"--target", dest='target_groups',type='string',action='append',
                        help="Residues to calculate dpKa values for. Specify multiple times. Default: ALL",default=[])
    parser.add_option('-i','--recalc_intpka',dest='recalc_intpka',action='store_true',
                      help='Recalculate intrinsic pKa for target residues. Default= %default',default=True)
    parser.add_option('-j','--no_recalc_intpka',dest='recalc_intpka',action='store_false',
                      help='Do not recalculate intrinsic pKa value for target')
    parser.add_option('-e','--recalc_intpka_dist',dest='recalc_intpka_dist',action='store',
                      help='Mutations closer than this distance to the target group will force a recalculation of the intrinsic pKa of the target. Default: %default A',
                      default=10)
    parser.add_option('-m','--ionic_strength',dest='ion',type='float',action='store',
                      help='ionic strength to use in the calculations. Default= %default',default=0.144)
    parser.add_option('-a','--all_mutations',dest='allmutations',action='store_true',
                      help='Ignore the mutation file and try all possible mutations. Default= %default',default=False)
    
    (options, args) = parser.parse_args()
    if len(args)!=2 and not (options.allmutations and len(args)==1):
        parser.error('You must specify a PDB file and a mutation file, or -a and a PDB file.')
    #
    # Call main
    #
    main(options,args)   
