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

#
"""This script performs a full Alanine scan on a protein and records the change in pKa value
for one or more target residues"""



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
    defaults['pHstart'][0]=0.1
    defaults['pHstop'][0]=12.0
    defaults['pHstep'][0]=0.25
    defaults['pKMCsteps'][0]=200000
    #
    # Design settings
    #
    # Target
    #
    target_residues=target_residues.split(',')
    target_text=''
    for target in target_residues:
        target_text=target_text+target+'=0.0,' # Dummy pKa value
    defaults['pKas'][0]=target_text[:-1]
    #
    # Method
    #
    defaults['MC'][0]=1
    defaults['tabulated'][0]=1
    defaults['MCsteps'][0]=0
    #
    # Be not-so-noisy
    #
    defaults['verbose'][0]=3
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
    defaults['recalc_intpka'][0]=recalc_intpka
    defaults['recalc_intpka_dist'][0]=10.0
    defaults['use_titration_curves'][0]=1
    defaults['calc_dpka'][0]=1
    defaults['generate_mutations'][0]=False
    defaults['ligands'][0]=[]
    defaults['allow_unknown_atoms'][0]=1
    defaults['unknown_rad'][0]=0.0
    return defaults


def main():
    """Load the PDB file, make all Ala mutations and calculate dpKa for the target residues"""
    try:
        import sys
        pdbfile=sys.argv[1]
        target_residues=sys.argv[2]
    except:
        print
        print 'Usage: pKa_alascan.py <pdbfile> <target residues>'
        print 'Example: pKa_alascan.py 2lzt.pdb :0035:GLU,:0052:ASP'
        print 'This command will perform a full Alanine scan and report the effect of each mutation on the pKa values of Glu 35 and Asp 52'
        print 'If ALL is supplied instead of a list of target residues, then dpKas will be calculated for all residues'
        print
        raise Exception,'Incorrect usage'
    #
    # Start the work
    #
    import Protool
    P=Protool.structureIO()
    P.readpdb(pdbfile)
    residues=P.residues.keys()
    residues.sort()
    #
    # All titgroups?
    #
    if target_residues=='ALL':
        titgroups=P.get_titratable_groups()
        import string
        target_residues=string.join(titgroups,',')
    #
    # Start looping
    #
    results={}
    import pickle, os
    for residue in residues:
        #
        # Define result filename
        #
        resultfile=os.path.join(os.getcwd(),'alascan_%s.result' %residue)
        if os.path.isfile(resultfile):
            fd=open(resultfile)
            results[residue]=pickle.load(fd)
            fd.close()
        else:
            if P.resname(residue)=='ALA' or P.resname(residue)=='GLY' or not P.isaa(residue):
                print 'Skipping',residue,P.resname(residue)
                continue
            print 'Calculating for residue',residue,P.resname(residue)
            recalc_intpka=1
            defaults=local_defaults(pdbfile,target_residues,recalc_intpka)
            #
            # Set the mutations
            #
            defaults['mutations'][0]='%s:%s:%s' %(residue,P.resname(residue),'ALA')
            import pKa.Design_pKa as Design_pKa
            #
            # Calculate the dpKas
            #
            solutions,pKd_dict=Design_pKa.run_opt(defaults)
            results[residue]=solutions.copy()
            #
            # Save this result
            #
            fd=open(resultfile,'w')
            pickle.dump(results[residue],fd)
            fd.close()
    #
    # Save all
    #
    name='%s.alascan.pickle' %pdbfile
    fd=open(name,'w')
    import pickle
    pickle.dump(results,fd)
    fd.close()

if __name__=="__main__":
    main()
