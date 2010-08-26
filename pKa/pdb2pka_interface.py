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

pdb2pka_path='/home/people/nielsen/lib/pdb2pqr/branches/yong-dev/pdb2pqr/pdb2pka'
import sys
sys.path.append(pdb2pka_path)
import pkanew as pdb2pka

class P2pI:

    def __init__(self,topdir,pdbfile,params,parent=None):
        """Init"""
        self.topdir=topdir
        self.pdbfile=pdbfile
        self.params=params
        self.parent=parent
        return

    #
    # -----
    #

    def pdb2pka_sugelm(self):
        """Explore all possible mutations and calculate a phimap for each using pdb2pka (APBS)"""
        import Protool
        P=Protool.structureIO()
        P.readpdb(self.pdbfile)
        P.RemoveALT()
        #import Protool.mutate
        #MUT=Protool.mutate.Mutate(P)
        #
        # Construct arrays
        #
        import pKD_dict
        self.data=pKD_dict.pKD_dict()
        self.atom_data=pKD_dict.pKD_dict()
        #
        # Create dir for mutant PDB files
        #
        import os
        mutdir=os.path.join(self.topdir,self.pdbfile+'.pdbs')
        if not os.path.isdir(mutdir):
            os.mkdir(mutdir)
        #
        # Loop over all residues
        #
        residues=P.residues.keys()
        residues.sort()
        for residue in residues:
            orgres=P.resname(residue)
            print 'Calculating for %s %s' %(residue,P.resname(residue))
            #
            # If neutral mutate to Asp, Glu, Lys, Arg, His
            #
            targets=[]
            for res in ['ARG','LYS','HIS','ASP','GLU']:
                if P.resname(residue)!=res:
                    targets.append(res)
            #if orgres=='GLU':
            #    targets.append('GLN')
            #elif orgres=='ASP':
            #    targets.append('ASN')
            #elif orgres=='HIS':
            #    targets.append('PHE')
            #elif orgres=='ARG' or P.resname(residue)=='LYS':
            #    targets.append('MET')
            #
            # Target identified. Now model each
            #
            for target in targets:
                import pKD_tools
                resid=pKD_tools.get_resid_from_res(residue)
                orgres=P.resname(residue)
                filename=os.path.join(mutdir,'%s:%s:%s.pdb' %(residue,orgres,target))
                mutation='%s:%s:%s' %(residue,orgres,target)
                if not os.path.isfile(filename):
                    import Design_pKa_help
                    Design_pKa_help.make_mutation(self.pdbfile,mutation)
                NP=Protool.structureIO()
                NP.readpdb(filename)
                NP.writepdb(filename,TER=None)
                #
                # Calculate the interaction energies
                #
                protein,routines,forcefield,apbs_setup,lig_titgrps = pdb2pka.pre_init(pdbfilename=filename,
                                                                                      ff='parse',
                                                                                      ligand=None,
                                                                                      verbose=1)
                mypkaRoutines = pdb2pka.pKaRoutines(protein, routines, forcefield,apbs_setup)
                #
                # Find our group
                #
                sp=residue.split(':')
                chainid=sp[0]
                resnum=int(sp[1])
                mypkaRoutines.findTitratableGroups()
                this_pKa=None
                for pKa in mypkaRoutines.pKas:
                    print pKa.residue.resSeq,resnum
                    print pKa.residue.chainID,chainid
                    print pKa.residue.name,target
                    print pKa.pKaGroup.name,target
                    print '--------------'
                    print 'ChainID',pKa.residue.chainID
                    if pKa.residue.resSeq==resnum and pKa.residue.chainID==chainid and pKa.residue.name==target and pKa.pKaGroup.name==target:
                        #print 'Found group',pKa.residue.resSeq,pKa.pKaGroup.name
                        this_pKa=pKa
                        break
                if not this_pKa:
                    raise Exception,'Could not find inserted titratable group'
                mypkaRoutines.get_interaction_energies_setup(this_pKa,mode='pKD')
                matrix=mypkaRoutines.matrix
                #
                # Dig the interaction energies out of the pdb2pka array
                #
                for titration1 in matrix[this_pKa].keys():
                    for state1 in matrix[this_pKa][titration1].keys():
                        grp_sub=matrix[this_pKa][titration1][state1]
                        if mypkaRoutines.is_charged(this_pKa,titration1,state1):
                            for pKa2 in grp_sub.keys(): 
                                import string
                                chainID2=pKa.residue.chainID
                                resid2='%s:%s' %(chainID2,string.zfill(pKa2.residue.resSeq,4))
                                for titration2 in grp_sub[pKa2].keys():
                                    for state2 in grp_sub[pKa2][titration2].keys():
                                        if mypkaRoutines.is_charged(pKa2,titration2,state2):
                                            #
                                            # Both states are charged, so now we can pull the
                                            # interaction energies out
                                            #
                                            if not self.data.has_key(mutation):
                                                self.data[mutation]={}
                                            self.data[mutation][resid2]=grp_sub[pKa2][titration2][state2]
                                            #
                                            # Get the potentials at all atoms too
                                            #
                                            all_pots=mypkaRoutines.all_potentials[this_pKa][titration1][state1]
                                            sub_all_pots=all_pots[pKa2][titration2][state2]
                                            for atom in sub_all_pots.keys():
                                                resid=mutation
                                                import pKD_tools
                                                resid2=pKD_tools.get_resid_from_res(atom)
                                                atomname=atom.split(':')[-1] #atom.name
                                                if atomname[0]=='H' or atomname in ['N','C','O']:
                                                    continue # Skip all H atoms and all non-CA backbone atoms to save memory
                                                if not self.atom_data.has_key(resid):
                                                    self.atom_data[resid]={}
                                                if not self.atom_data[resid].has_key(resid2):
                                                    self.atom_data[resid][resid2]={}
                                                self.atom_data[resid][resid2][atomname]=abs(sub_all_pots[atom])
        return self.data,self.atom_data

    #
    # ----
    #

    def pdb2pka_desolv_backgr(self,residue):
        """Calculate the desolvation and background interaction energy for a residue"""
        protein,routines,forcefield,apbs_setup,lig_titgrps = pdb2pka.pre_init(pdbfilename=self.pdbfile,
                                                                              ff='parse',
                                                                              ligand=None,
                                                                              verbose=1)
        mypkaRoutines = pdb2pka.pKaRoutines(protein, routines, forcefield,apbs_setup)
        #
        # Find our group
        #
        sp=residue.split(':')
        chainid=sp[0]
        #if chainid!='':
        #    raise Exception,'pKD cannot handle PDB files with ChainIDs!'
        resnum=int(sp[1])
        target=sp[2]
        mypkaRoutines.findTitratableGroups()
        this_pKa=None
        for pKa in mypkaRoutines.pKas:
            print pKa.residue
            print pKa.uniqueid
            print pKa.residue.resSeq,resnum, pKa.residue.chainID,'CID',chainid,pKa.residue.name,target,pKa.pKaGroup.name,target
            #if pKa.residue.resSeq==resnum and pKa.residue.chainID==chainid and pKa.residue.name==target and pKa.pKaGroup.name==target:
            if pKa.residue.resSeq==resnum and pKa.residue.chainID==chainid and pKa.residue.name==target and pKa.pKaGroup.name==target:
                this_pKa=pKa
                break
        if not this_pKa:
            raise Exception,'Could not find inserted titratable group'
        mypkaRoutines.calculateBackground(onlypKa=pKa)
        mypkaRoutines.calculateDesolvation(onlypKa=pKa)
        #
        # Get the intrinsic pKa
        #
        pKaGroup=pKa.pKaGroup
        Gtype=pKa.pKaGroup.type
        #
        # We measure intrinsic pKa values against a single reference state
        #
        desolv=[]
        backgr=[]
        import math
        ln10=math.log(10)
        for titration in pKaGroup.DefTitrations:
            #
            # Find an uncharged reference state
            #
            ref_state=mypkaRoutines.neutral_ref_state[pKa][titration]
            all_states=titration.allstates
            all_states.sort()
            for state in all_states:
                if mypkaRoutines.is_charged(pKa,titration,state)==1:
                    dpKa_desolv=(pKa.desolvation[state]-pKa.desolvation[ref_state])/ln10
                    dpKa_backgr=(pKa.background[state]-pKa.background[ref_state])/ln10
                    #
                    # Make acid and base modifications
                    #
                    if Gtype=='base':
                        dpKa_desolv=-dpKa_desolv
                        dpKa_backgr=-dpKa_backgr
                    #
                    # Now calculate intrinsic pKa
                    #
                    backgr.append(dpKa_backgr)
                    desolv.append(dpKa_desolv)
                    intpKa=titration.modelpKa+dpKa_desolv+dpKa_backgr
                    #print 'Energy difference for %s -> %s [reference state] is %5.2f pKa units' %(state,ref_state,intpKa)
                    pKa.intrinsic_pKa[state]=intpKa
                else:
                    #
                    # Neutral states - why do we treat them differently?
                    #
                    dpKa_desolv=(pKa.desolvation[state]-pKa.desolvation[ref_state])/ln10
                    dpKa_backgr=(pKa.background[state]-pKa.background[ref_state])/ln10
                    #
                    # Make acid and base modifications
                    #
                    if Gtype=='base':
                        dpKa_desolv=-dpKa_desolv
                        dpKa_backgr=-dpKa_backgr
                    backgr.append(dpKa_backgr)
                    desolv.append(dpKa_desolv)
                    dpKa=dpKa_desolv+dpKa_backgr
                    #print 'Energy difference for %s -> %s [reference state] is %5.2f kT' %(state,ref_state,dpKa)
                    pKa.intrinsic_pKa[state]=dpKa
        #print '-----------------'
        #
        # One value is zero, so just get the avg of the rest
        #
        des_sum=0.0
        for val in desolv:
            des_sum=des_sum+val
        des_sum=des_sum/float(len(desolv)-1)
        #
        bac_sum=0.0
        for val in backgr:
            bac_sum=bac_sum+val
        bac_sum=bac_sum/float(len(backgr)-1)
        return  {residue:des_sum},{residue:bac_sum}
