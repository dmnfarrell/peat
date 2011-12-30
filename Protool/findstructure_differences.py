#!/usr/bin/env python
#
# # Protool - Python class for manipulating protein structures
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

import Protool.PDBServices
PDB=Protool.PDBServices.PDBServices()

def readchain(chainnum):
    sp=chainnum.split(':')
    PDBID=sp[0]
    CIDnum=int(sp[1])
    chains=PDB.getChains(PDBID)
    CID=chains[CIDnum-1]
    return readpdb(PDBID,CID)

def readpdb(structure,CID):
    """Get a PDBID from the PDB website and read it into Protool"""
    print 'Getting PDBID:chain %s:%s from the PDB website' %(structure,CID)
    PDBID=structure
    pdblines=PDB.getPDB(PDBID)
    print 'Reading the PDB into Protool'
    import Protool
    X=Protool.structureIO()
    X.parsepdb(pdblines)
    # 
    # Delete all residues except the ones in the chain we need
    #
    chains=X.chains.keys()
    for chain in chains:
        if chain!=CID:
            print 'Deleting chain',chain
            for res in X.chains[chain]:
                X.Delete_residue(res,update=False)
    X.Update()
    return X

#
# ----
#

def alignsequences(pdb1,pdb2,name1='pdb1',name2='pdb2'):
    """Align the sequences in the two PDB files"""
    import PEATDB.sequence_alignment as SA
    seq1=pdb1.PirSeq()
    seq2=pdb2.PirSeq()
    ALIGN=SA.NW(seq1,seq2,gap=5.0)
    al1,al2,alres1,alres2=ALIGN.Align(verbose=False)
    #if ALIGN.sequence_identity<95.0:
    #    print 'Sequence identity too low: %f' %ALIGN.sequence_identity
    #    return False
    #
    # Now relate the PDB sequence numbers to each other
    #
    pdb1_residues=pdb1.residues.keys()
    pdb1_residues.sort()
    pdb2_residues=pdb2.residues.keys()
    pdb2_residues.sort()
    align_map={name1:{},name2:{}}
    count=0
    for res in alres1:
        if res!='-':
            align_map[name1][pdb1_residues[count]]=pdb2_residues[res]
        else:
            align_map[name1][pdb1_residues[count]]=None
        count=count+1
    count=0
    for res in alres2:
        if res!='-':
            align_map[name2][pdb2_residues[count]]=pdb1_residues[res]
        else:
            align_map[name2][pdb2_residues[count]]=None
        count=count+1
    return align_map

def superimpose(structure1,structure2):
    import PEATDB.sequence_alignment as SA
    pdb1=readchain(structure1)
    pdb2=readchain(structure2)
    print alignsequences(pdb1,pdb2)
    stop
    
def main():
    # Get the query structure
    PDBID='2lzt'
    CID='A'
    X=readpdb(PDBID,CID)

    print 'Getting similar structures'
    # Get its PIR sequence
    IDs=PDB.FASTa(PDBID,CID)

    print 'Found %d similar structures' %(len(IDs))
    #
    # Check that all of these structures are similar enough to the query
    #
    N=10
    ensemble=[]
    for structure in IDs:
        pdb=readchain(structure)
        seqid=alignsequences(X,pdb)
        if seqid>=95.0:
            ensemble.append(structure)
        if len(ensemble)>=10:
            break
    print 'The ensemble consists of %d structures' %len(ensemble)
    print ensemble
    print 'Now calculating per residue variability'
    #
    # Superimpose the first N structures and compile residue-by-residue variability distribution
    #

    Cdiffs={}
    for structure1 in ensemble[:N]:
        print 'Superimposing all on %s' %structure1
        for structure2 in ensemble[:N]:
            if structure1==structure2:
                continue
            superimpose(structure1,structure2)
    return
    
if __name__=='__main__':
    main()


