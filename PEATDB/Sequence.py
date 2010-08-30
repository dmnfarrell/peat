#!/usr/bin/env python
#
# Protein Engineering Analysis Tool DataBase (PEATDB)
# Copyright (C) 2010 Damien Farrell & Jens Erik Nielsen
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

import string

class SequenceOperations:

    @classmethod
    def makeOperations(self, sequence, operations):
        """Perform the specified operations on the sequence
        Sequence must be in the [[A:0001:ALA],[A:0002:GLU],['A:0003:THR'], ..] format
        Operations is a list of the following types:
        Mutations: A:0001:ALA:ASP
        Deletions: delete:A:0002:GLU
        Insertions: insert:1:A:0003:THR:ALA, insert:2:A:0003:THR:TRP (insert THR,TRP after A:0003:THR)
        Operations are always performed in sequence numbering order
        """
        if operations==[]:
            return sequence
        ops_sorted={}
        insertions=[]
        for operation in operations:
            s_op=operation.split(':')
            if s_op[0]=='insert':
                resid='%s:%s' %(s_op[2],s_op[3])
                if ops_sorted.has_key(resid):
                    ok=False
                    if type(ops_sorted[resid]) is type(list):
                        if ops_sorted[resid][0]=='insert':
                            ok=True
                    if not ok:
                        raise Exception('More than one operation on the same residue: %s' %resid)
                else:
                    ops_sorted[resid]=['insert',{}]
                #
                # Add the residue to be inserted
                #
                ins_num=s_op[1]
                org_typ=s_op[4]
                ins_typ=s_op[5]
                ops_sorted[resid][1][ins_num]=[org_typ,ins_typ]
            elif s_op[0]=='delete':
                resid='%s:%s' %(s_op[1],s_op[2])
                if ops_sorted.has_key(resid):
                    raise Exception('More than one operation on the same residue: %s' %resid)
                restyp=s_op[3]
                ops_sorted[resid]=['delete',restyp]
            else:
                # Normal mutation
                import pKa.pKD_tools as pKD_tools
                resid=pKD_tools.get_resid_from_mut(operation)
                if ops_sorted.has_key(resid):
                    raise Exception('More than one operation on the same residue: %s' %resid)
                ops_sorted[resid]=['mutate',operation]
        #
        # Perform the operations, one after one
        #
        new_seq=[]
        new_count=None
        new_chain=None
        for resid,restyp in sequence:
            # Make sure that the chain hasn't changed or if we are at the beginning then init
            if resid.split(':')[0]!=new_chain:
                #Initialise
                sp_resid=resid.split(':')
                new_chain=sp_resid[0]
                new_count=int(sp_resid[1])
                newresid='%s:%s' %(new_chain,string.zfill(new_count,4))
            # Does this residue have an operation?
            if ops_sorted.has_key(resid):
                op=ops_sorted[resid]
                if op[0]=='delete':
                    # Deletion
                    if op[1]==restyp:
                        pass # This deletes the residue
                    else:
                        raise Exception('Incorrect org residue in deletion: %s' %op)
                elif op[0]=='insert':
                    # Insertion
                    inserts=op[1].keys()
                    inserts.sort()
                    for i in inserts:
                        if i[0]==restyp:
                            new_seq.append([newresid,i[1]])
                            new_count=new_count+1
                            newresid='%s:%s' %(new_chain,string.zfill(new_count,4))
                elif op[0]=='mutate':
                    # Mutation
                    import pKa.pKD_tools as pKD_tools
                    orgres=pKD_tools.get_oldrestyp_from_mut(op[1])
                    if orgres==restyp:
                        new_seq.append([newresid,pKD_tools.get_newrestyp_from_mut(op[1])])
                        new_count=new_count+1
                        newresid='%s:%s' %(new_chain,string.zfill(new_count,4))
                    pass
                else:
                    raise Exception('Unknown mutations spec: %s' %op)
            else:
                new_seq.append([resid,restyp])
                new_count=new_count+1
                newresid='%s:%s' %(new_chain,string.zfill(new_count,4))
        return new_seq

    @classmethod
    def findSequenceDifferences(self, record_sequence, parent_sequence,
                                full_parent_sequence, PDBaln,recordALN):
        """
        # Find all amino acid differences
        """
        
        # Loop over the sequences - changed are record from parent -> child
        
        import string
        operations=[]
        import DNAtool.mutation
        Cterm_add=0
        insert_num=0
        full_parent_count=0
        for count in range(len(record_sequence)):
            parent_aa=parent_sequence[count]
            child_aa=record_sequence[count]
            if parent_aa!=child_aa:
                if len(full_parent_sequence)<=count:
                    Cterm_add=Cterm_add+1
                    aa_identifier=full_parent_sequence[-1][0]#+'+%d' %Cterm_add
                else:
                    # Find the PDB file number
                    aa_identifier=full_parent_sequence[full_parent_count][0]
                if aa_identifier[-1]==':':
                    aa_identifier=aa_identifier[:-1]
                #
                # Convert to 3letter format
                #
                if parent_aa!='-':
                    full_parent_count=full_parent_count+1
                    parent_aa=DNAtool.mutation.one_to_three[parent_aa]
                if child_aa!='-':
                    child_aa=DNAtool.mutation.one_to_three[child_aa]
                if parent_aa=='-':
                    operations.append('insert%d:%s:%s' %(insert_num,aa_identifier,child_aa))
                    insert_num=insert_num+1
                elif child_aa=='-':
                    insert_num=0
                    operations.append('delete:%s:%s' %(aa_identifier,parent_aa))
                else:
                    insert_num=0
                    operations.append('%s:%s:%s' %(aa_identifier,parent_aa,child_aa))
            else:
                full_parent_count=full_parent_count+1
        return operations

    @classmethod    
    def shortenOperations(self, operations):
        """Provide the more conventional short form of the mutations"""
        new_ops=[]
        import DNAtool.mutation
        import pKa.pKD_tools as pKD_tools
        for operation in operations:
            if operation.find('insert')!=-1:
                text=operation.replace('insert','').split(':')
                insertnum=int(text[0])
                resnum='%s:%s' %(text[1],text[2])
                resnum=int(pKD_tools.get_intresnum_from_res(resnum))
                chainID=text[1]
                if len(chainID)!=0:
                    chainID='(%s)' %chainID
                insertchar='abcdefghijklmnopqrstuvw'
                new=text[-1]
                new_ops.append('*%s%d%s%s' %(chainID,resnum,insertchar[insertnum],DNAtool.mutation.three_to_one[new]))
            elif operation.find('delete:')!=-1:
                text=operation.replace('delete:','')
                restext=text.split(':')
                chainID=restext[0]
                if len(chainID)!=0:
                    chainID='(%s)' %chainID
                resnum=int(pKD_tools.get_intresnum_from_res(text))
                old=restext[-1]
                new_ops.append('%s%s%d*' %(DNAtool.mutation.three_to_one[old],chainID,resnum))
            else:
                new=pKD_tools.get_newrestyp_from_mut(operation)
                old=pKD_tools.get_oldrestyp_from_mut(operation)
                chainID=operation.split(':')[0]
                resnum=pKD_tools.get_intresnum_from_res(operation)
                if len(chainID)!=0:
                    chainID='(%s)' %chainID
                new_ops.append('%s%s%d%s' %(DNAtool.mutation.three_to_one[old],
                                            chainID,
                                            resnum,
                                            DNAtool.mutation.three_to_one[new]))
        return new_ops


    def getAllRecordMutations(self):
        """Save a file specifying the mutations of all records"""
        refprot=self.get_refprotein()
        info=[]
        for protein in self.proteins:
            if refprot is None or not self.DB.has_key(refprot):
                info.append('%s,%s\n' %(protein,'NA'))
            #
            # Get the mutations
            #
            is_parent,operations=self.isparent(protein,refprot)
            if is_parent:
                import string
                mut_text=string.join(operations,'+')
            else:
                mut_text='NA'
            info.append('%s,%s\n' %(protein,mut_text))
        return info
        
    def getMutatedResidues(self,record):
        """Get all the mutated residues for a single protein"""
        refprot=self.get_refprotein()
        if refprot is None or not self.DB.has_key(refprot):
            info.append('%s,%s\n' %(protein,'NA'))
        #
        # Get the mutations
        #
        is_parent,operations=self.isparent(record,refprot)
        residues=[]
        if not operations:
            return []
        for op in operations:
            import pKa.pKD_tools
            residues.append(pKa.pKD_tools.get_resid_from_mut(op))
        return residues
