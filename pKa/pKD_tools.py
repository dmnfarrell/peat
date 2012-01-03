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


def filter_deletion_insertion(mutation):
    """Process insertions and mutations so we get the right residue numbers"""
    if mutation.split(':')[0]=='delete':
        # Right now we model deletions as glycines
        mutation=mutation.replace('delete:','')
        mutation=mutation+':GLY'
    elif mutation.find('insert')!=-1:
        sp=mutation.split(':')
        mutation='%s:%s:NON:%s' %(sp[1],sp[2],sp[3])
    elif mutation.find('+')!=-1:
        raise Exception('Something went wrong. I should have gotten a single mutation and got %s' %mutation)
    return mutation

def get_intresnum_from_mut(mutation):
    #
    # Given a mutation (A:0112:ASN:ARG) return the int residue number 112
    #
    return int(get_resnum_from_mut(mutation))

def get_resnum_from_mut(mutation):
    #
    # Given a mutation (A:0112:ASN:ARG) return
    # the residue number as a string '0112'
    #
    mutation=filter_deletion_insertion(mutation)
    import string    
    return string.split(mutation,':')[1]

def get_resid_from_mut(mutation):
    """Given a mutation (A:0112:ASN:ARG) return the residue id (A:0112)"""
    mutation=filter_deletion_insertion(mutation)
    return '%s:%s' %(mutation.split(':')[0],mutation.split(':')[1])

def get_newrestyp_from_mut(mutation):
    #
    # Given a mutation (A:0112:ASN:ARG) return the new residue (ARG)
    #
    import string
    mutation=filter_deletion_insertion(mutation)
    return string.split(mutation,':')[-1]

def get_newres_from_mut(mutation):
    #
    # Given a mutation (A:0112:ASN:ARG) return the new residue (A:0112:ARG)
    #
    mutation=filter_deletion_insertion(mutation)
    return get_resid_from_mut(mutation)+':'+get_newrestyp_from_mut(mutation)

def get_oldres_from_mut(mutation):
    #
    # Given a mutation (A:0112:ASN:ARG) return the old residue (A:0112:ASN)
    #
    mutation=filter_deletion_insertion(mutation)
    return get_resid_from_mut(mutation)+':'+get_oldrestyp_from_mut(mutation)

def get_oldrestyp_from_mut(mutation):
    #
    # Given a mutation (A:0112:ASN:ARG) return the old residue (ASN)
    #
    import string
    mutation=filter_deletion_insertion(mutation)
    return string.split(mutation,':')[-2]

def get_intresnum_from_res(residue):
    #
    # Given a residue (A:0012:HIS), return the residue number
    #
    import string
    return string.atoi(string.split(residue,':')[1])

def get_resnum_from_res(residue):
    #
    # Given a residue (:0012:HIS), return the residue number (0012)
    #
    import string
    return string.split(residue,':')[1]

def get_resid_from_mut(mutation):
    mutation=filter_deletion_insertion(mutation)
    return get_resid_from_res(mutation)

def get_resid_from_res(residue):
    return '%s:%s' %(residue.split(':')[0],residue.split(':')[1])

def get_resid_from_titgroup(group):
    return get_resid_from_res(group)

def get_titgroup_type_from_titgroup(group):
    """Given a titratable group unique id e.g. (A:0112:CTERM), return the titgroup type (CTERM)"""
    return group.split(':')[-1]

def get_restype_from_titgroup(group):
    """Given a titratable group unique id e.g. (A:0112:CTERM or A:0111:ASP), return the residue type (ASP)"""
    ptype=group.split(':')[-1].upper()
    if ptype=='NTERM' or ptype=='CTERM':
        return None
    else:
        return ptype

def check_mutation_syntax(operations,sequence=None):
    """Check the syntax of an operation string and check that the residue is present in the sequence"""
    single=get_single_operations(operations)
    for op in single:
        resid=get_resid_from_mut(op)
        new=get_newrestyp_from_mut(op)
        old=get_oldrestyp_from_mut(op)
        import Protool
        X=Protool.structureIO()
        if not X.aminoacids.has_key(old) or not X.aminoacids.has_key(new):
            raise Exception('New or old residue is not an amino acid: %s' %op)
        if sequence:
            if not [resid,old] in sequence:
                raise Exception('Original sequence does not contain this residue: %s:%s' %(resid,old))
    return combine_operations(single)

def get_single_operations(operations):
    """Get the single mutations/operations from a combined expression"""
    single=[]
    sp=operations.split('+')
    for s in sp:
        single.append(s.strip())
    return single

def interpret_mutations(operations):
    """Combine a string of operations to a set of PEAT operations.
    This function detects whether we have a classic or a PEAT format
    and cleans the operations afterwards"""
    if operations.find(':')!=-1:
        pass
    else:
        operations=convert_classic_to_PEAT(operations)
    if not operations:
        return False
    return check_mutation_syntax(operations)

def combine_operations(operations):
    """Combine single operations into a single operation string"""
    operations.sort()
    import string
    return string.join(operations,'+')

def convert_classic_to_PEAT(operations):
    """Convert a set of classic mutations to a set of PEAT operations
    The classic operations are in the format: A12G+R45V etc."""
    #
    # Do a quick sanity check
    #
    for junk in ['?','unknown','empty']:
        if operations.lower().find(junk)!=-1:
            return False
    #
    # Deal with the operations
    #
    sp=operations.split('+')
    import Protool, string
    P=Protool.structureIO()
    POP=[]
    for op in sp:
        if op=='wt':
            continue
        old=op[0]
        new=op[-1]
        number=int(op[1:-1])
        try:
            POP.append('%s:%s:%s:%s' %('',string.zfill(number,P.length_of_residue_numbers),P.one_to_three[old],P.one_to_three[new]))
        except KeyError:
            return False
    return string.join(POP,'+')

def convert_PEAT_to_classic(operations):
    operations=get_single_operations(operations)
    newops=[]
    import string, Protool
    PI=Protool.structure()
    for operation in operations:
        resid=get_resid_from_mut(operation)
        old=PI.three_to_one[get_oldrestyp_from_mut(operation)]
        new=PI.three_to_one[get_newrestyp_from_mut(operation)]
        resnum=get_resnum_from_mut(operation)
        resnum=string.lstrip(resnum,'0')
        classic='%s%s%s' %(old,resnum.strip(),new)
        newops.append(classic)

    classicops=string.join(newops,'+')
    return classicops

#
# -----
#

def get_mutations(parent,child,ignoreCterm=False):
    """Given two Protool instances, extract sequences, align and extract mutations
    PDB residue identifier information"""
    import PEATDB.sequence_alignment as SA
    seq1=parent.PirSeq()
    seq2=child.PirSeq()
    #
    ALIGN=SA.NW(seq1,seq2,gap=15.0)
    al1,al2,alres1,alres2=ALIGN.Align(verbose=False)
    if ALIGN.sequence_identity<95.0:
        print 'Sequence identity too low: %f' %ALIGN.sequence_identity
        return False
    #
    operations=findSequenceDifferences(al2,al1,parent.sequence,ignoreCterm)
    return operations

def findSequenceDifferences(child_sequence, parent_sequence,full_parent_sequence,ignoreCterm=False):
    """
    # Find all amino acid differences between child_sequence and parent_sequence
    Child sequence and parent sequence must be aligned and in 1-letter format:
    child_sequence, parent_sequence: AAADEFFG
    full parent sequence is a Protool.sequence object
    """
    #
    # Loop over the sequences - changes are record from parent -> child
    #
    import string
    operations=[]
    import Protool
    PI=Protool.structureIO()
    #
    Cterm_add=0
    insert_num=0
    full_parent_count=0
    #for count in range(len(record_sequence)):
    #    parent_aa=parent_sequence[count]
    #    child_aa=record_sequence[count]
    for parent_aa,child_aa in zip(parent_sequence,child_sequence):
        #
        #print parent_aa,child_aa
        if parent_aa!=child_aa:
            
            # Find the PDB file residue number
            if full_parent_count>=len(full_parent_sequence):
                # If we have an insertion at the Cterm
                aa_identifier=full_parent_sequence[-1][0]
                if ignoreCterm:
                    continue
            else:
                aa_identifier=full_parent_sequence[full_parent_count][0]
            #if aa_identifier[-1]==':':
            #    aa_identifier=aa_identifier[:-1]
            #
            # Convert to 3letter format
            #
            if parent_aa!='-':
                full_parent_count=full_parent_count+1
                parent_aa=PI.one_to_three[parent_aa]
            if child_aa!='-':
                child_aa=PI.one_to_three[child_aa]
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

#
# -----
#


def shortenOperations(operations):
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
