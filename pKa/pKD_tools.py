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
    elif mutation.find('insert:')!=-1:
        raise Exception('We cannot model insertions')
    elif mutation.find('+')!=-1:
        raise Exception('Something went wrong. I should have gotten a single mutation and got %s' %mutation)
    return mutation

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
    # Given a mutation (A:0112:ASN:ARG) return the new residue (A:0112:ASN)
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

def convert_classic_to_PEAT(operations):
    """Convert a set of classic mutations to a set of PEAT operations
    The classic operations are in the format: A12G+R45V etc."""
    sp=operations.split('+')
    import Protool, string
    P=Protool.structureIO()
    POP=[]
    for op in sp:
        old=op[0]
        new=op[-1]
        number=int(op[1:-1])
        POP.append('%s:%s:%s:%s' %('',string.zfill(number,P.length_of_residue_numbers),P.one_to_three[old],P.one_to_three[new]))
    return string.join(POP,'+')
        
