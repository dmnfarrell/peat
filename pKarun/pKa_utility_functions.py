#!/usr/bin/env python
#
# pKarun - scripts for running pKa calculations
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
# Definitions for all pKa rout

#
# Routines for getting WHAT IF residue identifiers
#

def getWI_resid(line):
    """
    # Take a WHAT IF line of the form:
    # Residue:     1 LYS  (   1  )      (Prp= 1.00)
    #
    # and return a unique resid
    """
    import string
    split=string.split(line)
    residuename=split[2]
    last=line.split('(')[1]
    number=last.split(')')[0]
    if len(string.strip(string.split(last,')')[1]))>0:
        chainID=string.split(string.split(last,')')[1])[0]
    else:
        chainID=''
    number=number.strip().replace(' ','')
    resid='%s:%4s:%3s' %(chainID,string.zfill(number,4),residuename)
    return resid

def getWI_resid2(line,format=None):
    """
     Convert a WHAT IF line of the form:
     43 ASP  (  43  )         xxxxxx
     or
     182 SER  ( 182  ) A   
     and return a unique resid
    

    If format is pdb2pka, then just return the line"""
    import string
    if format=='pdb2pka':
        return string.strip(line.split()[0])
    #
    # WHAT IF format
    #

    split=string.split(line)[:-1]
    line=string.join(split)
    #
    # Line reformatted
    #
    residuename=split[1]
    last=string.split(line,'(')[1]
    number=string.split(last,')')[0].replace(' ','')
    last_part=string.strip(string.split(last,')')[1])
    last_split=string.split(last_part)
    if len(last_split)==1:
        chainID=last_split[0]
    else:
        chainID=''
    number=string.strip(number)
    resid='%s:%4s:%3s' %(chainID,string.zfill(number,4),residuename)
    return resid

def getWI_resid3(line):
    #
    # Take a WHAT IF line of the form:
    # Residue:     1 THR      1    A     Prp= 0.00
    # and return a unique resid
    #
    import string
    split=string.split(line)
    residuename=split[2]
    number=split[3]
    cid=split[4]
    chainID=''
    if len(cid)==1:
        chainID=cid
    else:
        chainID=''
    number=string.strip(number)
    resid='%s:%4s:%3s' %(chainID,string.zfill(number,4),residuename)
    return resid

def getWI_resid4(line):
    #
    # Take a WHAT IF line of the form:
    #   1 THR  (   1  )       N   xxxxxxx
    # and return a unique resid
    #
    import string
    #
    # First trim the string
    #
    line=string.join(string.split(string.strip(line))[:-1])
    #
    # Now the line looks like: 1 THR ( 1 ) N
    #
    split=string.split(line)
    atomname=split[-1]
    residuename=split[1]
    last=string.split(line,'(')[1]
    number=string.strip(string.split(last,')')[0])
    # Chain ID
    laststrip=string.strip(last)
    lastsplit=string.split(last,')')[1]
    split_lastsplit=string.split(lastsplit)
    if len(split_lastsplit)>1:
        chainID=split_lastsplit[0]
    else:
        chainID=''
    resid='%s:%4s:%3s:%s' %(chainID,string.zfill(number,4),residuename,atomname)
    return resid

#
# -----
#

def get_resid(uniqueid):
    import string
    return string.join(string.split(uniqueid,':')[:2],':')

def get_resnum(uniqueid):
    import string
    # Given a uniqueid this function returns the residue number
    return string.split(uniqueid,':')[1]

def get_resname(uniqueid):
    import string
    return string.split(uniqueid,':')[2]

def get_chainid(self,uniqueid):
    import string
    return string.split(uniqueid,':')[0]

def is_terminal(uniqueid):
    #
    # Is this residue a terminal titratable group?
    #
    if string.split(uniqueid,':')[-1]=='TERM':
        return 1
    return None

charged=['ARG','LYS','HIS','ASP','GLU']

def is_normally_titratable(uniqueid):
    """Does this group have a pKa value in the range 2-12"""
    type=uniqueid.split(':')[-1]
    if istitratable(uniqueid) and type!='SER' and type!='THR':
        return 1
    return None

def is_titratable(uniqueid):
    return istitratable(uniqueid)

def istitratable(uniqueid):
    import string
    type=string.split(uniqueid,':')[-1]
    import pKaTool.pKadata
    if pKaTool.pKadata.acidbase.has_key(type):
        return 1
    else:
        return None

def is_charged(uniqueid):
    import string
    type=string.split(uniqueid,':')[-1]
    if type in charged:
        return 1
    return None

def isacid(uniqueid):
    #
    # Is the residue a base or an acid?
    #
    import string
    type=string.split(uniqueid,':')[-1]
    import pKaTool.pKadata
    val=pKaTool.pKadata.acidbase[type.upper()]
    if val==-1:
        return 1
    return None

def acibas(uniqueid):
    #
    # Return -1 (acid) or 1 (base)
    #
    if isacid(uniqueid):
        return -1
    return 1
#
# ----
#

def reformat_name(oldname,Nterm=None,format='WHAT IF'):
    """
    # Reformat the residue names 
    """
    if format=='WHAT IF':
        import copy
        residue=copy.copy(oldname)
        newname=''
        Tflag=None
        if residue[0]=='T':
            Tflag=1
            residue=residue[1:]
        number=residue[:4]
        name=residue[4:7]
        chainID=''
        if len(residue)==8:
            chainID=residue[7]
        newname='%s:%4s:%3s' %(chainID,number,name)
        if Tflag:
            #
            # Nterm is 1 for an Nterm and 0 for a Cterm
            #
            if Nterm==1:
                newname=newname.replace(name,'NTERM')
            else:
                newname=newname.replace(name,'CTERM')
    else:
        import string
        return string.strip(oldname)
    return newname
