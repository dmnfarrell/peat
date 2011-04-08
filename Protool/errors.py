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

# Updated by Michael Johnston: 29 May 2009
##
# This file contains error messages and a few flags that control
# the overall behaviour of Protool
#
#



# Invalid atoms - f.ex. side chain atoms for a phi/psi evaluation
#
class ProtoolError(Exception):    
    def __init__(self,message='Protool error'):
        Exception.__init__(self, message)

class TypeCheckError(ProtoolError):
	pass

class ParseLineError(ProtoolError):
    def __init__(self,message='Error parsing line'):
        ProtoolError.__init__(self, message)

class ParseAtomError(ProtoolError):
    def __init__(self,message='Error parsing atom'):
        ProtoolError.__init__(self, message)
        
class InvalidAtomError(ProtoolError):
    def __init__(self,message='Invalid atoms for this operation'):
        ProtoolError.__init__(self, message)

class HydrogenInTorsionError(ProtoolError):
    def __init__(self,message='Hydrogen atom in torsion angle'):
        ProtoolError.__init__(self, message)
                        
class NotAnAminoAcidError(ProtoolError):
    def __init__(self,residue=None):
        if residue is None:
            message = 'Non-amino acid residue'
        else:
            message='%s is not an amino acid' %residue
        ProtoolError.__init__(self, message)

# 'Not Found' Errors..
class AtomNotFoundError(ProtoolError):

    def __init__(self,uniqueid=None):
        if uniqueid is None:
            message = 'Atom not found - id not supplied'
        else:
            message = 'Atom not found: %s' % uniqueid
        ProtoolError.__init__(self, message)
        
class ResidueNotFoundError(ProtoolError):

    def __init__(self,uniqueid=None):
        if uniqueid is None:
            message = 'Residue not found - id not supplied'
        else:
            message = 'Residue not found: %s' % uniqueid

        print message
        ProtoolError.__init__(self, message)

# N- and C-terminals
class Nterm(ProtoolError):
    def __init__(self,message='N-terminal Residue'):
        ProtoolError.__init__(self, message)
        
class Cterm(ProtoolError):
    def __init__(self,message='C-terminal Residue'):
        ProtoolError.__init__(self, message)
        
# Incomplete Errors
#

class IncompletePositionError(ProtoolError):
    def __init__(self,message='Not all coordinates were found this atom'):
        ProtoolError.__init__(self, message)

#
# I/O errors
EndOfFileError='End Of File'
class FileNotFoundError(ProtoolError):
    def __init__(self,filename):
        ProtoolError.__init__(self,'File Not Found: %s' %filename)


#
# Flags
#
silent=1
