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

import Protool
X=Protool.structureIO()
print 'Reading PDB file'
X.readpdb('1arq.pdb')

#
# X.residues.keys() gives you residue numbers
#
residue_numbers=X.residues.keys()
#
# X.resname gives you the residue name
#
print 'Residue ten in each chain is:'
print 'Chain A:',X.resname('A:0010')
print 'Chain B:', X.resname('B:0010')
#
# The mutate module will mutate anything for you
#
print 
import Protool.mutate
MUT=Protool.mutate.Mutate(X)
MUT.mutate('A:0010','TYR')
MUT.mutate('B:0010','ALA')

#
# Finally save the PDB file
#
X.writepdb('mutated_pdb.pdb')
#
# ----
# New instance to check
#
X2=Protool.structureIO()
print 'Reading mutated PDB file'
X2.readpdb('mutated_pdb.pdb')
print 'Residue ten in each chain is:'
print 'Chain A:',X2.resname('A:0010')
print 'Chain B:', X2.resname('B:0010')
