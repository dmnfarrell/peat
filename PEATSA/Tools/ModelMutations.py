#! /bin/env python
#
# Protein Engineering Analysis Tool Structure Analysis (PEATSA)
# Copyright (C) 2010 Michael Johnston & Jens Erik Nielsen
#
# Author: Michael Johnston
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
'''Example script showing how to model a set of mutations on a pdb

Replace filename with an actual pdb'''
import sys, Protool
from PEATSA import Core

#Read in the pdb
filename = sys.argv[1]
pdb = Protool.structureIO()
pdb.readpdb(filename)

#Create the set of mutations to perform
mutationSet = Core.Data.MutationSet(name='MyClone1', code="A52H")

print '\n%s' % mutationSet
#To see the mutations relative to a sequence pass it a pdb
print 'Mutation I will perform is %s\n' % mutationSet.codeString(pdb)

#The MutantCollection class creates a set of mutants given a list of mutationSets
#You can 'load' any previously created set into another script aswell
#It can also be passed as an argument to ProteinDesignTool on the command line.
collection = Core.Data.MutantCollection(pdbFile=filename, name='ACollection', maxOverlap=1.4)
collection.createMutants(mutationList=[mutationSet]) 

print collection.mutantFiles()
