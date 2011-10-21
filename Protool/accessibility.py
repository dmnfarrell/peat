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

from numpy import *

class access:

    def __init__(self,protein):
        """Store the protool instance and calculate boxes"""
        self.PI=protein
        self.cutoff=6.0
        #
        # Count number of atoms within 6A of side chain atoms. Only heavy atoms are counted
        #
        self.residues=self.PI.residues.keys()
        self.residues.sort()
        return

    #
    # ----
    #

    def get_access(self,residue):
        """Get the number of atoms close"""
        acc=0.0
        count=0
        for res2 in self.PI.residues.keys():
            if res2==residue:
                continue
            if self.PI.dist(residue+':CA',res2+':CA')<25:
                counted={}
                for atom1 in self.PI.residues[residue]:
                    if self.PI.is_backbone(atom1) and atom1.split(':')[-1]!='CA':
                        continue
                    for atom2 in self.PI.residues[res2]:
                        if not counted.has_key(atom2):
                            if self.PI.dist(atom1,atom2)<self.cutoff:
                                acc=acc+1.0
                                counted[atom2]=1
        #
        #
        numatoms=float(len(self.PI.residues[residue])-4) #Subtract backbone atoms
        #
        # Correct for compactness of aromatics
        #
        if self.PI.resname(residue) in ['PHE','TYR','HIS']:
            numatoms=numatoms-3
        if self.PI.resname(residue) in ['TRP']:
            numatoms=numatoms-6
        if numatoms>0.0:
            acc=acc/numatoms
        else:
            acc=0
        return acc
                        
                
    
