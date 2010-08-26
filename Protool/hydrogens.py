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

try:
    import numpy
except:
    import Numeric as numpy

import Protool.mutate
#
# Import the error definitions
#
from errors import *

class Hydrogens(Protool.mutate.Mutate):

    def __init__(self,Protool):
        """Init the hydrogen building function"""
        self.PI=Protool
        #
        # Read definition files
        #
        self.read_aa_defs()
        return

    #
    # ----
    #

    def repair_all_atoms(self):
        """Build all Hs by simple translation from std. coordinates"""
        residues=self.PI.residues.keys()
        for residue in residues:
            if self.aadefs.has_key(self.PI.resname(residue)):
                #
                # Building hydrogens
                #
                for atom,coord,number in self.aadefs[self.PI.resname(residue)]['atoms']:
                    atom_name='%s:%s' %(residue,atom)
                    if not self.PI.atoms.has_key(atom_name):
                        self.build_atom(atom,residue)
                    else:
                        print 'Atom is present',atom_name
            else:
                print 'Unknown residue',residue
        return

    #
    # ----
    #

    def build_all_HNs(self):
        """Build all HNs"""
        residues=self.PI.residues.keys()
        residues.sort()
        for residue in residues:
            if self.PI.is_residue(residue):
                self.build_HN(residue)
        self.PI.Update()
        return

    #
    # ----
    #

    def build_HN(self,residue):
        """Build the HN"""
        if self.PI.Nterminal(residue):
            return None
        if self.PI.resname(residue)=='PRO':
            return None
        #
        # Definition of position
        #
        coords={'C':numpy.array([8.8,4.4,14.1]),
                'N':numpy.array([8.9,3.1,14.4]),
                'CA':numpy.array([8.9,2.1,13.3]),
                'H':numpy.array([9.0,2.8,15.3])}
        fit_coords=[]
        for atom in ['C','N','CA']:
            fit_coords.append(coords[atom])
        #
        # Now find the corresponding atoms in the structure
        #
        try:
            prevres=self.PI.PreviousResidue(residue)
            C_name='%s:C' %prevres
        except ResidueNotFoundError:
            print 'prev not found'
            return None
        ref_coords=[]
        ref_coords.append(self.PI.GetPosition(C_name))
        ref_coords.append(self.PI.GetPosition('%s:N' %residue))
        ref_coords.append(self.PI.GetPosition('%s:CA' %residue))
        #
        # Get rotation and translation for bb -> bb
        #
        import quatfit
        refcenter,fitcenter,rotation=quatfit.qfit(len(ref_coords),ref_coords,fit_coords)
        #
        # Get the atoms that should be transformed
        #
        trans_atoms=['H','H2']
        trans_coords=[coords['H'],coords['H']]
        #
        # Apply rotation and translation to trans_coords
        #
        newcoords = quatfit.qtransform(len(trans_coords), trans_coords, refcenter, fitcenter, rotation)
        #
        # Add the atom
        #
        struct_atom=self.PI.atoms['%s:CA' %residue]
        resnumber=struct_atom['RESNUM']
        restype=struct_atom['RESNAME']
        chainid=struct_atom['CHAINID']
        atomnumber=struct_atom['NUMBER']
        uniqueid='%s:%s:%s' %(chainid,resnumber,'H')
        count=0
        if not self.PI.atoms.has_key(uniqueid):
            self.PI.add_atom(uniqueid,atomnumber=atomnumber,
                             atomname='H',chainid=chainid,residuename=restype,
                             residuenumber=resnumber,
                             xcoord=newcoords[count][0],ycoord=newcoords[count][1],zcoord=newcoords[count][2],update=None)
        #
        # Fix atom numbers
        #
        self.PI.renumber_atoms()
        return
    
    #
    # ----
    #

    def build_atom(self,name,resnumber):
        """Build an atom"""
        #
        # Get the bound_to atoms
        #
        atom_name='%s:%s' %(resnumber,name)
        bound=self.getbound(atom_name)
        while len(bound)<3:
            for atom in bound[:]:
                bound=bound+self.getbound(atom)
                print bound
            #
            # Remove duplicates
            #
            b={}
            for at in bound:
                b[at]=1
            bound=b.keys()
        print 'Found three atoms',bound
        stop
