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

    def build_all_HAs(self):
        """Build all HAs"""
        residues=self.PI.residues.keys()
        residues.sort()
        for residue in residues:
            if self.PI.is_residue(residue):
                self.build_HA(residue)
        self.PI.Update()
        return

    #
    # ----
    #

    def build_HA(self,residue):
        """Build all or some hydrogens for this residue"""
        coords={'N':numpy.array([2.967,4.770,13.995]),
                'CA':numpy.array([2.755,5.653,12.837]),
                'C':numpy.array([3.345,7.017,13.183]),
                'HA':numpy.array([3.210,5.226,12.056])}
        buildname='HA'
        return self.build_hydrogen(residue,buildname,coords)

    def build_HN(self,residue):
        """Build the HN"""
        if self.PI.Nterminal(residue):
            return None
        if self.PI.resname(residue)=='PRO':
            return None
        #
        # Definition of position
        #
        coords={'C-1':numpy.array([8.8,4.4,14.1]),
                'N':numpy.array([8.9,3.1,14.4]),
                'CA':numpy.array([8.9,2.1,13.3]),
                'H':numpy.array([9.0,2.8,15.3])}
        buildname='H'
        return self.build_hydrogen(residue,buildname,coords)
        

    def build_hydrogen(self,residue,buildH,coords):
        """Build a hydrogen in the structure when given coordinates, the Hname and the residue"""
        fit_coords=[]
        findatoms=sorted(coords.keys())
        findatoms.remove(buildH)
        ref_coords=[]
        for atom in findatoms:
            fit_coords.append(coords[atom])
            #
            # Now find the corresponding atom in the structure
            #
            struct_res=residue
            if atom[-2:]=='-1':
                try:
                    struct_res=self.PI.PreviousResidue(residue)
                    atom=atom[:-2]
                except ResidueNotFoundError:
                    print 'prev not found'
                    return None
            ref_coords.append(self.PI.GetPosition('%s:%s' %(struct_res,atom)))
        #
        # Get rotation and translation for bb -> bb
        #
        import quatfit
        refcenter,fitcenter,rotation=quatfit.qfit(len(ref_coords),ref_coords,fit_coords)
        #
        # Get the atoms that should be transformed
        #
        trans_coords=[coords[buildH],coords[buildH]]
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
        uniqueid='%s:%s:%s' %(chainid,resnumber,buildH)
        count=0
        if not self.PI.atoms.has_key(uniqueid):
            self.PI.add_atom(uniqueid,atomnumber=atomnumber,
                             atomname=buildH,chainid=chainid,residuename=restype,
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
