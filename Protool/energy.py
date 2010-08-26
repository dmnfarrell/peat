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


class energy:

    def VdWr(self,atom):
        """Get the Van der Waals radius of an atom"""
        if self.is_carbon(atom):
            return self.vdwr['C']
        elif self.is_nitrogen(atom):
            return self.vdwr['N']
        elif self.is_oxygen(atom):
            return self.vdwr['O']
        elif self.is_phosphor(atom):
            return self.vdwr['P']
        elif self.is_sulphur(atom):
            return self.vdwr['S']
        else:
            return None

    #
    # ----
    #

    def Van_der_Waals_sum(self):
        """Get the sum of all Van der Waals radii"""
        sum=0.0
        print " vdwr ",self.vdwr
        for atom in self.atoms.keys():
            print "VW of ",atom," = ",self.VdWr(atom)
            sum_i = self.VdWr(atom)
            if sum_i !=None:
                sum=sum+sum_i
        return sum

    #
    # ----
    #

    def get_hbonds(self):
        """Return a list of all hydrogen bonds"""
        hbonds=[]
        for atom1 in self.atoms.keys():
            for atom2 in self.get_neighbour_atoms(atom1,box='small'):
                if atom2==atom1:
                    continue
                if self.is_hydrogen_bond(atom1,atom2):
                    hbond=[atom1,atom2]
                    if not hbond in hbonds:
                        hbonds.append(hbond)
        return hbonds

    #
    # ----
    #

    def construct_distance_cubes(self):
        """Divide atoms into cubes"""
        cube_size=6


        return

    #
    # ----
    #

    def find_coordinates_extremes(self):
        """Find the extremes of coordinates in all three directions"""
        minx=99999.9
        maxx=-99999
        miny=99999
        maxy=-999999
        minz=999999
        maxz=-999999
        for atom in self.atoms.keys():
            X=self.atoms[atom]['X']
            Y=self.atoms[atom]['Y']
            Z=self.atoms[atom]['Z']
            minx=min(minx,X)
            maxx=max(maxx,X)
            miny=min(miny,Y)
            maxy=max(maxy,Y)
            minz=min(minz,Z)
            maxz=max(maxz,Z)
        return minx,maxx,miny,maxy,minz,maxz
    
