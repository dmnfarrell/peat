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

from errors import *

import math, numpy
from numpy import *
from numpy.linalg import *
linear_least_squares=lstsq
inverse=inv
Float=float
eigenvectors=eig


def length(vector):
    # This function returns the length of vector
    import math
    sum=0.0
    for value in vector:
        sum=sum+math.pow(value,2)
    return math.sqrt(sum)

def normalise (vector):
    # normalise vector
    newvector=[]
    l=length(vector)
    if l!=0:
        for item in vector:
            newvector.append(float(item)/l)
    else:
        newvector=vector
    return array(newvector)



class geometry:
    """Geometry class"""

    def distance(self,atom1,atom2):
        """Return the distance between two atoms"""
        return self.dist(atom1,atom2)

    #
    # ----
    #

    def dist(self,atom1,atom2):
        """
        # Calculates the distance between two atoms
        """
        vector=self.GetPosition(atom1)-self.GetPosition(atom2)
        return length(vector)

    #
    # -----
    #
    
    def CA_positions(self):
        """Return a numpy array of all CA distances"""
        CAs=[]
        atoms=self.atoms.keys()
        atoms.sort()
        for atom in atoms:
            if atom[-3:]==':CA':
                CAs.append(self.GetPosition(atom))
        return numpy.array(CAs)

    #
    # ----
    #
    
    def CA_distance_matrix(self):
        """Calculate the CA distance matrix"""
        CAs=[]
        atoms=self.atoms.keys()
        atoms.sort()
        for atom in atoms:
            if atom[-3:]==':CA':
                CAs.append(atom)
        #
        # Calculate the matrix
        #
        CAs.sort()
        matrix=numpy.zeros([len(CAs),len(CAs)])
        matrix=[]
        for atom1 in CAs:
            row=[]
            for atom2 in CAs:
                row.append(self.dist(atom1,atom2))
            matrix.append(row)
        return numpy.array(matrix,dtype=float)

    #
    # ---
    #

    def superpose(self,reference_coords,fit_coords,return_pdbs=False):
        """Superpose the fit_coords on the reference_coords"""
        import quatfit
        self.refcenter,self.fitcenter,self.rotation=quatfit.qfit(len(reference_coords),reference_coords,fit_coords)
        self.newcoords = quatfit.qtransform(len(reference_coords), fit_coords, self.refcenter, self.fitcenter, self.rotation)
        #
        # Calc rmsd
        #
        rmsd=0.0
        self.rmsd_map=[]
        for i in range(len(self.newcoords)):
            disp = self.newcoords[i] - reference_coords[i]
            self.rmsd_map.append(length(disp))
            dist_squared = dot(disp,disp)
            rmsd = rmsd + dist_squared
        rmsd=rmsd/float(len(self.newcoords))
        import math
        rmsd=math.sqrt(rmsd)
        #
        # If asked then return pdblines for the org and aligned atoms (as CA atoms)
        #
        if return_pdbs:
            pdblines=[]
            for coords in [reference_coords,self.newcoords]:
                theselines=[]
                count=1
                for coord in coords:
                    line='ATOM  %5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f %5s %5s\n' %(count,' CA ','ALA','A',count,float(coord[0]),float(coord[1]),float(coord[2]),'1.0','10.0')
                    count=count+1
                    theselines.append(line)
                pdblines.append(theselines)
            return rmsd,pdblines
        return rmsd

    #
    # ----
    #
    

    def translate(self,ref,fit):
        #
        # Return the vector that translates the center of
        # gravity of ref into the center of gravity of ref
        #
        # Center of geometry of reference
        cogref = array([ 0., 0., 0.])
        for coord in ref:
            cogref = cogref + coord
        cogref = cogref/len(ref)

        # Center of geometry of fit
        cogfit = array([ 0., 0., 0.])
        for coord in fit:
            cogfit = cogfit + coord
        cogfit = cogfit/len(fit)

        # Find the vector to be added to fit to translate
        # COG to ref's COG
        tcogfit = cogref - cogfit
        return tcogfit

 
    #
    # ----
    #
    
    def center(self,ref,fit):
      """
      # Center both coordinate sets at 0,0,0
      """
      #
      # First check that we have some coordinates
      #
      if len(ref)==0 or len(fit)==0:
        raise ProtoolError, "Either reference of fit coordinate set contained no coordinates"
      #
      # Center of coordinates of reference
      #
      cogref = array([0.,0.,0.])
      for coord in ref:
        cogref = cogref + coord
      cogref = cogref/len(ref)
      #
      # Produce new coordinates
      #
      center_ref=[]
      vect_ref=array([0.,0.,0.])-cogref
      for coord in ref:
        center_ref.append(coord+vect_ref)


      # Center of coordinates of fit
      cogfit = array([ 0., 0., 0.])
      for coord in fit:
        cogfit = cogfit + coord
      cogfit = cogfit/len(fit)
      #
      # Produce new coordinates
      #
      center_fit=[]
      vect_fit=array([0.,0.,0.])-cogfit
      for coord in fit:
        center_fit.append(coord+vect_fit)

      return center_ref,center_fit,vect_ref,vect_fit
      
    #
    # ----
    #
    
    def get_center_of_coords(self):
        """Get the center of coordinates for this protein"""
        coor_sum=numpy.array([0.0,0.0,0.0])
        xmax=-99999 
        ymax=-99999
        zmax=-99999
        xmin=9999
        ymin=99999
        zmin=99999
        for atom in self.atoms:
            coords=self.GetPosition(atom)
            xmax=max(xmax,coords[0])
            xmin=min(xmin,coords[0])
            ymax=max(ymax,coords[1])
            ymin=min(ymin,coords[1])
            zmax=max(zmax,coords[2])
            zmin=min(zmin,coords[2])
            coor_sum=coor_sum+coords
        coor_sum=coor_sum/float(len(self.atoms))
        return coor_sum,array([abs(xmax-xmin),abs(ymax-ymin),abs(zmax-zmin)])


#
# Validate the superpos function every time this module is imported
#


import copy
rv=array([[1,1,1],[0,-1,1],[1,-2,1],[1,-3,1],[-1,-3,2]])
fv=array([[1,1,0],[-1,2,0],[-2,1,0],[-3,1,0],[-3,3,1]])

import quatfit
refcenter,fitcenter,rotation=quatfit.qfit(len(rv),rv,fv)
newcoords = quatfit.qtransform(len(rv), fv, refcenter, fitcenter, rotation)
#
# Calc rmsd
#
rmsd=0.0
for i in range(len(newcoords)):
    disp = newcoords[i] - rv[i]
    dist = dot(disp,disp)
    rmsd = rmsd + dist
if rmsd>0.0000001:
    print 'Superpositioning test failed!!! in geometry.py' 
    print 'RMSD',rmsd
    print newcoords
    import sys
    sys.exit(1)
else:
    None
    #print 'superpositioning algorithm validated'
