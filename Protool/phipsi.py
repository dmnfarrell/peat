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
from geometry import *


try:
    import numpy as Numeric
except:
    import Numeric
    
import math

#Conversion factor between degrees and radians
dr=360.0/(2.0*math.pi)


def angle(vector1,vector2,refvector):
    """This function returns the angle between vector1 and vector2
     By using refvector, the rotation is determined to be positive (counter-clockwise)
    # or negative (clockwise)"""
    dot=Numeric.dot(vector1,vector2)
    dot=dot/(length(vector1)*length(vector2))
    vectorprod=length(cross(vector1,vector2))/(length(vector1)*length(vector2))
    ang=dr*angleFromSineAndCosine(vectorprod,dot)
    if Numeric.dot(cross(vector1,vector2),refvector)<0.0:
        ang=-ang
    return ang

def cross (vector1,vector2):
    # Evaluates the crossproduct between the 3D vectors: vector1 and vector2
    vectorproduct=Numeric.array([vector1[1]*vector2[2]-vector1[2]*vector2[1],vector1[2]*vector2[0]-vector1[0]*vector2[2],vector1[0]*vector2[1]-vector1[1]*vector2[0]])
    return vectorproduct

    

def angleFromSineAndCosine(sine, cosine):
    # Copied from Scientific Python
    sine = min(1., max(-1., sine))
    angle = Numeric.arcsin(abs(sine))
    if sine*cosine < 0.:
	angle = -angle
    if cosine < 0.:
	angle = Numeric.pi + angle
    if angle < 0.:
	angle = 2.*Numeric.pi + angle
    return angle


def GetPhiPsi(molecule,resid):
    #
    # This function calculates phi, psi and omega for a residue
    #
    #
    # Phi
    #
    try:
        phi=GetPhi(molecule,resid)
    except (NotAnAminoAcidError,AtomNotFoundError,InvalidAtomError):
        phi=None
    #
    # Psi
    #
    try:
        psi=GetPsi(molecule,resid)
    except (NotAnAminoAcidError,AtomNotFoundError,InvalidAtomError):
        psi=None
    #
    # Omega
    #
    try:
        omega=GetOmega(molecule,resid)
    except (NotAnAminoAcidError,AtomNotFoundError,InvalidAtomError):
        omega=None
    return phi, psi, omega
    


def GetPhi(molecule,resid):
    #
    # Calculate Phi for resid
    # Is resid OK?
    #
    if not molecule.isaa(resid):
        raise NotAnAminoAcidError,resid
    # Get the previous residue
    try:
        prevres=molecule.PreviousResidue(resid)
    except Nterm:
        return None
    # is prevres OK?
    if not molecule.isaa(prevres):
        raise NotAnAminoAcidError,prevres
    C1=prevres+':C'
    N2=resid+':N'
    CA2=resid+':CA'
    C2=resid+':C'
    return GetDihedral(molecule,C1,N2,CA2,C2)

def GetPsi(molecule,resid):
    #
    # Calculate Psi for resid
    # is resid OK?
    #
    if not molecule.isaa(resid):
        raise NotAnAminoAcidError, resid
    # Get the next residue
    try:
        nextres=molecule.NextResidue(resid)
    except Cterm:
        return None
    # Is nextres OK?
    if not molecule.isaa(nextres):
        raise NotAnAminoAcidError, nextres
    N1=resid+':N'
    CA1=resid+':CA'
    C1=resid+':C'
    N2=nextres+':N'
    return GetDihedral(molecule,N1,CA1,C1,N2)

def GetDefAngle(molecule,resid,defangle_span=1):
    #
    # Calculate the deformation angle for resid
    # is resid OK?
    #
    if not molecule.isaa(resid):
        raise NotAnAminoAcidError, resid
    #
    # Get the next residue, defangle_span times
    #
    resid2=resid
    nextres2=None

    for x in range(defangle_span):
        try:
            nextres2=molecule.NextResidue(resid2,None)
        except:
            return None
        resid2=nextres2
    #print 'Nextres:',nextres
    #
    # Is nextres OK?
    #
    if not nextres2:
        return
    if not molecule.isaa(nextres2):
        raise NotAnAminoAcidError, nextres2
    #
    # Get the previous residue
    #
    prevres=resid
    for x in range(defangle_span):
        try:
            prevres=molecule.PreviousResidue(prevres,None)
        except Nterm:
            return None
    #
    # is prevres OK?
    #
    if not molecule.isaa(prevres):
        raise NotAnAminoAcidError,prevres
    #
    # Get the atoms
    #
    #N1=prevres+':CA'
    CA0=prevres+':CA'
    CA1=resid+':CA'
    CA2=nextres2+':CA'
    #
    # Get the two vectors
    #
    v1=molecule.GetPosition(CA1)-molecule.GetPosition(CA0)
    v2=molecule.GetPosition(CA2)-molecule.GetPosition(CA1)
    return angle(v1,v2,v1)

def GetOmega(molecule,resid):
    #
    # Calculate Omega for resid
    # is resid OK?
    #
    if not molecule.isaa(resid):
        raise NotAnAminoAcidError, resid
    #
    # Get the next residue
    #
    try:
        nextres=molecule.NextResidue(resid)
    except Cterm:
        return None
    #
    # Is nextres OK?
    #
    if not molecule.isaa(nextres):
        raise NotAnAminoAcidError, nextres
    CA1=resid+':CA'
    C1=resid+':C'
    N2=nextres+':N'
    CA2=nextres+':CA'
    return GetDihedral(molecule,CA1,C1,N2,CA2)

def GetDihedral(molecule,atom1,atom2,atom3,atom4):
    #
    # Calculates the dihedral of the bond from atom2 to atom3
    # by measuring the angle between atom1 and atom4.
    #
    # Check that the atom is present
    #
    if not molecule.ispresent(atom1):
        raise AtomNotFoundError,atom1
    if not molecule.ispresent(atom2):
        raise AtomNotFoundError,atom2
    if not molecule.ispresent(atom3):
        raise AtomNotFoundError,atom3
    if not molecule.ispresent(atom4):
        raise AtomNotFoundError,atom4
    #
    # Check that they are bound correctly
    #
    #if not molecule.bound(atom1,atom2) or not molecule.bound(atom2,atom3) or not molecule.bound(atom3,atom4):
        #print molecule.dist(atom1,atom2)
        #print molecule.dist(atom2,atom3)
        #print molecule.dist(atom3,atom4)
    #    raise InvalidAtomError, [atom1,atom2,atom3,atom4]
    if molecule.bound(atom1,atom3) or molecule.bound(atom1,atom4) or molecule.bound(atom2,atom4):
        raise InvalidAtomError, ['Atoms are incorrectly connected',atom1,atom2,atom3,atom4] 
    #
    # Get the positions
    #
    atom1=molecule.GetPosition(atom1)
    atom2=molecule.GetPosition(atom2)
    atom3=molecule.GetPosition(atom3)
    atom4=molecule.GetPosition(atom4)
    vector1=atom2-atom1
    vector2=atom3-atom2
    vector3=atom4-atom3
    cross1=cross(vector1,vector2)
    cross2=cross(vector2,vector3)
    refvector=vector3
    return angle(cross1,cross2,refvector)

class Plot:
    #
    # This class is for plotting phi-psi pairs with Gnuplot
    #
    def __init__(self):
        self.array=[]
        return

    def plot_phi_psi(self,phi,psi):
        self.array.append([phi,psi])
        return

    def show(self):
        import Gnuplot
        g=Gnuplot.Gnuplot()
        g('set pointsize 1')
        #data=Gnuplot.Data(self.array,title='Ramachandran',with='points 3')
        g.title('Ramachandran Plot')
        g.xlabel('Phi')
        g.ylabel('Psi')
        g('set yrange [-180:180]')
        g('set xrange [-180:180]')
        #g.plot(data,Gnuplot.Data([[180,0],[-180,0]],with='lines 1'),Gnuplot.Data([[0,-180],[0,180]],with='lines 1'))
        raw=raw_input('Press Enter to continue')
        return
    

def example():
    # This gives an angle of -135 degrees
    vector2=Numeric.array([-1,-1,0])
    vector1=Numeric.array([1,0,0])
    refvector=Numeric.array([0,0,1])
    print angle(vector1,vector2,refvector)
    
if __name__=="__main__":
    example()
    
