#!/usr/bin/env python
#
# pKaTool - analysis of systems of titratable groups
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
    import Numeric
except:
    import numpy


def length(vector):
    # This function returns the length of vector
    import math
    sumt=0.0
    for value in vector:
        sumt=sumt+math.pow(value,2)
    return math.sqrt(sumt)

class MolDyn:

    def __init__(self,grps):
        self.atoms={}
        count=1
        import random
        for grp in grps.keys():
            x=random.randint(-10,10)
            y=random.randint(-10,10)
            z=random.randint(-10,10)
            grp_name=grps[grp]
            pos=Numeric.array([200+x,200+y,200+z])
            vel=Numeric.array([0,0,0])
            acc=Numeric.array([0,0,0])
            self.atoms[grp_name]={'pos':pos,'vel':vel,'acc':acc}
            count=count+1
            #pass
        #
        # For each pair of atoms set the energy func
        #
        self.enes={}
        #for grp1 in grps.keys():
        #    self.enes[grp1]={}
        #    for grp2 in grps.keys():
        #        self.enes[grp1][grp2]=None
        #print self.enes
        return

    def set_eqdists(self,matrix):
        #
        # Set new equilibrium distances
        #
        for grp1 in matrix.keys():
            self.enes[grp1]={}
            for grp2 in matrix[grp1].keys():
                self.enes[grp1][grp2]=spring(50.0*abs(matrix[grp1][grp2]))
        return

    def EM(self,steps,timestep=0.01):
        #
        # Do step steps of EM
        #
        import math
        for step in range(steps):
            diff=0.0
            #
            # Noise
            #
            dists={}
            for grp1 in self.atoms.keys():
                #
                # Apply random noise
                #
                import random
                nlvl=5
                x=random.randint(-nlvl,nlvl)
                y=random.randint(-nlvl,nlvl)
                z=random.randint(-nlvl,nlvl)
                noise=Numeric.array([x,y,z])
                self.atoms[grp1]['pos']=self.atoms[grp1]['pos']+noise/10.0
            #
            # Calculate distances
            #
            for grp1 in self.atoms.keys():
                dists[grp1]={}
                for grp2 in self.atoms.keys():
                    if grp1!=grp2:
                        dists[grp1][grp2]=self.atoms[grp2]['pos']-self.atoms[grp1]['pos']
            #
            # Calculate forces
            #
            for grp1 in self.atoms.keys():
                force=0.0
                for grp2 in self.atoms.keys():
                    if grp1!=grp2:
                        f_contrib,power=self.enes[grp1][grp2].get_force(dists[grp1][grp2])
                        force=force+f_contrib
                        diff=diff+power
                #
                # As set all masses to 1
                #
                self.atoms[grp1]['acc']=force/1.0
            #
            # Update positions
            #
            for grp in self.atoms.keys():
                movement=self.atoms[grp]['acc']
                #print movement
                self.atoms[grp]['pos']=self.atoms[grp]['pos']+movement*math.pow(timestep,2)
        return diff

#
# ------
#

class Dist_geom_EM(MolDyn):
    
    def __init__(self,grps):
        #
        # Just store the groups
        #
        self.atoms={}
        for grp in grps.keys():
            name=grps[grp]
            null = Numeric.array([0,0,0])
            pos=Numeric.array([0,0,0]) # just init to zero for now
            self.atoms[name]={'pos':pos,'vel':null,'acc':null}
        return

    #
    # ---------
    #

    def set_eqdists(self,matrix):
        #
        # Get the distances that we should try to match
        #
        grps_ordered=matrix.keys()
        grps_ordered.sort()
        #
        self.dists=[]
        for grp1 in grps_ordered:
            dist_row=[]
            self.enes[grp1]={}
            for grp2 in grps_ordered:
                self.enes[grp1][grp2]=spring(50.0*abs(matrix[grp1][grp2]))
                dist_row.append(matrix[grp1][grp2])
            self.dists.append(dist_row)
        #
        # That's all ok, now do the distance geometry
        #
        self.dists=Numeric.array(self.dists)
        positions=self.distance_geometry(self.dists)
        count=0
        for grp1 in grps_ordered:
            pos=positions[count]
            self.atoms[grp1]['pos']=pos
        #
        # Now we have a starting state. Do 10000 steps of EM 
        #
        diff_coarse=self.EM(steps=1000,random=None,timestep=0.01)
        #
        # Fine-tune
        #
        diff_fine=self.EM(steps=1000,random=None,timestep=0.001)
        return diff_fine

    #
    # -------
    #

    def distance_geometry(self,dists,error=0.5):
        #
        """From a set of distances produce coordinates that are consistent with them"""
        """The distances can vary +/- error A"""
        #
        #
        # Do triangle smoothing
        #
        atoms=self.dists.keys()
        atoms.sort()
        not_converged=None
        round=0
        while not_converged>0:
            #
            # Maybe we will converge this time?
            #
            not_converged=0
            round=round+1
            #
            # Loop over all triangles
            #
            for atom1 in atoms:
                for atom2 in atoms:
                    for atom3 in atoms:
                        #
                        # Is dist atom1-atom2 <= atom1-atom3 + atom3-atom2?
                        #
                        if dists[atom1][atom2]>(dists[atom1][atom3]+dists[atom3][atom2])+3*error:
                            #
                            # Decrease the atom1,atom2 distance by error
                            #
                            dists[atom1][atom2]=dists[atom1][atom2]-error
                            dists[atom2][atom1]=dists[atom1][atom2]
                            not_converged=not_converged+1
                        #
                        # Is dist atom1-atom2 < dist atom1-atom3 - atom3-atom2 ?
                        #
                        if dists[atom1][atom2] < dists[atom1][atom3]-dists[atom3][atom2] -3*error:
                            #
                            # Increase the atom1,atom2 distance by error
                            #
                            dists[atom1][atom2]=dists[atom1][atom2]+error
                            dists[atom2][atom1]=dists[atom1][atom2]
                            not_converged=not_converged+1
            #
            # Converged?
            #
            print 'Round: %4d, Incorrect triangles: %4d ' %(round,not_converged)
        return

    #
    # -------
    #


class spring:
    #
    # Simple spring force
    #

    def __init__(self,eq_dist):
        #
        # Define the energy function
        #
        self.eq_dist=eq_dist
        return

    def get_force(self,vector):
        power=length(vector)-self.eq_dist
        #.print power,vector,power*vector
        return 2*power*vector,abs(power)
