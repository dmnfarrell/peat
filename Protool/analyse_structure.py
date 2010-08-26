#!/usr/bin/env python
#
# # Protool - Python class for manipulating protein structuress
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
import string

class find_pattern:

    def acid_amide_Gly(self):
        residues=self.residues.keys()
        residues.sort()
        prevtype=''
        sites=[]
        for residue in residues:
            resname=self.resname(residue)
            if prevtype!='':
                if string.strip(string.upper(resname))=='GLY':
                    if string.strip(string.upper(prevtype))=='ASP' or string.strip(string.upper(prevtype))=='ASN':
                        print 'Asp/Asn-Gly site found:',residue
                        sites.append(residue)
            prevtype=resname
        return sites

    #
    # =====================================================
    #

    def cis_Pro(self):
        import phipsi
        residues=self.residues.keys()
        residues.sort()
        for residue in residues:
            #
            # Get the next residue
            #
            skip=None
            try:
                nextres=self.NextResidue(residue)
            except Cterm:
                skip=1
            #
            # Is nextres OK?
            #
            if not self.isaa(nextres):
                skip=1
            if not skip:
                if string.strip(string.upper(self.resname(nextres)))=='PRO':
                    if phipsi.GetPhiPsi(self,residue)[2]:
                        if abs(phipsi.GetPhiPsi(self,residue)[2])<30.0:
                            print 'Cis proline peptide bond: ',residue,nextres
                    #print phipsi.GetPhiPsi(self,residue)[2],residue,self.resname(residue)
        return

    #
    # ==========================================================
    #

    def find_sequence(self,sequence):
        #
        # Find a sequence of amino acids in the present PDB
        # sequence holds a one-letter sequence
        #
        import string
        sequence=string.upper(sequence)
        residues=self.residues.keys()
        residues.sort()
        results=[]
        for residue in residues:
            if self.three_to_one.has_key(self.resname(residue)):
                resnam=self.three_to_one[self.resname(residue)]
            else:
                print 'Residue not found: %s' %residue
                #resnam='?'
            if resnam==sequence[0]:
                found=1
                lastres=residue
                stretch=[residue]
                for testres in sequence[1:]:
                    try:
                        nextres=self.NextResidue(lastres)
                    except Cterm:
                        found=None
                        break
                    if self.three_to_one.has_key(self.resname(nextres)):
                        resnam=self.three_to_one[self.resname(nextres)]
                    lastres=self.NextResidue(lastres)
                    if resnam!=testres:
                        found=None
                        break
                    stretch.append(lastres)
                if found:
                    #
                    # Add this stretch to the results
                    #
                    results.append(stretch)
        return results

#
# ---------------------------------------
#
                    
class distances:

    def calc_distance_matrix(self):
        CA={}
        epsilon=0.01
        count=1
        for atom in self.atoms.keys():
            if self.atname(atom)=="CA":
                CA[atom]=count
                count=count+1
        dists={}
        CAs=CA.keys()
        for pos1 in range(len(CAs)):
            for atom2 in CAs[pos1+1:]:
                atom1=CAs[pos1]
                dist=self.dist(atom1,atom2)
                if dist<920.0:
                    ID1=CA[atom1]
                    ID2=CA[atom2]
                    if not dists.has_key(ID1):
                        dists[ID1]={}
                    dists[ID1][ID2]=[self.dist(atom1,atom2)*(1-epsilon),self.dist(atom1,atom2)*(1+epsilon)]
        #
        # Write file
        #
        done={}
        fd=open('dg.data','w')
        distsum=0.0
        count=0
        max=-9.9
        for ID1 in dists.keys():
            for ID2 in dists[ID1].keys():
                fd.write('%4d\t%4d\t%7.3f\t%7.3f\n' %(ID1,ID2,dists[ID1][ID2][0],dists[ID1][ID2][1]))
                distsum=distsum+dists[ID1][ID2][1]
                count=count+1
                if dists[ID1][ID2][1]>max:
                    max=dists[ID1][ID2][1]
        fd.close()
        print 'Average CA-CA distance: %5.3f. Max CA-CA dist: %5.3f ' %(distsum/float(count),max)
        return
                         
    #
    # ----
    #

    def get_min_distance(self,residue1,residue2):
        """Get the minimum distance between residue1 and residue2"""
        min_dist=99999.9
        for atom1 in self.residues[residue1]:
            for atom2 in self.residues[residue2]:
                dist=self.dist(atom1,atom2)
                if dist<min_dist:
                    min_dist=dist
        return min_dist
