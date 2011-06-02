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
import  math


#
# Physical constants
#
e=1.60217646E-19
pi=3.14159265
k=1.3806503E-23
A=1E-10
e0=8.8542E-12
T=298.15
Na=6.02E23

class electrostatics:

    def get_titratable_groups(self):
        #
        # Defines self.titrateable_groups
        #
        self.titgroups={'ASP':{'charge':-1,'atoms':['OD1','OD2','CG'],'name':'ASP'},
                        'GLU':{'charge':-1,'atoms':['OE1','OE2','CD'],'name':'GLU'},
                        'ARG':{'charge':1,'atoms':['CZ','NE','NH1','NH2'],'name':'ARG'},
                        'LYS':{'charge':1,'atoms':['NZ'],'name':'LYS'},
                        'HIS':{'charge':1,'atoms':['ND1','NE2'],'name':'HIS'},
                        'CYS':{'charge':-1,'atoms':['SG'],'name':'CYS'},
                        'TYR':{'charge':-1,'atoms':['OH'],'name':'TYR'},
                        'NTERM':{'charge':1,'atoms':['N'],'name':'NTERM'},
                        'CTERM':{'charge':-1,'atoms':['C',"""O''""","""O'"""],'name':'CTERM'}}
        self.titratable_groups={}
        #
        # Find all the titratable groups
        #
        for residue in self.residues.keys():
            rname=self.resname(residue)
            if self.isaa(residue):
                if self.titgroups.has_key(rname) or self.Cterminal(residue) or self.Nterminal(residue):
                    self.titratable_groups[residue]=[]
                    
                    if self.titgroups.has_key(rname):
                        if not self.is_SSbonded(residue): #Check for Cys
                            self.titratable_groups[residue].append(self.titgroups[rname])
                    if self.Cterminal(residue):
                        self.titratable_groups[residue].append(self.titgroups['CTERM'])
                    if self.Nterminal(residue):
                        self.titratable_groups[residue].append(self.titgroups['NTERM'])
                    #
                    # If nothing was added (SS bond) then delete
                    #
                    if self.titratable_groups[residue]==[]:
                        del self.titratable_groups[residue]
                    
        #
        # Make sure that we have all the atoms
        #
        grps=self.titratable_groups.keys()
        grps.sort()

        for residue in grps:
            count=0
            for group in self.titratable_groups[residue]:
                okatoms=[]
                for atom in group['atoms']:
                    atomname='%s:%s' %(residue,atom)
                    if self.atoms.has_key(atomname):
                        okatoms.append(atom)
                    else:
                        error=0
                        if atom=="O'":
                            if self.atoms.has_key('%s:O' %residue):
                                okatoms.append('O')
                            else:
                                error=1
                        elif atom=="O''":
                            if self.atoms.has_key('%s:OXT' %residue):
                                okatoms.append('OXT')
                            elif self.atoms.has_key('%s:O2' %residue):
                                okatoms.append('O2')
                            else:
                                error=1
                        else:
                            error=1
                        if error==1:
                            print 'Could not find atom: %s\nIgnoring...' %atomname
                if okatoms==[]:
                    print 'Cannot find any atoms for this group'
                self.titratable_groups[residue][count]['atoms']=okatoms[:]
                count=count+1
        #
        # Return a list of the names
        #
        group_names=[]
        residues=self.titratable_groups.keys()
        residues.sort()
        for residue in residues:
            for group in self.titratable_groups[residue]:
                group_name='%s:%s' %(residue,group['name'])
                group_names.append(group_name)
        group_names.sort()
        return group_names
        
    #
    # ------
    #
    
    def calculate_titgroup_center(self,residue):
        """Given a residue, check if it's a titgroup and return its centers (there might be more than one)"""
        if not hasattr(self,'titratable_groups'):
            self.get_titratable_groups()
        centers={}
        if self.titratable_groups.has_key(residue):
            centers[residue]={}
            for titgroup in self.titratable_groups[residue]:
                centre=self.calculate_group_center(residue,titgroup['atoms'])
                centers[residue][titgroup['name']]=centre
        return centers
    #
    # ----
    #

    def calculate_group_center(self,residue,atoms):
        """Calculate the center of the atoms"""
        ratoms=[]
        for atom in atoms:
            ratoms.append('%s:%s' %(residue,atom))
        return self.calculate_center(ratoms)

    #
    # -----
    #

    def calculate_center(self,atoms):
        """Calculate the geometric center of the atoms in atoms"""
        x=0.0
        y=0.0
        z=0.0
        for atom in atoms:
            x=x+self.atoms[atom]['X']/float(len(atoms))
            y=y+self.atoms[atom]['Y']/float(len(atoms))
            z=z+self.atoms[atom]['Z']/float(len(atoms))
        return numpy.array([x,y,z])

    #
    # ----
    #

    def calculate_matrix_PBE(self,filename=None,eps_map=[],exdi=80,updater=None):
        """Calculate the interaction energy matrix using Chresten and Kaspers
        PBE solver"""
        #
        # And now we calculate the matrix
        #
        import string
        if not getattr(self,'titratable_groups',None):
            self.get_titratable_groups()
        self.matrix={}
        titgrps=self.titratable_groups
        residues=titgrps.keys()
        residues.sort()
        #
        # Get the total number of calcs
        #
        ncalcs=len(residues)*len(residues)
        count=0
        #
        # Loop over all residues
        #
        for residue in residues:
            for group in titgrps[residue]:
                name1=self.get_titgroup_name(residue,group) 
                if not self.matrix.has_key(name1):
                    self.matrix[name1]={}
                #
                # Calculate center of atoms for this group
                #
                center1=self.calculate_group_center(residue,group['atoms'])
                #
                # Place charge at this center and solve PBE
                #
                atomlist=self.writepdb('test.pqr',nowrite=1,pqr=1)
                #
                # Now for the second residue
                #
                for residue2 in titgrps.keys():
                    for group2 in titgrps[residue2]:
                        if updater:
                            updater('%d of %d (%4.1f%% completed)' %(count,ncalcs,float(count)/float(ncalcs)*100.0))
                        count=count+1
                        name2=self.get_titgroup_name(residue2,group2)
                        if self.matrix[name1].has_key(name2):
                            continue
                        #
                        # Get center of atoms for this res
                        #
                        center2=self.calculate_group_center(residue2,group2['atoms'])

        return

    #
    # ----
    #

    def calculate_matrix(self,filename=None,eps=8,exdi=80,updater=None):
        """
        # This function calculates a simplified version of the charge-charge interaction
        # energy matrix
        #
        #
        # Calculate the factor that converts formal charges and distances in A into kT energies
        # eps is the dielectric constant
        """
        factor=(1.0/(4.0*pi*e0))*(e*e/A)*1.0/eps*1.0/(k*T)
        #
        # Get the fast scaled accessibility
        #
        acc=self.get_atoms_around(updater=updater)
        #
        # And now we calculate the matrix
        #
        import string
        if not getattr(self,'titratable_groups',None):
            self.get_titratable_groups()
        self.matrix={}
        self.distmatrix={}
        titgrps=self.titratable_groups
        residues=titgrps.keys()
        residues.sort()
        #
        # Get the total number of calcs
        #
        ncalcs=len(residues)*len(residues)
        count=0
        for residue in residues:
            for group in titgrps[residue]:
                name1=self.get_titgroup_name(residue,group) 
                if not self.matrix.has_key(name1):
                    self.matrix[name1]={}
                if not self.distmatrix.has_key(name1):
                    self.distmatrix[name1]={}
                #
                # Calculate distance from center of atoms
                #
                x=0.0
                y=0.0
                z=0.0
                for atom in group['atoms']:
                    x=x+self.atoms[residue+':'+atom]['X']/float(len(group['atoms']))
                    y=y+self.atoms[residue+':'+atom]['Y']/float(len(group['atoms']))
                    z=z+self.atoms[residue+':'+atom]['Z']/float(len(group['atoms']))
                #
                # Now for the second residue
                #
                for residue2 in titgrps.keys():
                    for group2 in titgrps[residue2]:
                        if not updater is None:
                            updater('%d of %d (%4.1f%% completed)' %(count,ncalcs,float(count)/float(ncalcs)*100.0))
                        count=count+1
                        name2=self.get_titgroup_name(residue2,group2)
                        if self.matrix[name1].has_key(name2):
                            continue
                        #
                        # Get center of atoms for this res
                        #
                        x2=0.0
                        y2=0.0
                        z2=0.0
                        for atom in group2['atoms']:
                            x2=x2+self.atoms[residue2+':'+atom]['X']/float(len(group2['atoms']))
                            y2=y2+self.atoms[residue2+':'+atom]['Y']/float(len(group2['atoms']))
                            z2=z2+self.atoms[residue2+':'+atom]['Z']/float(len(group2['atoms']))
                        #
                        # Get the distance
                        #
                        v1=numpy.array([x,y,z])
                        v2=numpy.array([x2,y2,z2])
                        import geometry
                        dist=geometry.length(v1-v2)
 
                        if dist!=0:
                            crg1=group['charge']
                            crg2=group2['charge']
                            ene=factor*float(crg1*crg2)/(math.pow(dist,1.5))
                            #
                            # Get the influence of the acc
                            #
                            acc_frac=math.sqrt(acc[name1]*acc[name2])
                            acc_frac=math.pow(acc_frac,2.0)
                            ene=ene*1.99*acc_frac+(1.0-acc_frac)*float(eps)/float(exdi)*ene
                        else:
                            ene=0.0
                        #if ene:
                        #    print '%20s - %20s: %7.3f' %(name1,name2,ene)
                        self.matrix[name1][name2]=[ene,0.0,0.0,0.0]
                        self.distmatrix[name1][name2]=dist 
                        if not self.matrix.has_key(name2):
                            self.matrix[name2]={}
                            self.matrix[name2]={}
                        self.matrix[name2][name1]=[ene,0.0,0.0,0.0]
                        self.distmatrix[name2][name1]=dist
        #
        # Write the matrix file
        #
        if filename:
            import pKaIO
            X=pKaIO.pKaIO()
            X.matrix=self.matrix
            X.write_matrix(filename+'.MATRIX.DAT')
        return self.matrix
        
    #
    # -----
    #
    
    def Gaussian_chain_matrix(self,ss=False):
        """Calculate electrostatic interactions between residues based on a Gaussian-chain
        defined unfolded state"""
        #
        # Check if we have the titratable groups
        #
        if not getattr(self,'titratable_groups',None):
            self.get_titratable_groups()
        #
        self.matrix={}
        self.distmatrix={}
        titgrps=self.titratable_groups
        residues=titgrps.keys()
        residues.sort()
        #
        for grp1 in titgrps:
            self.matrix[grp1]={}
            for grp2 in titgrps:
                res_dist=self.res_reparation(self.resid(grp1),self.resid(grp2),ss=ss)
                effective_dist=7.5*math.sqrt(float(res_dist))+5
                d=effective_dist
                import numpy as np
                import math
                for probe_dist in np.arange(0,20,0.1):
                    p=math.pow(4*math.pi,2)*probe_dist*(3/(2*math.pi*d**2))**(3.0/2.0)*math.exp(-3*probe_dist**2/(2.0*d**2))                
        return matrix

    #
    # -----
    #

    def get_titgroup_name(self,residue,group):
        name=residue+':'+self.resname(residue)
        if group['name'] not in ['ASP','GLU','ARG','LYS','HIS','CYS','TYR']:
            name='%s:%s' %(name,group['name'])
        return name

    #
    # -----
    #
                        
    def get_net_charge_from_Bfact(self,uniqueid):
        #
        # This option simply sums the charges in the B-factor column
        #
        if self.is_residue(uniqueid):
            crg=0.0
            for atom in self.residues[uniqueid]:
                crg=crg+self.atoms[atom]['B-FACTOR']
            return crg
        elif self.is_atom(uniqueid):
            return self.atoms[uniqueid]['B-FACTOR']
        return None

    #
    # ----------------------------
    # 
    
    def get_atoms_around(self,cutoff_dist=8.0,updater=None):
        """Count the number of atoms close to a titratable group"""
        #
        # First we need the titratable groups
        #
        if not getattr(self,'titratable_groups',None):
            self.get_titratable_groups()
        titgrps=self.titratable_groups
        #
        # Start calculating
        #
        self.atom_close={}
        self.close_atoms_list={}
        residues=titgrps.keys()
        residues.sort()
        count=0
        for residue in residues:
            if updater:
                updater('Analyzing protein geometry.. %5.1f%% done' %(float(count)/len(residues)*100.0))
            count=count+1
            for group in titgrps[residue]:
                name1=self.get_titgroup_name(residue,group)
                #
                # Find atoms close to these atoms
                #
                close_atoms={}
                for residue2 in self.residues.keys():
                    test_atom=self.residues[residue2][0]
                    if self.dist(test_atom,residue+':'+group['atoms'][0])<cutoff_dist+25.0:
                        for atom2 in self.residues[residue2]:
                            for atom in group['atoms']:
                                if self.dist(residue+':'+atom,atom2)<cutoff_dist:
                                    if not atom2 in self.residues[residue]: #Do not count atoms in same residue
                                        close_atoms[atom2]=1
                #
                # Store the value
                #
                self.close_atoms_list[name1]=close_atoms.keys()
                self.atom_close[name1]=len(close_atoms.keys())
        #
        # Scale from 0 to 1. 0 contacts is the real minimum
        # Maximum is around 160
        #
        min=9999
        max=-9999
        for name in self.atom_close.keys():
            value=self.atom_close[name]
            if value<min:
                min=value
            if value>max:
                max=value
        #
        # Scale
        #
        span=160-0
        for name in self.atom_close.keys():
            self.atom_close[name]=float((self.atom_close[name]-min))/float(span)
        return self.atom_close

    #
    # -----
    #

    def divide_atoms(self):
        """Divide the atoms in the molecule into bins - 8A cubed for easy neighbour finding"""
        #
        # find maximum and minimum of x,y and z
        #
        xmax=-999
        xmin=999
        for atom in self.atoms.keys():
            pass
        return

    #
    # ------
    #

    def calculate_desolvation(self,filename=None,updater=None):
        """Calculate the desolvation energies for all titratable groups"""
        if not getattr(self,'atom_close',None):
            self.get_atoms_around()
        #
        #
        titgrps=self.titratable_groups
        residues=titgrps.keys()
        residues.sort()
        self.desolv={}
        for residue in residues:
            if updater:
                updater('Calculating for %s' %residue)
            for group in titgrps[residue]:
                name1=self.get_titgroup_name(residue,group)
                #
                # Calculate the desolvation
                #
                self.desolv[name1]=12.5*math.pow(self.atom_close[name1],1.7)
        #
        # Write the file
        #
        if filename:
            import pKaIO
            X=pKaIO.pKaIO()
            X.desolv=self.desolv
            X.write_desolv(filename+'.DESOLV.DAT')
        return self.desolv

    #
    # -----
    #

    def calculate_background(self,filename=None,updater=None):
        """Calculate the background interaction energy"""
        if not getattr(self,'atom_close',None):
            self.get_atoms_around()
        #
        # Calc the background ene
        #
        titgrps=self.titratable_groups
        residues=titgrps.keys()
        residues.sort()
        self.background={}
        count=0
        for residue in residues:
            if updater:
                updater('Calculating for %s' %residue)
            for group in titgrps[residue]:
                name1=self.get_titgroup_name(residue,group)
                #
                # Calculate the background interaction energy
                #
                # First approximation: we count potential hydrogen bonds
                #
                hbonds=0
                tot_bonds=0
                #
                # Get the sign right
                #
                counted={}
                sign=group['charge']
                for atom in self.close_atoms_list[name1]:
                    
                        for atom2 in group['atoms']:
                            #
                            # First for full-strength hydrogen bonds
                            #
                            hb=self.is_hydrogen_bond(atom,residue+':'+atom2,2.2,3.5)
                            
                            if hb==1:
                                if not counted.has_key(atom):
                                    #print 'Accepting from',atom,residue+':'+atom2
                                    hbonds=hbonds-1*sign
                            elif hb==2:
                                if not counted.has_key(atom):
                                    #print 'Donating to',atom,residue+':'+atom2
                                    hbonds=hbonds+1*sign
                            #
                            if hb:
                                tot_bonds=tot_bonds+1
                                #counted[atom]=1
                            #
                            # Then for weaker - purely electrostatic interactions
                            # These are only counted half
                            #
                            hb=self.is_hydrogen_bond(atom,residue+':'+atom2,3.5,7.0)
                            weak_hb=0.3
                            if hb==1:
                                if not counted.has_key(atom):
                                    #print 'Accepting from',atom,residue+':'+atom2
                                    hbonds=hbonds-weak_hb*sign
                            elif hb==2:
                                if not counted.has_key(atom):
                                    #print 'Donating to',atom,residue+':'+atom2
                                    hbonds=hbonds+weak_hb*sign
                            #
                            if hb:
                                tot_bonds=tot_bonds+weak_hb
                                #counted[atom]=1
                
                #
                #
                #print 'Hbonds: %3d, tot_bonds: %3d' %(hbonds,tot_bonds)
                self.background[name1]=0.69*hbonds*math.pow(self.atom_close[name1],1)
        #
        # Write the file
        #
        if filename:
            import pKaIO
            X=pKaIO.pKaIO()
            X.backgr=self.background
            X.write_backgr(filename+'.BACKGR.DAT')
        return self.background
        
    #
    # --------------
    #
    
    def find_salt_bridges(self,cutoff=8.0):
        """Find all salt bridges in the protein
        """
        SBs=[]
        self.get_titratable_groups()
        for res1 in self.titratable_groups.keys():
            r1_cens=self.calculate_titgroup_center(res1)
            for group1 in r1_cens[res1].keys():
                group1_ID='%s:%s' %(res1,group1)
                for res2 in self.titratable_groups.keys():
                    r2_cens=self.calculate_titgroup_center(res2)
                    for group2 in r2_cens[res2].keys():
                        group2_ID='%s:%s' %(res2,group2)
                        if group1_ID!=group2_ID:
                            if self.titgroups[group1]['charge']*self.titgroups[group2]['charge']==-1:
                                import geometry
                                dist=geometry.length(r1_cens[res1][group1]-r2_cens[res2][group2])
                                if dist<=cutoff:
                                    SBs.append([group1_ID,group2_ID,dist])
        #
        # Done
        #
        return SBs
                    
                    
                    
                
                
