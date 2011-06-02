#!/usr/bin/env python
#
# Protool - Python class for manipulating protein structures
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
class flags:
    
    def bound(self,atom1,atom2):
        #
        # Two atoms are covalently bound if they are closer than
        # 2.0 A. This is not true for hydrogens....
        #
        #print self.dist(atom1,atom2)
        if self.dist(atom1,atom2)<2.0:
            return True
        else:
            return False


    #
    # ----
    #

    def is_hydrogen_bond(self,atom1,atom2,low_cutoff=2.2,high_cutoff=3.5):
        """Returns 1 or 2 if there is a potential hydrogen bond between these two atoms.
        Returns None otherwise. Returns 1 if the first atom is the donor, and 2
        if the second atom is the donor"""
        #
        # INCLUDE ANGLE REQUIREMENTS!
        #
        if self.is_hbond_donor(atom1) and self.is_hbond_acceptor(atom2):
            if self.dist(atom1,atom2)>low_cutoff and self.dist(atom1,atom2)<=high_cutoff:
                return 1
        if self.is_hbond_acceptor(atom1) and self.is_hbond_donor(atom2):
            if self.dist(atom1,atom2)>low_cutoff and self.dist(atom1,atom2)<=high_cutoff:
                return 2
        return None

    #
    # -----
    #

    def is_hbond_donor(self,atom):
        """Returns 1 if the atom can be a hydrogen donor"""
        atomname=atom.split(':')[-1]
        residuename=self.resname(atom)
        donors={'LYS':['NZ'],
                'ARG':['NE','NH1','NH2'],
                'SER':['OG'],
                'THR':['OG1'],
                'TYR':['OH'],
                'HIS':['ND1','NE2'],
                'CYS':['SG'],
                'ASN':['ND2'],
                'GLN':['NE2'],
                'TRP':['NE1']}
        #
        # Backbone
        #
        if self.is_mainchain(atom) and atomname=='N' and residuename!='PRO':
            return 1
        #
        # Other residues
        #
        if donors.has_key(residuename):
            if atomname in donors[residuename]:
                return 1
        return None


    def is_hbond_acceptor(self,atom):
        """Returns 1 if the atom can be a hydrogen acceptor"""
        atomname=atom.split(':')[-1]
        residuename=self.resname(atom)
        acceptors={'ASP':['OD1','OD2'],
                   'GLU':['OE1','OE2'],
                   'SER':['OG'],
                   'THR':['OG1'],
                   'TYR':['OH'],
                   'HIS':['ND1','NE2'],
                   'ASN':['OD1'],
                   'GLN':['OE1']}
        #
        # Backbone
        #
        if self.is_mainchain(atom) and atomname=='O':
            return 1
        if self.is_mainchain(atom) and atomname=='N' and residuename=='PRO':
            return 1
        #
        # Other residues
        #
        if acceptors.has_key(residuename):
            if atomname in acceptors[residuename]:
                return 1
        return None

    #
    # ---------
    #

    def is_hydrophilic(self,atom):
        """Returns 1 if the atom is hydrophilic"""
        if self.is_oxygen(atom) or self.is_nitrogen(atom):
            return 1
        return None

    #
    # ----------
    #
        

    def Nterminal(self,residue):
        #
        # Returns 1 if <residue> is an N-terminal residue
        #
        try:
            x=self.PreviousResidue(residue)
        except Nterm:
            return 1
        except NotAnAminoAcidError:
            return 1
        return None
    
    def Cterminal(self,residue):
        #
        # Returns 1 if residue is an C-terminal residue
        #
        try:
            x=self.NextResidue(residue)
        except Cterm:
            return 1
        except NotAnAminoAcidError:
            return 1
        return None

    def is_backbone(self,atom):
        return self.is_mainchain(atom)

    #
    # ---
    #

    def is_mainchain(self,atom):
        #
        # Returns 1 if the atom is a mainchain (backbone) atom
        #
        import string
        atomname=string.split(atom,':')[-1]
        if atomname=='CA' or atomname=='N' or atomname=='C' \
           or atomname=='O' or atomname=='H':
            return 1
        else:
            return None

    #
    # ----
    #

    def is_CB(self,atom):
        """Returns 1 if the atom is CB"""
        import string
        atomname=string.split(atom,':')[-1]
        if atomname=='CB':
            return 1
        else:
            return None

    def phosphorylated(self, residue):
        #
        #Returns 1 if the residue is an phosphorylated residue
        #
        import string
        atomname=string.split(atom,':')[-1]
        if atomname=='P' or atomname=='O1P' or atomname=='O2P' or atomname=='O3P':
            return 1
        else:
            return None

    def is_residue(self,uniqueid):
        # Given a uniqueid, returns 1 if uniqueid is a residue, and returns None if uniqueid
        # is something else (e.g. a chain or an atom)
        if self.residues.has_key(uniqueid):
            return 1
        return None
    
    # --------------------------------

    def is_atom(self,uniqueid):
        # Given a uniqueid, returns 1 if uniqueid is an atom, and returns None if uniqueid
        # is something else (e.g. a chain or a residue)
        if self.atoms.has_key(uniqueid):
            return 1
        return None        

    def isaa(self,uniqueid):
        # Given an uniqueid, this function returns 1 if the residue/atom is part of an amino acid, and
        # None if it is not
        import string
        if self.aminoacids.has_key(string.upper(self.resname(uniqueid))):
            return 1
        return None

    def is_hydrogen(self,uniqueid):
        # Get the atom name
        import string
        atomname=self.atoms[uniqueid]['ATNAME']
        return self.atomname_is_hydrogen(atomname)

    #
    # ----
    #

    def atomname_is_hydrogen(self,atomname):
        import string
        if atomname[0]=='H':
            return 1
        if len(atomname)>2:
            if atomname[0] in string.digits and atomname[1]=='H':
                return 1
        return None

    #
    # -----
    #

    def is_oxygen(self,atom):
        """Returns 1 if the current atom is an oxygen"""
        atomname=atom.split(':')[-1]
        if atomname.find('O')!=-1:
            return 1
        return None

    def is_nitrogen(self,atom):
        """Returns 1 if the current atom is a nitrogen"""
        atomname=atom.split(':')[-1]
        if atomname.find('N')!=-1:
            return 1
        return None

    def is_carbon(self,atom):
        """Returns 1 if the current atom is a carbon"""
        atomname=atom.split(':')[-1]
        if atomname.find('C')!=-1:
            return 1
        return None

    def is_sulphur(self,atom):
        """Returns 1 if the current atom is a carbon"""
        atomname=atom.split(':')[-1]
        if atomname.find('S')!=-1:
            return 1
        return None

    def is_phosphor(self,atom):
        """Returns 1 if the current atom is a carbon"""
        atomname=atom.split(':')[-1]
        if atomname.find('P')!=-1:
            return 1
        return None

    def get_atom_element(self,atom):
        """Return the element type of an atom """
        #
        # Get rid of any ALT records
        #
        atom=atom.replace(',ALT','')
        atom=atom.replace(':ALT','')
        
        atomname=atom.split(':')[-1]
        residue_name=self.resname(atom)
        if self.isaa(atom):
            #
            # Amino acid, then return first letter
            #
            import string
            for char in atomname:
                if char in ['H','C','N','O','S','P']:
                    return char
            #
            # Unknown amino acid atom??
            #
            print atom
            return None
        else:
            #
            # Not a standard amino acid
            # Look for ions first
            #
            ions=['CU','CA','ZN','CO','NA','MN','MG']
            for ion in ions:
                if atomname.find(ion)!=-1 and residue_name in ['ION',ion]:
                    return ion
            #
            # No ions found
            # Simply return fist letter
            #
            import string
            for char in atomname:
                if char in string.letters:
                    return char
        #
        # No resolution
        #
        return None

    
    def is_backbone(self,uniqueid):
        # Get the atom name
        atomname=self.atoms[uniqueid]['ATNAME']
        backbone=['H','CA','N','O','C','HA']
        if atomname in backbone:
            return 1
        return None


    def is_CA(self,uniqueid):
        """Is this a CA atom?"""
        atomname=self.atoms[uniqueid]['ATNAME']
        if atomname=='CA':
            return 1
        return None

    #
    # ---
    #

    def is_SSbonded(self,uniqueid):
        """Is this residue or atom in a residue that's SS-bonded"""
        if self.resname(uniqueid)=='CYS':
            SG1='%s:%s' %(self.resid(uniqueid),'SG')
            for residue in self.residues.keys():
                if residue==self.resid(uniqueid):
                    continue
                if self.resname(residue)=='CYS':
                    SG2='%s:%s' %(residue,'SG')
                    if self.dist(SG1,SG2)<2.5:
                        return residue
        return None
    

