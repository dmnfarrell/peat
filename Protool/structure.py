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


#
from errors import *
import electrostatics 
from geometry import *
from paths import *
from Modify import *
import flags, analyse_structure, energy

try:
    import numpy as Numeric
except:
    import Numeric
    
try:
    from PEATDB.DNA_sequence import *
    #print 'Using PEATDB sequence module'
except:
    from dummy_sequence import *
    #print 'Cannot find PEATDB'
                             
##############


class display:

    #
    # Class for displaying structures
    #
    def showtot(self,colour=None,BF1=10,BF2=50):
        #
        # Show it all. Colour can hold a number of uniqueids (atoms or residues) that
        # will displayed in a different colour than the rest. colour is a dictionary
        #
        import string
        colcmd=''
        if colour:
            for residue in self.residues.keys():
                resname=string.split(residue,':')[1]
                if colour.has_key(residue):
                    for atom in self.residues[residue]:
                        self.atoms[atom]['B-FACTOR']=BF2
                else:
                    for atom in self.residues[residue]:
                        if colour.has_key(atom):
                            self.atoms[atom]['B-FACTOR']=BF2
                        else:
                            self.atoms[atom]['B-FACTOR']=BF1
            colcmd='colour atoms temperature\n'
        import tempfile, os
        tmp=tempfile.mktemp()
        self.writepdb(tmp,'load pdb inline \n'+colcmd+'exit\n')
        os.system (RASMOL+' -script '+tmp+' >/dev/null')
        os.unlink(tmp)
        return


#
# --------------------------------------------------------------------------------------
#
    
class structure(geometry,flags.flags,sequence,analyse_structure.find_pattern,
                display,modify,analyse_structure.distances,
                electrostatics.electrostatics,energy.energy):

    def __init__(self):
        self.aminoacids={'ALA':[],'CYS':[],'ASP':[],'GLU':[],'PHE':[],'GLY':[],'HIS':[],'ILE':[],'LYS':[],'LEU':[],'MET':[],'ASN':[],'PRO':[],'GLN':[],'ARG':[],'SER':[],'THR':[],'VAL':[],'TRP':[],'TYR':[],'TPO':[],'SEP':[]}
        self.trueaminoacids={'ALA':[],'CYS':[],'ASP':[],'GLU':[],'PHE':[],'GLY':[],'HIS':[],'ILE':[],'LYS':[],'LEU':[],'MET':[],'ASN':[],'PRO':[],'GLN':[],'ARG':[],'SER':[],'THR':[],'VAL':[],'TRP':[],'TYR':[]}
        self.three_to_one={'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','TPO':'Z','SEP':'X'}
        self.one_to_three={'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR','TPO':'Z','SEP':'X'}
        #
        # The attribute dictionary
        #
        self.attribute={'G':{'size':0,'class':{'flexible':1,'small':1}},
                        'A':{'size':1,'class':{'small':1,'aliphatic':1,'hydrophobic':1}},
                        'V':{'size':3,'class':{'hydrophobic':1,'betasub':1,'aliphatic':1}},
                        'L':{'size':4,'class':{'hydrophobic':1,'aliphatic':1}},
                        'I':{'size':4,'class':{'hydrophobic':1,'aliphatic':1}},
                        'S':{'size':2,'class':{'small':1,' polar':1},'pKa':15.0},
                        'C':{'size':2,'class':{'small':1,'polar':1,'negative':1,'charged':1},'pKa':8.7},
                        'T':{'size':3,'class':{'polar':1,'small':1,'betasub':1},'pKa':15.0},
                        'M':{'size':4,'class':{'aliphatic':1,'hydrophobic':1}},
                        'F':{'size':7,'class':{'aromatic':1,'hydrophobic':1}},
                        'Y':{'size':8,'class':{'aromatic':1,'polar':1,'hydrophobic':1},'pKa':9.6},
                        'W':{'size':10,'class':{'hydrophobic':1,'polar':1,'aromatic':1}},
                        'P':{'size':3,'class':{'small':1,'rigid':1,'hydrophobic':1}},
                        'H':{'size':6,'class':{'positive':1,'polar':1,'charged':1},'pKa':6.3},
                        'K':{'size':5,'class':{'positive':1,'polar':1,'charged':1},'pKa':10.4},
                        'R':{'size':7,'class':{'polar':1,'positive':1,'charged':1},'pKa':13.0},
                        'D':{'size':4,'class':{'negative':1,'charged':1,'polar':1},'pKa':4.0},
                        'E':{'size':5,'class':{'negative':1,'charged':1,'polar':1},'pKa':4.4},
                        'N':{'size':4,'class':{'polar':1}},
                        'Q':{'size':5,'class':{'polar':1}},
                        'Z':{'size':7,'class':{'phosphorylated':1}},
                        'X':{'size':6,'class':{'phosphorylated':1}}
                        }
        #
        # Van der Waals radii
        #
        self.vdwr={'C':1.7,
                  'H':1.2,
                  'N':1.55,
                  'O':1.52,
                  'S':1.8,
                  'P':1.9,
                  'CA':2.0,
                  'K':2.75,
                  'NA':2.27,
                   'ZN':1.39,
                   'MG':2.0,
                   'MN':2.0}
        #
        # Alternative atom names for use with AutoDock
        #
        self.autodock_names={}
        self.autodock_names['H']='HN  '
        self.autodock_names['2H']='HN2 '
        self.autodock_names['3H']='HN3 '
        #
        self.autodock_names['1HH1']='HH11'
        self.autodock_names['1HH2']='HH12'
        self.autodock_names['2HH1']='HH21'
        self.autodock_names['2HH2']='HH22'
        # LYS
        self.autodock_names['3HZ']='HZ3 '
        self.autodock_names['2HZ']='HZ2 '
        self.autodock_names['1HZ']='HZ1 '
        # ASN
        self.autodock_names['1HD2']='HD12'
        self.autodock_names['2HD2']='HD22'
        # GLN
        self.autodock_names['1HE2']='HE21'
        self.autodock_names['2HE2']='HE22'

        self.autodock_names["O'"]='O  '
        self.autodock_names["O''"]='OXT'
        
        
        #
        # Length of the residue numbers
        #
        self.length_of_residue_numbers=4
        #
        # Flag for checking for TER records and assigning new chain IDs
        #
        self.parse_terms=1
        #
        # Init vars
        #
        self.atoms={}
        self.chains={}
        self.residues={}
        
        return
        
    def resid(self,atomname):
        """Given an atom name this function return a resid"""
        sp=atomname.split(':')
        return '%s:%s' %(sp[0],sp[1])

    def resnum(self,uniqueid):
        """Given a uniqueid this function returns the residue number"""
        import string
        if len(string.split(uniqueid,','))>1:
            return string.split(uniqueid,':')[0]+':'+string.zfill(string.split(uniqueid,':')[1],self.length_of_residue_numbers)+','+string.split(uniqueid,',')[1]
        else:
            return string.split(uniqueid,':')[0]+':'+string.zfill(string.split(uniqueid,':')[1],self.length_of_residue_numbers)
        


    def resname(self,uniqueid):
        """ Given an uniqueid or a residue number, this function returns the residue name"""
        import string
        if self.atoms.has_key(uniqueid):
            return self.atoms[uniqueid]['RESNAME']
        elif self.residues.has_key(uniqueid):
            #
            # Get the first atom of the residue and return the resname there
            #
            atomname=self.residues[uniqueid][0]
            if not self.atoms.has_key(atomname):
                raise ProtoolError, 'Inconsistent self.atoms and self.residues: %s is not found in atoms.' %atomname
            return self.atoms[atomname]['RESNAME']
        elif self.residues.has_key(':'+uniqueid):
            return self.atoms[self.residues[':'+uniqueid][0]]['RESNAME']
        #
        # Check for non-existing atom, but existing residue
        #
        elif self.residues.has_key(string.join(uniqueid.split(':')[:-1],':')):
            uniqueid=string.join(uniqueid.split(':')[:-1],':')
            atomname=self.residues[uniqueid][0]
            if not self.atoms.has_key(atomname):
                raise ProtoolError, 'Inconsistent self.atoms and self.residues: %s is not found in atoms.' %atomname
            return self.atoms[atomname]['RESNAME']
        elif self.residues.has_key(':'+uniqueid.split(':')[1]):
            res=':'+uniqueid.split(':')[1]
            at0=self.residues[res][0]
            resname=self.atoms[at0]['RESNAME']
            if resname==uniqueid.split(':')[2]:
                return self.atoms[at0]['RESNAME']
            else:
                raise ResidueNotFoundErrror,'Residue not found for %s' %uniqueid
        else:
            raise ResidueNotFoundError,'Residue not found for %s' %uniqueid

    def pKa_ID(self,uniqueid):
        #
        # Given a residue id, this function returns a uniqueID
        # for use with the pKa routines
        #
        resname=self.resname(uniqueid)
        resnumber=self.resnum(uniqueid)
        #chainID=self.chainid(uniqueid)
        return resnumber+':'+resname

    def chainid(self,uniqueid):
        # Given an uniqueid this function returns the chain id#
        if self.atoms.has_key(uniqueid):
            return uniqueid.split(':')[0]
        else:
            if self.residues.has_key(uniqueid):
                atom0=self.residues[uniqueid][0]
                return atom0.split(':')[0]
        raise Exception('Cannot find uniquid: %s' %uniqueid)

    def atname(self,uniqueid):
        # Given an uniqueid this function returns the atom name
        import string
        return string.join(string.split(uniqueid,':')[2:],':')

    def isnormal(self,uniqueid):
        # Given an uniqueid or a residue number, this function returns 1 if the residue/atom is
        # a non-ALT entity, and None if it is an 'ALT' entity
        if uniqueid[-3:]=='ALT':
            return None
        return 1

    def ispresent(self,uniqueid):
        # Given an uniqueid or a residue number, this function returns 1 if the residue/atom is
        # present, and None if it is not
        if self.atoms.has_key(uniqueid) or self.residues.has_key(uniqueid):
            return 1
        return None

    # --------------------------------


    #
    # ----------
    #
    
    def PreviousResidue(self,resid,checkforbound=1):
        #
        # Both this function and the Next Residue should be rewritten!!!!
        #
        #
        # This function returns the residue one before resid in the sequence in the same molecule
        # No checking for that the two are connected!!
        #
        if not self.isaa(resid):
            raise NotAnAminoAcidError
        residues=self.residues.keys()
        residues.sort()
        previous=None
        for residue in residues:
            if residue==resid:
                if previous:
                    if self.chainid(resid)==self.chainid(previous) and self.isnormal(resid)==self.isnormal(previous):
                        #
                        # Check that the N of resid is bound to the C of previous
                        #
                        N=resid+':N'
                        C=previous+':C'
                        if not self.bound(N,C) and checkforbound:
                            raise Nterm()
                        return previous
                    else:
                        raise Nterm()
                else:
                    raise Nterm()
            previous=residue
        raise ResidueNotFoundError(resid)

    #
    # --------------
    #

    def NextResidue(self,resid,checkforbound=1):
        #
        # This function returns the residue one after resid in the sequence in the same molecule
        # No checking for that the two are connected!!
        #
        if not self.isaa(resid):
            raise NotAnAminoAcidError('%s %s' %(resid,self.resname(resid)))
        residues=self.residues.keys()
        residues.sort()
        for x in range(len(residues)):
            if residues[x]==resid:
                if x!=len(residues)-1:
                    if self.chainid(resid)==self.chainid(residues[x+1]) and self.isnormal(resid)==self.isnormal(residues[x+1]):
                        # Check that the C of resid is bound to the N of the residues[x+1]
                        C=resid+':C'
                        N=residues[x+1]+':N'
                        if not self.bound(N,C) and checkforbound:
                            raise Cterm
                        return residues[x+1]
                    else:
                        raise Cterm
                else:
                    raise Cterm
        raise 'Not Found','The residue was not found in this molecule'
    
        


    def GetPosition(self,uniqueid):
        """
        # This function loads the X, Y and Z coordinates of the atom
        # in an Numeric Python vector
        """
        if not self.atoms.has_key(uniqueid):
            raise AtomNotFoundError,uniqueid
        if not self.atoms[uniqueid].has_key('X') or not self.atoms[uniqueid].has_key('Y') or not self.atoms[uniqueid].has_key('Z'):
            print self.atoms[uniqueid]
            raise IncompletePositionError, uniqieid
        return Numeric.array([self.atoms[uniqueid]['X'],self.atoms[uniqueid]['Y'],self.atoms[uniqueid]['Z']])

    #
    # ----
    #
    
    def add_atom(self,uniqueid,
                 atomnumber=0,atomname='X',
                 chainid='X',residuename='DUM',residuenumber='000',
                 xcoord=0.0,ycoord=0.0,zcoord=0.0,update=1,BFACTOR=None,OCCUPANCY=None,CHARGE=None,RADIUS=None,tag=None,accept_duplicate=False):
        """Add an atom to the arrays"""
        #
        # Correct for duplicate atom names
        #
        Fail=False
        if accept_duplicate:
            import string
            orguid=uniqueid
            
            if self.atoms.has_key(uniqueid):
                Fail=True
                for letter in string.letters:
                    if not self.atoms.has_key(uniqueid):
                        Fail=False
                        atomname=uniqueid.split(':')[-1]
                        print 'Changed uniqueid and atom name from: %s to %s' %(orguid,uniqueid)
                        break
                    uniqueid=orguid+letter
        #
        # Add the atom
        #
        if not self.atoms.has_key(uniqueid) and not Fail:
            self.atoms[uniqueid]={}
            self.atoms[uniqueid]['NUMBER']=atomnumber
            self.atoms[uniqueid]['ATNAME']=atomname
            self.atoms[uniqueid]['CHAINID']=chainid
            self.atoms[uniqueid]['RESNAME']=residuename
            self.atoms[uniqueid]['RESNUM']=residuenumber
            self.atoms[uniqueid]['X']=xcoord
            self.atoms[uniqueid]['Y']=ycoord
            self.atoms[uniqueid]['Z']=zcoord
            if not OCCUPANCY:
                OCCUPANCY=1.0
            self.atoms[uniqueid]['OCCUPANCY']=OCCUPANCY
            if not BFACTOR:
                BFACTOR=0.0
            self.atoms[uniqueid]['B-FACTOR']=BFACTOR
            if CHARGE:
                self.atoms[uniqueid]['CHARGE']=CHARGE
            if RADIUS:
                self.atoms[uniqueid]['RADIUS']=RADIUS
            if tag:
                self.atoms[uniqueid]['tag']=tag
            if update:
                self.Update()
            return 1
        else:
            raise Exception('ProTool: Trying to add an atom that is already present!: %s %s' %(uniqueid,str(self.atoms[uniqueid])))
        return None
        
    def renumber_residues(self,residues,start_num=1):
        """Renumber the given residues"""
        residues.sort()
        number=start_num
        for residue in residues:
            for atom in self.residues[residue]:
                new_number=string.zfill(number,self.length_of_residue_numbers)
                oldname=atom.split(':')
                new_name='%s:%s:%s' %(oldname[0],new_number,oldname[2])
                self.atoms[new_name]=self.atoms[atom].copy()
                self.atoms[new_name]['RESNUM']=new_number
                del self.atoms[atom]
            number=number+1
        self.Update()
        return

    #
    # ----
    #
    
    def remove_residue(self,unique_identifier,update=True):
        """Remove a residue. This includes deleted all atoms in the residue"""
        if not self.residues.has_key(unique_identifier):
            raise ResidueNotFoundError(unique_identifier)
        for atom in self.residues[unique_identifier]:
            self.remove_atom(atom,update=False)
        if update:
            self.Update()
        return

    def remove_atom(self,unique_identifier,update=True):
        """
        # Remove an atom
        """
        if self.atoms.has_key(unique_identifier):
            del self.atoms[unique_identifier]
            if update:
                self.Update()
        else:
            raise AtomNotFoundError(unique_identifier)
        return

    def RemoveALT(self):
        #
        # This function removes all alternative atoms from the structure
        #
        for atom in self.atoms.keys():
            if not self.isnormal(atom):
                del self.atoms[atom]
        self.Update()
        return
        
    def Remove_All_NonAminoAcids(self):
        """Remove all non amino acid residues"""
        for residue in self.residues.keys():
            if not self.isaa(residue):
                for atom in self.residues[residue]:
                    del self.atoms[atom]
        self.Update()
        return

    def Delete_residue(self,residue,update=True):
        """Delete a residue"""
        for atom in self.residues[residue]:
            del self.atoms[atom]
        if update:
            self.Update()
        return
        
    def remove_atoms_with_tag(self,tag=''):
        """Remove all atoms with the specified tag"""
        for atom in self.atoms.keys():
            if self.atoms[atom].has_key('tag'):
                if self.atoms[atom]['tag']==tag:
                    self.remove_atom(atom)
        return
        
    def remove_all_hydrogens(self):
        """Remove all hydrogen atoms"""
        for atom in self.atoms.keys():
            if self.is_hydrogen(atom):
                self.remove_atom(atom,update=False)
        self.Update()
        return

    def Update(self):
        """
        This function reconstructs self.residues and self.sequence from self.atoms
        and additionally performs some checks
        """
        self.makeresidues()
        self.makesequence()
        self.makechains()
        self.make_boxes()
        self.name_checks()
        return

    #
    # ----
    #

    def makeresidues(self):
        #
        # Construct self.residues
        #
        self.residues={}
        for atom in self.atoms.keys():
            resnum=self.resnum(atom)
            if self.residues.has_key(resnum):
                self.residues[resnum].append(atom)
            else:
                self.residues[resnum]=[atom]
        return 

    #
    # ----
    #

    def makechains(self):
        """Construct self.chains"""
        self.chains={}
        for residue in self.residues:
            cid=residue.split(':')[0]
            if not self.chains.has_key(cid):
                self.chains[cid]=[]
            self.chains[cid].append(residue)
        return

    #
    # ----
    #


    def makesequence(self):
        #
        # Construct self.sequence [Does not take into account that there can be different mols]
        #
        self.sequence=[]
        residues=self.residues.keys()
        residues.sort()
        for resnum in residues:
            #
            # We do not want to include 'ALT' residues in the sequence
            #
            if self.isnormal(resnum):
                #
                # We don't want HETATM only residues either
                #
                onlyhet=1
                for atom in self.residues[resnum]:
                    if not self.atoms[atom].has_key('HETATM'):
                        onlyhet=None
                        break
                if not onlyhet:
                    resname=self.resname(resnum)
                    self.sequence.append([resnum,resname])

        return
        
    #
    # ----
    #
    
    def make_boxes(self):
        """Divide all atoms into 6A boxes"""
        self.boxes={}
        self.small_boxes={}
        for atom in self.atoms:
            x_index=int(self.atoms[atom]['X']/6.0)
            y_index=int(self.atoms[atom]['Y']/6.0)
            z_index=int(self.atoms[atom]['Z']/6.0)
            index='%d_%d_%d' %(x_index,y_index,z_index)
            if not self.boxes.has_key(index):
                self.boxes[index]=[]
            self.boxes[index].append(atom)
            #
            self.atoms[atom]['box']=[x_index,y_index,z_index]
            #
            # --- Now the smaller box
            #
            x_index=int(self.atoms[atom]['X']/3.5)
            y_index=int(self.atoms[atom]['Y']/3.5)
            z_index=int(self.atoms[atom]['Z']/3.5)
            index='%d_%d_%d' %(x_index,y_index,z_index)
            if not self.small_boxes.has_key(index):
                self.small_boxes[index]=[]
            self.small_boxes[index].append(atom)
            #
            self.atoms[atom]['smallbox']=[x_index,y_index,z_index]
        return
            
    #
    # -----
    #
    
    def get_neighbour_atoms(self,atom,box='large'):
        """Get all atoms in this box and in neighbouring boxes"""
        atoms=[]
        if box=='large':
            name='box'
        else:
            name='smallbox'
        #
        index=self.atoms[atom][name]
        x=index[0]
        y=index[1]
        z=index[2]
        for xmod in range(-1,2):
            for ymod in range(-1,2):
                for zmod in range(-1,2):
                    box_index='%d_%d_%d' %(x+xmod,y+ymod,z+zmod)
                    if name=='box':
                        if self.boxes.has_key(box_index):
                            atoms=atoms+self.boxes[box_index]
                    else:
                        if self.small_boxes.has_key(box_index):
                            atoms=atoms+self.small_boxes[box_index]
        return atoms

    #
    # ----
    #
        
    def name_checks(self):
        """Do a simple name check"""
        residues=self.residues.keys()
        residues.sort()
        for residue in residues:
            changed=None
            if self.Cterminal(residue) and self.isaa(residue):
                for atom in self.residues[residue]:
                    if self.atoms.has_key(atom):
                        if (self.atname(atom)=="O''" or self.atname(atom)=='O2'):
                            import copy
                            cid,resid,name=copy.copy(atom).split(':')
                            name=name.replace(self.atname(atom),"OXT")
                            uniqueid=cid+':'+resid+':'+name
                            self.atoms[uniqueid]=self.atoms[atom].copy()
                            self.atoms[uniqueid]['ATNAME']=name
                            del self.atoms[atom]
                            changed=1
                        elif (self.atname(atom)=="O'" or self.atname(atom)=='O1'):
                            import copy
                            cid,resid,name=copy.copy(atom).split(':')
                            name=name.replace(self.atname(atom),"O")
                            uniqueid=cid+':'+resid+':'+name
                            self.atoms[uniqueid]=self.atoms[atom].copy()
                            self.atoms[uniqueid]['ATNAME']=name
                            del self.atoms[atom]
                            changed=1
            #
            # If something was changed then we have to update self.residues
            #
            if changed:
                self.makeresidues()
        return

    #
    # ------------------
    #

    def identical_positions(self):
        #
        # Checks that all atom positions are differ by at least 0.5 A in one
        # coordinate
        #
        for atom in self.atoms:
            pass
        return

    #
    # =========================================================================
    #

    def select(self,list):
        #
        # Selects the atoms or residues in list. All atoms of the selection are
        # stored in a new instance of self 
        #
        selection={}
        for id in list:
            if self.atoms.has_key(id):
                selection[id]=self.atoms[id]
            elif self.residues.has_key(id):
                for atom in self.residues[id]:
                    selection[atom]=self.atoms[atom]
            else:
                print id
                raise ResidueNotFoundError,'Residue or atom not found'
        #
        # Make new instance of self
        #
        N=structureIO()
        N.atoms=selection
        N.Update()
        return N

    #
    # ---------------
    #

    def hasMissingMainChainAtoms(self):
        """Check that each amino acid residue has all heavy main chain atoms"""
        missing=[]
        for residue in self.residues.keys():
            if not self.isaa(residue):
                continue
            #print self.residues[residue]
            for atom in ['O','C','CA','N']:
                if not '%s:%s' %(residue,atom) in self.residues[residue]:
                    missing.append('%s:%s' %(residue,atom))
        if len(missing)==0:
            return False
        else:
            return missing

    #
    # ----
    #

    def hasChainBreak(self):
        """Check if the soup contains chains with chain breaks"""
        for chain in self.chains.keys():
            residues=sorted(self.chains[chain])
            Cterms=0
            Nterms=0
            for residue in residues:
                if not self.isaa(residue):
                    continue
                try:
                    nextres=self.NextResidue(residue)
                except Cterm:
                    Cterms=Cterms+1
                try:
                    prevres=self.PreviousResidue(residue)
                except Nterm:
                    Nterms=Nterms+1
            if Nterms>1 or Cterms>1:
                return True
        return False

                
    
#--------------------------------------------------------------------------------------------
#--------------- This class can (sometimes) parse PDB files ---------------------------------
#--------------------------------------------------------------------------------------------
class PDBparser:
    
    def __init__(self,lines,parse_terms):
        #
        # Initialise atom arrays
        #
        import string
        self.lines=[]
        for line in lines:
            self.lines.append(string.strip(line))
        self.sequence=[]
        self.residues={}
        self.atoms={}
        self.attribute={}
        #
        # The remaining arrays are for the header information
        #
        self.header=''
        self.compnd=''
        self.source=''
        self.author=''
        self.revdat=''
        self.remark=''
        self.seqres=[]
        self.helix=[]
        self.sheet=[]
        self.turn=[]
        self.ssbond=[]
        self.cryst=[]
        self.orig=[]
        self.scale=[]
        self.NMRmodel=False
        self.readmodel=True
        self.readmodels=1
        #
        # Parse TER records and assign new chain IDs?
        #
        self.parse_terms=parse_terms
        #
        # Initialise some help variables
        #
        self.linenumber=0
        self.atomlinedefs={}
        self.chains={}
        self.useids=['J','K','L','M','O','P','Q','R','S','T','U','V','W','X','Y','Z']
        for x in string.letters:
            self.useids.append(string.lower(x))
        for x in string.digits:
            self.useids.append(string.upper(x))
        self.oldchainid=None
        self.freeresnum=9000
        #
        # Defining the junk headers
        #
        self.junkheaders={'TITLE':'REMARK','KEYWDS':'REMARK','DBREF':'REMARK','SEQADV':'REMARK','HET':'REMARK','HETNAM':'REMARK','CISPEP':'REMARK','LINK':'REMARK','MTRIX1':'REMARK','MTRIX2':'REMARK','MTRIX3':'REMARK','SITE':'REMARK','HETSYN':'REMARK','ANISOU':'JUNK','SPRSDE':'REMARK','MODRES':'REMARK','HYDBND':'REMARK','TVECT':'REMARK','SLTBRG':'REMARK','SIGUIJ':'REMARK','SIGATM':'REMARK','CAVEAT':'REMARK','JNRL':'JRNL'}
        return

    def parse(self):
        #
        # Loop through all lines and extract all readable information
        # (for now only ATOM lines are parsed)
        #
        for line in self.lines:
            try:
                self.parseline(line)
            except ParseLineError:
                ok=None
        #
        # Check for name clashes and atoms with identical positons.
        #
        self.atoms=self.nameclash()
        #
        # Renumber
        #
        self.renumber()
        return

    def nameclash(self):
        #
        # This routine checks for atom/residue name clashes in self.atomlinedefs
        # and gives later entries the key 'ALT'
        # So 'A:231:CD' would become 'A:231:CD,ALT' if there was a 'A:231:CD' earlier
        # in the pdb file. Identical atom names in different molecules (separated by a 'TER')
        # are not affected since Protool does not allow identical chain identifiers.
        # If there is already an '*,ALT' then the name becomes '*,ALT:ALT' and so forth...
        #
        clashdict={}
        lines=self.atomlinedefs.keys()
        lines.sort()
        for line in lines:
            atom=self.atomlinedefs[line]
            uniqueid=atom['CHAINID']+':'+atom['RESNUM']+':'+str(atom['ATNAME'])
            clashdict,newuniqueid=self.enternewname(clashdict,uniqueid,atom)
        return clashdict

    def enternewname(self,clashdict,uniqueid,atom):
        #
        # inserts uniqueid:atom in clashdict. If clashdict already contains
        # uniqueid, then uniqueid gets ',ALT' appended
        #
        import string
        if clashdict.has_key(uniqueid):
            #
            # If the residue names are identical, or if the last three letters are identical
            # [like in AGLN], then the add the ALT tag
            # If the resnames are different, then we assume that it is an insertion, and add
            # an 'A' or a 'B' to the residue number
            # This should be written more robustly so that insertations that are identical in
            # type to the previous residue are handled correctly...
            #
            if atom['RESNAME']==clashdict[uniqueid]['RESNAME'] or (len(atom['RESNAME'])==4 and atom['RESNAME'][1:]==clashdict[uniqueid]['RESNAME'][1:]):
                if uniqueid[-4:]==',ALT':
                    uniqueid=uniqueid+':ALT'
                else:
                    uniqueid=uniqueid+',ALT'
                clashdict,uniqueid=self.enternewname(clashdict,uniqueid,atom)
                return clashdict,uniqueid
            else:
                resnum=string.split(uniqueid,':')[1]
                if resnum==resnum[:-1]+'A':
                    resnum=resnum[:-1]+'B'
                    uniqueid=string.split(uniqueid,':')[0]+':'+resnum+':'+string.join(string.split(uniqueid,':')[2:],':')
                else:
                    resnum=resnum+'A'
                    uniqueid=string.split(uniqueid,':')[0]+':'+resnum+':'+string.join(string.split(uniqueid,':')[2:],':')
                clashdict,uniqueid=self.enternewname(clashdict,uniqueid,atom)
                return clashdict,uniqueid
        else:
            clashdict[uniqueid]=atom
            return clashdict,uniqueid
    

    def parseline(self,line):
        import string
        #
        # This function calls other routines to parse the line
        #
        type=string.strip(string.upper(line[:6]))
        #
        # There are so many junk headers, so we simplify a bit
        #
        if self.junkheaders.has_key(type):
            type=self.junkheaders[type]
        #
        # Should we ignore this model?
        #
        if self.readmodel is False and type!='ENDMDL':
            return
        #
        # Now call the appropritate function
        #
        if type=='ATOM':
            self.parseatom(line)
        elif type=='HETATM':
            self.parsehet(line)
        elif type=='HEADER':
            self.parseheader(line)
        elif type=='COMPND':
            self.parsecompnd(line)
        elif type=='SOURCE':
            self.parsesource(line)
        elif type=='AUTHOR':
            self.parseauthor(line)
        elif type=='REVDAT':
            self.parserevdat(line)
        elif type=='REMARK':
            self.parseremark(line)
        elif type=='SEQRES':
            self.parseseqres(line)
        elif type=='HELIX':
            self.parsehelix(line)
        elif type=='SHEET':
            self.parsesheet(line)
        elif type=='TURN':
            self.parseturn(line)
        elif type=='SSBOND':
            self.parsessbond(line)
        elif type=='CRYST1':
            self.parsecryst(line)
        elif type=='ORIGX1' or type=='ORIGX2' or type=='ORIGX3':
            self.parseorig(line)
        elif type=='SCALE1' or type=='SCALE2' or type=='SCALE3':
            self.parsescale(line)
        elif type=='CONECT':
            self.parseconnect(line)
        elif type=='MASTER':
            self.parsemaster(line)
        elif type=='TER':
            self.parseter(line)
        elif type=='JRNL':
            self.parsejournal(line)
        elif type=='FTNOTE':
            self.parsefootnote(line)
        elif type=='HET':
            self.parsehetheader(line)
        elif type=='FORMUL':
            self.parseformul(line)
        elif type=='MODEL':
            self.parsemodel(line)
        elif type=='ENDMDL':
            self.parseendmodel(line)
        elif type=='END':
            self.parseend(line)
        elif type=='EXPDTA':
            self.parseexpdta(line)
        elif type=='JUNK':
            pass
        else:
            raise ParseLineError, 'Unknown Header type: %s' %type
        return

    #
    # ----
    #
        
    def parseatom(self,line,hetatm=None):
        #
        # Parse a single ATOM line
        #
        import typecheck, string
        G=structure()
        #
        # This routine parses every atom line and stores it in self.atomlinedefs
        #
        self.linenumber=self.linenumber+1
        #
        # Check the atomnumber
        #
        atomnumber=line[6:11]
        if typecheck.isint(atomnumber):
            atomnumber=string.atoi(atomnumber)
        else:
            print '\nAtom number:',atomnumber
            raise ParseAtomError,'Non-int in atonnumber'
        #
        # All we now about the atomname is that it must contain a nondigit character
        #
        atomname=string.strip(line[12:16])
        if not typecheck.containsletter(atomname):
            raise ParseAtomError,'Atom name does not contain a character: %s,\nLINE: %s' %(atomname,line)
        #
        # Also the residuename must contain a letter
        #
        residuename=string.strip(line[17:21])
        if not typecheck.containsletter(residuename):
            print '\nResidue name:',residuename
            raise ParseAtomError,'Residuename does not contain a character'
        # ---------------------/------------------------------------------
        #
        # The ChainID must be a letter, a blank or a digit
        #
        if self.NMRmodel is False or self.readmodels==1:
            chainid=line[21]
        else:
            chainid='%d_%s' %(self.NMRmodel,line[21])
        #
        # Check the chainID
        #
        if not typecheck.containsletter(chainid) and chainid!=' ' and not typecheck.isint(chainid):
            print '\nChainID',chainid
            raise ParseAtomError,'Chain ID must be a letter,blank or digit'
        else:
            chainid=string.strip(chainid)
            #print 'Chainid: %s' %chainid
        newchainid=chainid
        #
        # Check that we are still building the same chain
        #
        if not self.chains.has_key(chainid):
            for cid in self.chains.keys():
                if self.chains[cid]=='building' or (self.oldchainid!=newchainid and self.chains[cid]=='building - dummy'):
                    self.chains[cid]='terminated'
            self.chains[chainid]='building'
        else:
            if self.chains[chainid]!='building':
                fail=1
                if chainid=='':
                    #
                    # Find a new chainid
                    #
                    for cid in self.chains.keys():
                        if self.chains[cid]=='building - dummy':
                            chainid=cid
                            fail=None
                            break
                    if chainid=='':
                        for cid in self.useids:
                            if not self.chains.has_key(cid):
                                self.chains[cid]='building - dummy'
                                chainid=cid
                                fail=None
                                break
                if fail:
                    for cid in self.chains.keys():
                        if self.chains[cid]=='building - dummy' and self.oldchainid==newchainid:
                            chainid=cid
                            fail=None
                            break
                    if fail:
                        for cid in self.useids:
                            if not self.chains.has_key(cid):
                                self.chains[cid]='building - dummy'
                                chainid=cid
                                fail=None
                                break
                    if fail:
                        print self.chains.keys()
                        raise ProtoolError, 'I ran out of alternative chain ID names'
        self.oldchainid=newchainid
        # -------------------------------------------------------------------------------------
        #
        # Check the residuenumber
        #
        residuenumber=line[22:27]
        done=False
        suffix=''
        while not done:
            if typecheck.isint(residuenumber):
                residuenumber=string.zfill(string.atoi(residuenumber),G.length_of_residue_numbers)
                residuenumber=residuenumber+suffix
                done=True
            else:
                suffix=residuenumber[-1]+suffix
                residuenumber=residuenumber[:-1]
        #
        # Get and check the coordinates
        #
        if len(line)<54:
            print line
            xcoord=None
            ycoord=None
            zcoord=None
        else:
            xcoord=line[30:38]
            ycoord=line[38:46]
            zcoord=line[46:54]
            if not typecheck.isnumber(xcoord) or not typecheck.isnumber(ycoord) or not typecheck.isnumber(zcoord):
                print '\nX,Y,Z',xcoord,ycoord,zcoord
                raise ParseAtomError,'All coordinates must be numbers'
            else:
                xcoord=string.atof(xcoord)
                ycoord=string.atof(ycoord)
                zcoord=string.atof(zcoord)
        #
        # Get the occupancy
        #
        if len(line)<60:
            print '\nNo occupancy for this line'
            occupancy=None
        else:
            occupancy=line[55:60]
            if not typecheck.isnumber(occupancy):
                print '\nOccupancy:',occupancy
                raise ParseAtomError,'Occupancy is not a number'
            else:
                occupancy=string.atof(occupancy)
        #
        # Get the B-factor
        #
        if len(line)<66:
            print '\nNo B-factor for this atom',line
            bfactor=None
        else:
            bfactor=line[60:66]
            if not typecheck.isnumber(bfactor):
                print 'B-factor:',bfactor
                bfactor=None
            else:
                bfactor=string.atof(bfactor)
        #
        # Put the rest of the line in a junk array
        #
        if len(line)>66:
            junk=line[66:]
        else:
            junk=None
        #
        # Finally put everything in self.atomlinedefs
        #
        self.atomlinedefs[self.linenumber]={'NUMBER':atomnumber}
        self.atomlinedefs[self.linenumber]['ATNAME']=atomname
        self.atomlinedefs[self.linenumber]['CHAINID']=chainid
        self.atomlinedefs[self.linenumber]['RESNAME']=residuename
        self.atomlinedefs[self.linenumber]['RESNUM']=residuenumber
        if xcoord!=None:
            self.atomlinedefs[self.linenumber]['X']=xcoord
            self.atomlinedefs[self.linenumber]['Y']=ycoord
            self.atomlinedefs[self.linenumber]['Z']=zcoord
        if occupancy or occupancy==0.0:
            self.atomlinedefs[self.linenumber]['OCCUPANCY']=occupancy
        if bfactor or bfactor==0.0:
            self.atomlinedefs[self.linenumber]['B-FACTOR']=bfactor
        if hetatm:
            self.atomlinedefs[self.linenumber]['HETATM']=1
        if junk:
            #
            # This is a tag
            #
            import string
            self.atomlinedefs[self.linenumber]['tag']=string.strip(junk)
                
        #
        # What is this atom bound to?
        #
        self.atomlinedefs[self.linenumber]['BOUND_TO']=[]
        return

    def parsehet(self,line):
        self.parseatom(line,1)
        return

    def parseheader(self,line):
        return

    def parsecompnd(self,line):
        return
    
    def parsesource(self,line):
        return

    def parseauthor(self,line):
        return

    def parserevdat(self,line):
        return

    def parseremark(self,line):
        return

    def parsescale(self,line):
        return

    def parsejournal(self,line):
        return

    def parseseqres(self,line):
        return

    def parsehelix(self,line):
        return
    
    def parsesheet(self,line):
        return

    def parseturn(self,line):
        return

    def parsessbond(self,line):
        return

    def parsecryst(self,line):
        return

    def parseorig(self,line):
        return

    def parseconnect(self,line):
        return

    def parsemaster(self,line):
        return

    def parseter(self,line):
        #
        # Molecule stops here
        #
        for cid in self.chains.keys() :
            if (self.chains[cid]=='building' or self.chains[cid]=='building - dummy') and self.parse_terms:
                self.chains[cid]='terminated'
                break
        return

    def parsejournal(self,line):
        return

    def parsefootnote(self,line):
        return

    def parsehetheader(self,line):
        return

    def parseformul(self,line):
        return

    def parseend(self,line):
        return

    def parseexpdta(self,line):
        return

    def parsemodel(self,line):
        self.NMRmodel=int(line.strip().split()[1])
        if self.NMRmodel>self.readmodels:
            self.readmodel=False
        else:
            self.readmodel=True
        return

    def parseendmodel(self,line):
        self.NMRmodel=False
        self.readmodel=True
        return

    def renumber(self):
        #
        # Sort all atoms and renumber them
        #
        atoms=self.atoms.keys()
        atoms.sort()
        count=1
        for atom in atoms:
            self.atoms[atom]['NUMBER']=count
            count=count+1
        return
    

#---------------------------------------------------------------------------------------------------
#--------------------- General class for structure I/O ---------------------------------------------
#---------------------------------------------------------------------------------------------------
class structureIO(structure):


    def readpdb(self,filename=None,data=None,parse=1):
        """
        # Reads a pdb file and lets PDBparser.parse do the parsing
        """
        import os
        if data != None:
            import StringIO
            stream = StringIO.StringIO(data)            
            self.lines = stream.readlines()            
        elif os.path.isfile(filename):
            if not silent:
                print 'Reading: ',filename
            fd=open(filename)
            self.lines=fd.readlines()
            fd.close()
        else:
            raise FileNotFoundError(filename)            
        if parse:           
            Y=PDBparser(self.lines,self.parse_terms)
            Y.parse()
            self.atoms=Y.atoms
            self.Update()
            self.attribute=Y.attribute
        return    
        
    #
    # ----
    #

    def renumber_atoms(self):
        """Renumber all atoms based on their names"""
        count=0
        atoms=self.atoms.keys()
        atoms.sort()
        for atom in atoms:
            self.atoms[atom]['NUMBER']=count+1
            count=count+1
        return

    #
    # -------------------------------------------------------------------
    #

    def writepqr(self,filename):
        """Write a PQR file"""
        self.writepdb(filename,pqr=True)
        return

    #
    # ---
    #

    def writepdb(self,filename,header=None,nowrite=False,dont_write_HETATMS=None,TER=1,pqr=None,style=None):
        """
        # Writes a PDB file from the informations in self.atoms
        """
        #
        # Renumber all atoms first
        #
        self.renumber_atoms()
        #
        # Write the file
        #
        class virtual_file:

            def __init__(self):
                self.lines=[]

            def write(self,text):
                self.lines.append(text)

            def close(self):
                return
        #
        #
        #
        import string, os
        if nowrite:
            fd=virtual_file()
        else:
            fd=open(filename,'w')
        #
        # Continue
        #
        if header:
            fd.write(header)
        atoms={}
        for atom in self.atoms.keys():
            atoms[self.atoms[atom]['NUMBER']]=atom
        atomnums=atoms.keys()
        atomnums.sort()
        lastchainid=None
        for atomno in atomnums:
            atom=atoms[atomno]
            if self.atoms[atom].has_key('HETATM'):
                if self.atoms[atom]['HETATM']==1 and dont_write_HETATMS:
                    continue
            atomname=string.ljust(self.atoms[atom]['ATNAME'],4)
            if atomname[0] not in string.digits:
                atomname=' '+atomname[:3]
            if self.atoms[atom].has_key('OCCUPANCY'):
                occupancy='%5.2f' %self.atoms[atom]['OCCUPANCY']
            else:
                occupancy='     '
            #
            # Should we use Autodock or WHAT IF names
            #
            import string
            atomname=string.strip(atomname)
            if style:
                if style.lower()=='autodock':
                    if self.autodock_names.has_key(string.strip(atomname)):
                        atomname=self.autodock_names[string.strip(atomname)]
                elif style.lower()=='whatif':
                    if self.Cterminal(self.resnum(atom)):
                        if atomname=='O':
                            atomname="O'"
                        elif atomname=='OXT':
                            atomname="O''"
            atomname=string.ljust(atomname,4)
            if atomname[0] not in string.digits:
                atomname=' '+atomname[:3]
            #
            # Continue
            #
            if self.atoms[atom].has_key('B-FACTOR'):
                bfactor='%5.2f' %self.atoms[atom]['B-FACTOR']
            else:
                bfactor='     '
            #
            # New Chain?
            #
            if self.atoms[atom]['CHAINID']!=lastchainid and lastchainid:
                fd.write('TER\n')
            lastchainid=self.atoms[atom]['CHAINID']
            if not pqr:
                try:
                    resnum=int(self.atoms[atom]['RESNUM'])
                except:
                    resnum=self.atoms[atom]['RESNUM']
                if type(resnum)==type(3):
                    line='ATOM  %5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f %5s %5s' %(self.atoms[atom]['NUMBER'],atomname,self.atoms[atom]['RESNAME'],self.atoms[atom]['CHAINID'],resnum,self.atoms[atom]['X'],self.atoms[atom]['Y'],self.atoms[atom]['Z'],occupancy,bfactor)
                else:
                    suffix=resnum[-1]
                    resnum=int(resnum[:-1])
                    line='ATOM  %5i %4s %3s %1s%4i%1s   %8.3f%8.3f%8.3f %5s %5s' %(self.atoms[atom]['NUMBER'],atomname,self.atoms[atom]['RESNAME'],self.atoms[atom]['CHAINID'],resnum,suffix,self.atoms[atom]['X'],self.atoms[atom]['Y'],self.atoms[atom]['Z'],occupancy,bfactor)
                
            else:
                charge=self.atoms[atom]['CHARGE']
                radius=self.atoms[atom]['RADIUS']
                line='ATOM  %5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f %5s %5s' %(self.atoms[atom]['NUMBER'],atomname,self.atoms[atom]['RESNAME'],self.atoms[atom]['CHAINID'],string.atoi(self.atoms[atom]['RESNUM']),self.atoms[atom]['X'],self.atoms[atom]['Y'],self.atoms[atom]['Z'],charge,radius)
            #
            # Is there a tag?
            #
            if self.atoms[atom].has_key('tag'):
                line=line+'\t %s' %self.atoms[atom]['tag']
            fd.write(line+'\n')
        if TER:
            fd.write('TER\n')
        fd.close()
        #
        # If no write then we simply return the lines
        #
        if nowrite:
            lines=fd.lines
            return lines
        else:
            return
        return

    #
    # -----
    #

    def parsepdb(self,lines):
        """
         Parses a PDB file, which is contained in lines (a list of lines)
        """
        Y=PDBparser(lines,self.parse_terms)
        Y.parse()
        self.atoms=Y.atoms
        self.Update()
        self.attribute=Y.attribute        
        return

    #
    # ----
    #

    def addpdblines(self,lines):
        """Add the atoms in lines to the current set of atoms"""
        Y=PDBparser(lines,self.parse_terms)
        Y.parse()
        error=False
        for atom in Y.atoms.keys():
            if not self.atoms.has_key(atom):
                self.atoms[atom]=Y.atoms[atom]
            else:
                print 'Name clash for',atom
                error=True
        if not error:
            self.Update()
            self.attribute=Y.attribute        
            return
        else:
            print 'Protool does not allow name clashes'
            import os
            os._exit(0)


    #
    # -----
    #
        
    def Read_DSSP(self,filename):
        import string
        residues={}
        file=open(filename)
	while 1:
            line=file.readline()
            #print line
            if not line:
                file.close()
                raise EOFError,' not enough lines'
            if line[2:18]=='#  RESIDUE AA ST':
                break
        while 1:
            line=file.readline()
            if not line:
                break 
            x=string.find(string.strip(line[95:101]),'-')
            if x==-1 or x==0:
                if string.strip(line[95:101])!='':
                    phi=string.atof(line[95:101])
                    x=string.find(string.strip(line[101:107]),'-')
                    if x==-1 or x==0:
                        psi=string.atof(line[101:107])
                        chainid=string.strip(line[11])
                        restyp=line[13]
                        resnum=string.strip(line[6:10])
                        if resnum:
                            resnum=string.zfill(resnum,self.length_of_residue_numbers)
                            sec=line[16]
                            acc=string.atoi(line[35:38])
                            hyd1=line[40:75]
                            resid=chainid+':'+resnum
                            residues[resid]={'restyp':restyp,'phi':phi,'psi':psi,'SS':sec,'H-bond':hyd1,'ACC':acc}
        return residues
        
    def rundssp(self):
        #
        # Runs DSSP for the PDB file in this instance
        #
        if not silent:
            print 'Running DSSP ....',
        import tempfile, os
        dir=tempfile.mktemp()
        os.mkdir(dir)
        pdbname=os.path.join(dir,'temp.pdb')
        dsspname=os.path.join(dir,'temp.dssp')
        self.writepdb(pdbname)
        status=os.system(DSSP+' '+pdbname+' '+dsspname+' >/dev/null')
        dsspresidues=self.Read_DSSP(dsspname)
        for residue in self.residues.keys():
            if dsspresidues.has_key(residue):
                for uniqueid in self.residues[residue]:
                    self.atoms[uniqueid]['SS']=dsspresidues[residue]['SS']
                    self.atoms[uniqueid]['ACC']=dsspresidues[residue]['ACC']
                    #print dsspresidues[residue]['ACC']
            else:
                if not silent:
                    print 'DSSP file and PDB file do not agree...??'
                    print 'This residue was not found in the DSSP output: ',residue
                for uniqueid in self.residues[residue]:
                    self.atoms[uniqueid]['SS']=None
                    self.atoms[uniqueid]['ACC']=None
        #
        # Remove the temporary directory
        #
        import ostools
        ostools.rmr(dir)
        if not silent:
            print 'done.'
        return

#
# Fast PDB parsing for Kristin
#
import typecheck, string
class structureIO_fast(structure):

    def readpdb(self,filename,parse=1):
        #
        # Reads a pdb file and lets PDBparser.parse do the parsing
        #
        import os
        if os.path.isfile(filename):
            fd=open(filename)
            self.lines=fd.readlines()
            fd.close()
            if parse:
                Y=fast_PDBparser(self.lines,self.parse_terms)
                Y.parse()
                self.atoms=Y.atoms
                self.Update()
                self.attribute=Y.attribute
            return
        else:
            raise FileNotFoundError,filename

    def parsepdb(self,lines):
        Y=fast_PDBparser(lines,self.parse_terms)
        Y.parse()
        self.atoms=Y.atoms
        self.Update()
        self.attribute=Y.attribute
        return

    #
    # -------------------------------------------------------------------
    #

    def writepdb(self,filename,header=None,nowrite=None):
        #
        # Writes a PDB file from the informations in self.atoms
        #
        import string, os
        fd=open(filename,'w')
        if header:
            fd.write(header)
        atoms={}
        for atom in self.atoms.keys():
            atoms[self.atoms[atom]['NUMBER']]=atom
        atomnums=atoms.keys()
        atomnums.sort()
        lastchainid=None
        for atomno in atomnums:
            atom=atoms[atomno]
            atomname=string.ljust(self.atoms[atom]['ATNAME'],4)
            if atomname[0] not in string.digits:
                atomname=' '+atomname[:3]
            if self.atoms[atom].has_key('OCCUPANCY'):
                occupancy='%5.2f' %self.atoms[atom]['OCCUPANCY']
            else:
                occupancy='     '
            if self.atoms[atom].has_key('B-FACTOR'):
                bfactor='%5.2f' %self.atoms[atom]['B-FACTOR']
            else:
                bfactor='     '
            #
            # New Chain?
            #
            if self.atoms[atom]['CHAINID']!=lastchainid and lastchainid:
                fd.write('TER\n')
            lastchainid=self.atoms[atom]['CHAINID']
            line='ATOM  %5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f %5s %5s' %(self.atoms[atom]['NUMBER'],atomname,self.atoms[atom]['RESNAME'],self.atoms[atom]['CHAINID'],string.atoi(self.atoms[atom]['RESNUM']),self.atoms[atom]['X'],self.atoms[atom]['Y'],self.atoms[atom]['Z'],occupancy,bfactor)
            fd.write(line+'\n')
        fd.write('TER\n')
        fd.close()
        if nowrite:
            fd=open(filename)
            lines=fd.readlines()
            fd.close()
            os.unlink(filename)
            return lines
        else:
            return
        return


import string
class fast_PDBparser(PDBparser):

    def parseline(self,line):
        type=string.strip(string.upper(line[:6]))
        if type=='ATOM':
            self.parseatom(line)
        elif type=='HETATM':
            self.parsehet(line)
        else:
            pass
        return    

    #
    # ---------------------------------------
    #

    def parseatom(self,line,hetatm=None):
        self.linenumber=self.linenumber+1
        self.length_of_residue_numbers=4
        atomnumber=line[6:11]
        atomname=string.strip(line[12:16])
        residuename=string.strip(line[16:21])
        chainid=string.strip(line[21])
        newchainid=chainid
        if not self.chains.has_key(chainid):
            for cid in self.chains.keys():
                if self.chains[cid]=='building' or (self.oldchainid!=newchainid and self.chains[cid]=='building - dummy'):
                    self.chains[cid]='terminated'
            self.chains[chainid]='building'
        else:
            if self.chains[chainid]!='building':
                fail=1
                if chainid=='':
                    #
                    # Find a new chainid
                    #
                    for cid in self.chains.keys():
                        if self.chains[cid]=='building - dummy':
                            chainid=cid
                            fail=None
                            break
                    if chainid=='':
                        for cid in self.useids:
                            if not self.chains.has_key(cid):
                                self.chains[cid]='building - dummy'
                                chainid=cid
                                fail=None
                                break
                if fail:
                    for cid in self.chains.keys():
                        if self.chains[cid]=='building - dummy' and self.oldchainid==newchainid:
                            chainid=cid
                            fail=None
                            break
                    if fail:
                        for cid in self.useids:
                            if not self.chains.has_key(cid):
                                self.chains[cid]='building - dummy'
                                chainid=cid
                                fail=None
                                break
                    if fail:
                        raise ProtoolError, 'I ran out of alternative chain ID names'
        #
        # -----
        #
        self.oldchainid=newchainid
        residuenumber=string.zfill(string.atoi(line[22:27]),self.length_of_residue_numbers)
        xcoord=string.atof(line[30:38])
        ycoord=string.atof(line[38:46])
        zcoord=string.atof(line[46:54])
        self.atomlinedefs[self.linenumber]={'NUMBER':atomnumber}
        self.atomlinedefs[self.linenumber]['ATNAME']=atomname
        self.atomlinedefs[self.linenumber]['CHAINID']=chainid
        self.atomlinedefs[self.linenumber]['RESNAME']=residuename
        self.atomlinedefs[self.linenumber]['RESNUM']=residuenumber
        self.atomlinedefs[self.linenumber]['X']=xcoord
        self.atomlinedefs[self.linenumber]['Y']=ycoord
        self.atomlinedefs[self.linenumber]['Z']=zcoord
        if hetatm:
            self.atomlinedefs[self.linenumber]['HETATM']=1
        return
