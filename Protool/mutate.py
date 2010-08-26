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

import numpy
#
# Import the error definitions
#
from errors import *


class Mutate:

    def __init__(self,Protool=None,onlydefs=None,max_bump=0.5):
        """Init the mutate function"""
        if Protool:
            self.PI=Protool
        #
        # Read definition files
        #
        print 'Initialising Protool mutate class ...'
        self.max_tolerated_bump=max_bump
        print 'Setting max tolerated bump to %5.3f' %self.max_tolerated_bump
        self.read_aa_defs()
        if onlydefs:
            print 'done'
            return
        import rotamer_lib
        self.rots=rotamer_lib.rots
        print 'Done'
        return

    #
    # ----
    #

    def new_PDB(self,Protool):
        """Set a new Protool instance"""
        self.PI=Protool
        return

    #
    # -----
    #

    def read_aa_defs(self):
        """Read the amino acid defininition file"""
        import Protool, os
        location=os.path.split(Protool.__file__)[0]
        location=os.path.join(location,'AA.DAT')
        fd=open(location)
        line=fd.readline()
        #
        #
        #
        self.aadefs={}
        while line:
            import string
            line=string.strip(line)
            if line[:2]=='//':
                line=fd.readline()
                continue
            #
            # Start of record
            #
            if string.strip(line)=='*END':
                break
            elif line[0]=='*':
                line=string.strip(fd.readline()).split()
                name=line[0]
                number=line[1]
                line=string.strip(fd.readline()).split()
                atoms=[]
                count=0
                while line[0][0] not in string.digits:
                    count=count+1
                    atname=line[0]
                    x=float(line[1])
                    y=float(line[2])
                    z=float(line[3])
                    atoms.append([atname,numpy.array([x,y,z]),count])
                    line=string.strip(fd.readline()).split()
                self.aadefs[name]={'atoms':atoms[:]}
                bonds=[]
                for count in range(0,len(line),2):
                    bonds.append([int(line[count]),int(line[count+1])])
                self.aadefs[name]['bonds']=bonds
                #
                # Process torsions
                #
                line=string.strip(fd.readline()).split()
                torsions=[]
                if len(line)>1:
                    for count in range(0,len(line)-1,4):
                        torsions.append([int(line[count+1]),int(line[count+2]),int(line[count+3]),int(line[count+4])])
                self.aadefs[name]['torsions']=torsions[:]
            line=fd.readline()
        fd.close()
        return

    #
    # ----
    #

    def get_nonH_torsions(self,residue):
        """Get the number of non-hydrogen torsion angles"""
        restype=self.PI.resname(residue)
        torsions=self.aadefs[restype]['torsions']
        count=0
        for torsion in torsions:
            fail=False
            for atom_num in torsion:
                atom_name=self.aadefs[restype]['atoms'][atom_num-1][0]
                if self.PI.atomname_is_hydrogen(atom_name):
                    fail=True
                    break
            if fail:
                break
            count=count+1
        return count

    #
    # -----
    #

    def get_torsion_atoms(self,residue,chinum):
        """Get the atoms that define a given torsional angle
        The first torsion angle is always 1"""
        restype=self.PI.resname(residue)
        #
        # Search for the atoms
        #
        chinum=chinum-1
        torsions=self.aadefs[restype]['torsions']
        this_torsion=torsions[chinum]
        atoms=[]
        for atom_num in this_torsion:
            atom_name=self.aadefs[restype]['atoms'][atom_num-1][0]
            if self.PI.atomname_is_hydrogen(atom_name):
                raise HydrogenInTorsionError
            for atom in self.PI.residues[residue]:
                if self.PI.atname(atom)==atom_name:
                    atoms.append(atom)
                    break
        return atoms


    #
    # -----
    #

    def get_chi(self,residue,chinum):
        """Get the torsion angle for a specific residue"""
        atoms=self.get_torsion_atoms(residue,chinum)
        import phipsi
        return phipsi.GetDihedral(self.PI,atoms[0],atoms[1],atoms[2],atoms[3])

    #
    # -----
    #

    def set_chi(self,residue,chinum,angle):
        """Set the chi angle of a residue to angle
        The first chiangle in a residue is number 1
        """
        oldchi=self.get_chi(residue,chinum)
        difchi=angle-oldchi
        #
        # Get the aa def
        #
        tors_atoms=self.get_torsion_atoms(residue,chinum)
        vector=self.PI.GetPosition(tors_atoms[2])-self.PI.GetPosition(tors_atoms[1])
        #
        # Get the atoms that should be moved
        #
        moveatoms=[tors_atoms[3]]
        stop_atoms=[tors_atoms[2],tors_atoms[1]]
        added=1 # To get the loop started
        while added:
            added=None
            for atom in moveatoms+[tors_atoms[2]]:
                for bound_atom in self.getbound(atom):
                    if not bound_atom in moveatoms and not bound_atom in stop_atoms:
                        moveatoms.append(bound_atom)
                        added=1
        movecoords=[]
        #print 'Moving these atoms',moveatoms
        for atom in moveatoms:
            movecoords.append(self.PI.GetPosition(atom)-self.PI.GetPosition(tors_atoms[1]))
        import quatfit
        new_coords=quatfit.qchichange(vector,movecoords,difchi)
        count=0
        for array in new_coords:
            atom=moveatoms[count]
            self.PI.atoms[atom]['X']=array[0]+self.PI.GetPosition(tors_atoms[1])[0]
            self.PI.atoms[atom]['Y']=array[1]+self.PI.GetPosition(tors_atoms[1])[1]
            self.PI.atoms[atom]['Z']=array[2]+self.PI.GetPosition(tors_atoms[1])[2]
            count=count+1
        return

    #
    # ----
    #

    def getbound(self,atom):
        """Get all atom names that are bound to this atom"""
        chainid=self.PI.chainid(atom)
        restyp=self.PI.resname(atom)
        resnum=self.PI.resnum(atom)
        atname=self.PI.atname(atom)
        bonds=self.aadefs[restyp]['bonds']
        atoms=self.aadefs[restyp]['atoms']
        atnumber=None
        for local_atom in atoms:
            if local_atom[0]==atname:
                atnumber=local_atom[2]
        if not atnumber:
            if self.PI.atomname_is_hydrogen(local_atom[0]):
                pass
            else:
                print 'No atnumber',local_atom
                raise Exception,"huh?"
        bonded=[]
        for b1,b2 in bonds:
            if b1==atnumber:
                bonded.append(atoms[b2-1][0])
            elif b2==atnumber:
                bonded.append(atoms[b1-1][0])
        #
        # If the atom is an N or a C, then add the appropriate atom from the previous or next residue
        #
        if atname=='N':
            try:
                prevres=self.PI.PreviousResidue(resnum)
                C_name='%s:C' %prevres
                if self.PI.atoms.has_key(C_name):
                    bonded.append(C_name)
            except Nterm:
                pass
        elif atname=='C':
            try:
                nextres=self.PI.NextResidue(resnum)
                N_name='%s:N' %nextres
                if self.PI.atoms.has_key(N_name):
                    bonded.append(N_name)
            except Cterm:
                pass
                #print 'Cterm',resnum
        #
        # Convert names
        #
        rbonded=[]
        for atom in bonded:
            if not self.PI.atoms.has_key(atom):
                for real_atom in self.PI.residues[resnum]:
                    if self.PI.atname(real_atom)==atom:
                        rbonded.append(real_atom)
                        break
            else:
                rbonded.append(atom)
        return rbonded
    
    #
    # -----
    #   

    def get_bb_rotamers(self,residue):
        """Find the subset of rotamers that we can use for this bb conformation"""
        resname=self.PI.resname(residue)
        if not self.rots.has_key(resname):
            return []
        #
        # Yes, found
        #
        use_rots=[]
        rots=self.rots[resname]
        import phipsi
        phi,psi,omega=phipsi.GetPhiPsi(self.PI,residue)
        crit=0
        for comp in [phi,psi]:
            if comp:
                crit=crit+1

        for rot in rots:
            ok=0
            for lib,real in [[rot['phi'],phi],[rot['psi'],psi]]:
                if lib and real:
                    if abs(lib-real)<=10.0:
                        ok=ok+1
                if ok>=crit:
                    print crit,len(use_rots)
                    use_rots.append(rot)
        use_rots=use_rots[:300]
        for rot in use_rots:
            print rot
        return use_rots
    #
    # ----
    #
    def new_mutation(self):
        """Reinitialise self.mutate_operations"""
        print 'Initialising mutate_operations'
        self.mutate_operations={}
        return


    def mutate(self,residue,newtype,orgtype=None):
        """Mutate the residue to the new type
        Normal substitutions: A:0035,GLN,GLU
        Deletions: A:0035,delete,GLU"""
        #
        # Reset score
        #
        score=None
        #
        # Store the operations
        #
        if not getattr(self,'mutate_operations',None):
            self.mutate_operations=[]
        this_operation=[]
        local_operation=[]
        #
        # Get phi and psi for this residue
        #
        #import phipsi
        #phi,psi,omega=phipsi.GetPhiPsi(self.PI,residue)
        newtype=newtype.upper()
        aas=self.rots.keys()
        aas.sort()
        if not self.aadefs.has_key(newtype):
            raise Exception,"invalid aa type"
        #
        if not self.PI.residues.has_key(residue):
            res=self.PI.residues.keys()
            res.sort()
            raise Exception,"Invalid residue number: %s\nKnown residues are:\n%s" %(residue,str(res))
        #
        # Special cases
        #
        oldtype=self.PI.resname(residue)
        #
        # If orgtype was specified then check that there's agreement with the pdb file
        #
        if orgtype:
            if orgtype!=oldtype:
                raise Exception('Incorrect original residue type: %s %s. You said: %s' %(residue,oldtype,orgtype))
        #
        #
        #
        if newtype=='ALA' and (oldtype!='GLY' and oldtype!='PRO'):
            this_operation=self.fast_mutatetoALA(residue)
            score=0.0
        elif newtype=='GLY':
            this_operation=self.fast_mutatetoGLY(residue)
            score=0.0
        elif oldtype=='ASP' and newtype=='ASN':
            #
            # Just rename
            #
            for atom in self.PI.residues[residue]:
                this_operation.append(['delete',atom,self.PI.atoms[atom]])
                self.PI.atoms[atom]['RESNAME']='ASN'
                if self.PI.atname(atom)=='OD2':
                    self.PI.atoms[atom]['ATNAME']='ND2'
                this_operation.append(['add',atom,self.PI.atoms[atom]])
            score=0.0
        elif oldtype=='ASN' and newtype=='ASP':
            #
            # Just rename
            #
            for atom in self.PI.residues[residue]:
                this_operation.append(['delete',atom,self.PI.atoms[atom]])
                self.PI.atoms[atom]['RESNAME']='ASP'
                if self.PI.atname(atom)=='ND2':
                    self.PI.atoms[atom]['ATNAME']='OD2'
                this_operation.append(['add',atom,self.PI.atoms[atom]])
            score=0.0
        elif oldtype=='GLU' and newtype=='GLN':
            #
            # Just rename
            #
            for atom in self.PI.residues[residue]:
                this_operation.append(['delete',atom,self.PI.atoms[atom]])
                self.PI.atoms[atom]['RESNAME']='GLN'
                if self.PI.atname(atom)=='OE2':
                    self.PI.atoms[atom]['ATNAME']='NE2'
                this_operation.append(['add',atom,self.PI.atoms[atom]])
            score=0.0
        elif oldtype=='GLN' and newtype=='GLU':
            #
            # Just rename
            #
            for atom in self.PI.residues[residue]:
                this_operation.append(['delete',atom,self.PI.atoms[atom]])
                self.PI.atoms[atom]['RESNAME']='GLU'
                if self.PI.atname(atom)=='NE2':
                    self.PI.atoms[atom]['ATNAME']='OE2'
                this_operation.append(['add',atom,self.PI.atoms[atom]])
            score=0.0
        elif oldtype=='TYR' and newtype=='PHE':
            #
            # Delete OH
            #
            for atom in self.PI.residues[residue]:
                if self.PI.atname(atom)=='OH':
                    this_operation.append(['delete',atom,self.PI.atoms[atom]])
            score=0.0
        else:
            #
            # Use generic code
            #
            #
            # Get the standard conformation
            #
            std_conf=self.aadefs[newtype]
            #
            # See which kind of mutation we do
            # 1. Just put standard conformation
            # 2. Standard conf + rotamer search
            # 3. Use wt rotamer (if possible, otherwise use best rotamer)
            # All can be used with and without debumping.
            #
            #
            # Get the torsion angles for this residue
            #
            org_torsangles=[]
            number_of_chis=self.get_nonH_torsions(residue)
            #print 'Number of chi values for %s is %d' %(residue,number_of_chis)
            for chinum in range(1,number_of_chis+1):
                org_torsangles.append(self.get_chi(residue,chinum))
            #
            # Get the backbone atoms from the structure
            #
            bb_struct=[]
            ref_coords=[]
            for atom in self.PI.residues[residue]:
                if (self.PI.is_backbone(atom) or self.PI.is_CB(atom)) and not self.PI.is_hydrogen(atom):
                    bb_struct.append(atom)
                    ref_coords.append(self.PI.GetPosition(atom))
            #
            # Get the corresponding atoms from the template
            #
            bb_template=[]
            fit_coords=[]
            for s_atom in bb_struct:
                s_atname=self.PI.atname(s_atom)
                for atom in std_conf['atoms']:
                    if atom[0]==s_atname:
                        bb_template.append(atom[0])
                        fit_coords.append(atom[1])
                        break
            #
            # Get rotation and translation for bb -> bb
            #
            import quatfit
            refcenter,fitcenter,rotation=quatfit.qfit(len(ref_coords),ref_coords,fit_coords)
            #
            # Get the atoms that should be transformed
            #
            trans_atoms=[]
            trans_coords=[]
            for atom in std_conf['atoms']:
                if not atom in bb_template and atom[0]!='O':
                    trans_atoms.append(atom[0])
                    trans_coords.append(atom[1])
            #
            # Apply rotation and translation to trans_coords
            #
            newcoords = quatfit.qtransform(len(trans_coords), trans_coords, refcenter, fitcenter, rotation)
            #
            # Delete the old side chain atoms and rename the ones that stay
            #
            for atom in self.PI.residues[residue]:
                if atom.split(':')[-1] in ['OXT','O2',"O''"]:
                    # Rename the OXT
                    local_operation.append(['delete',atom,self.PI.atoms[atom]])
                    self.PI.atoms[atom]['RESNAME']=newtype
                    local_operation.append(['add',atom,self.PI.atoms[atom]])
                elif not atom in bb_struct:
                    this_operation.append(['delete',atom,self.PI.atoms[atom]])
                    self.PI.remove_atom(atom)
                else:
                    local_operation.append(['delete',atom,self.PI.atoms[atom]])
                    self.PI.atoms[atom]['RESNAME']=newtype
                    local_operation.append(['add',atom,self.PI.atoms[atom]])
            #
            # Copy coordinates to Protool array
            #
            struct_atom=self.PI.atoms[bb_struct[0]]
            resnumber=struct_atom['RESNUM']
            chainid=struct_atom['CHAINID']
            atomnumber=struct_atom['NUMBER']
            count=0
            for atom in trans_atoms:
                if atom[0]!='H':
                    uniqueid='%s:%s:%s' %(chainid,resnumber,atom)
                    if not self.PI.atoms.has_key(uniqueid):
                        self.PI.add_atom(uniqueid,atomnumber=atomnumber,
                                         atomname=atom,chainid=chainid,residuename=newtype,
                                         residuenumber=resnumber,
                                         xcoord=newcoords[count][0],ycoord=newcoords[count][1],zcoord=newcoords[count][2])
                        this_operation.append(['add',uniqueid,self.PI.atoms[uniqueid]])
                count=count+1
            #
            # Fix atom numbers
            #
            self.PI.renumber_atoms()
            #
            # Set the chiangles just as the parent residue
            #
            new_number_of_chis=self.get_nonH_torsions(residue)
            count=1
            for chi in org_torsangles:
                if count>new_number_of_chis:
                    break
                self.set_chi(residue,count,chi)
                count=count+1
            #
            # Score this rotamer, if bump is zero then keep it
            #
            score=self.get_bump_value(residue)
            if (score is None or score>0.0) and newtype!='ALA':
                #
                # Do rotamer search
                #
                score=self.get_best_rotamer(residue)
                if score is None:
                    return None
        # ----------------------------------------------------
        #
        # Keep track of the operations
        #
        updated_op=[]
        for op,uniqueid,atom in this_operation:
            if op=='delete':
                updated_op.append([op,uniqueid,{}])
            else:
                added_atom=self.PI.atoms[uniqueid].copy()
                not_needed=['CHAINID','tag','BOUND_TO','RESNUM']
                for delete in not_needed:
                    if added_atom.has_key(delete):
                        del added_atom[delete]
                updated_op.append([op,uniqueid,self.PI.atoms[uniqueid]])
        self.mutate_operations=self.mutate_operations+updated_op
        #
        # Score is bump score
        #
        return score
            
    #
    # ----
    #

    def fast_mutatetoALA(self,residue):
        """Mutate to Ala by removing all atoms that are not backbone or CB"""
        this_operation=[]
        for atom in self.PI.residues[residue]:
            this_operation.append(['delete',atom,self.PI.atoms[atom]])
            if not self.PI.is_backbone(atom) and not self.PI.is_CB(atom):
                self.PI.remove_atom(atom)
            else:
                self.PI.atoms[atom]['RESNAME']='ALA'
                this_operation.append(['add',atom,self.PI.atoms[atom]])
        return this_operation

    #
    # ---
    #

    def fast_mutatetoGLY(self,residue):
        """Mutate to Gly by removing all side chain atoms"""
        this_operation=[]
        for atom in self.PI.residues[residue]:
            this_operation.append(['delete',atom,self.PI.atoms[atom]])
            if not self.PI.is_backbone(atom):
                self.PI.remove_atom(atom)
            else:
                self.PI.atoms[atom]['RESNAME']='GLY'
                this_operation.append(['add',atom,self.PI.atoms[atom]])
        return this_operation

    #
    # -----
    #

    def undo_last_mutate_operation(self):
        """Undo the last mutate operation"""
        raise Exception('Not implemented yet')
        return

    #
    # ----
    #

    def get_number_of_chis(self,residue):
        """Get the number of chis for this residue"""
        pass
        return
    

    #
    # ----
    #

    def get_best_rotamer(self,residue):
        """Find the best rotamer in the library"""
        step=0
        close_atoms=self.find_debump_atoms(residue,30.0)
        scores=[]
        #
        # Get the list of rotamers for this bb conformation
        #
        rotamer_list=self.get_bb_rotamers(residue)
        #
        # Score all rotamers
        #
        number_of_chis=len(self.aadefs[self.PI.resname(residue)]['torsions'])
        scores=[]
        count=0
        text='%5d' %len(rotamer_list)
        print text,
        import sys
        sys.stdout.flush()
        for rotamer in rotamer_list:
            text='\b\b\b\b\b\b%5d' %(len(rotamer_list)-count)
            print text,
            sys.stdout.flush()
            #
            # Set this rotamer
            #
            for chi in rotamer.keys():
                if chi[:3]!='chi':
                    continue
                number=int(chi[3])
                if number>number_of_chis:
                    continue
                try:
                    self.set_chi(residue,number,rotamer[chi])
                except AtomNotFoundError:
                    print 'Wrong chi?',rotamer[chi],chi,residue                    
                except HydrogenInTorsionError:
                    #
                    # We don't set torsions with hydrogens
                    #
                    pass
            #
            # All angles set
            #
            score=self.find_bumps(residue,close_atoms,use_backbone=None) #dont debump backbone atoms of this residue
            scores.append([score,rotamer]) 
            if score==0.0:
                # At the moment we cannot get better than 0.0, so we just pick the first rotamer with 0.0
                new_score=0.0
                break
            count=count+1
        scores.sort()
        found_rotamer=False
        for score,rotamer in scores:
            if not score is None:
                self.use_rotamer(residue,rotamer,number_of_chis)
                new_score=self.get_bump_value(residue,use_backbone=None,selection_radius=30.0)
                found_rotamer=True
                break
        if not found_rotamer:
            new_score=False
        return new_score

    #
    # ----
    #

    def use_rotamer(self,residue,rotamer,number_of_chis):
        """Set this rotamer"""
        for chi in rotamer.keys():
            if chi[:3]!='chi':
                continue
            number=int(chi[3])
            if number>number_of_chis:
                continue
            try:
                self.set_chi(residue,number,rotamer[chi])
            except AtomNotFoundError:
                print 'Wrong chi?',rotamer[chi],chi,residue
            except HydrogenInTorsionError:
                #
                # We don't set torsions with hydrogens
                #
                pass
        return

    #
    # ----
    #

    def get_bump_value(self,residue,use_backbone=1,selection_radius=10.0):
        """Get the bump value for the residue"""
        close_atoms=self.find_debump_atoms(residue,selection_radius)
        #print 'close atoms',close_atoms
        bump_val=self.find_bumps(residue,close_atoms,use_backbone=use_backbone)
        return bump_val

    #
    # ----
    #

    def find_debump_atoms(self,residue,radius):
        """Get atoms close to the residue"""
        #close_atoms=[]
        #for atom in self.PI.residues[residue]:
        #    for catom in self.PI.get_neighbour_atoms(atom):
        #        if not catom in close_atoms:
        #            close_atoms.append(catom)
        #return close_atoms
        close_atoms=[]
        if self.PI.resname(residue)!='GLY':
            CB=residue+':CB'
        else:
            # Use CA instead
            CB=residue+':CA'
        for atom in self.PI.atoms.keys():
            if self.PI.dist(CB,atom)<radius:
                if not atom in self.PI.residues[residue]:
                    close_atoms.append(atom)
        #print 'I am using %d debump atoms' %(len(close_atoms))
        return close_atoms

    #
    # ----
    #

    def show_all_bumps(self):
        """Loop over all residues and find all bumps in the PDB file"""
        residues=self.PI.residues.keys()
        residues.sort()
        for residue in residues:
            close_atoms=self.find_debump_atoms(residue,10)
            self.find_bumps(residue,close_atoms)
        return

    #
    # ----
    #

    def find_bumps(self,residue,close_atoms,use_backbone=1):
        """Find bumps between the residue and the close atoms

        If use_backbone is 1, then also count the bumps that the backbone atoms of residue make.
        If use_backbone is None then only debump the side chain
        """
        bump_val=0.0
        for atom1 in self.PI.residues[residue]:
            if not use_backbone:
                if self.PI.is_backbone(atom1):
                    continue
            bound_atoms=self.getbound(atom1)
            #print 'Bound to',atom1,bound_atoms
            for atom2 in close_atoms:
                if not atom2 in bound_atoms:
                    bump_val=bump_val+self.get_overlap(atom1,atom2)
                    #print 'Overlap',atom1,atom2,self.get_overlap(atom1,atom2)
                else:
                    bump_val=bump_val
            if bump_val>self.max_tolerated_bump:
                return None
        return bump_val

    #
    # ----
    #

    def are_bound(self,atom1,atom2):
        """Find out if two atoms are covalently bound"""
        if atom2 in self.getbound(atom1):
            return True
        return False

    #
    # ----
    #

    def are_1_3_bound(self,atom1,atom2):
        """1-3 bound atoms?"""
        for atom_next in self.getbound(atom1):
            if self.are_bound(atom_next,atom2):
                return True
        return False

    #
    # ----
    #

    def are_1_4_bound(self,atom1,atom2):
        """Find out if two atoms are 1-4 bound"""
        for atom_next in self.getbound(atom1):
            if self.are_1_3_bound(atom_next,atom2):
                return True
        return False

    #
    # ----
    #

    def get_overlap(self,atom1,atom2):
        """Find out if two atoms overlap"""
        vdw1=self.PI.vdwr[self.PI.get_atom_element(atom1)]
        vdw2=self.PI.vdwr[self.PI.get_atom_element(atom2)]
        dist=self.PI.dist(atom1,atom2)
        sub_dist=dist-vdw1-vdw2
        tolerance=0.30 #Changed from 0.6
        
        if sub_dist+tolerance<0.0:
            #
            # don't count S-S bridges
            #
            if self.PI.atname(atom1)=='SG' and self.PI.atname(atom2)=='SG' and self.PI.resname(atom1)=='CYS' and self.PI.resname(atom2)=='CYS':
                return 0.0
            #
            # Check for 1-3 and 1-4 bound atoms
            #
            if self.are_1_3_bound(atom1,atom2):
                tolerance=1.3
            elif self.are_1_4_bound(atom1,atom2):
                tolerance=0.8
            if abs(sub_dist+tolerance)<0.1:
                print 'BUMP',atom1,atom2,dist,vdw1,vdw2,abs(sub_dist+tolerance)
                return abs(sub_dist+tolerance)
            else:
                return 0.0
        else:
            return 0.0



#
# -----
# 
# Convenient mutate functions


#
# -----
#

def Model_Mutations(pdbfile,mol2files,mutations,max_overlap=0.5,max_totalbump=1.0,return_score=False):
    """Model a number of mutations in a pdbfile when one or more ligands are present"""
    #
    # Check for stupidity
    #
    if max_overlap>max_totalbump:
        max_totalbump=max_overlap
        print 'Adjusted total bump cutoff to %5.2f' %max_totalbump
    #
    # Read PDB file
    #
    import Protool
    P=Protool.structureIO()
    P.readpdb(pdbfile)
    P.remove_all_hydrogens()
    #
    # Read mol2 file
    #
    L=Protool.ligand(P)
    for mol2file in mol2files:
        print "Added %s with tag 'LIGAND'" %mol2file
        L.readmol2(mol2file,tag='LIGAND')
    #
    # Get the pdb lines
    #
    pdblines=P.writepdb('junk.pdb',nowrite=True)
    #
    # Pass the lines to FFF
    #
    import FFF.FFFcontrol
    myFFF=FFF.FFFcontrol.FFF()
    myFFF.parse_lines(pdblines)
    #myFFF.soup_stat()
    Model=FFF.FFFcontrol.model_class(myFFF,Rotamerlib,FFFaadef_dir)
    #
    import pKa.pKD_tools as pKD_tools
    total_bump=0.0
    for mutation in mutations:
        resid=pKD_tools.get_resid_from_mut(mutation)
        chainid=resid.split(':')[0]
        resid=resid.split(':')[1]
        #
        # Get rid of the leading zeros
        #
        done=False
        while not done:
            if resid[0]=='0' and len(resid)>1:
                resid=resid[1:]
            else:
                done=True
        #
        newres=pKD_tools.get_newrestyp_from_mut(mutation)
        oldres=pKD_tools.get_oldrestyp_from_mut(mutation)
        opttype=3 # Rotamer library
        energies=Model.Mutate(chainid,resid,newres,opttype,max_overlap)
        bump_score=energies[0]
        Sum=energies[1]
        Coulomb=energies[2]
        #
        total_bump=total_bump+bump_score
        print 'Bump score: %5.2f, total bump: %5.2f' %(bump_score,total_bump)
        if bump_score>max_overlap or total_bump>max_totalbump:
            print 'Cannot model this set of mutations - too many bumps'
            if return_score:
                return False, total_bump
            else:
                return False
        print 'Bump score for %s: %5.3f' %(mutation,bump_score)
    print 'Total bump score for all mutations: %5.3f' %(total_bump)
    class FFF_fix:

        def __init__(self,FFF):
            self.PI=FFF
            return
    #
    # Return
    #
    if (return_score):
        return FFF_fix(myFFF),total_bump
    return FFF_fix(myFFF)

#
# ----
#

def Model_Mutations_old(pdbfile,mol2files,mutations,max_overlap=0.5,return_score=False):
    """Model a number of mutations in a pdbfile when one or more ligands are present"""
    #
    # Initialise mutate routines
    #
    MUT=Mutate(max_bump=max_overlap)
    #
    # Read PDB file
    #
    import Protool
    P=Protool.structureIO()
    P.readpdb(pdbfile)
    P.remove_all_hydrogens()
    #
    # Read mol2 file
    #
    L=Protool.ligand(P)
    for mol2file in mol2files:
        print "Added %s with tag 'LIGAND'" %mol2file
        L.readmol2(mol2file,tag='LIGAND')
    #
    # Pass combined pdb file to mutate routines and mutate
    #
    MUT.new_PDB(P)
    import pKa.pKD_tools as pKD_tools
    total_bump=0.0
    #
    # Model
    #
    for mutation in mutations:
        #
        # Get info
        #
        resid=pKD_tools.get_resid_from_mut(mutation)
        newres=pKD_tools.get_newrestyp_from_mut(mutation)
        oldres=pKD_tools.get_oldrestyp_from_mut(mutation)
        bump_score=MUT.mutate(resid,newres,orgtype=oldres)
        if bump_score is None or bump_score is False or bump_score>max_overlap:
            print 'Cannot model this set of mutations - too many bumps'
            return False,20.0
        print 'Bump score for %s: %5.3f' %(mutation,bump_score)
        total_bump=total_bump+bump_score
    print 'Total bump score for all mutations: %5.3f' %(bump_score)
    if return_score:
        return MUT,bump_score
    return MUT

#
# -----
#

def test_function(options):
    #X=Mutate()
    import os
    mutations=options.mutations
    pdbfile=options.pdbfile
    if not os.path.isfile(pdbfile):
        raise Exception('PDB file not found: %s' %pdbfile)
    MUT,bs=Model_Mutations(pdbfile,[],mutations,return_score=True,max_overlap=options.bump,max_totalbump=options.totalbump)
    if MUT:
        MUT.PI.writepdb(options.outfile)
    else:
        print 'Cannot model mutations'
    return
    #
    # Test new score
    #
    import Protool
    X=Protool.structureIO()
    X.readpdb('1atp.pdb')
    res=X.residues.keys()
    res.sort()
    aas=X.trueaminoacids.keys()
    x=[]
    y=[]
    import random
    bettermodel=[]
    while len(x)<2000:
        try:
            resi=random.choice(res)
            aa=random.choice(aas)
            mutation='%s:%s:%s' %(resi,X.resname(resi),aa)
            print 'Calculating for %s' %mutation
            oldscore=None
            newscore=None
            for function in ['Model_Mutations','Model_Mutations_old']:
                resultfile=os.path.join(os.getcwd(),'scores/'+mutation+str(function))
                if not os.path.isfile(resultfile):
                    M,score=eval(function)('1atp.pdb',[],[mutation],return_score=True)
                    fd=open(resultfile,'w')
                    import pickle
                    A=score
                    pickle.dump(A,fd)
                    fd.close()
                else:
                    import pickle
                    fd=open(resultfile)
                    score=pickle.load(fd)
                    fd.close()
                #
                # Assign scores to the right variables
                #
                if resultfile[-4:]=='_old':
                    oldscore=score
                else:
                    newscore=score
            if oldscore>10:
                bettermodel.append(newscore)
                continue
            x.append(oldscore)
            y.append(newscore)
        except:
            pass
    print 'In %d cases FFF was able to construct a model where Protool was not' %len(bettermodel)
    print bettermodel
    import pylab
    pylab.scatter(x,y)
    pylab.show()
    return

#
# Read the rotamer library for FFF
#
import Protool, os, cPickle
#
# Get the FFF aadef file
#
try:
    import FFF.FFFcontrol
    FFFaadef_file=os.path.split(FFF.FFFcontrol.__file__)[0]
    FFFaadef_dir=os.path.join(FFFaadef_file,'parameters')
    #
    location=os.path.join(FFFaadef_dir,'bbdep02.May.sortlib')
    text='Importing rotamer library from %s.....' %location
    Rotamerlib=FFF.FFFcontrol.Rotamer_class(location)
except:
    print 'FFF not available'

if __name__=='__main__':
    print
    print 'Modelling of a single mutant protein containing 1 or more mutations'
    print
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]',version='%prog 1.0')
    parser.add_option('-p','--pdb',dest='pdbfile',action='store',type='string',default='2lzt.pka.pdb',
                      help='The PDB file. Default: %default')
    parser.add_option('-o','--outfile',dest='outfile',action='store',type='string',default='2lzt.mut.pdb',
                      help='The output (mutated) PDB file. Default: %default')
    
    parser.add_option('-m','--mutations',type='string',dest='mutations',action='append',
                      help="Mutations to use model. Format: A:0001:ASP:ALA",default=[])
    parser.add_option('-b','--bump',type='float',dest='bump',action='store',
                      help='Maximum bump value for a single mutation. Default: %default',default=0.5)
    parser.add_option('-t','--totalbump',type='float',dest='totalbump',action='store',
                      help='Maximum total bump value for all mutations. Default: %default',default=1.0)
    (options, args) = parser.parse_args()
    test_function(options)
    #print "Don't run this module directly"
    
    
