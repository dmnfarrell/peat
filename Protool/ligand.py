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

import string, copy

class ligand:
    """Read a mol2 file"""
    
    def __init__(self,PI):
        self.PI=PI
        return

    #
    # ----
    #
    
    def readmol2(self,filename,tag=''):
        """
        Routines for reading MOL2 file
        """
        self.tag=tag
        self.filename = filename
        data = open(self.filename).read()
        # ATOM section
        start = data.find("@<TRIPOS>ATOM")
        stop = data.find("@<TRIPOS>BOND")
        atoms = data[start+14:stop-2].split("\n")
        # BOND section
        start = data.find("@<TRIPOS>BOND")
        stop = data.find("@<TRIPOS>SUBSTRUCTURE")
        bonds = data[start+14:stop-1].split("\n")
        self.parse_mol2lines(atoms)
        #self.parseBonds(bonds)
        #self.createlBondedAtoms()
        return

    #
    # ----
    #

    def parse_mol2lines(self,lines):
        """This parses a single ligand and inserts it as a new chain in Protool"""
        #
        # Find the chain ID
        #
        fail=True
        for cid in ['L','I','J','K','L','M','O','P','Q','R','S','T','U','V','W','X','Y','Z']:
            if not self.PI.chains.has_key(cid):
                chainid=cid
                fail=False
                break
        if fail:
            print 'I ran out of alternative chain ID names'
            raise Exception()
        #
        # We got the chain ID
        # Residue number is 1
        #
        shortened_message={}
        resnumber='0001'
        for atomline in lines:
            sline=atomline.split()
            if sline==[]:
                continue
            atomnumber=0
            atomname=sline[1]
            x=float(sline[2])
            y=float(sline[3])
            z=float(sline[4])
            resname=sline[7]
            if len(resname)>3:
                if not shortened_message.has_key(resname):
                    print 'Shortened residue name from "%s" to "%s"' %(resname,resname[:3])
                    shortened_message[resname]=1
                resname=resname[:3]

            sybyltype=sline[5]
            charge=sline[-1]
            uniqueid='%s:%s:%s' %(cid,resnumber,atomname)
            self.PI.add_atom(uniqueid,
                             atomnumber=atomnumber,atomname=atomname,
                             chainid=chainid,residuename=resname,residuenumber=resnumber,
                             xcoord=x,ycoord=y,zcoord=z,update=None,
                             BFACTOR=0.0,OCCUPANCY=1.0,CHARGE=charge,RADIUS=None,tag=self.tag,accept_duplicate=True)
        self.PI.Update()
        return
    
        
    def parseBonds(self,BondList):
        """
        for parsing @<TRIPOS>BOND
        """
        for BondLine in BondList:
            SeparatedBondLine = BondLine.split()
            thisBond = MOL2BOND(
                int(SeparatedBondLine[1]), # bond frm
                int(SeparatedBondLine[2]), # bond to
                SeparatedBondLine[3],      # bond type
                int(SeparatedBondLine[0])  # bond id
                )
            self.lBonds.append(thisBond)

    #
    # ----
    #

    def createlBondedAtoms(self):
        """
        Creates for each atom a list of the bonded Atoms
        
        This becomes one attribute of MOL2ATOM!
        """
        for bond in self.lBonds:
            self.lAtoms[bond.frm-1].lBondedAtoms.append(
                self.lAtoms[bond.to-1])

            self.lAtoms[bond.to-1].lBondedAtoms.append(
                self.lAtoms[bond.frm-1])

            atbond = copy.deepcopy(bond)
            atbond.other_atom=self.lAtoms[bond.to-1]
            self.lAtoms[bond.frm-1].lBonds.append(atbond)

            atbond = copy.deepcopy(bond)
            atbond.other_atom=self.lAtoms[bond.frm-1]
            self.lAtoms[bond.to-1].lBonds.append(atbond)
        return

    #
    # ----
    #
    
    def createPDBlineFromMOL2(self):
        FakeType = "HETATM"
        return ('%s%5i%5s%4s%2s%5s   %8.3f%8.3f%8.3f\n' %
                (FakeType,            self.serial,   self.name,
                 self.resName, ' L',          self.resSeq, 
                 self.x,self.y, self.z)) 


