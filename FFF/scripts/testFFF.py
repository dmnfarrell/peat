#!/usr/bin/env python
#
# FFF - Flexible Force Field
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
class FFF:

    def __init__(self):
        import Protool
        self.PI=Protool.structureIO()
        self.aas=self.PI.trueaminoacids.keys()
        self.PI.readpdb('test.pdb')
        #
        import FFF.FFFcontrol as FFFC
        import os, sys
        scriptdir=os.path.split(os.path.abspath(__file__))[0]
        FFFdir=os.path.split(scriptdir)[0]
        Rotamerlib=FFFC.Rotamer_class(os.path.join(FFFdir,'parameters/small_lib'))
        self.FFF=FFFC.FFF()
        self.FFF.read_pdb('test.pdb')
        #self.Model=FFFC.pKa_class(self.FFF,Rotamerlib,os.path.join(FFFdir,'parameters'))
        self.Model=FFFC.model_class(self.FFF,Rotamerlib,os.path.join(FFFdir,'parameters'))
        #
        # Test mutations
        #
        self.mutate_test()
        #self.Model.repair_all()
        #
        # Build all hydrogens - standard protonation state
        #
        #self.Model.build_hydrogens()
        #self.FFF.write_pqr('2lzt.pqr.pdb')
        return
        


    def mutate_test(self):
        """Full mutation scan test"""
        energies={}
        import os
        if not os.path.isdir('alascan'):
            os.mkdir('alascan')
        for resid in (self.PI.residues.keys()):
            if not energies.has_key(resid):
                energies[resid]={}
            for aa in self.aas:
                print 'Mutating residue %s to %s' %(resid,aa)
                max_clash=1.0
                energy=self.Model.Mutate(resid,aa,3,max_clash)
                energies[resid][aa]=energy
                print 'Acc for %s is %5.3f' %(resid,self.Model.get_accessibility(resid))
                print 'Mutate done'
                self.FFF.write_pdb('alascan/2lzt_%d_%s.pdb' %(count,aa))
                print 'PDB file written'
                self.Model.undo_mutation()
                print 'Undid mutation'
        #
        # Analyse the energies
        #
        count=0
        for count in energies.keys():
            for aa in energies[count].keys():
                if energies[count][aa]<0.1:
                    count=count+1
        print '%d mutations modelled with energy of 0.1 or less' %count

        return

if __name__=='__main__':
    X=FFF()
