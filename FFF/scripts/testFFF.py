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
        X=Protool.structureIO()
        aas=X.trueaminoacids.keys()
        #
        import FFF.FFFcontrol as FFFC
            
        Rotamerlib=FFFC.Rotamer_class('parameters/bbdep02.May.sortlib')
        self.FFF=FFFC.FFF()
        self.FFF.read_pdb('2oh4.pdb')
        self.FFF.write_pdb('teset1000.pdb')
        stop
        #lines=X.writepdb('junk',nowrite=True)
        #self.FFF=FFFcontrol.FFF()
        #self.FFF.parse_lines(lines)
        self.Model=FFFC.model_class(self.FFF,Rotamerlib,'parameters')
        #self.Model.repair_all()
        energies={}
        for count in range(15,350):
            if not energies.has_key(count):
                energies[count]={}
            for aa in aas:
                print 'Mutating residue %d to %s' %(count,aa)
                max_clash=1.0
                energy=self.Model.Mutate('E',str(count),aa,3,max_clash)
                energies[count][aa]=energy
                #print 'Mutate done'
                self.FFF.write_pdb('alascan/1atp_%d_%s.pdb' %(count+15,aa))
                #print 'PDB file written'
                self.Model.undo_mutation()
                #print 'Undid mutation'
        #
        # Analyse the energies
        #
        count=0
        for count in energies.keys():
            for aa in energies[count].keys():
                if energies[count][aa]<0.1:
                    count=count+1
        print '%d mutations modelled with energy of 0.1 or less' %count

if __name__=='__main__':
    X=FFF()
