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


class APBS_calc:

    def __init__(self, Protool_instance,pdie=8.0,sdie=80.0,pbetype='intene'):
        """Run apbs for this protein"""
        protein=Protool_instance
        #
        import inputgen_APBS
        igen = inputgen_APBS.inputGen(protein)
        #
        # For convenience
        #
        igen.pdie=pdie
        print 'Setting protein dielectric constant to ',igen.pdie
        igen.sdie=sdie
        #
        # Set the type
        #
        igen.set_type(pbetype)
        #
        # Save the input file
        #
        apbs_inputfile=igen.printInput()
        #
        # Run APBS
        #
        import apbs
        self.APBS=apbs.runAPBS()
        self.potentials = self.APBS.runAPBS(protein, apbs_inputfile)
        self.APBS.cleanup()
        self.APBS=None
        return


if __name__=='__main__':
    import Protool
    P=Protool.structureIO()
    P.readpdb('2lzt.pka.pdb')
    # Build hydrogens
    import Protool.hydrogens
    H=Protool.hydrogens.Hydrogens(P)
    H.build_all_HNs()
    #
    # Assign charges
    #
    import Protool.assign_parameters
    PQR=Protool.assign_parameters.assign_parameters(P)
    PQR.clear_all_charges_and_radii()
    PQR.set_all_radii()
    PQR.assign_charge(':0035','GLU-')
    X=APBS_calc(P)
