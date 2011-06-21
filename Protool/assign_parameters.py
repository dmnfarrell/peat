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


class assign_parameters:

    def __init__(self,Protool_instance):
        self.PI=Protool_instance
        #
        # Read charge file
        #
        import os
        thisfile=__file__
        thisdir=os.path.split(thisfile)[0]
        filename=os.path.join(thisdir,'OPLS.CRG')
        #
        # Read the file
        #
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        self.charges={}
        count=0
        line=lines[count]
        while line:
            import string
            line=string.strip(line)
            if line=='*':
                count=count+1
                nline=string.strip(lines[count])
                new_entry=nline.split()[0]
                self.charges[new_entry]={}
                count=count+1
            elif line=='':
                count=count+1
            elif line=='#END':
                break
            else:
                sp=line.split()
                self.charges[new_entry][sp[0]]=float(sp[1])
                count=count+1
            #
            line=lines[count]
        return

    #
    # ----
    #

    def clear_all_charges_and_radii(self):
        """Set all charges to 0.0"""
        for atom in self.PI.atoms.keys():
            self.PI.atoms[atom]['CHARGE']=0.0
            self.PI.atoms[atom]['RADIUS']=0.0
        return

    #
    # ----
    #

    def set_all_radii(self):
        """Assign all radii"""
        for atom in self.PI.atoms.keys():
            element=self.PI.get_atom_element(atom)
            if self.PI.vdwr.has_key(element):
                self.PI.atoms[atom]['RADIUS']=self.PI.vdwr[element]
            else:
                print 'Cannot find',atom,element
                raise Exception('Cannot assign parameters. Values missing for %s %s' %(atom,element))
        return

    #
    # ----
    #

    def set_radius(self,atom,radius=1.0):
        """Set the radius for this atom"""
        self.PI.atoms[atom]['RADIUS']=radius
        return

    def set_charge(self,atom,charge=0.0):
        """Set the charge for this atom"""
        self.PI.atoms[atom]['CHARGE']=charge
        return

    #
    # ----
    #

    def assign_charge(self,residue,template):
        """Assign charges to a residue based on the template"""
        atoms=self.PI.residues[residue]
        template=self.charges[template]
        for atom in template.keys():
            PIatom='%s:%s' %(residue,atom)
            if self.PI.atoms.has_key(PIatom):
                self.PI.atoms[PIatom]['CHARGE']=template[atom]
            elif template[atom]==0.0:
                pass # We don't care about zero charges
            else:
                print 'Could not find',PIatom
                raise Exception, 'Could not assign charges'
                
