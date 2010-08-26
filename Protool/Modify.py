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

from errors import *

class modify:

    def delete(self,id):
        #
        # This function deletes the residue or atom in id
        #
        if self.ispresent(id):
            if self.isaa(id):
                for atom in self.residues[id]:
                    del self.atoms[atom]
                del self.residues[id]
            else:
                del self.atoms[atom]
        else:
            raise NotFoundError,id
        return
    
