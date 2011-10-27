/*
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
 */
#include "Access.h"


double Access::get_energy(int chainnumber,int resnumber) {
  // Count the number of atoms within 6A
  double interactions=0;
  for (unsigned int atom=0;atom<_P.chains[chainnumber].residues[resnumber].atoms.size();atom++) {
    vector<atom_class*> close_atoms=(*(_P.BOXLJ)).get_close_atoms(_P.chains[chainnumber].residues[resnumber].atoms[atom]);
    for (unsigned int count2=0;count2<close_atoms.size();count2++) {
      interactions=interactions+1.0;
    }
  }
  return interactions;
}



