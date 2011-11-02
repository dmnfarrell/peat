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
  // Count the number of atoms within  10A
  double interactions=0.0;
  // Set the counter to False everywhere
   for (unsigned int atom1=0;atom1<_P.all_atoms.size();atom1++) {
     (*_P.all_atoms[atom1]).counted=false;
   }
  //
   double countatoms=0;
  if (check_indexes(chainnumber,resnumber)) {
    for (unsigned int atom=0;atom<_P.chains[chainnumber].residues[resnumber].atoms.size();atom++) {
      atom_class* atom1=&(_P.chains[chainnumber].residues[resnumber].atoms[atom]);
      //if ((*atom1).is_backbone() || (*atom1).is_hydrogen()) {
      if ((*atom1).is_hydrogen()) {
	continue;
      }
      countatoms++;
      vector<atom_class*> close_atoms=(*(_P.BOX10A)).get_close_atoms(atom1);
      for (unsigned int count2=0;count2<close_atoms.size();count2++) {
	atom_class* atom2=close_atoms[count2];
	if ((*atom2).inresidue!=resnumber || (*atom2).inchain!=chainnumber) {
	  if (not (*atom2).is_hydrogen() && not (*atom2).counted) {
	    if (Dist(*atom1,*atom2)<=6.0 && Dist(*atom1,*atom2)>1.5) {
	      interactions=interactions+1.0;
	      (*atom2).counted=true;
	    }
	  }
	}
      }
    }
    interactions=max(0.0,80-interactions)*42.0/55.0;
    return interactions;
  }
  else {
    throw exception();
  }
  return 0.0;
}



