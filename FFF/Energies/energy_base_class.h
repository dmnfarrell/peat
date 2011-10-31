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
#ifndef ENERGY_BASE_CLASS_H
#define ENERGY_BASE_CLASS_H
#include <math.h>
#include "matrix.h"
#include "atom_class.h"
#include "fff.h"

class energy_base_class {
public:
    energy_base_class(FFF& P) : _P(P) {
        printf ("Energy base initialized\n");
    };
    FFF& _P;
    //
    virtual double calculate_interaction(const atom_class& atom1,const atom_class& atom2)=0;
    virtual double get_energy() {
    return 1.0;
  }
    virtual double get_energy(int chainnumber, int residuenumber) {
    return 1.0;
    }; // Get energy for one residue

    double get_external_energy(int chainnumber, int resnumber) {
        // Get the energy of all interactions between the side chain atoms and everything else nearby
        double energy=0.0;
        for (unsigned int atom=0;atom<_P.chains[chainnumber].residues[resnumber].atoms.size();atom++) {
            if (not _P.chains[chainnumber].residues[resnumber].atoms[atom].is_backbone()) {
                vector<atom_class*> close_atoms=(*(_P.BOXLJ)).get_close_atoms(_P.chains[chainnumber].residues[resnumber].atoms[atom]);
                for (unsigned int count2=0;count2<close_atoms.size();count2++) {
                    atom_class* atom2=close_atoms[count2];
                    if (not both_in_same_residue(_P.chains[chainnumber].residues[resnumber].atoms[atom],*atom2)) {
                        energy=energy+calculate_interaction(_P.chains[chainnumber].residues[resnumber].atoms[atom],*atom2);
                    }
                }
            }
        }
        return energy;
        
    };

    bool check_indexes(int chainnumber,int resnumber) {
      // Check that we're not outside the number of chain and residues
      if (chainnumber>=static_cast<int>(_P.chains.size())) {
	printf ("Unknown chain: %d\n",chainnumber);
	return false;
      }
      if (resnumber>=static_cast<int>(_P.chains[chainnumber].residues.size())) {
	printf ("Unknown residue number %d in chain %d\n",resnumber,chainnumber);
	return false;
      }
      return true;
    }
    //
    //
    virtual double get_external_energy() {
      return 0.0;
    };
    //
    // Constants
    //
    //double kB;
};

//
// Distances
//
class Distances : public energy_base_class {
public:
    Distances(FFF& P) : energy_base_class(P) {
    };
    //
    double calculate_interaction(const atom_class& atom1,const atom_class& atom2) {
    //
    // Calculate the distance
    //
    return sqrt((atom1.x-atom2.x)*(atom1.x-atom2.x)+
		(atom1.y-atom2.y)*(atom1.y-atom2.y)+
		(atom1.z-atom2.z)*(atom1.z-atom2.z));
  }
};


#endif
