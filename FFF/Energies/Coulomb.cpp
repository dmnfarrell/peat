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
#include "Coulomb.h"

double Coulomb::calculate_interaction(const atom_class& atom1,const atom_class& atom2) {
    // Quick screen for a steric clash. Note that tolerance values are set high here
    // to quickly identify positions that should have their energy evaluated in more detail (with UFFBAPBS or
    // the internal FFF energy function)
    // If the two atoms are bonded then return 0.0
    if (atom1==atom2) {
        return 0.0;
    }
    if (_P.are_bonded(atom1,atom2)) { 
        return 0.0; 
    }
    if (_P.could_be_Hbonded(atom1,atom2)) {
        return 0.0;
    }
    double factor=1.0; //(1.0/6.242E18)*(1.0/6.242E18)/(4.0*pi)/1E-10;
    double distance=Dist(atom1,atom2);
    double coulomb_energy=factor*atom1.charge*atom2.charge/(distance);
    //atom1.print_detail();
    //atom2.print_detail();
    //printf ("Coulomb energy: %6.2f\n",coulomb_energy);
    return coulomb_energy;
}

double Coulomb::get_energy() {
    // Get the energy for all atoms in _P
    //Boxes BOXLJ(_P.all_atoms,6.0); // Make 6A boxes
    double energy=0.0;
    for (unsigned int atom1=0;atom1<_P.all_atoms.size();atom1++) {
        vector<atom_class*> close_atoms=(*(_P.BOXLJ)).get_close_atoms(_P.all_atoms[atom1]);
        for (unsigned int count2=0;count2<close_atoms.size();count2++) {
            atom_class* atom2=close_atoms[count2];
            energy=energy+calculate_interaction(*(_P.all_atoms[atom1]),*atom2);
        }
    }
    //printf ("SC energy: %5.2f\n",energy);
    return energy;
}

double Coulomb::get_energy(int chainnumber,int resnumber) {
    //Boxes BOXLJ(_P.all_atoms,6.0); // Make 6A boxes
    double energy=0.0;
    for (unsigned int atom=0;atom<_P.chains[chainnumber].residues[resnumber].atoms.size();atom++) {
        //_P.chains[chainnumber].residues[resnumber].atoms[atom].print_detail();
        vector<atom_class*> close_atoms=(*(_P.BOXLJ)).get_close_atoms(_P.chains[chainnumber].residues[resnumber].atoms[atom]);
        for (unsigned int count2=0;count2<close_atoms.size();count2++) {
            atom_class* atom2=close_atoms[count2];
            energy=energy+calculate_interaction(_P.chains[chainnumber].residues[resnumber].atoms[atom],*atom2);
        }
    }
    return energy;
}

/*
double Coulomb::get_external_energy(int chainnumber, int resnumber) {
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
}
*/


