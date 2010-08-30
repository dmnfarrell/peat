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

#include "Titration.h"

double Titration_Energy::calculate_interaction(const atom_class& atom1,const atom_class& atom2) {
    // Not defined here
    printf ("Titration Energy has no atom-atom interaction");
    throw Titration_Error();
}

double Titration_Energy::get_titgroup_energy(TitGroup& thisgroup) {
    //
    // Get the energy of this titgroup at this pH
    //
    return 2.3*(_P._pH-thisgroup.intpKa)*thisgroup.acidbase;
}

double Titration_Energy::get_energy() {
    // Get the energy for all atoms in _P
    double energy=0.0;
    for (unsigned int group=0;group<_P._titratable_groups.size();group++) {
        energy=energy+get_titgroup_energy(_P._titratable_groups[group]);
    }
    return energy;
}






