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

#include "Model.h"

using namespace std;

double model_class::resolve_clashes() {
    //
    // Resolve all clashes in the soup
    //
    // First find all clashing residues
    //
    _P.update_BOXLJ(); // We must update the LJ boxes before calculating energies
    //
    // List clashing residues
    //
    vector<double> ene;
    for (unsigned int chain=0;chain<_P.chains.size();chain++) {
        for (unsigned int res=0;res<_P.chains[chain].residues.size();res++) {
            ene=_ENERGY.get_external_energy(chain,res);
            if (ene[0]>0.0) {
                printf ("Chain: %d, residue: %d, clash: %6.2f\n",chain,res,ene[0]);
            }
        }
    }
    return 0.0;
}
