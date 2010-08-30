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

#include "energy.h"

void energy_class::init_all() {
  printf ("-----------------------------\n");
  printf ("Entering energy initialisation\n");
  //
  // Get the number of atoms
  //
  int num_atoms=_P.all_atoms.size();
  //
  // Construct a matrix
  //
 // matrix_class<double> M(num_atoms,num_atoms);
  //
  // Construct all the distances
  //
  //Distances dists;
  //energy_function distance_matrix(_P,_dists,0);
  //printf ("Distances calculated\n");
  //
  // Now do all the real energy functions
  //
  //
  // Instantiate the LJ function
  //
  _energy_components.push_back(new Steric_Clash (_P));
  printf ("Steric clash function initialised\n");
  //
  // Coulomb energies
  //
  _energy_components.push_back(new Coulomb (_P));
  // 
  // Titration energies (intpKa based)
  //
  //_energy_components.push_back(new Titration_Energy (_P));
  //
  // Desolvation energies 
  //
  //_energy_components.push_back(new Desolvation (_P));
  //
  printf("Exiting energy initialisation\n");
}

//
// ----
//

vector<double> energy_class::get_energy() {
  vector<double> energies;
  for (unsigned int i=1;i<_energy_components.size();i++) {
    energies.push_back((*(_energy_components[i])).get_energy());
  }
    return return_energies(energies);
}

vector<double> energy_class::get_energy(int chainnumber,int resnumber) {
    vector<double> energies;
    for (unsigned int i=1;i<_energy_components.size();i++) {
    energies.push_back((*(_energy_components[i])).get_energy(chainnumber,resnumber));
    }
    return return_energies(energies);
}

vector<double> energy_class::get_external_energy(int chainnumber,int resnumber) {
   vector<double> energies;
    for (unsigned int i=0;i<_energy_components.size();i++) {
        energies.push_back((*(_energy_components[i])).get_external_energy(chainnumber,resnumber));
    }
    return return_energies(energies);
}

vector<double> energy_class::get_external_energy() {
    vector<double> energies;
    for (unsigned int i=0;i<_energy_components.size();i++) {
        energies.push_back((*(_energy_components[i])).get_external_energy());
    }
    return return_energies(energies);
}


vector<double> energy_class::return_energies(vector<double> energies) {
  // Return all the energies in the right form
  // [0]: Steric Clash
  // [1]: Sum of all energies 2- (exluding Steric Clash)
  // [2]: Coulomb
  vector<double> results;
  results.push_back(energies[0]);
  results.push_back(0.0);
  double total=0.0;
  for (unsigned int i=1;i<energies.size();i++) {
    total=total+energies[i];
    results.push_back(energies[i]);
  }
  results[1]=total;
  return results;
}
    


