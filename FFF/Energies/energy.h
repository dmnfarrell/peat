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
#ifndef ENERGY_CLASS_H
#define ENERGY_CLASS_H
#include "fff.h"
#include "matrix.h"
#include "Lennard_Jones.h"
#include "Coulomb.h"
#include "energy_base_class.h"
#include "Titration.h"

class energy_class {
 public:
  energy_class(FFF &P): _P(P) {
    init_all();
  }
  // Functions
  void init_all();
  vector<double> get_energy();
    vector<double> get_energy(int chainnumber, int resnumber);
    vector<double> get_external_energy(int chainnumber,int resnumber);
  vector<double> get_external_energy();
    vector<double> return_energies(vector<double> energies);
  // Variables
  FFF& _P;
    //Distances _dists(_P);
  vector<energy_base_class*> _energy_components;
};

#endif
