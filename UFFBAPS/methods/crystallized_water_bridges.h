/*
 #
 # UFFBAPS - Unified Force Field for Binding And Protein Stability
 # Copyright (C) 2010 Jens Erik Nielsen & Chresten Soendergaard
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

#ifndef CRYSTALLIZED_WATER_BRIDGES_H
#define CRYSTALLIZED_WATER_BRIDGES_H
#include "method.h"
#include "hbonds.h"
#include "distances.h"

using namespace std;

class Crystallized_Water_Bridges:  Method{

 public:
  Crystallized_Water_Bridges();
  Crystallized_Water_Bridges(FILE * reporter);
  Crystallized_Water_Bridges(FILE * resultsFile, FILE * reporter, vector<Soup *> A);
  Crystallized_Water_Bridges(FILE * resultsFile, FILE * reporter, vector<Soup *> A, vector<Soup *> B);

  ~Crystallized_Water_Bridges();

  float calculate(Soup * A);
  float calculate(Soup * A, Soup * B);

 private:
  float check_for_interal_water_bridge(Atom * oxygen);
  float check_for_water_bridge(Atom * oxygen);
  
  vector<Atom*> internal_dry_atoms, dry_atoms_A, dry_atoms_B;

  float same_residue_factor,neighbour_residue_factor,water_bridge_energy;

  
  Hbonds HB;
  Distances dist;

  void set_parameters();
  void writeDescription(FILE * reporter);
};

#endif
