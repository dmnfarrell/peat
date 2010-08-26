#ifndef ATOM_H
#define ATOM_H
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


#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>
#include <ctype.h>
#include "stringtools.h"
#include "algebra.h"

using namespace std;

class Atom {

 public:
  Atom(string line);
  Atom();
  Atom(Vector v, int residue_number);

  double x, y, z;
  
  string name,element,residue,adv_type,sybyl_type,generic_key;
  double occupancy, tempFactor;
  int residue_no;

  //force field parameters
  double vdw_dist, vdw_energy, charge, fragmental_volume, group_fragmental_volume, desolvation_factor, occ;
  bool is_hydrophobic, is_hydrogen, hydrogen_bond_could_be_replaced_by_water;

  vector<Atom*> bonded_atoms, hydrogen_bonds;

  void display();
  Vector get_pos();

  // Box mode
  int bx,by,bz;

 private:


    
};
#endif
