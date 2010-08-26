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

#ifndef GENERALPARAMETERS_H
#define GENERALPARAMETERS_H
#include "method.h"
#include <map>
#include "soupmanip.h"
#include "atom.h"
#include "prepinreader.h"
#include "molecule.h"
#include "correct_atomnames.h"

using namespace std;

class Generalparameters:  Method{

 public:
  Generalparameters(FILE * reporter, vector<Soup *> A);
  Generalparameters(FILE * reporter);
  Generalparameters();
  ~Generalparameters();

  void set_vdw_radii(vector<Atom *> all_atoms);

 private:

  void set_protein_charges(vector<ProteinChain *> P);
  void set_ion_charges(vector<Molecule *> all_atoms);
  void print_all(vector<Soup *> A);
  void check_parameters(vector<Atom *> all_atoms);


  void set_parameters();
  void apply_parameters(vector<Soup *> A);
  void writeDescription(FILE * reporter);

  FILE * rep;
  map<string,float> vdw_radius;
  map<string,float> ion_charge;

};

#endif
