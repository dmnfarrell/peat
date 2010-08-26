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
#ifndef SIMPLE_PARAMETERS_H
#define SIMPLE_PARAMETERS_H
#include "method.h"
#include "generalparameters.h"
#include <vector>

using namespace std;

class Simple_Parameters:  Method{

 public:
  Simple_Parameters();
  Simple_Parameters(FILE * reporter);
  Simple_Parameters(FILE * reporter, vector<Soup*> A);
  ~Simple_Parameters();

  
  void prepare(vector<Soup*> A);
  void prepare(Soup * A);



 private:
  void set_electrostatic_parameters(Soup * A);
  void set_dispersive_parameters(Soup * A);

  void set_electrostatic_parameters_for_protein_chain(ProteinChain * p);
  void set_electrostatic_parameters_for_molecule(Molecule * m);
  void set_electrostatic_parameters_for_water(Water * w);


  void set_parameters();
  void writeDescription(FILE * reporter);

  map<string, float> atomic_charges;
  
  Generalparameters GP;


};

#endif
