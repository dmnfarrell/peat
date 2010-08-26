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
#ifndef NN_DISCRIMINATION_PREPARATION_H
#define NN_DISCRIMINATION_PREPARATION_H

#include "ligand_entropy.h"

using namespace std;

/*
Preparation of inputs for nn_discrimination
-------------------------------------------

* Inputs:

Hbonds: Fraction of possible ligand donors/acceptors forming a hydrogen bond

Ligand entropy: fraction of rotatable ligand bonds blocked

Hydrophobicity: fraction of ligand hydrophobic atoms matched by protein hydrophobic atoms




*/





class NN_DISCRIMINATION_PREPARATION{

 public:
  NN_DISCRIMINATION_PREPARATION(FILE * resultsFile,FILE * reporter,, vector<Soup*> A, vector<Soup*> B);
  ~NN_DISCRIMINATION_PREPARATION();

  void go();

 private:
  ofstream out_file;


};


#endif
