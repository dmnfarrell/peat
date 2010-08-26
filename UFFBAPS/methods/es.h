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

#ifndef ES_H
#define ES_H
#include "method.h"
#include <math.h>
#include "prepinreader.h"
#include "distances.h"
#include "generalparameters.h"
#include "bondmatrix.h"
#include "simple_parameters.h"

using namespace std;

class Es:  Method{

 public:
  Es();
  Es(FILE * reporter);
  Es(FILE * resultsFile, FILE * reporter, vector<Soup *> A, vector<Soup *> B);
  Es(FILE * resultsFile, FILE * reporter, vector<Soup *> A);
  ~Es();

  double calculate(Soup * A, Soup * B, double cutoff);
  double calculate(Soup * A, double cutoff);

 private:
  double calculate_es(Atom * atom1, Atom * atom2, double dist);
  double permittivity(Atom * atom1, Atom * atom2, double r);
  void set_parameters();
  void writeDescription(FILE * reporter);

  string parameter_type,permittivity_type;
  double constant, D, res, result, same_residue_factor, neighbour_residue_factor;
  int residue1, residue2,include_surf_factor,include_debye_huckel;

  Distances dist;
  BondMatrix BM;

  FILE * rep;

  Sigmoidal S;

};

#endif
