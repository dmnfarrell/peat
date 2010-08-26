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

#ifndef CORRECT_ATOMNAMES_H
#define CORRECT_ATOMNAMES_H
#include "method.h"
#include "distances.h"


using namespace std;

class Correct_Atomnames:  Method{

 public:
  Correct_Atomnames(FILE * reporter, vector<Soup*> A);
  Correct_Atomnames();
  ~Correct_Atomnames();

  void go(vector<Soup*> A);

 private:
  void correct_heavy_atoms();
  void correct_hydrogen_atoms();
  vector<int> find_bounded_hydrogens(int a);
  int get_number_of_bounded_heavy_atoms(int a);
  void set_parameters();
  float get_distance(int i, int j); //to save memory instead of doing a distance matrix

  vector<Atom*> atoms;
  double SSdist,Hdist,defaultDist,squareSSdist,squareHdist,squaredefaultDist;

  void writeDescription(FILE * reporter);
};

#endif
