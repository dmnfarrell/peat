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

#ifndef FRAGMENTAL_VOLUME_H
#define FRAGMENTAL_VOLUME_H
#include "method.h"
#include "bondmatrix.h"
#include <cstdlib> 
#include <map>
#include <sstream>
#include <math.h>

using namespace std;

class Fragmental_Volume:  Method{

 public:
  Fragmental_Volume();
  Fragmental_Volume(FILE * reporter, string mode);
  Fragmental_Volume(FILE * reporter, vector<Soup*> A, string mode);
  Fragmental_Volume(FILE * reporter, SoupObject* A, string mode);
  ~Fragmental_Volume();

  void set_mode(string mode);

  void calculate(vector<Soup*> A);
  void calculate(Soup * A);
  void calculate(SoupObject* A);

  void set_precalculated_fragmental_volumes(Soup * A);

  void calculate_protein_fraction_in_ellipse(Soup* A, int a, int b);
  float calculate_molecular_volume(Soup * A);


 private:
  void find_touching_atoms(int i);
  float atom_fragmental_volume(int i);
  float calculate_molecular_volume_by_integration();
  inline bool inside(float x, float y, float z, int i);

  void apply_precalculated_volumes(SoupObject *A);
  void set_bonding_atoms(int i);

  float atomic_radius(string element);

  void writeDescription(FILE * reporter);
  void set_parameters();

  // maps for storing radii
  map<string,float> radii;
  map<Atom*,float> atom_radii;

  BondMatrix BM; 

  int points_pr_cubic_A;
  float total_volume,accumulated_sphere_volume,integrated_total_volume,accumulated_group_fragmental_volume,sphere_volume;

  vector<Atom*> touching_atoms,bonding_atoms,atoms;
  vector<vector<int> > bonds;

  string operation_mode;

  map<string,float> precalculated_volumes;


  // ellipse mode stuff
  vector<int> focal_point_atoms1,focal_point_atoms2;
  float ellipse_radius_factor;

  FILE * out;

};

#endif
