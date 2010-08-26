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

#ifndef HBONDS_H
#define HBONDS_H
#include <math.h>
#include "method.h"
#include "distances.h"
#include "topology.h"
#include "bondmatrix.h"
#include "bfactor.h"
#include "Boxes.h"

using namespace std;

class Hbonds:  Method{
  friend class Water_Bridges;

 public:
  Hbonds();
  Hbonds(FILE * reporter);
  Hbonds(FILE * resultsFile,FILE * reporter, vector<Soup*> A, vector<Soup*> B);
  Hbonds(FILE * resultsFile,FILE * reporter, vector<Soup*> A);
  ~Hbonds();
 
  float calculate(Soup * A, Soup * B);
  float calculate(Soup * A);
  float hbond_factor(Atom * acceptor,int ai, Atom * donor,int bi, double dist);
  bool is_hbond(Atom * acceptro, Atom * donor, double dist);
  bool is_potential_hbond(Atom * acceptor, Atom * donor);
  bool water_bridge_hbond(Atom * acceptor, Atom * donor, double dist); //used by the water bridge method
  bool is_potential_acceptor(Atom * acceptor);
  bool is_potential_donor(Atom * donor);


  
  void set_parameters();
  
 protected:
  vector<Atom *> atomsA, atomsB;  
  vector<vector<int> > bondsA,bondsB;



 private:
  void writeDescription(FILE * reporter);

  bool number_of_hbonds_ok(Atom * acceptor, Atom * donor);

  float distance_factor(float distance);
  float angle_factors(int a, int b);
  float angle_factors_internal(int a, int b);

  float angle_factor(float angle);
  float find_angle(vector<float> p,vector<float> q,vector<float> r);
  float find_lenght(vector<float> a);
  vector<float> find_root_A(int a);
  vector<float> find_root_B(int b);
  int number_bounded_hydrogens(Atom * a);
  void register_bond(Atom * acceptor, Atom * donor);
  void make_list_of_number_of_bonded_hydrogens();


  float min_dist, max_dist,piConv, replacement_by_water_factor;
  int include_surf_factor,include_dist_factor,include_angle_factor;

  map<Atom*,Atom*> bonds; //list of all found hbonds bonds[donor] = acceptor
  map<Atom*,int> donor_count, acceptor_count, no_bounded_hydrogens;

  Distances dist;

  map<string,int> no_lone_pairs;
  vector <string> possible_acceptors;
  vector <string> possible_donors;
  vector<vector<float> > distances;  
  vector<float> ligand_surf_dists,protein_surf_dists;
  FILE * rep;

};

#endif
