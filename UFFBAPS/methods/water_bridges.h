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
#ifndef WATER_BRIDGES_H
#define WATER_BRIDGES_H
#include "method.h"
#include "hbonds.h"
#include "distances.h"
#include "bondmatrix.h"
#include "write_file.h"
#include <math.h>
#include <deque>
//#include "algebra.h"


#include <sstream>

using namespace std;

typedef vector<int> index_list;


struct potential_water_struct
{
  Atom water_atom;
  bool active;
  vector<Atom*> contacts;

  potential_water_struct(Atom water_in)
  {
    water_atom = water_in;
    active = true;
  }

  void add_contact(Atom * contact)
  {
    contacts.push_back(contact);
  }

};


class Water_Bridges:  Method{

 public:
  Water_Bridges();
  Water_Bridges(FILE * reporter);
  Water_Bridges(FILE * resultsFile, FILE * reporter, vector<Soup *> A);
  Water_Bridges(FILE * resultsFile, FILE * reporter, vector<Soup *> A, vector<Soup *> B);

  ~Water_Bridges();

  float calculate(Soup * A);
  float calculate(Soup * A, Soup * B);
  void find_water_contacts(Soup * A);

 private:

  void place_water_contacts();
  void place_potential_water_contacts(Atom* atom);
  void make_boxes();

  Vector find_root(Atom * atom);
  void write_water_atoms(string name);

  void remove_overlaps_with_protein(vector<int> indices, vector<int> dry_neighbour_indices);
  void remove_water_bumps(vector<int> indices, vector<int> neighbour_indices);
  void merge_waters(vector<int> indices, vector<int> neighbour_indices);
  void merge(potential_water_struct * potential_water1, potential_water_struct * potential_water2);

  bool no_heavy_atoms_bonded(Atom * atom);

  vector<Atom*> dry_atoms, dry_atoms_A, dry_atoms_B;

  float same_residue_factor,neighbour_residue_factor,water_bridge_energy,
    optimal_hb_dist,optimal_hb_dist_sq,protein_bump_dist,
    protein_bump_dist_sq,water_join_dist,water_join_dist_sq,min_hb_dist_sq,
    min_hb_dist, water_desolvation_penalty, double_optimal_hb_dist_sq;

  vector<float> optimal_hb_angles,d_along,d_perp, steps;

  Hbonds HB;
  Distances dist;
  BondMatrix BM;
  
  int include_surf_factor;
  
  void set_parameters();
  void writeDescription(FILE * reporter);

  vector<potential_water_struct> potential_waters;
  vector<index_list> potential_water_boxes, potential_water_boxes_incl_neighbours, atom_boxes, atom_boxes_incl_neighbours;

  


};

#endif
