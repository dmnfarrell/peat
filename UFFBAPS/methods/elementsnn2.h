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

#ifndef ELEMENTSNN2_H
#define ELEMENTSNN2_H
#include "nnmodel.h"
#include "bondmatrix.h"
#include "topology.h"
#include "hbonds.h"

using namespace std;

struct pose;

class Elementsnn2: public Nnmodel{

 public:
  Elementsnn2(FILE * reporter);
  ~Elementsnn2();

  void create_network();
  void build_inputs(pose * A);
  void define_input();
  void test_pose(pose * A);
  void train_pose(pose * A);
  void update_weights();
  void write_network(string name);
  void read_network(string name);

  int get_no_epochs();
  int get_max_epochs_since_record();
  int get_smoothing();
  float get_target_factor();
  string get_dataset();

  float get_train_pct();
  float get_validate_pct();
  float get_test_pct();

  bool get_include_co_factors();
  bool get_include_water();


 private:
  vector<float> get_ligand_freqs(int p);
  vector<float> get_protein_freqs(int p);

  int no_protein_types, no_ligand_types,no_inputs,no_protein_bonds_incl,no_ligand_bonds_incl, smoothing;

  float learning_rate, threshold, linear_coef;

  NN network;


  vector<int> topology;
  vector<int> tf;

  map<string, vector<float> > protein_types;
  map<string, vector<float> > ligand_types;

  vector<vector<int> > protein_bonds;
  vector<vector<int> > ligand_bonds;

  vector<float> ligand_surf_dists;
  vector<float> protein_surf_dists;

  vector<Atom *> protein_atoms;
  vector<Atom *> ligand_atoms;

  vector<string> include_elements_p;
  vector<string> ignore_elements_p;

  vector<string> include_elements_l;
  vector<string> ignore_elements_l;

  Hbonds HB;

  int no_epochs, max_epochs_since_record, include_co_factors, include_water,include_surf_dist,include_hb_info;
  float target_factor, train_pct, validate_pct, test_pct;
  string dataset;

  vector<vector<int> >  first_build(vector<Atom*> protein_atoms, vector<Atom*> ligand_atoms);
  void set_parameters();

  FILE * rep;
  void writeDescription(FILE * reporter);
  string description, version;


};

#endif
