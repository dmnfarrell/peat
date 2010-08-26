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
#ifndef NN_FITTER_H
#define NN_FITTER_H
#include "targets.h"
#include "ligand_entropy.h"
#include "../tools/NN/NN.h"

using namespace std;

struct pose;

class NN_FITTER{

 public:
  NN_FITTER(FILE * resultsFile,FILE * reporter);
  ~NN_FITTER();

  void organise();  

 private:
  int no_inputs;
  vector<string> input_lines;
  
  float learning_rate, threshold, linear_coef;
  vector<int> tf;
  vector<int> topology;

  NN network;
  
  int no_epochs, max_epochs_since_record;
  float target_factor, train_pct, validate_pct, test_pct;
  string dataset;

  FILE * rep;
  FILE * results;

  void set_parameters();
  void set_targets(string dataset);
  void create_network();
  void build_inputs();

  void go();
  float test();
  float train();
  float validate();
  void write_network(string name);
  void read_network(string name);
  void writeDescription(FILE * reporter, FILE * resultsFile );
  float get_contacts(string pdb, string dist);
  string description, version;

  map<string, vector<float> > inputs;
  map<string, float> targets;

  vector<vector<float> > training_set,test_set,validation_set;
  vector<vector<float> > training_set_targets,test_set_targets,validation_set_targets;
  vector<string> test_names;
  vector<string> tested_names;

  float average_test_error;

  vector<string> include_contacts;

};

#endif
