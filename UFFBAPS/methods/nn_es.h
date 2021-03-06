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
#ifndef NN_ES_H
#define NN_ES_H
#include "nnmodel.h"

using namespace std;

struct pose;

class NN_ES: public Nnmodel{

 public:
  NN_ES(FILE * reporter);
  ~NN_ES();

  void create_network();
  void build_inputs(pose * A);
  void define_input();
  void test_pose(pose * A);
  void train_pose(pose * A);
  void write_network(string name);
  void read_network(string name);

  int get_no_epochs();
  int get_max_epochs_since_record();
  float get_target_factor();
  string get_dataset();

  float get_train_pct();
  float get_validate_pct();
  float get_test_pct();

  bool get_include_co_factors();
  bool get_include_water();

 private:
  int no_types;
  int no_inputs;


  float learning_rate, threshold, linear_coef;
  vector<int> tf;
  vector<int> topology;

  vector<vector<NN> > networks;
  
  map<string, int> element_indexes;
  vector<int> powers;
  vector<string> include_elements;
  vector<string> ignore_elements;

  int no_epochs, max_epochs_since_record, include_co_factors, include_water;
  float target_factor, train_pct, validate_pct, test_pct;
  string dataset;


  void set_parameters();
  void writeDescription(FILE * reporter);
  string description, version;


};

#endif
