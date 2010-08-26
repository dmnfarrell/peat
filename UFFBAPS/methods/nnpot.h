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
#ifndef NNPOT_H
#define NNPOT_H
#include "method.h"
#include "../tools/NN/NN.h" //should probaly go
#include "distances.h"
#include "targets.h"
#include "soupmanip.h"
#include <map>
#include <iostream>
#include <algorithm>
#include "nnmodel.h"
#include "elementsnn.h"
#include "elementsnn2.h"
#include "distancenn.h"
#include "pose.h"
#include "generalparameters.h"
#include "nn_es.h"
#include "nn_vdw.h"
#include "nn_ls.h"
#include "nn_hbond.h"
#include "genericnn.h"

using namespace std;

class Nnpot:  Method{

 public:
  Nnpot(FILE * resultsFile, FILE * reporter, vector<Soup *>  protein_set, vector<Soup *>  ligand_set, int model_id, bool train);
  ~Nnpot();


 private:
  
  void organise(bool do_training, vector<Soup *>  protein_set, vector<Soup *>  ligand_set);
  void train();
  void test();
  void write_poses();
  void writeDescription(FILE * reporter);
  void prepare_pose(Soup * protein, Soup * ligand);
  void prepare_decoy_pose(vector<Soup *> proteins, Soup * ligand); 
  void set_targets(string dataset);
  void write_status();
  void write_train_status(int epoch);
  void reorganise_lists(vector<Soup *>  protein_set, vector<Soup *>  ligand_set);

  float train_pct,test_pct,validate_pct,target_factor;

  bool incl_co_factors, incl_water, inputs_on_the_fly;
  
  ofstream stats_file;
  ofstream error_file;

  Nnmodel * model;
  
  vector<pose> all_poses;
  vector<pose *> train_poses;
  vector<pose *> validate_poses;
  vector<pose *> test_poses;
  vector<pose> decoy_poses;

  vector<Soup *> new_ligands, decoys;

  string dataset;
  int no_epochs,max_epochs_since_record;

  FILE * report;
  FILE * results_f;
};

#endif
