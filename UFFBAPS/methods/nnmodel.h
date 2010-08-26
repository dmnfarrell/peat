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
#ifndef NNMODEL_H
#define NNMODEL_H
//#include <stdlib.h>
//#include "nnpot.h"
#include "../tools/NN/NN.h"
#include "../soupobject.h"
#include <map>
#include "soupmanip.h"
#include "pose.h"
#include <iostream>
#include <fstream>

using namespace std;

class Nnmodel{

 public:
  Nnmodel();
  virtual ~Nnmodel();

  virtual void create_network(){}
  virtual void define_input(){}
  virtual void test_pose(pose * A){}
  virtual void train_pose(pose * A){}
  virtual void update_weights(){cout<<"Virtual function - you should not see this!\n";}
  virtual void build_inputs(pose * A){}
  virtual void write_network(string name){}
  virtual void read_network(string name){}

  virtual int get_no_epochs(){return 0;}
  virtual int get_max_epochs_since_record(){return 0;}
  virtual float get_target_factor(){return 0.0;}
  virtual string get_dataset(){return "";}
  virtual float get_train_pct(){return 0.0;}
  virtual float get_validate_pct(){return 0.0;}
  virtual float get_test_pct(){return 0.0;}
  virtual int get_smoothing(){return 0;}
  virtual bool get_include_co_factors(){return false;}
  virtual bool get_include_water(){return false;}
  /*
  void create_network();
  void define_input();
  void test_pose(pose * A);
  void train_pose(pose * A);
  void prepare_pose(pose * A);
  void build_inputs(pose * A);
  void write_network(string name);
  void read_network(string name);
  */
};

#endif
