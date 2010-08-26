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
#ifndef NN_DISCRIMINATION_H
#define NN_DISCRIMINATION_H
#include "targets.h"
#include <stdio.h>
#include <string.h>
#include "ligand_entropy.h"
#include "../tools/NN/NN.h"

using namespace std;

struct pose;

struct DecoyCmp 
{
  bool operator()(string d1, string d2 ) 
    const {
    return (atof(d1.c_str()) <  atof(d2.c_str()));
  }
};


struct PDBDecoyCmp 
{
  bool operator()(string d1, string d2 ) 
    const {
    string pdb1 = d1.substr(0,4);
    string pdb2 = d2.substr(0,4);

    string decoy1 = d1.substr(4);
    string decoy2 = d2.substr(4);


    /*    cout<<"pdb1 "<<pdb1
	<<" pdb2 "<<pdb2
	<<" decoy1 "<<decoy1
	<<" decoy2 "<<decoy2<<endl;
    */

    if( pdb1 != pdb2 )
      return (strcmp(pdb1.c_str(),pdb2.c_str())<0);
    else
      return (atof(decoy1.c_str()) <  atof(decoy2.c_str()));
  }
};



struct DecoyDecoyCmp 
{
  bool operator()(vector<string> d1, vector<string> d2 ) 
    const {

    /*    cout<<" d1[1] "<<d1[1]
	<<" d1[2] "<<d1[2]
	<<" d2[1] "<<d2[1]
	<<" d2[2] "<<d2[2]
	<<endl;
    */

    if (d1[0] !=  d2[0])
      return (strcmp(d1[0].c_str(),d2[0].c_str())<0);
    else if (atof(d1[1].c_str()) !=  atof(d2[1].c_str()))
      return (atof(d1[1].c_str()) <  atof(d2[1].c_str()));
    else
      return (atof(d1[2].c_str()) <  atof(d2[2].c_str()));
  }
};


class NN_DISCRIMINATION{

 public:
  NN_DISCRIMINATION(FILE * resultsFile,FILE * reporter);
  ~NN_DISCRIMINATION();

  void organise();  

 private:
  int no_inputs;
  
  float learning_rate, threshold;
  vector<int> tf;
  vector<int> topology;

  NN network;
  
  int no_epochs, max_epochs_since_record;
  float target_factor, train_pct, validate_pct, test_pct;
  string dataset;

  FILE * rep;
  FILE * resultsF;

  void set_parameters();

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

  //  pdb  decoy1-decoy2  inputs         
  map<string, map<string, vector<float>, DecoyCmp> > inputs;
  
  //  pdb  decoy1-decoy2  targets
  map<string, map<string, vector<float>, DecoyCmp> > targets;


  //  pdb         decoy   rmsd
  map<string, map<string, float> > rmsds;

  vector<vector<float> > training_set, test_set, validation_set;
  vector<vector<float> > training_set_targets, test_set_targets, validation_set_targets;
  vector<vector<string> > test_names;
  vector<string> tested_names;
  
  float average_test_error;

  vector<string> include_contacts;


 private:
  void prepare_set(int set, string pdb);
  vector<float> get_target(string pdb, string decoy1, string decoy2);
  void rank_results(map<vector<string>, vector<float>,DecoyDecoyCmp > results);

  ofstream error_file;


};


#endif
