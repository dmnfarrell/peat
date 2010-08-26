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
#ifndef PREPINREADER_H
#define PREPINREADER_H
#include "method.h"
#include "soup.h"
#include <map>
#include <cctype>

using namespace std;

class PrepinReader:  Method{

 public:
  PrepinReader(FILE * reporter, vector<Soup *> A);
  PrepinReader(FILE * reporter);
  PrepinReader();
  
  ~PrepinReader();
 
  void go(vector<Soup *> A);
  void find_protein_chain_charges(vector<ProteinChain*> P);


 private:
  void writeDescription(FILE * reporter);

  //mapping partial charge, adv types and masses
  map<string,double> pc;
  map<string,string> advtype;
  map<string,float> nonbondEnergy;
  map<string,float> nonbondR;

  vector<string> keys;
  vector<string> nonbondKeys;

  void read_prep(string fn);
  void read_parm(string fn);
  void compare(vector<Atom *> atoms);
  void set_parameters(vector<Atom *> atoms);

  string make_search_key(string residue, string name);
  bool look_up_key(string key);
  bool look_up_non_bound_key(string key);

  void do_proteinChain(ProteinChain * p);
  void do_molecule(Molecule * m, string soup_name);
  void do_water(Water * w);

  void set_protonation_state_HIS(vector<Atom *> atoms);
  void set_protonation_state_LYS(vector<Atom *> atoms);
  void set_protonation_state_GLU(vector<Atom *> atoms);
  void set_protonation_state_ASP(vector<Atom *> atoms);
  void set_protonation_state_CYS(vector<Atom *> atoms);

  string correct_1_3_error(string name, string residue);
  string correct_hydrogen_name(string name);
  string correct_WI_names(string name, string residue);

  vector<string> words;
  string amino;

  ifstream file;

  



};

#endif
