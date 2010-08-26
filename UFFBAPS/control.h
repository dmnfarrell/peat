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

#ifndef CONTROL_H
#define CONTROL_H
#include <stdio.h>
#include <string>
#include <iostream> 
#include <vector>
#include <sstream>
#include <exception>
#include <algorithm>
#include <cctype>
#include "soup.h"
#include "stringtools.h"
#include "taskmaster.h"
#include "methods/method.h"
#include "methods/mutant_stability_change.h"
#include "methods/energy.h"
#include "methods/matrix_io.h"
#include "methods/ligand_entropy.h"

//#include "atom.h"

using namespace std;

class Control{

 public:
  Control();
  ~Control();  
  
  void command(string line);
  void command();
  void commands(vector<string> lines);
  void commands(string file);
  
 private:

  void interpretator(string line);
  vector<string> lines;
  vector<Soup> soups;
  TaskMaster taskMaster;
  vector<Soup*> list1;
  vector<Soup*> list2;

  //script related stuff
  bool script_initiated;
  MSC Stability;
  Energy Binding;
  LigandEntropy Ligand_entropy;

  vector<vector<double> > A;
  vector<double> b;
  vector<string> names;
  void do_script(vector<string> words);
  void initiate_scripting();

  void doTask(vector<string> words);
  void load(string pdbfile, string name);
  void loadToSoup(int offset, string pdbfile, string name);
  void stat();
  void list();
  void clear(); 
  vector<string> readFile(string filename);
  void addSoupToList(string name, int list);
  vector<vector<string> > divide_mol2(vector<string> lines);
  vector<vector<string> > divide_pdb(vector<string> lines);
};

#endif
