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

#ifndef LIGAND_ENTROPY_H
#define LIGAND_ENTROPY_H
#include <map>
#include "method.h"
#include "soupmanip.h"
#include "distances.h"
#include "bondmatrix.h"
#include "desolvation.h"
#include "fragmental_volume.h"

using namespace std;



class LigandEntropy:  Method{

 public:
  LigandEntropy();
  LigandEntropy(FILE * reporter);
  LigandEntropy(FILE * resultsFile, FILE * reporter, vector<Soup*> A);
  LigandEntropy(FILE * resultsFile, FILE * reporter, vector<Soup*> A, vector<Soup*> B);
  ~LigandEntropy();


  float calculate_using_desolvation(Soup * protein, Soup * ligand);
  float calculate(Soup * A);
  float calculate(vector<Atom *> A);
  float no_rot_bonds();

 private:

  void set_parameters();
  bool is_terminal(int i);
  void determine_ring_memberships(int no_atoms);
  bool is_ring_bond(int i, int j);
  bool is_ring_member(int current, int last, int target, int no);

  vector<Atom*> atoms;
  vector<vector<int> > ring_bonds;
  vector<vector<int> > bonds;

  float temperature;
  string method;
  FILE * rep;
  void writeDescription(FILE * reporter);

  Desolvation DE;
  Fragmental_Volume FV;

};

#endif
