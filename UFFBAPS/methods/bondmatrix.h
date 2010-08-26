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

#ifndef BONDMATRIX_H
#define BONDMATRIX_H
#include "distances.h"
#include "method.h"
#include "Boxes.h"

using namespace std;

class BondMatrix:  Method{

 public:
  BondMatrix(FILE * reporter, SoupObject * A);
  BondMatrix(FILE * reporter, vector<Atom*> A);
  BondMatrix(FILE * reporter);
  BondMatrix();
  ~BondMatrix();

  void assign_bonds(Soup * A);

  void calculate(vector<Atom*> atoms);
  void calculate_list(vector<Atom*> atoms);
  void calculate(SoupObject * A);

  vector<vector<int> > getResult();

  bool connection(vector<vector<int> > * bonds, int i, int j, int no_bonds, int cur_no_bonds,int last);
  void store_bonding_information_in_atoms(vector<Atom*> atoms);

 private:
  FILE * rep;
  void writeDescription(FILE * reporter);
  vector<vector<int> > result;
  int noBondsCompExp(int a, vector<Atom*> atoms);

  int no_bonds(int i);
  int no_bonds_matrix(int i);
  void only_bond(int i,int j,int type);
  void make_list_entry(int i,int j,int type);
  bool is_in_list(int i,int j);
  void change_type(int i,int j,int newtype);

  Distances d;


  vector<int> find_bonding_atoms(vector<vector<int> > * bonds, int a);


  vector<vector<int> > _bondmatrix;

};

#endif
