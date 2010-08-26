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

#ifndef CHEMICAL_SHIFT_CHANGES_H
#define CHEMICAL_SHIFT_CHANGES_H
#include "method.h"
#include "soupmanip.h"
#include "bondmatrix.h"
#include "distances.h"
#include "pdbwriter.h"
#include <sstream>

using namespace std;

class Chemical_Shift_Changes:  Method{

 public:
  Chemical_Shift_Changes();
  Chemical_Shift_Changes(FILE * reporter);
  Chemical_Shift_Changes(FILE * reporter, vector<Soup*> A);
  ~Chemical_Shift_Changes();


  void amide_shifts(Soup * A);
  void amide_shifts(vector<ProteinChain *> protein_chains);

  map<string,int> titratable_residues;
  map<string,string> names3to1,names1to3;

  float induced_shift(Atom * charged_atom, Atom * N, Atom * H, Atom * C);

 private:
  float angle(vector<float> v1, vector<float> v2);
  void generate_list_of_amides(vector<Atom*> atoms);
  vector<Atom*>::iterator has_element_bonded(Atom * atom, string element);
  vector<Atom*>::iterator has_atom_name_bonded(Atom * atom, string name);

  int no_hydrogen_bounded(Atom * atom);
    

  void write_out_chemical_shifts(string fn, vector<Atom*> atoms);


  void set_parameters();
  void writeDescription(FILE * reporter);

  BondMatrix BM;
  Distances D;
  Pdbwriter PDB_writer;
  
  float alpha,pi,eps_0,e,a0,Eh,eps_rel;


  vector<vector<Atom*> > amides;


};

#endif
