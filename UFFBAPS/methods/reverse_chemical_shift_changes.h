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
#ifndef REVERSE_CHEMICAL_SHIFT_CHANGES_H
#define REVERSE_CHEMICAL_SHIFT_CHANGES_H
#include "method.h"
#include "chemical_shift_changes.h"
#include "distances.h"

using namespace std;

struct amide{

  Atom* N,*H,*C;
  float experimental_chemical_shift,dielectric,electrostatic_field,electrostatic_potential;

};


class Reverse_Chemical_Shift_Changes:  Method{

 public:
  Reverse_Chemical_Shift_Changes(FILE * reporter);
  Reverse_Chemical_Shift_Changes(FILE * reporter,vector<Soup*> A);

  ~Reverse_Chemical_Shift_Changes();



  void calculate_dielectrics(vector<Soup*> A);
  void calculate_dielectrics(Soup * A);
  void calculate_dielectrics(vector<Atom*> atoms_in);


 private:
  void read_in_chemical_shifts(string filename);
  void set_affector_atom(string residue);
  void insert_amide(int residue, float chem_shift);

  void set_parameters();
  void writeDescription(FILE * reporter);

  Atom * affector_atom;
  float affector_atom_charge;

  vector<amide> amides;
  vector<Atom*> atoms;

  float alpha,pi,eps_0,e,a0,Eh;
  
  Chemical_Shift_Changes CSC;
  Distances D;
  
  string data_file;

};




#endif
