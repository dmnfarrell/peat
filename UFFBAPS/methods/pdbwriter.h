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
#ifndef PDBWRITER_H
#define PDBWRITER_H
#include "method.h"

using namespace std;

class Pdbwriter:  Method{

 public:
  Pdbwriter(FILE * reporter, string fn, Soup * A);
  Pdbwriter(FILE * reporter);
  Pdbwriter();
  ~Pdbwriter();

  void write(string fn, Soup * A);

 private:
  void writeDescription(FILE * reporter);
  void make_atoms(string tag, vector<Atom *> a);

  FILE *pdb;
  int atom_no;
  int resi_no;

};

#endif
