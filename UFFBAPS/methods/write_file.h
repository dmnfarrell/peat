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
#ifndef WRITE_FILE_H
#define WRITE_FILE_H
#include "method.h"

using namespace std;

class Write_File:  Method{

 public:
  Write_File();
  Write_File(FILE * reporter);
  Write_File(FILE * reporter, vector<Soup*> A, string format);
  ~Write_File();

  void write(vector<Soup*> A, string format, vector<string> remarks, string filename);
  void write_atoms(vector<Atom*> atoms, string format, string prefix, string filename);

 private:


  void write_objects(vector<SoupObject*> objects, string format, string prefix);
  void write_atoms(vector<Atom*> atoms, string format, string prefix);


  void writeDescription(FILE * reporter);

  FILE * out;
  int atom_serial, residue_serial;

};



#endif
