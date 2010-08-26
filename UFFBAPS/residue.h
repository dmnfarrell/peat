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
#ifndef RESIDUE_H
#define RESIDUE_H
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include "atom.h"
#include "soupobject.h"

using namespace std;

class Residue: public SoupObject{
 
 public:
  ~Residue();
  

  Residue(vector<string> lines);
  vector<Atom*> getAtoms();
  vector<Atom> atoms;
  void listAtoms();
  void delete_backbone_atoms();
  void find_atom_named(string name, Atom * a);

  string residueType;
  int residueNumber; //number in read pdb file

  float phi, psi;

};

#endif
