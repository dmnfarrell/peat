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
#ifndef PROTEINCHAIN_H
#define PROTEINCHAIN_H
#include <stdio.h>
#include <string>
#include <vector>
#include "soupobject.h"
#include "residue.h"


using namespace std;

class ProteinChain: public SoupObject{
 
 public:
  ProteinChain(vector<string> lines);
  ProteinChain();
  ~ProteinChain();

  int numberOfAtoms();

  void listAtoms();
  vector<Atom*> getAtoms();
  vector<Residue*> getResidues();

  vector<Residue> copyResidues();
  vector<Residue>* get_pointer_to_vector_of_residues();

  void add_residue(Residue res);
  vector<Residue>::iterator delete_residue(int i);
  vector<Residue>::iterator delete_residue(vector<Residue>::iterator it);

 private:
  vector<Residue> residues;

};

#endif
