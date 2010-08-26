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
#ifndef SOUP_H
#define SOUP_H
#include <stdlib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include "proteinchain.h"
#include "water.h"
#include "molecule.h"
//#include "soupmanip.h"

using namespace std;

class Soup {
public:
  Soup();
  ~Soup();

  void addToSoup(string filename);
  void addToSoup(vector<string> lines, string filename);
  void soupStat();
  void listAtoms();

  Soup operator+(Soup &other);
  

  vector<Atom *> getAtoms();
  vector<Atom *> get_dry_atoms();
  vector<ProteinChain *> getProteinChains();
  vector<Water *> getWaters();
  vector<Molecule *> getMolecules();

  vector<ProteinChain> copyProteinChains();

  vector<vector<SoupObject *> > getSoupObjects();

  int getTotalNumberOfAtoms();
  void clearSoup();

  string name,assigned_parameters;
  bool isActive,bonds_assigned,sybyl_types_assigned;
  int totalNumberOfAtoms;

  
  void add_protein_chain(ProteinChain object);
  void add_molecule(Molecule object);
  void add_water(Water object);

  void delete_protein_chain(int i);

  void set_sybyl_atom_types();


  /*  void add_soup_object(SoupObject object);
  void delete_soup_object(string name);
  */
  
 private:
  void interpretatePDBLines(vector<string> lines, string filename);
  void interpretateMOL2Lines(vector<string> lines, string filename);

  vector<ProteinChain> proteinChains;
  vector<Water> waters;
  vector<Molecule> molecules;

  map<string,string> sybyl_types;

};

#endif
