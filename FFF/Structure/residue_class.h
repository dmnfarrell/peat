/*
 #
 # FFF - Flexible Force Field
 # Copyright (C) 2010 Jens Erik Nielsen
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
#ifndef RESIDUE_CLASS_H
#define RESIDUE_CLASS_H
#include <vector>
#include <string>
#include <stdio.h>
//
// Local header files
//
#include "atom_class.h"
// Forward declaration
class chain_class;
//
class residue_class {
public:
  //
  // First initialisation is when reading the pdb file
  //
  residue_class(chain_class* C,vector<string> lines);
    residue_class(string resname,vector<string> charge_lines) {
        //
        // Initialise charges and radii
        //
        name=resname;
        for (unsigned int count=0;count<charge_lines.size();count++) {
            vector<string> sp=split(charge_lines[count]);
            //printf ("Atom: %s charge: %6s radius: %6s\n",sp[0].c_str(),sp[1].c_str(),sp[2].c_str());
            atoms.push_back(atom_class(sp[0],sp[1],sp[2]));
        }
    }
  //
  // This initialisation is for putting in instances of atom_class
  //
  //residue_class(residue_class* C,vector<atom_class>);
  //
  // Update 
  //
  void update();
    void print();
     //bool is_aa();
  //
  //Variables
  //
  int number, inchain;
  vector<atom_class> atoms;
  string name, pdbnum, chainname;
  double phi, psi, omega;
  string chain;
  atom_class *N,*CA,*C;
    chain_class* _C;
   
} ;
#endif
