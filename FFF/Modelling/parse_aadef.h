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
#ifndef PARSE_AADEF_H
#define PARSE_AADEF_H
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
//
// Local header files
//
#include "string_tools.h"
#include "vector3d.h"
#include "atom_class.h"

//
// ---------------------------------------------------------------
//

class parse_aadef {
  //
  // Parses a single record of the amino acid definition file
  //
public:
  //
  // Dummy constructor - used when parsing rotamers
  //
  parse_aadef();
  //
  // Real constructor - used when parsing the amino acid definition file
  //
  parse_aadef(string name,vector<string> lines);
  //
  // Array that holds the atoms in each unit (rotamer or aa definition)
  //
  vector <atom_class> atoms;
  //
  // Holds vectors that each hold four atoms - the atoms defining the dihedral angle
  //
  vector <atom_class> dihedral_atoms;
  int numchi;
  string name;
  //
  // Array for storing the # of bonds that an atom is from the CA
  //
  map<string,int> bonds_from_CA;
  //
  // Inline functions for getting atoms. 
  // The functions simply return the correct instance of the
  // atom_class class
  //
  atom_class inline CA() {
    if (c_alpha!=-1) {
      return atoms[c_alpha];
    }
    cout << "CA not defined\n";
    exit(0);
  }
  atom_class inline N() {
    if (n_mainchain!=-1) {return atoms[n_mainchain];}
    cout << "N not defined\n";
    exit(0);
  }
  atom_class inline C() {
    if (c_mainchain!=-1) {return atoms[c_mainchain];}
    cout << "C not defined\n";
    exit(0);
  }
    //
    // Function for getting the atoms that are bonded
    //
    vector<string> get_bonded_atoms(string atomname) {
        for (unsigned int latom=0;latom<atoms.size();latom++) {
            if (atoms[latom].name==atomname) {
                vector<string> b_atoms;
                for (unsigned int bond=0;bond<atoms[latom].bonds.size();bond++) {
                    int bond_atom=atoms[latom].bonds[bond];
                    b_atoms.push_back(atoms[bond_atom].name);
                }
                return b_atoms;
            }
        }
        vector<string> b_atoms;
        b_atoms.resize(0);
        return b_atoms;
    }
 private:
  int c_alpha,n_mainchain,c_mainchain;
};

//
// ------------------------------------------------------------
//

class parse_aadefs  {
  //
  // Parses the amino acid definition file 
  // 
 public:
  parse_aadefs(vector<string> aadeflines);
  vector<parse_aadef> aadefs;
};

#endif
