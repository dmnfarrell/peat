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
#ifndef PARSEHCLASS_H
#define PARSEHCLASS_H
#include <string>
#include <vector>
#include "atom_class.h"
#include "string_tools.h"
using namespace std;

class hydrogen_conformation_class {
    //
    // This class holds the information on a single hydrogen conformation
    //
public: 
    hydrogen_conformation_class(string name,vector<string> lines);
    string name;  // Name of the conformation
    string Hname; // Name of hydrogen to be added/removed
    vector<atom_class> coordinates;
};

//
// ----------
//

class titgroup_conformations {
    //
    // This class holds the information on all conformations for a given titgroup
    //
public:
    titgroup_conformations(string groupname, string acidbase, vector<string> lines); // Parse the lines for one titgroup
    string name; // Name of the titgroup
    bool acidbase;
    vector<hydrogen_conformation_class> conformations;
};

//
// ----------
//

class hydrogens_defclass {
    //
    // This class holds all information in HYDROGENS.DAT
    //
public:
    hydrogens_defclass(vector<string> lines);
    vector<titgroup_conformations> Hbuild_titgroups;
};

#endif