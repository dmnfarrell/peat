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
#ifndef CHAIN_H
#define CHAIN_H
#include <vector>
#include <string>
//
// Local header files
//
#include "residue_class.h"
#include "atom_class.h"
//
// Forward declaration
//
class FFF;

class chain_class {
public:
       //
       // Functions
       //
       // Initialise with a pdb file
       //
       chain_class(FFF* P,vector<string> lines);
       //
       // Intilialise with a list of resiudes
       //
       //chain_class(FFF& P,vector<residue_class> residues);
       //
       // Update pointers
       //
       void update();
       //
       // Variables
       //
       vector<residue_class> residues;
       vector<atom_class*> atoms;
       int number;
       string name;
    FFF* _P;
};
#endif
