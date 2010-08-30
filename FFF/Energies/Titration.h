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

#ifndef TITRATION_H
#define TITRATION_H
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <exception>

#include "matrix.h"
#include "energy_base_class.h"
#include "Boxes.h"
#include "fff.h"
#include "atom_class.h"
#include "vector3d.h"

//
// Energies related to the titration of acids and bases
//
class Titration_Energy : public energy_base_class {
public:
    //
    // Constructor
    //
    Titration_Energy(FFF& P) : energy_base_class(P) {};
    //
    // Function for calculating energy
    //
    virtual double calculate_interaction(const atom_class& atom1,const atom_class& atom2);  
    double get_energy();
    double get_titgroup_energy(TitGroup& titgroup_number);
};

class Titration_Error {};

#endif
