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
#ifndef ASSIGN_H
#define ASSIGN_H
#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <unistd.h>
#include "residue_class.h"

using namespace std;

class crgrad_class {
public:
    crgrad_class (const std::string parmdir); 
    void assign_all(residue_class residue);
    vector<residue_class> types;
};

#endif