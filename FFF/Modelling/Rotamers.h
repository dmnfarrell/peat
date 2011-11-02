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
#ifndef ROTAMER_H
#define ROTAMER_H
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <exception>
#include <fstream>
#include "string_tools.h"

using namespace std;

class Rotamer_class {
 public:
  //
  // Public functions
  //
  //
  // Constructor
  //
    Rotamer_class(const std::string filename);   //
    vector<vector<float> > get_rotamers(string resname,double phi, double psi);
  // Constructor end
  //
    vector<string> _names;
 private:
  vector<vector<vector<float> > >_rotamers;
  
};

#endif
