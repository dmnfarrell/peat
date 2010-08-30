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

#ifndef BOXES_H
#define BOXES_H
#include <math.h>
#include <vector>
#include <string>
#include "atom_class.h"
//#include "fff.h"

using namespace std;

class Boxes {
 public:
   // Boxes(FFF &P,float boxsize) : _P(P); // Boxes of a molecule
  Boxes(vector<atom_class *> atoms,float boxsize); //Boxes of a vector of atoms
  void make_boxes(vector<atom_class *> atoms,float boxsize); // Routine that does the work
  vector<atom_class *> get_close_atoms(atom_class * atom); // Get the atoms that are close  
    vector<atom_class *> get_close_atoms(atom_class& atom);

  int xnum,ynum,znum;
  double xmax,xmin,ymax,ymin,zmax,zmin;
 protected:
  vector<atom_class *> atomsA, atomsB; 
  vector<vector<atom_class *> > boxes;
  float _boxsize;
    //FFF& _P;
};

#endif
