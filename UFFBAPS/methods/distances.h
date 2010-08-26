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

#ifndef DISTANCES_H
#define DISTANCES_H
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include "soupobject.h"
#include "soupmanip.h"
#include "method.h"

using namespace std;

class Distances:  Method{

 public:
  Distances(FILE * reporter);
  Distances();
  ~Distances();

  void calculate(SoupObject * A, bool sqrt);
  void calculate(vector<Atom*> atoms, bool doSqrt);
  void calculate(SoupObject * A, SoupObject * B, bool sqrt);
  void calculate(vector<SoupObject *> A, vector<SoupObject *> B, bool doSqrt);
  void calculate(vector<Atom *> atomsA, vector<Atom *> atomsB, bool doSqrt);
  void calculate(Soup * A, Soup * B, bool doSqrt);
  void calculate(Atom * atomA, vector<Atom *> atomsB, bool doSqrt);

  float calculate(Atom * atomA, Atom * atomB, bool doSqrt);

  vector<vector<float> > getResult();
  vector<vector<float> > * getResultPointer();
  void print_result();


 private:
  vector<vector<float> > result;
  void writeDescription(FILE * reporter);

  float res,dx,dy,dz;


};


#endif
