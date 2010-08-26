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
#ifndef SIMPLELJ_H
#define SIMPLELJ_H
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include "soupobject.h"
#include "atom.h"
#include "method.h"

using namespace std;

class SimpleLJ: Method{

 public:
  SimpleLJ(FILE * resultsFile,FILE * reporter, vector<SoupObject *> A, vector<SoupObject *> B);
  ~SimpleLJ();
  double getResult();
 private:
  double result;
  map<string, int> conversion;
  float epsilon[6][6];
  float rmin[6][6];
  
  void writeDescription(FILE * reporter,vector<SoupObject *> A, vector<SoupObject *> B);

};

#endif
