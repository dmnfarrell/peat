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
#ifndef VDW_LINE_H
#define VDW_LINE_H
#include "method.h"
#include <math.h>
#include <algorithm>
#include "distances.h"
#include "soupmanip.h"
#include "prepinreader.h"
#include "generalparameters.h"
#include "topology.h"
#include "bfactor.h"
#include "bondmatrix.h"


using namespace std;

class Vdw_Line:  Method{

 public:
  Vdw_Line();
  Vdw_Line(FILE * reporter);
  Vdw_Line(FILE * resultsFile, FILE * reporter, vector<Soup *> A, vector<Soup *> B);
  Vdw_Line(FILE * resultsFile, FILE * reporter, vector<Soup *> A);
  ~Vdw_Line();

  double calculate(Soup * A, Soup * B, double cutoff);
  double calculate(Soup * A, double cutoff);
 private:
  void writeDescription(FILE * reporter);
  void set_parameters();
  double calculate_vdw_line(Atom * atom1, Atom * atom2, double dist);

  Distances dist;
  double width, top, deltaR, result, same_residue_factor, neighbour_residue_factor, R, E, res, fraction, fraction2, fraction3, fraction6, fraction12;
  string parameter_type;
  int include_surf_factor, residue1, residue2, LJ;

  BondMatrix BM;

  vector<float> ligand_surf_dists, protein_surf_dists;
  FILE * rep;
};

#endif
