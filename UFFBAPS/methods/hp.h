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

#ifndef HP_H
#define HP_H
#include "method.h"
#include "bondmatrix.h"
#include "generalparameters.h"
#include "bfactor.h"
#include "topology.h"

using namespace std;

class Hp:  Method{

 public:
  Hp(FILE * reporter);
  Hp(); //temporary constructor for use in water entropy
  Hp(FILE * resultsFile, FILE * reporter, vector<Soup *> A, vector<Soup *> B);
  Hp(FILE * resultsFile, FILE * reporter, vector<Soup *> A);

  ~Hp();

  float calculate(Soup * A, Soup * B);
  float calculate(Soup * A);

  float calculate(vector<Atom*> protein_atoms, 
		  vector<Atom*> ligand_atoms,
		  vector<vector<float> > distances, 
		  string assigned_parameters,
		  string assigned_parametersB);

  float calculate(vector<Atom*> atoms,
		  string assigned_parameters);

  void find_hydrophobic_atoms(vector<Atom*> atoms);

 private:

  Distances dist;
  BondMatrix BM;

  bool is_non_hetero_carbon(Atom * atom);
  bool is_non_hetero_hydrogen(Atom * atom);
  bool is_non_ionic_halogen(Atom * atom);

  float interaction(float dist, float sum_vdw_radii);

  void set_parameters();
  void writeDescription(FILE * reporter);

  int include_non_hetero_carbon,include_non_hetero_hydrogen,include_halogens,include_surf_factor;

  vector<float> ligand_surf_dists, protein_surf_dists;

  FILE * rep;

  vector<string> halogens;

};

#endif
