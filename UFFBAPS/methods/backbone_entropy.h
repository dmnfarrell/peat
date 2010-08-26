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

#ifndef BACKBONEENTROPY_H
#define BACKBONEENTROPY_H
#include "method.h"
#include <vector>
#include <map>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <list>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "proteinchain.h"
#include "residue.h"
using namespace std;

namespace la = boost::numeric::ublas;


struct propensity_map {
  
  vector<vector<int> > angles;
  vector<vector<float> > energies;
  vector<float> phis,psis;

};

class BackboneEntropy:  Method{

 public:
  BackboneEntropy();
  BackboneEntropy(FILE * reporter);
  BackboneEntropy(FILE * reporter, vector<Soup*> A, bool train);
  ~BackboneEntropy();

  void read_in_maps();
  void precalculate_energies();

  void calculate(vector<Soup*> A);
  float calculate(Soup* A);
  float calculate(ProteinChain * protein_chain);
  

  void do_statistics(vector<Soup*> A);
  bool set_phi_and_psi(Residue * prev_residue, Residue * residue, Residue * next_residue);
  la::vector<float> cross(la::vector<float> a, la::vector<float> b);

  float dot(la::vector<float> a, la::vector<float> b);
  float norm(la::vector<float> a);

  int make_index(float angle);

 private:
  void writeDescription(FILE * reporter);
  
  propensity_map new_propensity_map();
  void insert_angles(propensity_map * map, float phi, float psi);
  void print_map(propensity_map * map, string filename);
  void print_all_maps();


  map<string,propensity_map> propensity_maps;
  

};

#endif
