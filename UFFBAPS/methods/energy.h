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

#ifndef ENERGY_H
#define ENERGY_H
#include "method.h"
#include "vdw.h"
#include "vdw_flat.h"
#include "vdw_line.h"
#include "es.h"
#include "soupmanip.h"
#include "distances.h"
#include "bondmatrix.h"
#include "hbonds.h"
#include "hp.h"
#include "desolvation.h"
#include "protein_entropy.h"
#include "ligand_entropy.h"
#include "water_entropy.h"
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "matrix_io.h"
#include "targets.h"
#include "generalparameters.h"
#include "simple_parameters.h"
#include "water_bridges.h"

namespace la = boost::numeric::ublas;

using namespace std;

class Energy:  Method{

 public:
  Energy(FILE * reporter, vector<Soup *> proteins, vector<Soup *> ligands);
  Energy();
  ~Energy();

  void initiate();
  double calculate_binding(Soup * protein, Soup * ligand);
  int get_no_components();
  vector<double> get_components();


 private:
  void writeDescription(FILE * reporter);

  Distances dist;
  BondMatrix BM;

  Vdw_Line vdw_line;
  Es es;
  ProteinEntropy pe;
  Desolvation de;
  Hbonds hb;
  Hp hp;
  LigandEntropy le;
  Water_Entropy we;
  Water_Bridges wb;



  void test();
  void train();
  void set_coefficients(la::vector<double> c);
  void stream_results();
  double calc_total_binding_energy(int i);
  void read_parameters();
  void generate_target_vector();
  void calculate_result();

  la::matrix<double> A;
  la::vector<double> b,coefficients,components;

  vector<Soup *> proteins,ligands;
  
  Simple_Parameters simple_parameters;

  string mode,dataset;
  int no_components;
  double r_vdw, r_es, r_hb, r_hp, r_pe, r_le, r_we, de_holo, de_apo, r_wb;
  vector<vector<float> > intermolecular_distances;
  vector<vector<int> > ligand_bonds;
  vector<string> names;

  double 
    vdw_component, es_component, pe_component, 
    hb_component, de_component, le_component, bb_component,
    we_component, wb_component, total_binding;

  bool are_coefficients_set;

  FILE * rep;
};

#endif
