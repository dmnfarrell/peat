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

#ifndef MSC_H
#define MSC_H
#include "method.h"
#include "vdw.h"
#include "vdw_flat.h"
#include "vdw_line.h"
#include "es.h"
#include "hbonds.h"
#include "hp.h"
#include "protein_entropy.h"
#include "ligand_entropy.h"
#include "soupmanip.h"
#include "distances.h"
#include "desolvation.h"
#include "backbone_entropy.h"
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "matrix_io.h"
#include "simple_parameters.h"
#include "water_bridges.h"

using namespace std;
namespace la = boost::numeric::ublas;

class MSC:  Method{

 public:
  MSC(FILE * reporter, FILE * results, vector<Soup *> wts_in, vector<Soup *> mutants_in);
  MSC(); // Constructer for scripting

  ~MSC();

  void initiate();
  double calculate_stability(Soup * soup);
  int get_no_components();
  vector<double> get_components();

 private:

  vector<Soup *> wts, mts;

  double 
    wt_vdw, wt_es, wt_pe, wt_hb, wt_de, wt_le, wt_bb,
    mut_vdw, mut_es, mut_pe, mut_hb, mut_de, mut_le, mut_bb,
    D_vdw, D_es, D_pe, D_hb, D_de, D_le, D_bb;

  la::matrix<double> A,wt_comps,mut_comps;
  la::vector<double> b,coefficients,components;

  double 
    vdw_component, es_component, pe_component, 
    hb_component, de_component, le_component, bb_component,
    wb_component, total_stability;

  


  Distances dist;
  Vdw_Line vdw_line;
  Es es;
  Hbonds hb;
  //Hp hp ;
  ProteinEntropy pe;
  LigandEntropy le;
  Desolvation de;
  BackboneEntropy bb;
  Water_Bridges wb;

  map<char, int> aa_map;

  vector<int> residues;
  vector<string> names;

  string mode;
  int no_components,include_ff_components,include_aa_components,coefficients_are_set,count;

  PrepinReader prepin;
  Correct_Atomnames correct_names;
  Simple_Parameters simple_parameters;

  string print_matrix(la::vector<double> v);
  string print_matrix(la::matrix<double> m);
  void set_coefficients(la::vector<double> in);

  void train();
  void test();

  double calc_total_DDG(int i);

  void stream_results();
  void calculate_result();
    
  void read_parameters();
  void writeDescription(FILE * reporter);
  FILE * rep;
  FILE * output;
};

#endif
