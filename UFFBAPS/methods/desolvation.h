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

#ifndef DESOLVATION_H
#define DESOLVATION_H
#include "method.h"
#include <math.h>
#include <map>
#include "fragmental_volume.h"
#include "bondmatrix.h"
#include "hp.h"
#include "Boxes.h"

using namespace std;

class Desolvation:  Method{

 public:
  Desolvation(FILE * reporter);
  Desolvation();
  ~Desolvation();
  Desolvation(FILE * resultsFile,FILE * reporter, vector<Soup*> A, vector<Soup*> B);
  Desolvation(FILE * resultsFile,FILE * reporter, vector<Soup*> A);

  float calculate_contact_between_A_and_B(Soup * A, Soup * B);
  float calculate(Soup * A, Soup * B);
  float calculate(Soup * A);


 private:
  void writeDescription(FILE * reporter);

  float atomic_desolvation(Atom * atom);

  void prepare_soup(Soup * A);
  float sum_atoms(vector<Atom*> all);

  void set_desolvation_factors(Soup * A);
  void set_desolvation_factors(Soup * A, Soup * B);
  void set_desolvation_factors(vector<Atom*> atoms, bool include_hydrophilic);

  float get_min_desolvation_factor(string key);
  float get_max_desolvation_factor(string key);

  void set_parameters();

  Fragmental_Volume FV;
  BondMatrix BM;
  Hp HP;

  /*
  struct desolvation_parameter
  { 
    //    string aa;
    //    string atom;
    float volume;
    float Occ;
    float Occmax;
    float vdw; 
    float enerG;
    float Occ_short;
    float Occmax_short;
    float radius;
  };
  
  map<string,desolvation_parameter> desolvation_parameters;

  */


  map<string,float> desolvation_parameters, min_occupancies, max_occupancies;

  float desolvation_distance, avr_min_desolvation_factor, avr_max_desolvation_factor;

  FILE * out;


  
};

#endif
