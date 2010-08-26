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
#ifndef TARGETS_H
#define TARGETS_H
#include "method.h"
#include <iostream>



using namespace std;

class Targets:  Method{

 public:
  Targets(FILE * reporter, string dataset);
  ~Targets();
  
  float get_mlogKi(string key);
  float get_G(string key);
  float get_resolution(string key);
  float get_no_rotor(string key);
  int find_key(string key);

  float get_decoy_G(string key, float rmsd);
  float get_rmsd(string key, int decoy_no);


 private:

  vector<string> keys;
  vector<float> mlogKis;
  vector<float> Gs;
  vector<float> resolutions;
  vector<int> no_rotors;

  void read_file_wang(string filename);
  void read_file_pdbbind(string filename);
  void read_file_fictive(string filename);

  void writeDescription(FILE * reporter);



  float T,R,log_conversion;


};

#endif
