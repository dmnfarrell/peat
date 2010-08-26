/*!
 #
 # pKaTool - analysis of systems of titratable groups
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
#ifndef PKAMC_H
#define PKAMC_H

#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

class MC {
public:
  //
  // Functions
  //
 MC(vector<double> intpKas, vector<double> lin_matrix, vector<double> ab,int use_MC): _intpKas(intpKas),_lin_matrix(lin_matrix),_acid_base(ab) {
    reformat_arrays();
    // Set default value for MCsteps
    _MCsteps=20000;
    //Store Boltzmann/Monte Carlo flag
    _use_MC=use_MC;
  };
  //
  //void set_acid_base(vector<int> acid_base);
  //
  // Reformats the matrix
  void reformat_arrays();
  //
  void set_MCsteps(int MCsteps) {
    _MCsteps=MCsteps;
    return;
  }
  vector<float> calc_pKas(float pH_start,float pH_end, float pH_step);
  double calc_pKa(vector<float> charges, vector<double> pHs,double acid_base);
  vector<float> calc_charge(float pH);
  vector<float> calc_charge_Boltzmann(float pH);
  vector<int> get_state_specification(int state_number,int num_groups);
  double get_energy(float pH, vector<int> state);
  double get_energy_fast(float pH,vector<int> state,int change_group,double old_energy);

  void set_monitor_states(vector<int> mon_state);
  //
  // utility
  //
  void print_vector(vector<int> state);
  //
  //
  // Variables
  //
  vector<double> _intpKas, _lin_matrix, _acid_base;
  vector<vector<double> > _matrix, _charges;
  int _groups, _MCsteps, _use_MC, _do_monitor;
  vector<int> _monitor_state;
  vector<double> _CCPS_population;
  double lnten;
};
#endif
