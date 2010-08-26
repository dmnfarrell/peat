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
#include "pMC.h"
//
//
//

void MC::reformat_arrays() {
  //
  // Make the matrix
  //
  _groups=static_cast<int>(_intpKas.size());
  int count=0;
  for (int row=0;row<_groups;row++) {
    vector<double> row_vals;
    for (int column=0;column<_groups;column++) {
      row_vals.push_back(_lin_matrix[count]);
      count=count+1;
    }
    _matrix.push_back(row_vals);
  }
  //
  // Set natural log
  //
  lnten=log(10);
  return;
}

//
// ---------------------
//

vector<float> MC::calc_pKas(float pH_start,float pH_end, float pH_step) {
  //
  // Calculate pKa values for the system
  //
  // First get charges at all pH values
  //
  float max_pH=0.0;
  vector< vector<float> > charges;
  for (float pH=pH_start;pH<pH_end;pH=pH+pH_step) {
    if (_use_MC==1) {
      // Use Monte Carlo
      charges.push_back(calc_charge(pH));
    } else {
      // Evaluate the Boltzmann sum
      charges.push_back(calc_charge_Boltzmann(pH));
    }
    max_pH=pH;
  }
  // 
  // Now determine pKa values
  //
  int datapoints=5; // How many data points do we use for pKa determination?
  datapoints=(datapoints-1)/2;
  vector<float> pKas;
  for (int group=0;group<_groups;group++) {
    int count=0;
    float pKa=-999.9;
    float last_crg=charges[count][group];
    //
    for (float pH=pH_start;pH<pH_end;pH=pH+pH_step) {
      if ((pH-max_pH)>0.0) {
	continue;
      }
      float this_crg=charges[count][group];
      if (_acid_base[group]==1.0) {
	if (this_crg<=0.5 && last_crg>0.5) {
	  //pKa=(last_crg-0.5)/(last_crg-this_crg)*pH_step+(pH-pH_step);
	  //
	  // Get ph,charge sets and calc pKa from those
	  //
	  vector<double> pHs_pKadet;
	  vector<float> charges_pKadet;
	  int count2=count-static_cast<int>(datapoints);
	  if (count2<0) {count2=0;}
	  for (double pH2=max(pH_start,pH-datapoints*pH_step);pH2<min(pH_end,pH+datapoints*pH_step);pH2=pH2+pH_step) {
	    pHs_pKadet.push_back(pH2);
	    charges_pKadet.push_back(charges[count2][group]);
	    count2=count2+1;
	  }
	  pKa=calc_pKa(charges_pKadet,pHs_pKadet,_acid_base[group]);
	}
      } else {
	if (this_crg<=-0.5 && last_crg>-0.5) {
	  //pKa=(last_crg-(-0.5))/(last_crg-this_crg)*pH_step+(pH-pH_step);
	  //
	  // Get ph,charge sets and calc pKa from those
	  //
	  vector<double> pHs_pKadet;
	  vector<float> charges_pKadet;
	  int count2=count-static_cast<int>(datapoints);
	  if (count2<0) {count2=0;}
	  for (double pH2=max(pH_start,pH-datapoints*pH_step);pH2<min(pH_end,pH+datapoints*pH_step);pH2=pH2+pH_step) {
	    pHs_pKadet.push_back(pH2);
	    charges_pKadet.push_back(charges[count2][group]);
	    count2=count2+1;
	  }
	  pKa=calc_pKa(charges_pKadet,pHs_pKadet,_acid_base[group]);
	}
      }
      last_crg=this_crg;
      count=count+1;
    }
    pKas.push_back(pKa);
  }   
  //
  // We also return all the charges
  //
  // The first numbers we return is the start pH, the pH step, and the number of values
  //
  int num_pHs=0;
  for (float pH=pH_start;pH<pH_end;pH=pH+pH_step) {
    num_pHs++;
  }
  pKas.push_back(pH_start);
  pKas.push_back(pH_step);
  pKas.push_back(static_cast<float>(num_pHs));
  // 
  // Now add the charges
  //
  float this_crg;
  int count=0;
  for (int group=0;group<_groups;group++) {
    count=0;
    for (float pH=pH_start;pH<pH_end;pH=pH+pH_step) {
      pKas.push_back(pH);
      this_crg=charges[count][group];
      pKas.push_back(this_crg);
      count=count+1;
    }
    pKas.push_back(999.0);
    pKas.push_back(-999.0);
  }
  return pKas;
}

//
// ---------------------
//

double MC::calc_pKa(vector<float> charges,vector<double> pHs,double acid_base) {
  //
  // Calculate the pKa value from a selection of charges and pH values
  //
  // Assume perfect Henderson-Hasselbalch behaviour
  //
  // acid_base = -1.0 for acids
  // acid_base = 1.0 for bases
  double ratio=0.0;
  vector<double> pKas;
  double pKa=0.0;
  int points=static_cast<int>(charges.size());
  for (int count=0;count<points;count++) {
    if (acid_base!=1.0) {
      ratio=fabs(charges[count])/(1.0-fabs(charges[count]));
    } else {
      ratio=(1.0-fabs(charges[count]))/fabs(charges[count]);
    }
    pKas.push_back(pHs[count]-log10(ratio));
  }
  //
  // Find the average of the pKa values
  //
  double sum=0.0;
  for (int count=0;count<static_cast<int>(pKas.size());count++) {
    sum=sum+pKas[count];
  }
  pKa=sum/static_cast<double>(pKas.size());
  return pKa;
}
    

//
// ---------------------
//

vector<float> MC::calc_charge(float pH) {
  //
  // Calculate the fractional charges at this pH
  //
  // Initialise random number generator
  //
  srand(time(NULL));
  //
  // Get a random starting state
  //
  vector<int> current_state;
  vector<int> try_state;
  vector<int> sum_state;
  for (int group=0;group<_groups;group++) {
    current_state.push_back(static_cast<int>(rand()%2));
    if (current_state[group]==2) {
      current_state[group]=1;
    }
    // 
    // Dummy initialisation of try_state and sum_state
    //
    try_state.push_back(0);
    sum_state.push_back(0);
  }
  //
  // Get the energy of the starting state 
  //
  double current_energy=get_energy(pH,current_state);
  //
  // Start the MC loop
  //
  int eqsteps=_MCsteps/10;
  int keep=0;
  double tilf=0.0;
  double try_energy_new=0.0;
  double diff=0.0;
  int CCPS_count=0;
  for (int step=0;step<_MCsteps;step++) {
    //
    // Copy the current state to trystate
    //
    for (int count=0;count<_groups;count++) {
      try_state[count]=current_state[count];
    }
    //
    // Change a random group
    //
    int rand_group=static_cast<int>(rand()%_groups);
    try_state[rand_group]=abs(try_state[rand_group]-1);
    //
    // Get the energy of the new state
    //
    try_energy_new=get_energy_fast(pH,try_state,rand_group,current_energy);
    //
    // Keep or reject?
    //
    diff=try_energy_new-current_energy;
    if (diff<0.0) {
      //
      // Keep
      //
      for (int count=0;count<_groups;count++) {
	current_state[count]=try_state[count];
      }
      current_energy=try_energy_new;
    } else {
      if (diff<20.0) {
	tilf=static_cast<double>(rand()) / (static_cast<double>(RAND_MAX)+static_cast<double>(1));
	if (tilf<exp(-diff)) {
	  //
	  // Keep
	  //
	  for (int count=0;count<_groups;count++) {
	    current_state[count]=try_state[count];
	  }
	  current_energy=try_energy_new;
	}
      }
    }
    //
    // Record the state if we have equilibrated
    //
    if (step>eqsteps) {
      for (int count=0;count<_groups;count++) {
        sum_state[count]=sum_state[count]+current_state[count];
      }
      //
      // If we are monitoring a specific CCPS then record it
      //
      int found=1;
      if (_do_monitor==1) {
        for (int count=0;count<_groups;count++) {
            if (_monitor_state[count]!=current_state[count] && _monitor_state[count]!=-99) {
                found=0;
            }
        }
        if (found==1) {
            CCPS_count=CCPS_count+1;
        }
      }
    }
  }
  //
  // Calculate fractional charge
  //
  int sample_steps=_MCsteps-eqsteps;
  vector<float> charges_thispH;
  for (int count=0;count<_groups;count++){
    float charge=0.0;
    charge=static_cast<float>(_acid_base[count]);
    charge=charge*static_cast<float>(sum_state[count]);
    charge=charge/(static_cast<float>(sample_steps));
    charges_thispH.push_back(charge);
  }
  if (_do_monitor==1) {
    _CCPS_population.push_back(static_cast<double>(pH));
    _CCPS_population.push_back(static_cast<double>(CCPS_count)/static_cast<double>(sample_steps));
    }

  return charges_thispH;
}

vector<float> MC::calc_charge_Boltzmann(float pH) {
  //
  // Calculate the charge on each group at this pH by
  // evaluating the Boltzmann sum
  //
  // Loop over all states and store energy
  //
  double denominator=0.0;
  vector<float> energies;
  int numstates=static_cast<int>(pow(2,_groups));
  for (int state_number=0;state_number<numstates;state_number++) {
    vector<int> state_spec=get_state_specification(state_number,_groups);
    float energy=get_energy(pH,state_spec);
    energies.push_back(energy);
    // Update the denominator
    denominator=denominator+exp(-energy);
  }
  //
  // Calculate the fractional charge of each group
  //
  vector<float> charge;
  for (int group=0;group<_groups;group++) {
    charge.push_back(0.0);
  }
  for (int state_number=0;state_number<numstates;state_number++) {
    double population=exp(-energies[state_number])/denominator;
    vector<int> state_spec=get_state_specification(state_number,_groups);
    for (int group=0;group<_groups;group++) {
      if (state_spec[group]==1){
	charge[group]=charge[group]+population*static_cast<float>(_acid_base[group]);
      }
    }
  }
  return charge; 
}

void MC::print_vector(vector<int> state) {
  //
  // Print a state
  //
  int length=static_cast<int>(state.size());
  for (int pos=0;pos<length;pos++) {
    printf ("%d ",state[pos]);
  }
  printf ("\n");
  printf ("Returning frmo print_vector\n");
  return;
}

vector<int> MC::get_state_specification(int state_number,int num_groups) {
  //
  // given the state number return the state specification
  //
  //printf ("Int get_state_sepc\n");
  vector<int> state;
  for (int position=num_groups-1;position>-1;position--) {
    if (static_cast<int>(pow(2,position))>state_number) {
      state.push_back(0);
    } else {
      state.push_back(1);
      state_number=state_number-static_cast<int>(pow(2,position));
    }
  }
  //printf ("State number; %d",state_number);
  //print_vector(state);
  return state;
}

//
// --------------------
// 

double MC::get_energy(float pH,vector<int> state) {
  //
  // Calculate the energy of the present state
  //
  //printf ("Entering get_energy\n");
  double pH_value=static_cast<double>(pH);
  double energy=0.0;
  for (int group1=0;group1<_groups;group1++) {
    //
    // Add the energy from the intrinsic pKa
    //
    if (state[group1]==1) {
      energy=energy+_acid_base[group1]*lnten*(pH_value-_intpKas[group1]);
      //
      // Add the charged-charged energies
      //
      for (int group2=0;group2<_groups;group2++) {
      	if (state[group2]==1 and group2!=group1) {
       energy=energy+_matrix[group1][group2]/2.0;
      	}
      }
    }
  }
  //printf ("Exiting get_energy\n");
  return energy;
}
	
//
// --------------------
// 
  
double MC::get_energy_fast(float pH,vector<int> state,int change_group,double old_energy) {
  //
  // Calculate the energy of the present state
  //
  int startpointer;
  double energy=old_energy;
  double energy_diff=0.0;
  //
  // Add the energy from the intrinsic pKa
  //
  energy_diff=_acid_base[change_group]*lnten*(pH-_intpKas[change_group]);
  //
  // Add the charged-charged energies
  //
  for (int group2=0;group2<_groups;group2++) {
    if (state[group2]==1 && group2!=change_group) {
      energy_diff=energy_diff+_matrix[change_group][group2];
    }
  }
  //
  // Should we add or subtract the energy
  //
  if (state[change_group]==1) {
    //
    // Group became charged - add energy
    //
    energy=energy+energy_diff;
  } else {
    // 
    // Group became uncharged
    //
    energy=energy-energy_diff;
  }
  return energy;
}

//
// ---------------------------
//


void MC::set_monitor_states(vector<int> mon_state) {
  //
  // Set the state that we will monitor in the simulation
  //
  _do_monitor=1;
  _monitor_state=mon_state;
  //print_vector(_monitor_state);
  //printf ("Hello CCPS monitoring ok\n");
  return;
}




