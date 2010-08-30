/*
 #
 # FFF - Flexible Force Field
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

#include "pKa.h"


void pKa_class::find_titratable_groups() {
    printf ("Finding titratable groups\n");
    return;
}

vector<chargestate_class> pKa_class::calculate_titration_curves(unsigned int MCSTEPS) {
    //
    // Calculate pKa values 
    //
    _numtitgroups=_P._titratable_groups.size();
    printf ("Calculating pKas\n");
    //
    // Loop over all pH values
    //
    for (double pH=0.0;pH<=20.0;pH=pH+0.1) {
        _P.set_pH(pH); // Set the pH for the system
        //
        // Calculate the fractional charge for each group at this pH
        //
        titration_curve.push_back(calculate_fractional_charge(MCSTEPS));
    }
    return titration_curve;
}

//
// ---------
// 

chargestate_class pKa_class::calculate_fractional_charge(unsigned int MCSTEPS) {
    //
    // Calculate the fractional protonation state at this pH
    //
    // Clear all counters
    //
    _current_PSstate.clear();
    _try_PSstate.clear();
    _sum_PSstate.clear();
    //
    // Set starting protonation state
    //
    guess_best_protonation_state();
    //
    // Get the energy of the starting state 
    // 
    double current_energy=get_current_energy();
    //
    // First do 0.1 MCSTEPS equilibration steps and then sample for MCSTEPS steps
    //
    unsigned int record_point=static_cast<unsigned int>(0.1*static_cast<float>(MCSTEPS));
    unsigned int totsteps=MCSTEPS+record_point;
    for (unsigned int step=0;step<totsteps;step++) {
        //
        // Change something (in a clever way)
        //
        change_prot_or_H_state();
        try_energy_new=get_current_energy();
        //
        // Keep or reject?
        //
        diff=try_energy_new-current_energy;
        if (diff<0.0) {
            //
            // Keep
            //
            for (int count=0;count<_numtitgroups;count++) {
                _current_PSstate[count]=_try_PSstate[count];
            }
            current_energy=try_energy_new;
        } else {
            if (diff<20.0) {
                tilf=static_cast<double>(rand()) / (static_cast<double>(RAND_MAX)+static_cast<double>(1));
                if (tilf<exp(-diff)) {
                    //
                    // Keep
                    //
                    for (unsigned int count=0;count<_numtitgroups;count++) {
                        _current_PSstate[count]=_try_PSstate[count];
                    }
                    current_energy=try_energy_new;
                }
            }
        }
        //
        // Record if we have equilibrated
        //
        if (step>record_point) {
            record_protstate();
        }
    }
    //
    // Now calculate the fractional charge from the statistics
    //
    printf ("Write code for calculating fractional charge\n");
}

//
// --------------
//
    
double pKa_class::get_current_energy() {
    //
    // Get the energy of the system in the current configuration
    //
    _P.update_BOXLJ();  // Always update boxes when something has changed
    vector<double> energies=get_soup_energy(); // Get the full energy of the system
    return energies[0];
}

//
// -----
//

void pKa_class::guess_best_protonation_state() {
    //
    // Get a random starting state (for now)
    //
    for (unsigned int group=0;group<_numtitgroups;group++) {
        _current_PSstate.push_back(static_cast<int>(rand()%2));
        if (_current_PSstate[group]==2) {
            _current_PSstate[group]=1;
        }
        // 
        // Dummy initialisation of try_state and sum_state
        //
        _try_PSstate.push_back(0);
        _sum_PSstate.push_back(0);
    }
    return;    
}

//
// -----
//


void pKa_class::change_prot_or_H_state() {
    //
    // Find out what to change
    //
    return;
}

void pKa_class::change_Hstate() {
    // 
    // Change a hydrogen position or flip a His/Asn/Gln
    //
    return;
}

void pKa_class::change_protstate() {
    //
    // Change a protonation state
    //
    // Copy the current state to trystate
    //
    for (int count=0;count<_numtitgroups;count++) {
        _try_PSstate[count]=_current_PSstate[count];
    }
    //
    // Change a random group
    //
    int rand_group=static_cast<int>(rand()%_numtitgroups);
    _try_PSstate[rand_group]=abs(_try_PSstate[rand_group]-1);
    //
    // Set the state
    //
    model_current_state(_try_PSstate,_rotamer_state);
    return;
}
    
//
// ------------
//

void pKa_class::record_protstate() {
    //
    // Record the protonation state
    //
    for (int count=0;count<_numtitgroups;count++) {
        _sum_PSstate[count]=_sum_PSstate[count]+_current_PSstate[count];
    }
    return;
}

//
// -----------
//

void pKa_class::model_current_state(vector<int> PS_state,vector<int> Rotamer_state) {
    // 
    // Model the current protonation state/Hbond conformation and make changes to the structure/charges etc
    //
    return;
}

//
// ------------
//

void pKa_class::read_titgroup_definitions() {
}
