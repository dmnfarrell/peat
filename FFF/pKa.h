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
#ifndef pKa_H
#define pKa_H
#include <vector>
#include <string>
#include <iostream>
#include <exception>

#include "Model.h"
#include "fff.h"
#include "Boxes.h"

class chargestate_class;

class pKa_class : public model_class {
public:
    pKa_class(FFF &P,Rotamer_class &ROT,const std::string aadef_dir):model_class(P,ROT,aadef_dir) {
        printf ("Instantiating pKa class\n");
        read_titgroup_definitions();
        find_titratable_groups();
	set_standard_protonation_state(7.0);
        //
        // Initialise random number generator
        //
        srand(time(NULL));
    }
    //
    void find_titratable_groups();
    chargestate_class calculate_fractional_charge(unsigned int MCSTEPS);
    vector<chargestate_class> calculate_titration_curves(unsigned int MCSTEPS);
    //
    // 
private:
    // 
    // Functions
    //
    void change_prot_or_H_state();
    void change_Hstate();
    void change_protstate();
    double get_current_energy();
    void guess_best_protonation_state();
    void set_standard_protonation_state(float pH);
    void test_protstate();
    void record_protstate();
    void read_titgroup_definitions(); // Read the definition file for the titratable groups
    void model_current_state(vector<int> protonation_state,vector<int> rotamer_state); // Model the current protonation state (i.e. make changes to the structure and assign charges)
    
    //
    // Variables
    //
    unsigned int _numtitgroups; // Number of titratable groups in system
    vector<chargestate_class> titration_curve;
    vector<int> _current_PSstate;
    vector<int> _try_PSstate;
    vector<int> _sum_PSstate;
    vector<int> _rotamer_state; // Vector holding a specification of the current rotamer state
    //
    // Misc
    //
    double tilf, diff, try_energy, try_energy_new;
};

#endif

