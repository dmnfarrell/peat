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
#include "TitGroup_class.h"
#include "HydrogenAmb.h"


class pKa_class : public model_class {
public:
    pKa_class(FFF &P,Rotamer_class &ROT,const std::string aadef_dir):model_class(P,ROT,aadef_dir) {
        printf ("Instantiating pKa class\n");
        // Read TITRATION.DAT
        read_titgroup_definitions(aadef_dir);
        //
        // Read HYDROGENS.DAT - file specifying how to build alternative protonation states
        //
        read_hydrogen_definitions(aadef_dir);
        //
        // Check the connection between HYDROGENS.DAT, TITRATION.DAT and AA.CRGRAD
        // Find out which states are charged and neutral, and figure out 
        // If we have charge distributions for all states
        //
        printf ("Checking the parameter files are internally consistent\n");
        int typefound=0;
        for (unsigned int tgD=0;tgD<_TitGroupDefs.size();tgD++) {
            for (unsigned int count=0;count<_TitGroupDefs[tgD].allstates.size();count++) {
                // 
                // See if we can find this conformation in HYDROGENS.DAT
                //
                bool found_inHYD=false;
                for (unsigned int deftype=0;deftype<_HydrogenAmbDefs.size();deftype++) {
                    if (_HydrogenAmbDefs[deftype]._TitID==_TitGroupDefs[tgD]._TitID) {
                        for (unsigned int confnumber=0;confnumber<_HydrogenAmbDefs[deftype]._conformations.size();confnumber++) {
                            //printf ("'%s' - HDEF: '%s'\n",_HydrogenAmbDefs[deftype]._conformations[confnumber]._name.c_str(),_TitGroupDefs[tgD].allstates[count].c_str());
                            if (_HydrogenAmbDefs[deftype]._conformations[confnumber]._name==_TitGroupDefs[tgD].allstates[count]) {
                                found_inHYD=true;
                                _TitGroupDefs[tgD]._HydrogenAmbDef_number=deftype;
                            }
                        }
                    }
                }
                if (!found_inHYD) {
                    printf ("Could not find conformation info for state %s in TitID: %s in HYDROGENS.DAT. The state was specified in TITRATION.DAT.\nCheck TITRATION.DAT and HYDROGENS.DAT\n",
                            _TitGroupDefs[tgD].allstates[count].c_str(),_TitGroupDefs[tgD]._TitID.c_str());
                    exit(0);
                }
                //
                // See if we can find this conformation in AA.CRGRAD
                //
                bool found=false;
                for (unsigned int typeno=0;typeno<(*_crgrad).types.size();typeno++) {
                    if ((*_crgrad).types[typeno].name.c_str()==_TitGroupDefs[tgD].allstates[count]) { 
                        found=true;
                        typefound=typeno;
                        break;
                    }
                }
                if (!found) {
                    printf ("Could not find charge/radius entry for state %s specified in TITRATION.DAT.\nCheck TITRATION.DAT and AA.CRGRAD\n",_TitGroupDefs[tgD].allstates[count].c_str());
                    exit(0);
                }
                //
                // Calculate the net charge of this state
                //
                double sum_crg=0.0;
                for (unsigned int atom=0;atom<(*_crgrad).types[typefound].atoms.size();atom++) {
                    sum_crg=sum_crg+(*_crgrad).types[typefound].atoms[atom].charge;
                }
                //printf ("titgroup: %s, state: %s, charge: %5.2f \n",_TitGroupDefs[tgD]._TitID.c_str(),
                //       _TitGroupDefs[tgD].allstates[count].c_str(),
                //        sum_crg);
                if (sum_crg<0.01 && sum_crg>-0.01) {
                    _TitGroupDefs[tgD].neutral_states.push_back(_TitGroupDefs[tgD].allstates[count]);
                } else {
                    _TitGroupDefs[tgD].charged_states.push_back(_TitGroupDefs[tgD].allstates[count]);
                }
            }
        }
        printf ("Parameter files are ok\n");
        //
        // Find Titratable groups
        //
        find_titratable_groups();
        // Set the standard protonation state at pH 7.0
        set_standard_protonation_state(7.0);
        //
        // Initialise random number generator
        //
        srand(time(NULL));
    }
    //
    void find_titratable_groups();
    protonation_state_class calculate_fractional_charge(unsigned int MCSTEPS);
    vector<protonation_state_class> calculate_titration_curves(unsigned int MCSTEPS);
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
    void read_titgroup_definitions(string parmdir); // Read the definition file for the titratable groups
    void read_hydrogen_definitions(string parmdir);
    void model_current_state(vector<int> protonation_state,vector<int> rotamer_state); // Model the current protonation state (i.e. make changes to the structure and assign charges)
    void Build_Protonation_State(vector<string> ProtonationState);
    void Build_Group_ProtonationState(int TG,string ProtonationState);
    void Build_Group_Conformation(int Tg,string ProtonationState);
    void ReBuild_TitGroup_Hydrogens(); // Delete and rebuild all hydrogens of 
    //
    // Variables
    //
    unsigned int _numtitgroups; // Number of titratable groups in system
    vector<protonation_state_class> titration_curve;
    vector<string> _current_Protonation_state;
    vector<int> _current_PSstate;
    vector<int> _try_PSstate;
    vector<int> _sum_PSstate;
    vector<int> _rotamer_state; // Vector holding a specification of the current rotamer state
    vector<TitGroup_class> _TitGroupDefs;
    vector<TitGroup_class> _TitGroups;
    vector<HydrogenAmb_class> _HydrogenAmbDefs;
    //
    // Misc
    //
    double tilf, diff, try_energy, try_energy_new;
};

#endif

