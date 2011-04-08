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
    // 
    // find titratable groups as defined in TITRATION.DAT
    //
    printf ("Finding titratable groups\n");
    _TitGroups.resize(0);
    for (unsigned int chain=0; chain<_P.chains.size();chain++) {
        for (unsigned int residue=0;residue<_P.chains[chain].residues.size();residue++) {
            for (unsigned int tg=0;tg<_TitGroupDefs.size();tg++) {
                if (_P.chains[chain].residues[residue].name==_TitGroupDefs[tg]._residuename) {
                    _TitGroups.push_back(TitGroup_class(_TitGroupDefs[tg],chain,residue));
                    printf ("Found: %s in chain: %4d residue: %3d %4s\n",
                            _TitGroups.back()._TitID.c_str(),
                            chain,
                            residue,
                            _P.chains[chain].residues[residue].name.c_str());
                } else if (_TitGroupDefs[tg]._residuename=="CTERM" && _P.chains[chain].residues[residue].is_Cterm && _P.chains[chain].residues[residue].is_aa) {
                    _TitGroups.push_back(TitGroup_class(_TitGroupDefs[tg],chain,residue));
                    printf ("Found: %s in chain: %4d residue: %3d %4s\n",
                            _TitGroupDefs[tg]._TitID.c_str(),
                            chain,
                            residue,
                            _P.chains[chain].residues[residue].name.c_str());
                } else if (_TitGroupDefs[tg]._residuename=="NTERM" && _P.chains[chain].residues[residue].is_Nterm && _P.chains[chain].residues[residue].is_aa) {
                    _TitGroups.push_back(TitGroup_class(_TitGroupDefs[tg],chain,residue));
                    printf ("Found: %s in chain: %4d residue: %3d %4s\n",
                            _TitGroupDefs[tg]._TitID.c_str(),
                            chain,
                            residue,
                            _P.chains[chain].residues[residue].name.c_str());
                }
            }
        }
    }
    return;
}

//
// ----------
//

vector<protonation_state_class> pKa_class::calculate_titration_curves(unsigned int MCSTEPS) {
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

protonation_state_class pKa_class::calculate_fractional_charge(unsigned int MCSTEPS) {
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
            for (unsigned int count=0;count<_numtitgroups;count++) {
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
    protonation_state_class PSC;
    return PSC;
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
// --------------
//

void pKa_class::set_standard_protonation_state(float pH) {
  //
  // Set the standard protonation state for this pH
  //
    _current_Protonation_state.resize(0);
    for (unsigned int tg=0;tg<_TitGroups.size();tg++) {
        _TitGroups[tg].print();
        if (_TitGroups[tg]._acidbase==-1) {
            if (pH>_TitGroups[tg]._avg_modelpKa) {
                _current_Protonation_state.push_back(_TitGroups[tg].charged_states[0]); 
            } else {
                _current_Protonation_state.push_back(_TitGroups[tg].neutral_states[0]); 
            }
        } else if (_TitGroups[tg]._acidbase==1) {
            if (pH<=_TitGroups[tg]._avg_modelpKa) {
                _current_Protonation_state.push_back(_TitGroups[tg].charged_states[0]); 
            } else {
                _current_Protonation_state.push_back(_TitGroups[tg].neutral_states[0]); 
            }
        } else {
            printf ("acidbase is not 1 or -1.: %d\n",_TitGroups[tg]._acidbase);
        }
    }
    Build_Protonation_State(_current_Protonation_state);
    printf ("Set standard protonation state at pH: %5.2f \n",pH);
    return;
}



//
//-------------
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
    for (unsigned int count=0;count<_numtitgroups;count++) {
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
    for (unsigned int count=0;count<_numtitgroups;count++) {
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
// -----
//

void pKa_class::Build_Protonation_State(vector<string> ProtonationState) {
    //
    // Build the atoms required to construct the protonation state
    //
    ReBuild_TitGroup_Hydrogens();
    for (unsigned int TG=0;TG<_TitGroups.size();TG++) {
        printf ("==============\n");
        _TitGroups[TG].print();
        printf ("Protonation state: %s\n",ProtonationState[TG].c_str());
        Build_Group_ProtonationState(TG,ProtonationState[TG]);
    }
    return;
}

//
// ----------
//

void pKa_class::ReBuild_TitGroup_Hydrogens() {
    //
    // Delete all hydrogens in titratable groups
    //
    _use_hydrogens=true;
    int chain, residue;
    for (unsigned int TG=0;TG<_TitGroups.size();TG++) {
        chain=_TitGroups[TG]._chain;
        residue=_TitGroups[TG]._residue;
        delete_hydrogens(_P.chains[chain].residues[residue]); // Delete
        Mutate_2(chain,residue,_P.chains[chain].residues[residue].name,1,10.0); //Build
    }
    _P.update_all_atoms();
    return;
}

void pKa_class::Build_Group_ProtonationState(int TG,string ProtonationState) {
    //
    // Build the given protonation state (a string) for the TitGroup in question
    //
    // First find the definition
    //
    int HA_def=_TitGroups[TG]._HydrogenAmbDef_number;
    for (unsigned int conf=0;conf<_HydrogenAmbDefs[HA_def]._conformations.size();conf++) {
        if (_HydrogenAmbDefs[HA_def]._conformations[conf]._name==ProtonationState) {
            //int chain=_TitGroups[TG]._chain;
            //int residue=_TitGroups[TG]._residue;
            // Build others
            
            // Delete others
            //for (unsigned int count=0;count<delete_other_conformations.size()l;count++) {
                // Get the atom name in the other definition
                
              //  remove_atom(chain,residue,atomname);
            // Build this one
        }
            //vector<atom_class> superpos_atoms=_HydrogenAmbDefs[HA_def]._conformations[conf].get_superpos_atoms();
    }
return;
}

//
// ------------
//

void pKa_class::read_titgroup_definitions(string parmdir) {
    //
    // Read TITRATION.DAT
    //
    string titfile=parmdir+"/TITRATION.DAT";
    printf ("Reading %s\n",titfile.c_str());
    ifstream file(titfile.c_str());
    if (!file) {
        printf ("ERROR: File %s was not found!\n",titfile.c_str());
        exit(0);
    }
    string line;
    vector<string> lines;
    while(getline(file,line)) {
        lines.push_back(line);
    }
    file.close();
    printf ("%s has been read\n",titfile.c_str());
    //
    // Parse the lines
    //
    bool done=false;
    unsigned int count=0;
    string name;
    //
    // Define variables
    //
    bool record=false;
    vector<string> thisgroup;
    //
    // 
    while (!done) {
        line=strip(lines[count]);
        // Comments
        if (!line.substr(0,2).compare("//") || !line.substr(0,1).compare("#"))  {
            count++;
            continue;
        }
        // Start record
        if (line.substr(0,1).compare("*")==0 && count<=lines.size()-1) {
            thisgroup.resize(0);
            thisgroup.push_back(line);
            record=true;
            count++;
            continue;
        }
        // End record
        if (line.substr(0,3).compare("END")==0 && strip(line)!="END OF FILE") { 
            thisgroup.push_back(line);
            _TitGroupDefs.push_back(TitGroup_class(thisgroup));
            printf ("Parsing: %s\n",name.c_str());
            count++;
            record=false;
            continue;
        }
        // Normal record
        if (record) {
            thisgroup.push_back(line);
            count++;
            continue;
        }
        // End of file
        if (count>=lines.size() || strip(line)=="END OF FILE") {
            done=true;
        }
        count++;
    }
    printf ("Done parsing TITRATION.DAT\n");
    return;
}
    
//
// -----------
//

void pKa_class::read_hydrogen_definitions(string parmdir) {
    //
    // Read HYDROGENS.DAT
    //
    string titfile=parmdir+"/HYDROGENS.DAT";
    printf ("Reading %s\n",titfile.c_str());
    ifstream file(titfile.c_str());
    if (!file) {
        printf ("ERROR: File %s was not found!\n",titfile.c_str());
        exit(0);
    }
    string line;
    vector<string> lines;
    while(getline(file,line)) {
        lines.push_back(line);
    }
    file.close();
    printf ("%s has been read\n",titfile.c_str());
    //
    // Parse the lines
    //
    bool done=false;
    unsigned int count=0;
    string name;
    //
    // Define variables
    //
    bool record=false;
    vector<string> thisgroup;
    //
    // 
    while (!done) {
        line=strip(lines[count]);
        if (strip(line).size()==0) {
            count++;
            continue;
        }
        // Comments
        if (!line.substr(0,2).compare("//") || !line.substr(0,1).compare("#"))  {
            count++;
            continue;
        }
        // Start record
        if (line.substr(0,1).compare("*")==0 && count<=lines.size()-1) {
            thisgroup.resize(0);
            thisgroup.push_back(line);
            record=true;
            count++;
            continue;
        }
        // End record
        if (line.substr(0,3).compare("END")==0 && strip(line)!="END OF FILE") { 
            thisgroup.push_back(line);
            _HydrogenAmbDefs.push_back(HydrogenAmb_class(thisgroup));
            count++;
            record=false;
            continue;
        }
        // Normal record
        if (record) {
            thisgroup.push_back(line);
            count++;
            continue;
        }
        // End of file
        if (count>=lines.size() || strip(line)=="END OF FILE") {
            done=true;
        }
        count++;
    }
    printf ("Done parsing HYDROGENS.DAT\n");
    return;
}
