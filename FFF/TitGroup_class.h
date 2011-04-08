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
#ifndef TITGROUP_H
#define TITGROUP_H
#include <vector>
#include <string>
#include <iostream>
#include <exception>

// ==========================

class Titration_curve_class {
    //
    // Holds a single titration curve
    //
public:
    string name;
    vector<double> pH;
    vector<double> population;
};

// ==========================
class Transition_class {
    //
    // Transition_class holds information on a particular transition
    //
public:
    Transition_class(string start_state, string end_state, float modelpKa): _start_state(start_state),_end_state(end_state),_modelpKa(modelpKa) {
        // Instantiate
        print();
        return;
    }
    void print() {
        printf ("Start: %s, End: %s, model pKa: %6.3f\n",_start_state.c_str(),_end_state.c_str(),_modelpKa);
    }
    string _start_state, _end_state;
    double _modelpKa;
    double _intpKa;
    vector<Titration_curve_class> _population;
};

// ==========================

class TitGroup_class {
    //
    // TitGroup class holds information on a particular Titratable group
public:
    //
    TitGroup_class() {};
    //
    TitGroup_class(TitGroup_class &TGDEF,int chain,int residue) {
        // Copy most from Definition
        _TitID=TGDEF._TitID;
        _residuename=TGDEF._residuename;
        _chain=chain;
        _residue=residue;
        _acidbase=TGDEF._acidbase;
        _avg_modelpKa=TGDEF._avg_modelpKa;
        _HydrogenAmbDef_number=TGDEF._HydrogenAmbDef_number;
        //
        // Copy the transitions
        //
        for (unsigned int ts=0;ts<TGDEF._transitions.size();ts++) {
            _transitions.push_back(TGDEF._transitions[ts]);
        }
        // Copy neutral states
        for (unsigned int count=0;count<TGDEF.neutral_states.size();count++) {
            neutral_states.push_back(TGDEF.neutral_states[count]);
        }
        for (unsigned int count=0;count<TGDEF.charged_states.size();count++) {
            charged_states.push_back(TGDEF.charged_states[count]);
        }
    }
    //
    // ----------------------
    //
    TitGroup_class(vector<string> lines) {
        // Provide all lines for a TitGroup definition
        // Class will parse and instantiate sub-classes
        int count=0;
        string line=strip(lines[count]);
        string TitID(line.substr(2,5));
        _TitID=TitID;
        printf ("titid: %s\n",_TitID.c_str());
        //
        // Get the residue name
        //
        count++;
        line=lines[count];
        vector<string> words0=split(line,":");
        _residuename=strip(words0[1]);
        //
        // Get the acidbase definition
        //
        count++;
        line=lines[count];
        vector<string> words=split(line);
        _acidbase=0;
        if (words[1]=="Acid") {
            _acidbase=-1;
        } else if (words[1]=="Base") {
            _acidbase=1;
        };
        if (_acidbase==0) {
            printf ("Incorrect acidbase spec: %s\n",words[1].c_str());
            exit(0);
        }
        //
        // Now for the transitions
        //
        count++;
        line=lines[count];
        vector<string> words2=split(line,":");
        string startstates_s=split(words2[1],">")[0];
        vector<string> startstates=split(startstates_s,",");
        string endstates_s=split(words2[1],">")[1];
        vector<string> endstates=split(endstates_s,",");
        //
        // Model pKas
        //
        count++;
        line=lines[count];
        string modpkas=split(line,":")[1];
        vector<float> modelpka_values=getfloats(modpkas);
        //
        // Instantiate all the transitions
        //
        //_transitions.resize(0);
        float sum_modpka=0.0;
        unsigned int modpka_count=0;
        for (unsigned int s_state=0;s_state<startstates.size();s_state++) {
            for (unsigned int e_state=0;e_state<endstates.size();e_state++) {
                _transitions.push_back(Transition_class(strip(startstates[s_state]),
                                                        strip(endstates[e_state]),
                                                        modelpka_values[modpka_count]));
                sum_modpka=sum_modpka+modelpka_values[modpka_count];
                modpka_count++;
            }
        }
        _avg_modelpKa=sum_modpka/static_cast<float>(modpka_count);
        // 
        // Get a list of all the states
        //
        for (unsigned int count=0;count<startstates.size();count++) {
            allstates.push_back(strip(startstates[count]).c_str());
        }
        for (unsigned int count=0;count<endstates.size();count++) {
            allstates.push_back(strip(endstates[count]));
        }
        //
        // Check that we have an END statement
        //
        count++;
        line=strip(lines[count]);
        if (line!="END") {
            printf ("Incorrect or missing END statement: %s\n",line.c_str());
            exit(0);
        }
    }
    void print() {
        printf ("Titratable group: %s, chain: %d, residue: %d, residuename: %s\n",
                _TitID.c_str(),_chain,_residue,_residuename.c_str());
        return;
    }
    //
    // Variables
    //
    string _TitID; //TitID links to the corresponding entry in HYDROGENS.DAT
    string _residuename;
    int _chain; // the chain
    int _residue; // the residue
    vector<Transition_class> _transitions; // The pKas for the neutral->charged transitions as specified by start_states and end_states.
    vector<string> allstates, charged_states, neutral_states; // All states that this titgroup can occupy
    Titration_curve_class _Titration_curve;
    double _pKa, _avg_modelpKa;
    int _acidbase; // acidbase is -1 for acids and +1 for bases
    int _HydrogenAmbDef_number; // Pointer to the definition in _HydrogenAmbDef
};

//
// -------------
//

class protonation_state_class {
    //
    // Holds information on the fractional protonation state of all groups
    //
public:
    protonation_state_class () {
        // Nothing happening here
    };
    vector<string> TitID;
    vector<double> frac_charge;
};
    
#endif
