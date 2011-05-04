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
#ifndef HYDROGEN_AMB_H
#define HYDROGEN_AMB_H
#include <string>
#include <stdio.h>
#include <map>
#include <cstdlib>
//
// Local include files
//
#include "string_tools.h"
#include "vector3d.h"
#include "atom_class.h"

class Hydrogen_Conformation_class {
public:
    Hydrogen_Conformation_class(vector<string> lines) {
        _name=lines[0].substr(1,lines[0].size());
        vector<string> words;
        // First line is atom name, attachment point and bond length
        printf ("Line1: %s\n",lines[1].c_str());
        unsigned int count_s=1;
        words=split(strip(lines[count_s]));
        if (strip(lines[count_s])=="DEFAULT STATE") {
            //
            // If this is the default state then we don't build anything
            //
            _build=false;
            return;
        } else {
            if (words[0]=="BUILD:") {
                //
                // If there's a BUILD statement, then we build other hydrogens too
                // There might or might not be any hydrogens to be build in the definition itself
                //
                for (unsigned int os=1;os<words.size();os++) {
                    _build_other_conformations.push_back(words[os]);
                }
                count_s++;
                if (count_s>=lines.size()) {
                    _build=false;
                    return;
                }
                words=split(strip(lines[count_s]));
            }
            if (words[0]=="DELETE:") {
                //
                // If there's a DELETE statement, then we delete the conformations specified
                // There might or might not be any hydrogens to be build in the definition itself
                //
                for (unsigned int os=1;os<words.size();os++) {
                    _delete_other_conformations.push_back(words[os]);
                }
                count_s++;
                if (count_s>=lines.size()) {
                    _build=false;
                    return;
                }
                words=split(strip(lines[count_s]));
            }
        }
        if (strip(lines[count_s])=="END") {
            _build=false;
            return;
        }
        // More lines so there must be something to build
        count_s++;
        _build=true; // Should anything be built?
        printf ("LK: %s\n",lines[count_s-1].c_str());
        _add_atomname=words[0];
        _anchor_atom=words[1];
        _bondlength=atof(words[2].c_str());
        for (unsigned int count=count_s;count<lines.size();count++) {
            if (strip(lines[count])=="END") {
                break;
            }
            words=split(lines[count]);
            _atoms.push_back(atom_class(words[0],atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str())));
        }
        return;
    }
    string _name;
    string _add_atomname;
    string _anchor_atom;
    float _bondlength;
    bool _build;
    vector<string> _build_other_conformations, _delete_other_conformations;
    vector<atom_class> _atoms;
};

class HydrogenAmb_class {
public:
    HydrogenAmb_class(vector<string> lines) {
        //
        // Parse a single group entry for a hydrogen ambiguity definition
        //
        unsigned int count=0;
        string line=strip(lines[count]);
        if (line!="*") {
            printf ("Incorrect start line: %s\n",line.c_str());
            exit(0);
        }
        count++;
        line=strip(lines[count]);
        vector<string> words=split(line);
        _TitID=words[0];
        printf ("Parsing Hydrogen Ambiguity: %s\n",_TitID.c_str());
        // AcidBase
        if (strip(words[1])=="acid") {
            _acidbase=-1;
        } else if (strip(words[1])=="base") {
            _acidbase=1;
        } else if (strip(words[1])=="None") {
            _acidbase=0;
        } else {
            printf ("Incorrect acid/base specification: %s'n",words[1].c_str());
            exit(0);
        }
        // Number of conformations
        _num_conformations=atoi(words[3].c_str());
        //
        // Optcode
        //
        _optcode=atof(words[7].c_str());
        count++;
        //count++;
        //
        // Parse each of the conformations
        //
        vector<string> sublines;
        for (unsigned int count2=count;count2<lines.size();count2++) {
            line=strip(lines[count2]);
            if (line.substr(0,1)==">" && sublines.size()>1) {
                _conformations.push_back(Hydrogen_Conformation_class(sublines));
                sublines.resize(0);
                sublines.push_back(line);
            } else {
                sublines.push_back(line);
            }
        }
        _conformations.push_back(Hydrogen_Conformation_class(sublines));
        return;
    }
    string _TitID;
    int _acidbase;
    float _optcode;
    int _num_conformations;
    vector<Hydrogen_Conformation_class> _conformations;
    vector<string> _conformation_names;
};
        
#endif