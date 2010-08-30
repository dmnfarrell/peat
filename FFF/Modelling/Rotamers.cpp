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
#include "Rotamers.h"

using namespace std;

Rotamer_class::Rotamer_class (const std::string filename) {
    //
    // Read the rotamer library
    //
    string line;
    string name;
    vector<vector<float > > residue_temp;
    name="none";
    printf("Starting to read rotamer library: %s\n",filename.c_str());
    ifstream file(filename.c_str());
    if (!file) {
        printf("Rotamer library not found!: %s\n",filename.c_str());
        exit(0);
    }
    int jcount=0;
    // Put the entire file in lines
    while(getline(file,line)) {
        vector<string> words=split(line);
        // Get the name right and skip comment lines
        if (words[0]=="#") continue;
        if (words[0]!=name) {
            if (name!="none") {
                _rotamers.push_back(residue_temp);
                residue_temp.resize(0);
                _names.push_back(name);
            }
            name=words[0];
        }
        //
        vector<float> this_rotamer;
        for (unsigned int count=1;count<words.size();count++) {
            //float value;
            //from_string<float>(value,words[count],std::dec);
            //
            this_rotamer.push_back(atof(words[count].c_str()));
        }
        if (this_rotamer.size()==0) {
            printf ("storing empty rotamer\n");
        }
        residue_temp.push_back(this_rotamer);
        jcount++;
    }
    // Push back the last residue
    _rotamers.push_back(residue_temp);
    _names.push_back(name);
    // Close the file
    file.close();
    printf ("Read %d rotamers\n",jcount);
    for (unsigned int j=0;j<_names.size();j++) {
        printf ("Name: %s # rotamers: %d\n",_names[j].c_str(),static_cast<int>(_rotamers[j].size()));
    }
};

vector<vector<float> > Rotamer_class::get_rotamers(string resname,double phi, double psi) {
    //
    // Get the rotamers that fit the bill
    //
    unsigned int rot_residue=9999;
    for (unsigned int j=0;j<_names.size();j++) {
        if (resname.compare(_names[j])==0) {
            rot_residue=j;
            break;
        }
    }
    if (rot_residue>9998) {
        printf ("Request for invalid rotamer residue: %s\n",resname.c_str());
        exit(0);
    }
    vector<vector<float > > use_rotamers;
    for (unsigned int rotamer=0;rotamer<_rotamers[rot_residue].size();rotamer++) {
        //printf ("length: %d\n",_rotamers[rot_residue][rotamer].size());
        if (fabs(_rotamers[rot_residue][rotamer][0]-phi)<=10.0 && fabs(_rotamers[rot_residue][rotamer][1]-psi)<=10.0) {
            use_rotamers.push_back(_rotamers[rot_residue][rotamer]);
	    //for (unsigned int j=0;j<_rotamers[rot_residue][rotamer].size();j++) {
	    //  printf ("%6.2f ",_rotamers[rot_residue][rotamer][j]);
	    //}
	    //printf ("\n");
        }
    }
    return use_rotamers;
}

    
    
