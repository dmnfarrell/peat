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
#include "assign.h"

using namespace std;

crgrad_class::crgrad_class (const std::string parmdir) {
    // Read a charge or radius file
    string filename=parmdir+"/AA.CRGRAD";
    printf ("Reading charges and radii from: %s \n",filename.c_str());
    ifstream file(filename.c_str());
    if (!file) {
        printf ("ERROR: File %s was not found!\n",filename.c_str());
        exit(0);
    }
    string line;
    vector<string> lines;
    while(getline(file,line)) {
        lines.push_back(line);
    }
    file.close();
    //
    // Parse
    //
    unsigned int count=0;
    while (count<lines.size()) {
        if (strip(lines[count])=="#BEGIN") {
            while (true) {
                if (strip(lines[count])=="*") {
                    count++;
                    string name=strip(lines[count]);
                    count++;
                    // Get the chargelines for this name
                    vector<string> parm_lines;
                    while (strip(lines[count])!="") {
                        parm_lines.push_back(strip(lines[count]));
                        count++;
                    }
                    types.push_back(residue_class(name,parm_lines));
                }
                if (strip(lines[count])=="#END") {
                    break;
                }
                count++;
            }
        }
        count++;
    }
}


//
// -----------
//

void crgrad_class::assign_all(residue_class residue) {
    printf ("Assigning charge and radius\n");
    return;
}
