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
#include "chain_class.h"
chain_class::chain_class(FFF* P,vector<string> lines): _P(P) {
  //printf ("\nLines passed to chain_class:\n");
  string last_resnum="DUMMY";
  vector<string> residue_lines;
  //
  // Get the chain name
  //
  string chainid(lines[0],21,1);
    //printf ("ChainID in chain_class: '%s'\n",chainid.c_str());
  name=strip(chainid);
  for (unsigned int ii=0;ii<lines.size();ii++) {
    //
    // Get the residue name and number
    //
    string resname(lines[ii],17,3);
    string resnum(lines[ii],22,5);
    if (resnum!=last_resnum) {
      //
      // new residue
      //
      if (last_resnum!="DUMMY") {
          //printf ("Instantiating residue_class\n");
	residue_class X(this,residue_lines);
          //printf ("done resclass\n");
	residues.push_back(X);
	residue_lines.resize(0);
      }
      last_resnum=resnum;
    }
    residue_lines.push_back(lines[ii]);
  }
  //
  // And the last residue...
  //
  if (residue_lines.size()!=0) {
    residue_class X(this,residue_lines);
    residues.push_back(X);
  }
  //printf("Chain name: %s found\n",name.c_str());
  return;
}

//
// --------------------------------------------------
//

void chain_class::update() {
  //
  // Update the numbers of all residues
  //
  for (unsigned int residue=0;residue<residues.size();residue++) {
    residues[residue].number=residue;
    residues[residue].inchain=number;
      residues[residue].chainname=strip(name);
    residues[residue].update();
  }
  return;
}
