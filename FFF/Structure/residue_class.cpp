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
#include "residue_class.h"
#include "chain_class.h"
#include "fff.h"

//
//
//
using namespace std;
//
//
//
residue_class::residue_class(chain_class* C,vector<string> lines): _C(C) {
   string resname(lines[0],17,3);
   resname=strip(resname);
   string c_resnum(lines[0],22,5);
   //
   // Loop over all lines
   //
   for (unsigned int count=0;count<lines.size();count++) {
     string new_resname(lines[0],17,3);
     string new_c_resnum(lines[0],22,5);
     new_resname=strip(resname);
     //
     // If the residue number changes, then we have a problem
     //
     if (new_resname!=resname || c_resnum!=new_c_resnum) {
       printf ("Residue is corrupt!:\n");
       for (unsigned int ii=0;count<lines.size();count++) {
	 cout << lines[ii] << endl;
       }
       exit(0);
     }
     atom_class X(lines[count]);
     atoms.push_back(X);
   }
   //
   // Set the residue name
   //
   name=resname;
   //
   // and the number that the residue has in the pdb file
   //
   string pnum(lines[0],22,6);
   pdbnum=strip(pnum);
   //
   // Call the update function
   //
   //update();
   return;
};

//
// ----------------------------------------------------------------------
//

void residue_class::update() {
  //
  // See if we should delete any of the atoms
  //
  vector<atom_class> newatoms;
  unsigned int numatoms=atoms.size();
  for (unsigned int atom=0;atom<numatoms;atom++) {
    if (atoms[atom].name=="DELETE") {
        printf ("Deleting atom: %s\n",atoms[atom].name.c_str());
      for (unsigned int at2=atom;at2<(atoms.size()-1);at2++) {
          atoms[at2]=atoms[at2+1];
      }
      atoms.resize(atoms.size()-1);
      atom--;
      numatoms=atoms.size();
    }
  }
  //
  // Set a few pointers
  //
  for (unsigned int at=0;at<atoms.size();at++) {
    if (atoms[at].name=="N") {
      N=&(atoms[at]);
    }
    else if (atoms[at].name=="CA") {
      CA=&(atoms[at]);
    }
    else if (atoms[at].name=="C") {
      C=&(atoms[at]);
    }
    //
    // Set the inresidue pointer
    //
    atoms[at].inchain=inchain;
    atoms[at].inresidue=number;
    atoms[at].inresiduename=name;
    atoms[at].number=at;
      atoms[at].pdbresnum=pdbnum;
      atoms[at].pdbchainname=strip(chainname);
      //printf ("Assigning chain name: '%s'\n",strip(chainname).c_str());
      atoms[at].update();
  }
      
  
  return;
}

void residue_class::print() {
    printf ("Residue C: %d R: %d %s\n",inchain,number,name.c_str());
    for (unsigned int at=0;at<atoms.size();at++) {
        atoms[at].print_detail();
    }
}
