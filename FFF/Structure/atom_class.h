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
#ifndef ATOM_CLASS_H
#define ATOM_CLASS_H
#include <string>
#include <stdio.h>
#include <map>
#include <cstdlib>
//
// Local include files
//
#include "string_tools.h"
#include "vector3d.h"

class atom_class : public vector3d {
public:
  //
  // Functions
  //
  atom_class() { 
    name="uninitialized";
    tag="uninitialized";
  }
  //
  atom_class (string atomname,double xcoord, double ycoord, double zcoord) {
    name=atomname;
    x=xcoord;
    y=ycoord;
    z=zcoord;
    initialize();
    set_vdw();
  }
  atom_class (string atomname,double xcoord, double ycoord, double zcoord,float bf) {
    name=atomname;
    x=xcoord;
    y=ycoord;
    z=zcoord;
    bfactor=bf;
    initialize();
    set_vdw();
  }
    // 
    // For keeping charges and radii (parameter storage)
    //
    atom_class (string atomname,string crg,string rad) {
        name=atomname;
        charge=atof(crg.c_str());
        radius=atof(rad.c_str());
        return;
    }
  // Functions
  void initialize();
  void set_vdw();
  //
  atom_class (string line);
  //
  // Variables
  //
  float bfactor;
  int inresidue, inchain, number;
  double radius, charge, vdw;
  string name, inresiduename, pdbchainname, pdbresnum, _element;
  string tag; //for identifying ligands etc.  
  int index;
  bool is_donor, is_acceptor;
  map<string,double> attr;
  //
  // Bonds
  //
  vector<int> bonds;
  vector<atom_class*> covalent_bonds;
  //
  void update();
  //
  // Other functions
  //
  bool is_hydrogen() {  
    //
    // Is the atom a hydrogen?
    //
    name=strip(name);
    string first(name,0,1);
    if (first=="H") {
      return true;
    }
    else {
      if (first=="1" || first=="2" || first=="3") {
	string second(name,1,1);
	if (second=="H") {
	  return true;
	}
      }
    }
    return false;
  }
  //
  // Other flags
  //
  bool is_backbone() const;
  bool isacceptor() const;
  bool isdonor() const;
  bool is_metal() const;
  //
  // Print function
  //
  void print() const;
  void print_detail() const;
};

//inline atom_class operator==(const atom_class atom1,const atom_class atom2) {
    
inline bool operator==(const atom_class& atom1,const atom_class& atom2) {
  if (atom1.inresidue==atom2.inresidue &&
      atom1.inchain==atom2.inchain &&
      atom1.name==atom2.name) {
    //printf ("returning true\n");
    return true;
  }
  //printf ("Returning false\n");
  return false;
}

inline bool both_in_same_sidechain(const atom_class& atom1, const atom_class& atom2) {
  if (atom1.inresidue==atom2.inresidue &&
      atom1.inchain==atom2.inchain &&
      not atom1.is_backbone() && 
      not atom2.is_backbone()) {
    //printf ("returning true\n");
    return true;
  }
  //printf ("Returning false\n");
  return false;
}   

inline bool both_in_same_residue(const atom_class& atom1, const atom_class& atom2) {
  if (atom1.inresidue==atom2.inresidue &&
      atom1.inchain==atom2.inchain) {
    return true;
  }
  return false;
}   
#endif
