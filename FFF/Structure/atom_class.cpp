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
#include "atom_class.h"
atom_class::atom_class(string line) {
  //
  // Get the atom name
  //
  string N(line,12,4); 
    N=strip(N);
  name=N;
  //
  // Get X, Y, Z coordinates
  //
  string X(line,30,8); X=strip(X);
  x=atof(X.c_str());
  string Y(line,38,8); Y=strip(Y);
  y=atof(Y.c_str());
  string Z(line,46,8); Z=strip(Z);
  z=atof(Z.c_str());
  string BF(line,60,5);
  bfactor=atof(BF.c_str());
  set_vdw();
  initialize();
  // Do we have a tag? Always initialise with last word
  vector<string> sp=split(line);
  tag=strip(sp.back());
}

void atom_class::initialize() {
    // Initialize covalent_bonds
    covalent_bonds.resize(0);
    // Initialize the rest
    inchain=-999;
    inresidue=-999;
    tag="uninitialized";
    return;
}

void atom_class::update() {
    // Update by setting the vdw radii
    set_vdw();
    return;
}

void atom_class::set_vdw() {
    // This needs to be redone based on a parameter file
    string element(strip(name),0,1);
    if (name!="DUM") {
        if (element=="1" || element=="2" || element=="3") {
            string el2(strip(name),1,1);
            element=el2;
        }
        vdw=0.0;
        if (element=="C") vdw=1.7;
        if (element=="H") vdw=1.2;
        if (element=="N") vdw=1.55;
        if (element=="O") vdw=1.52;
        if (element=="S") vdw=1.8;
        if (element=="P") vdw=1.9;
        if (vdw==0.0) {
            if (name=="MN") {
                vdw=2.0;
                element="MN";
            }
            if (vdw==0.0) {
            printf ("unknown element: %s in: '%s'\n",element.c_str(),strip(name).c_str());
            }
        }
    }
    _element=element; // set the element
    // Set all the flags
    is_donor=isdonor();
    is_acceptor=isacceptor();
}


void atom_class::print_detail() const {
    printf ("Chain: '%1s', Residue: %4s %4s, Atom: %4s  %8.3f %8.3f %8.3f tag: %s\n",pdbchainname.c_str(),
                                                                    pdbresnum.c_str(),
                                                                    inresiduename.c_str(),
	    name.c_str(),x,y,z,
	    tag.c_str());
}
void atom_class::print() const {
  printf ("Atom: %4s  %7.3f %7.3f %7.3f\n",name.c_str(),x,y,z);
  return;
}

bool atom_class::is_backbone() const {
  //
  // Is this atom a backbone atom
  //
  string test(strip(name));
  if (test=="N" || test=="CA" || test=="C" || test=="O" 
      || test=="O2" || test=="HA" || test=="HN" || test=="H"
      || test=="tN") {
    return true;
  }
  return false;
}

//
// --------
//

bool atom_class::isacceptor() const {
    if (_element=="O" || 
        (name=="ND1" && inresiduename=="HIS") ||
        (name=="NE2" && inresiduename=="HIS")) {
        return true;
    }
    if (is_metal()) return true;
    if (name=="O" && inresiduename=="HOH") return true;
    return false;
}

bool atom_class::is_metal() const {
    if (_element=="MN") return true;
    return false;
}
//
// ---------
//

bool atom_class::isdonor() const {
    if (_element=="N") return true;
    if (_element=="O" and name!="O") return true; // This is only very approximately true
    return false;
}
