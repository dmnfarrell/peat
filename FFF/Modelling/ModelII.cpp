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
#include "Model.h"

using namespace std;

void model_class::setphi(int chainnumber, int residuenumber,double newangle) {
  //
  // Set the phi angle for a residue
  //
  // C-1, N, CA, C 
  // 
  int prev_res=_get_prev_residue(chainnumber,residuenumber);
  if (prev_res!=-1) {
    vector<atom_class> tors_atoms;
    tors_atoms.push_back(_get_atom(prev_res,"C"));
    tors_atoms.push_back(_get_atom(residuenumber,"N"));
    tors_atoms.push_back(_get_atom(residuenumber,"CA"));
    tors_atoms.push_back(_get_atom(residuenumber,"C"));
    double oldphi=Dihedral(tors_atoms[0],tors_atoms[1],
			   tors_atoms[2],tors_atoms[3]);
    vector3d V;
    int basis=1;
    int point=2;
    V=tors_atoms[point]-tors_atoms[basis];
    quatfit_class Q;
    double difchi=newangle-oldphi;
    Q.rotate_around_vector(V,difchi);
    //
    // Find all the atoms that should be moved
    //
    vector<atom_class> transform_coords;
    int numpoints=0;
    //
    // First fill in the atoms in the present residue
    // For phi this starts with CA + sidechain
    //
    vector<int> moveatoms;
    vector<int> first_last=_get_first_last(residuenumber);
    printf ("Atoms to be transformed in first residue\n");
    for (int atom=first_last[0];atom<=first_last[1];atom++) {
      string atomname(_P.pdb[atom].atomtype);
      atomname=strip(atomname);
      if (atomname!="N" && atomname!="tN") {
	atom_class X(atomname,
		     _P.pdb[atom].xcoord-tors_atoms[basis].x,
		     _P.pdb[atom].ycoord-tors_atoms[basis].y,
		     _P.pdb[atom].zcoord-tors_atoms[basis].z);
	moveatoms.push_back(atom);
	transform_coords.push_back(X);
	numpoints++;
      }
    }
    //
    // And then we fill in all atoms to the C-term of this one
    //
    if (residuenumber<_P.residuenum) {
      for (int residue=residuenumber+1;residue<=_P.residuenum;residue++) {
    	vector<int> first_last=_get_first_last(residue);
    	for (int atom=first_last[0];atom<=first_last[1];atom++) {
    	  string atomname(_P.pdb[atom].atomtype);
    	  atomname=strip(atomname);
    	  atom_class X(atomname,
    		       _P.pdb[atom].xcoord-tors_atoms[basis].x,
    		       _P.pdb[atom].ycoord-tors_atoms[basis].y,
    		       _P.pdb[atom].zcoord-tors_atoms[basis].z);
    	  transform_coords.push_back(X);
    	  moveatoms.push_back(atom);
    	  numpoints++;
    	}
      }
    }
    //
    // Transform the coordinates
    //
    Q.transform(numpoints,transform_coords);
    //
    // Put the new coordinates in the array
    //
    int counter=0;
    for (int n=0;n<static_cast<int>(moveatoms.size());n++) {
      _P.pdb[moveatoms[n]].xcoord=transform_coords[counter].x+tors_atoms[basis].x;
      _P.pdb[moveatoms[n]].ycoord=transform_coords[counter].y+tors_atoms[basis].y;
      _P.pdb[moveatoms[n]].zcoord=transform_coords[counter].z+tors_atoms[basis].z;
      counter++;
    }
  }
  //
  // Calculate the new phi angle
  //
  vector<atom_class> tors_atoms2;
  tors_atoms2.push_back(_get_atom(prev_res,"C"));
  tors_atoms2.push_back(_get_atom(residuenumber,"N"));
  tors_atoms2.push_back(_get_atom(residuenumber,"CA"));
  tors_atoms2.push_back(_get_atom(residuenumber,"C"));
  printf ("Torsion atoms 2:\n");
  for (unsigned int i=0;i<tors_atoms2.size();i++) {
    tors_atoms2[i].print();
  }
  double newphi=Dihedral(tors_atoms2[0],tors_atoms2[1],
  			 tors_atoms2[2],tors_atoms2[3]);
  //
  // Update the value in FFmeta
  //
  _P.sequence[residuenumber].phi=newphi;
  return;
}

//
// ------------------------------------------------
//


//
// ------------------------------------------------
//

void model_class::setpsi(int residuenumber,double newangle) {
  //
  // Set the psi angle for a residue
  //
  int next_res=_get_next_residue(residuenumber);
  //
  // Set the psi angle for a residue
  //
  // N, CA, C, N+1 
  // 
  if (next_res!=-1) {
    vector<atom_class> tors_atoms;
    tors_atoms.push_back(_get_atom(residuenumber,"N"));
    tors_atoms.push_back(_get_atom(residuenumber,"CA"));
    tors_atoms.push_back(_get_atom(residuenumber,"C"));
    tors_atoms.push_back(_get_atom(next_res,"N"));
    double oldpsi=Dihedral(tors_atoms[0],tors_atoms[1],
			   tors_atoms[2],tors_atoms[3]);
    vector3d V;
    int basis=1;
    int point=2;
    V=tors_atoms[point]-tors_atoms[basis];
    quatfit_class Q;
    double difchi=newangle-oldpsi;
    Q.rotate_around_vector(V,difchi);
    //
    // Find all the atoms that should be moved
    //
    vector<atom_class> transform_coords;
    int numpoints=0;
    //
    // First fill in the atoms in the present residue
    // For psi this is C and O
    //
    vector<int> moveatoms;
    vector<int> first_last=_get_first_last(residuenumber);
    for (int atom=first_last[0];atom<=first_last[1];atom++) {
      string atomname(_P.pdb[atom].atomtype);
      atomname=strip(atomname);
      if (atomname=="C" || atomname=="O" || atomname=="OXT") {
	atom_class X("DUM",
		     _P.pdb[atom].xcoord-tors_atoms[basis].x,
		     _P.pdb[atom].ycoord-tors_atoms[basis].y,
		     _P.pdb[atom].zcoord-tors_atoms[basis].z);
	moveatoms.push_back(atom);
	transform_coords.push_back(X);
	numpoints++;
      }
    }
    //
    // Fill in all atoms to the C-term of this one
    //
    if (residuenumber<_P.residuenum) {
      for (int residue=residuenumber+1;residue<=_P.residuenum;residue++) {
    	vector<int> first_last=_get_first_last(residue);
	for (int atom=first_last[0];atom<=first_last[1];atom++) {
	  atom_class X("DUM",
		       _P.pdb[atom].xcoord-tors_atoms[basis].x,
		       _P.pdb[atom].ycoord-tors_atoms[basis].y,
		       _P.pdb[atom].zcoord-tors_atoms[basis].z);
	  transform_coords.push_back(X);
	  moveatoms.push_back(atom);
	  numpoints++;
	}
      }
    }
    //
    // Transform the coordinates
    //
    Q.transform(numpoints,transform_coords);
    //
    // Put the new coordinates in the array
    //
    int counter=0;
    for (int n=0;n<static_cast<int>(moveatoms.size());n++) {
      _P.pdb[moveatoms[n]].xcoord=transform_coords[counter].x+tors_atoms[basis].x;
      _P.pdb[moveatoms[n]].ycoord=transform_coords[counter].y+tors_atoms[basis].y;
      _P.pdb[moveatoms[n]].zcoord=transform_coords[counter].z+tors_atoms[basis].z;
      counter++;
    }
  }
  vector<atom_class> tors_atoms2;
  tors_atoms2.push_back(_get_atom(residuenumber,"N"));
  tors_atoms2.push_back(_get_atom(residuenumber,"CA"));
  tors_atoms2.push_back(_get_atom(residuenumber,"C"));
  tors_atoms2.push_back(_get_atom(next_res,"N"));
  double newpsi=Dihedral(tors_atoms2[0],tors_atoms2[1],
  			 tors_atoms2[2],tors_atoms2[3]);
  //
  // Update the value in FFmeta
  //
  _P.sequence[residuenumber].psi=newpsi;
  return;
}  

//
// ------------------------------------------------
//

//
// ------------------------------------------------
//

void model_class::setomega(int residuenumber,double newangle) {
  //
  // Set the omega angle for a residue
  //
  int next_res=_get_next_residue(residuenumber);
  //
  // CA, C, N+1, CA+1 
  // 
  if (next_res!=-1) {
    vector<atom_class> tors_atoms;
    tors_atoms.push_back(_get_atom(residuenumber,"CA"));
    tors_atoms.push_back(_get_atom(residuenumber,"C"));
    tors_atoms.push_back(_get_atom(next_res,"N"));
    tors_atoms.push_back(_get_atom(next_res,"CA"));
    double oldpsi=Dihedral(tors_atoms[0],tors_atoms[1],
			   tors_atoms[2],tors_atoms[3]);
    vector3d V;
    int basis=1;
    int point=2;
    V=tors_atoms[point]-tors_atoms[basis];
    quatfit_class Q;
    double difchi=newangle-oldpsi;
    Q.rotate_around_vector(V,difchi);
    //
    // Find all the atoms that should be moved
    //
    vector<atom_class> transform_coords;
    int numpoints=0;
    vector<int> moveatoms;
    //
    // Fill in all atoms to the C-term of this one
    //
    if (residuenumber<_P.residuenum) {
      for (int residue=residuenumber+1;residue<=_P.residuenum;residue++) {
    	vector<int> first_last=_get_first_last(residue);
	for (int atom=first_last[0];atom<=first_last[1];atom++) {
	  atom_class X("DUM",
		       _P.pdb[atom].xcoord-tors_atoms[basis].x,
		       _P.pdb[atom].ycoord-tors_atoms[basis].y,
		       _P.pdb[atom].zcoord-tors_atoms[basis].z);
	  transform_coords.push_back(X);
	  moveatoms.push_back(atom);
	  numpoints++;
	}
      }
    }
    //
    // Transform the coordinates
    //
    Q.transform(numpoints,transform_coords);
    //
    // Put the new coordinates in the array
    //
    int counter=0;
    for (int n=0;n<static_cast<int>(moveatoms.size());n++) {
      _P.pdb[moveatoms[n]].xcoord=transform_coords[counter].x+tors_atoms[basis].x;
      _P.pdb[moveatoms[n]].ycoord=transform_coords[counter].y+tors_atoms[basis].y;
      _P.pdb[moveatoms[n]].zcoord=transform_coords[counter].z+tors_atoms[basis].z;
      counter++;
    }
  }
  //vector<atom_class> tors_atoms2;
  //tors_atoms2.push_back(_get_atom(residuenumber,"CA"));
  //tors_atoms2.push_back(_get_atom(residuenumber,"C"));
  //tors_atoms2.push_back(_get_atom(next_res,"N"));
  //tors_atoms2.push_back(_get_atom(next_res,"CA"));
  //double newpsi=Dihedral(tors_atoms2[0],tors_atoms2[1],
  //			 tors_atoms2[2],tors_atoms2[3]);
  return;
}  

//
// ------------------------------------------------
//

double model_class::getomega(int residuenumber) {
  //
  // Get the omega angle
  //
  int next_res=_get_next_residue(residuenumber);
  //
  // CA, C, N+1, CA+1 
  // 
  if (next_res!=-1) {
    vector<atom_class> tors_atoms;
    tors_atoms.push_back(_get_atom(residuenumber,"CA"));
    tors_atoms.push_back(_get_atom(residuenumber,"C"));
    tors_atoms.push_back(_get_atom(next_res,"N"));
    tors_atoms.push_back(_get_atom(next_res,"CA"));
    return Dihedral(tors_atoms[0],tors_atoms[1],
			   tors_atoms[2],tors_atoms[3]);
  }
  return -9999.9;
}

//
// ------------------------------------------------
//

vector<atom_class> model_class::init_build() {
  //
  // Build the first residue
  //
  string aatype("GLY");
  vector<atom_class> atoms;
  //
  // Find the definition of the residue
  //
  int defres=_get_AA_template(aatype);
  if (defres==-1) {
    //
    // Couldn't find the residue
    //
    printf ("I don't know of a residue called: %s\n",aatype.c_str());
    exit(0);
  }
  //
  // Fill the atoms into the vector
  //
  for (int atom=0;atom<static_cast<int>(AAsets[0].aadefs[defres].atoms.size());atom++) {
    if (strip(AAsets[0].aadefs[defres].atoms[atom].name)=="N" ||
	strip(AAsets[0].aadefs[defres].atoms[atom].name)=="CA" ||
	strip(AAsets[0].aadefs[defres].atoms[atom].name)=="C" ||
	strip(AAsets[0].aadefs[defres].atoms[atom].name)=="O") {
      atom_class X(strip(AAsets[0].aadefs[defres].atoms[atom].name),
		   AAsets[0].aadefs[defres].atoms[atom].x,
		   AAsets[0].aadefs[defres].atoms[atom].y,
		   AAsets[0].aadefs[defres].atoms[atom].z);
      atoms.push_back(X);
    }
  }
  return atoms;
}

vector<atom_class> model_class::add_C_residue(vector<atom_class> prev_residue) {
  //
  // Add a residue to the C-term of the atoms in prev_residue
  //
  // We superpos on the CA, C and O to get the backbone position right
  //
  //
  // Use this scaffold (this could probably be done in a smarter way)
  //
  vector<atom_class> scaffold;
  string SCA("CA");
  atom_class CA(SCA,0.001,0.001,0.001);
  string SC("C");
  atom_class C(SC,-1.249,0.882,0.001);
  string SO("O");
  atom_class O(SO,-2.184,0.661,-0.783);
  string SN("N");
  atom_class N1(SN,-1.294,1.910,0.867);
  atom_class CA1(SCA,-2.504,2.739,0.960);
  atom_class C1(SC,-2.732,3.466,-0.366);
  atom_class O1(SO,-1.943,4.102,-1.065);
  scaffold.push_back(N1);
  scaffold.push_back(CA1);
  scaffold.push_back(C1);
  scaffold.push_back(O1);
  //
  // Construct fitcoords 
  //
  vector<atom_class> fitcoords;
  fitcoords.push_back(CA);
  fitcoords.push_back(C);
  fitcoords.push_back(O);
  //
  // Construct refcoords
  //
  vector<atom_class> refcoords;
  // CA first
  for (unsigned int atom=0;atom<prev_residue.size();atom++) {
    if (strip(prev_residue[atom].name)=="CA") {
      refcoords.push_back(prev_residue[atom]);
      break;
    }
  }
  // C
  for (unsigned int atom=0;atom<prev_residue.size();atom++) {
    if (strip(prev_residue[atom].name)=="C") {
      refcoords.push_back(prev_residue[atom]);
      break;
    }
  }
  // O
  for (unsigned int atom=0;atom<prev_residue.size();atom++) {
    if (strip(prev_residue[atom].name)=="O") {
      refcoords.push_back(prev_residue[atom]);
      break;
    }
  }
  //
  // Find the rotation matrix and vector to put the scaffold CA, C and O on the
  // same atoms of the previous residue
  //
  quatfit_class Q;
  Q.fit(3,refcoords,fitcoords);
  //
  // Now transform all the scaffold coordinates
  //
  Q.transform(4,scaffold);
  //
  // And voila, we are done
  //
  return scaffold;
}
  
