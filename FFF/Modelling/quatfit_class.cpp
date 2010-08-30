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
#include "quatfit.h"

//
// Constructor
//
quatfit_class::quatfit_class() {
  vector<double> rotmatrix;
  center_on_transform=true;
}

//quatfit_class::~quatfit_class() {
  //
  // Destructor
  //
//}

void quatfit_class::fit(int numpoints,vector<atom_class> refpoints,vector<atom_class> fitpoints) {
  //
  // Calls quatfit and get the rotation matrix <rotmatrix> and the translation vector <tvector>
  //
  double* weight = new double[MAXPOINTS];
  double Q[4];
  double R[4][4];
  int maxsweeps=30;
  double drefpoints[4][MAXPOINTS];
  double dfitpoints[4][MAXPOINTS];
  //
  // Initialise
  //
  for (int i=0;i<4;i++) {
    fit_center[i]=0.0;
    ref_center[i]=0.0;
    Q[i]=0.0;
    R[i][0]=0.0;
    R[i][1]=0.0;
    R[i][2]=0.0;
    R[i][3]=0.0;
  }
  //
  // The Quatfit arrays all start at 1 for some reason...
  //
  for (int i=1;i<=numpoints;i++) {
    weight[i]=1.0;
    drefpoints[1][i]=static_cast<double>(refpoints[i-1].x);
    drefpoints[2][i]=static_cast<double>(refpoints[i-1].y);
    drefpoints[3][i]=static_cast<double>(refpoints[i-1].z);
    dfitpoints[1][i]=static_cast<double>(fitpoints[i-1].x);
    dfitpoints[2][i]=static_cast<double>(fitpoints[i-1].y);
    dfitpoints[3][i]=static_cast<double>(fitpoints[i-1].z);
  }
  //
  // Translate the center of refpoints to 0,0,0
  //
  center(numpoints,drefpoints,weight, 1,ref_center);
  //printf ("Refcenter returns: %7.3f %7.3f %7.3f\n",ref_center[1],ref_center[2],ref_center[3]);
  //
  // Translate center of fitpoints to 0,0,0
  //
  center(numpoints,dfitpoints,weight,1,fit_center);
  //
  // Calculate the Rotation matrix
  //
  qtrfit(numpoints,dfitpoints,drefpoints,weight,Q,R,maxsweeps);
  for (int p1=1;p1<=3;p1++) {
    for (int p2=1;p2<=3;p2++) {
      rotmatrix.push_back(R[p1][p2]);
    }
  }
  //
  // Rotate the fitted coordinates
  //
  rotmol(numpoints,dfitpoints,dfitpoints,R);
  //
  // Move the center of the fitted coordinates to the center of the reference
  //
  center(numpoints,dfitpoints,weight,3,ref_center);
    delete []weight;
  return;
}

//
// ---------------------------------------------------
//

void quatfit_class::transform_all(FFF &Protein) {
  //
  // Transform all atoms in the molecule
  //
  vector <atom_class> transform_coords;
  for (unsigned int chain=0;chain<Protein.chains.size();chain++) {
    for (unsigned int residue=0;residue<Protein.chains[chain].residues.size();residue++) {
      for (unsigned int atom=0;atom<Protein.chains[chain].residues[residue].atoms.size();atom++) {
	transform_coords.push_back(Protein.chains[chain].residues[residue].atoms[atom]);
      }
    }
  }
  _transform(transform_coords);
  //
  // Change the real protein coordinates
  //
  int count=0;
  for (unsigned int chain=0;chain<Protein.chains.size();chain++) {
    for (unsigned int residue=0;residue<Protein.chains[chain].residues.size();residue++) {
      for (unsigned int atom=0;atom<Protein.chains[chain].residues[residue].atoms.size();atom++) {
	Protein.chains[chain].residues[residue].atoms[atom].x=transform_coords[count].x;
	Protein.chains[chain].residues[residue].atoms[atom].y=transform_coords[count].y;
	Protein.chains[chain].residues[residue].atoms[atom].z=transform_coords[count].z;
	count++;
      }
    }
  }
  return;
}

//
// ----------------------------------------------------
//
void quatfit_class::_transform(vector<atom_class> & fit_coords) {
  //
  // Apply rotmatrix and vector3d to fit_coords
  // This routine actually modifies atom_class instances!!
  //
  double* weight = new double[MAXPOINTS];
  double R[4][4];
  double fitpoints[4][MAXPOINTS];
  //
  // Initialise
  //
  unsigned int numpoints=fit_coords.size();
  if (numpoints>=MAXPOINTS) {
    printf ("\nError: Too many points to transform: %d\n",numpoints);
    printf ("Current value of MAXPOINTS is: %d, increase MAXPOINTS in quatfit.h and recompile.\n",MAXPOINTS);
    printf ("Array overflow. Exiting\n\n");
    exit(0);
  }
  //
  // The Quatfit arrays all start at 1 for some reason...
  //
  for (unsigned int i=1;i<=numpoints;i++) {
    weight[i]=1.0;
    fitpoints[1][i]=static_cast<double>(fit_coords[i-1].x);
    fitpoints[2][i]=static_cast<double>(fit_coords[i-1].y);
    fitpoints[3][i]=static_cast<double>(fit_coords[i-1].z);
  }
  //
  // Fill R
  //
  int count=0;
  for (int p1=1;p1<=3;p1++) {
    for (int p2=1;p2<=3;p2++) {
      R[p1][p2]=rotmatrix[count];
      R[p1][p2]=rotmatrix[count];
      R[p1][p2]=rotmatrix[count];
      count++;
    }
  }
  //
  // Translate center of fitpoints to 0,0,0
  //
  center(numpoints,fitpoints,weight,1,fit_center);
  //
  // Rotate the fitted coordinates
  //
  rotmol(numpoints,fitpoints,fitpoints,R);
  //
  // Move the center of the fitted coordinates back to where it was
  //
  center(numpoints,fitpoints,weight,3,fit_center);
  //
  // Put the new coordinates in the result array
  //
  for (unsigned int i=1;i<=numpoints;i++) {
    fit_coords[i-1].x=static_cast<float>(fitpoints[1][i]);
    fit_coords[i-1].y=static_cast<float>(fitpoints[2][i]);
    fit_coords[i-1].z=static_cast<float>(fitpoints[3][i]);
  }
  //
  // Done
  //
    delete []weight;
  return;
}

//
// ----------------------------------------------------
//

void quatfit_class::transform(int numpoints,vector<atom_class> & fit_coords) {
  //
  // Apply rotmatrix and vector3d to fit_coords
  // This routine actually modifies atom_class instances!!
  //
  double* weight = new double[MAXPOINTS];
  double R[4][4];
  double fitpoints[4][MAXPOINTS];
  //
  // Initialise
  //
  if (numpoints>=MAXPOINTS) {
    printf ("\nError: Too many points to transform: %d\n",numpoints);
    printf ("Current value of MAXPOINTS is: %d, increase MAXPOINTS in quatfit.h and recompile.\n",MAXPOINTS);
    printf ("Array overflow. Exiting\n\n");
    exit(0);
  }
  //
  // The Quatfit arrays all start at 1 for some reason...
  //
  for (int i=1;i<=numpoints;i++) {
    weight[i]=1.0;
    fitpoints[1][i]=static_cast<double>(fit_coords[i-1].x);
    fitpoints[2][i]=static_cast<double>(fit_coords[i-1].y);
    fitpoints[3][i]=static_cast<double>(fit_coords[i-1].z);
  }
  //
  // Fill R
  //
  int count=0;
  for (int p1=1;p1<=3;p1++) {
    for (int p2=1;p2<=3;p2++) {
      R[p1][p2]=rotmatrix[count];
      R[p1][p2]=rotmatrix[count];
      R[p1][p2]=rotmatrix[count];
      count++;
    }
  }
  //
  // Translate center of fitpoints to 0,0,0
  //
  if (center_on_transform) {
    center(numpoints,fitpoints,weight,2,fit_center);
  }
  //
  // Rotate the fitted coordinates
  //
  rotmol(numpoints,fitpoints,fitpoints,R);
  //
  // Move the center of the fitted coordinates to the center of the reference
  //
  if (center_on_transform) {
    center(numpoints,fitpoints,weight,3,ref_center);
  }
  //
  // Put the new coordinates in the result array
  //
  for (int i=1;i<=numpoints;i++) {
    fit_coords[i-1].x=static_cast<float>(fitpoints[1][i]);
    fit_coords[i-1].y=static_cast<float>(fitpoints[2][i]);
    fit_coords[i-1].z=static_cast<float>(fitpoints[3][i]);
  }
  //
  // Done
  //
    delete []weight;
  return;
}

void quatfit_class::rotate_around_vector(vector3d V,double angle) {
  //
  // Create a rotation matrix that gives an <angle> degree rotation around the vector V 
  // V is specified by xin,yin,zin
  //
  double L[4];
  double R[4][4];
  double radangle;
  //
  // Initialise
  //
  for (int i=0;i<4;i++) {
    L[i]=0.0;
    R[i][0]=0.0;
    R[i][1]=0.0;
    R[i][2]=0.0;
    R[i][3]=0.0;
  }
  //
  // Convert angle to radians
  //
  radangle=angle;
  radangle=deg2rad(radangle);
  //
  // Normalise V and put it in L
  //
    V.Normalize();
  L[1]=static_cast<double>(V.x);
  L[2]=static_cast<double>(V.y);
  L[3]=static_cast<double>(V.z);
  //
  // Construct the rotation matrix
  //
  R[1][1]=cos(radangle) + L[1]*L[1] * (1.0-cos(radangle));
  R[2][2]=cos(radangle) + L[2]*L[2] * (1.0-cos(radangle));
  R[3][3]=cos(radangle) + L[3]*L[3] * (1.0-cos(radangle));
  R[1][2]=L[1]*L[2]*(1.0-cos(radangle))-L[3]*sin(radangle);
  R[1][3]=L[1]*L[3]*(1.0-cos(radangle))+L[2]*sin(radangle);
  R[2][1]=L[2]*L[1]*(1.0-cos(radangle))+L[3]*sin(radangle);
  R[2][3]=L[2]*L[3]*(1.0-cos(radangle))-L[1]*sin(radangle);
  R[3][1]=L[3]*L[1]*(1.0-cos(radangle))-L[2]*sin(radangle);
  R[3][2]=L[3]*L[2]*(1.0-cos(radangle))+L[1]*sin(radangle);
  //
  // Fill the class variables
  //
  for (int p1=1;p1<=3;p1++) {
    // No translation
    tvector.push_back(0.0);
    for (int p2=1;p2<=3;p2++) {
      rotmatrix.push_back(R[p1][p2]);
    }
  }
  //
  // Set the center_on_transform flag to false, since the user has to subtract the basis vector
  //
  center_on_transform=false;
  return;
}


