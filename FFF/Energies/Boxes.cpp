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
#include "Boxes.h"

//Boxes::Boxes(FFF &P,float boxsize) {
  //cout << "box begin"<<endl<<flush;
//  atomsA= A->all_atoms;
//  make_boxes(atomsA,boxsize);
//}

Boxes::Boxes(vector<atom_class*> atomsA,float boxsize) {
  //
  // Call make_boxes directly
  //
  make_boxes(atomsA,boxsize);
  return;
}
//
// ----
//

void Boxes::make_boxes(vector<atom_class*> atomsA,float boxsize) {
  //
  // Get extent of protein
  //
  _boxsize=boxsize;
  //_boxsize=100;
  xmax=-9999;
  xmin=9999;
  ymax=-9999;
  ymin=9999;
  zmax=-9999;
  zmin=9999;
    //printf ("Creating box");
  for (unsigned int atom=0;atom<atomsA.size();atom++) {
    double X=atomsA[atom]->x;
    double Y=atomsA[atom]->y;
    double Z=atomsA[atom]->z;
    xmax=max(xmax,X);
    ymax=max(ymax,Y);
    zmax=max(zmax,Z);
    xmin=min(xmin,X);
    ymin=min(ymin,Y);
    zmin=min(zmin,Z);
      //atomsA[atom]->print_detail();
  }
  //
  // Box size 
  //
  xnum=int((xmax-xmin)/_boxsize)+1;
  ynum=int((ymax-ymin)/_boxsize)+1;
  znum=int((zmax-zmin)/_boxsize)+1;
  //printf ("xnum: %d, ynum: %d, znum %d\n",xnum,ynum,znum);
  //
  vector<atom_class *> dummy;
  //int count=0;
  for (unsigned int x=0;x<static_cast<unsigned int>((xnum+1)*(ynum+1)*(znum+1));x++) {
    boxes.push_back(dummy);
    //count++;
  }
  //cout <<"Total number of boxes:"<<count<<endl;
  //
  for (unsigned int atom=0;atom<atomsA.size();atom++) {
    int Xn=int((atomsA[atom]->x-xmin)/_boxsize);
    int Yn=int((atomsA[atom]->y-ymin)/_boxsize);
    int Zn=int((atomsA[atom]->z-zmin)/_boxsize);
    int boxnum=xnum*ynum*Zn+xnum*Yn+Xn;
    //if (boxnum>=count||boxnum<0) {
    //  printf ("xnum: %d, ynum: %d, znum %d\n",Xn,Yn,Zn);
    //  cout << "HEEELP "<<boxnum<<endl;
    //}
    boxes[boxnum].push_back(atomsA[atom]);
  }
  //cout<<"box end"<<endl<<flush;
  return;
}

//
// --------------
//

vector<atom_class *> Boxes::get_close_atoms(atom_class& atom) {
    return get_close_atoms(&atom);
}

vector<atom_class *> Boxes::get_close_atoms(atom_class * atom) {
  //
  // Find atoms that are in the box + neighbouring boxes
  //
  vector<atom_class *> close_atoms;
  int Xn=int(((atom->x)-xmin)/_boxsize);
  int Yn=int(((atom->y)-ymin)/_boxsize);
  int Zn=int(((atom->z)-zmin)/_boxsize);
  //
  for (int x=max(Xn-1,0);x<=min(Xn+1,xnum);x++) {
    for (int y=max(Yn-1,0);y<=min(Yn+1,ynum);y++) {
      for (int z=max(Zn-1,0);z<=min(Zn+1,znum);z++) {
	int boxnum=xnum*ynum*z+xnum*y+x;	
	for (unsigned int count=0;count<boxes[boxnum].size();count++) {
	  close_atoms.push_back(boxes[boxnum][count]);
	}
      }
    }
  }
  //printf ("%d close atoms\n",close_atoms.size());
  return close_atoms;
}
  

