#include "Boxes.h"

Boxes::Boxes(Soup * A,float boxsize) {
  //cout << "box begin"<<endl<<flush;
  atomsA= A->get_dry_atoms();
  make_boxes(atomsA,boxsize);
}

Boxes::Boxes(vector<Atom*> atomsA,float boxsize) {
  //
  // Call make_boxes directly
  //
  make_boxes(atomsA,boxsize);
  return;
}
//
// ----
//

void Boxes::make_boxes(vector<Atom*> atomsA,float boxsize) {
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
  for (unsigned int atom=0;atom<atomsA.size();atom++) {
    float X=atomsA[atom]->x;
    float Y=atomsA[atom]->y;
    float Z=atomsA[atom]->z;
    xmax=max(xmax,X);
    ymax=max(ymax,Y);
    zmax=max(zmax,Z);
    xmin=min(xmin,X);
    ymin=min(ymin,Y);
    zmin=min(zmin,Z);
  }
  //
  // Box size 
  //
  xnum=int((xmax-xmin)/_boxsize)+1;
  ynum=int((ymax-ymin)/_boxsize)+1;
  znum=int((zmax-zmin)/_boxsize)+1;
  //printf ("xnum: %d, ynum: %d, znum %d\n",xnum,ynum,znum);
  //
  vector<int> dummy;
  //int count=0;
  for (unsigned int x=0;x<(xnum+1)*(ynum+1)*(znum+1);x++) {
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
    boxes[boxnum].push_back(atom);
  }
  //cout<<"box end"<<endl<<flush;
  return;
}

//
// --------------
//

vector<int> Boxes::get_close_atoms(Atom * atom) {
  //
  // Find atoms that are in the box + neighbouring boxes
  //
  vector<int> close_atoms;
  int Xn=int(((atom->x)-xmin)/_boxsize);
  int Yn=int(((atom->y)-ymin)/_boxsize);
  int Zn=int(((atom->z)-zmin)/_boxsize);
  //
  for (unsigned int x=max(Xn-1,0);x<=min(Xn+1,xnum);x++) {
    for (unsigned int y=max(Yn-1,0);y<=min(Yn+1,ynum);y++) {
      for (unsigned int z=max(Zn-1,0);z<=min(Zn+1,znum);z++) {
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
  

