#ifndef BOXES_H
#define BOXES_H
#include <math.h>
#include "method.h"
#include "distances.h"
#include "topology.h"
#include "bondmatrix.h"
#include "bfactor.h"

using namespace std;

class Boxes {

 public:
  Boxes(Soup * A,float boxsize); // Boxes of a soup
  Boxes(vector<Atom *> atomsA,float boxsize); //Boxes of a vector of atoms
  void make_boxes(vector<Atom *> atomsA,float boxsize); // Routine that does the work
  vector<int> get_close_atoms(Atom * atom); // Get the atoms that are close  

  int xnum,ynum,znum;
  float xmax,xmin,ymax,ymin,zmax,zmin;
 protected:
  vector<Atom *> atomsA, atomsB; 
  vector<vector<int> > boxes;
  float _boxsize;
};

#endif
