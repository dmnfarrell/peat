/*
 #
 # UFFBAPS - Unified Force Field for Binding And Protein Stability
 # Copyright (C) 2010 Jens Erik Nielsen & Chresten Soendergaard
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
#ifndef TOPOLOGY_H
#define TOPOLOGY_H
#include <math.h>
#include "method.h"
#include "atom.h"
#include "generalparameters.h"

using namespace std;

class Topology:  Method{

 public:
  Topology(FILE * reporter);
  Topology(FILE * resultsFile, FILE * reporter, vector<Soup *> A, vector<Soup *> B);
  ~Topology();

  void generate_topology(Soup * protein, Soup * ligand);
  void generate_topology(vector<Atom*> protein, vector<Atom*> ligand);
  void generate_topology(vector<Atom*> all_atoms);
  float find_surface_distace(Atom * atom);

  vector<float> get_protein_surface_distances();
  vector<float> get_ligand_surface_distances();


 private:

  vector<float> get_coordinates(int index_x, int index_y, int index_z);
  vector<int> get_indexes(float x, float y, float z);
  float search_box(Atom * atom);
  void set_parameters();
  void writeDescription(FILE * reporter);

  vector<Atom *> protein_atoms, ligand_atoms, all_atoms;

  //  vector<float> protein_surface_distances,ligand_surface_distances;
  vector<vector<vector<int> > > grid;

  float xmin,ymin,zmin,xmax,ymax,zmax;
  float spacing,buffer,alpha_factor;
  
  vector<int> no_points;
  vector<int> start;
  vector<int> end;


  int index_buffer;

  FILE * rep;

};

#endif
