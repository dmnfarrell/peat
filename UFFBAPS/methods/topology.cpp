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

#include "topology.h"


Topology::Topology(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  set_parameters();
}

Topology::Topology(FILE * resultsFile, FILE * reporter, vector<Soup *> A, vector<Soup *> B):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  set_parameters();
  
  generate_topology(A[0],B[0]);

}


Topology::~Topology(){}


void Topology::generate_topology(Soup * protein, Soup * ligand)
{

  vector<Atom*> p_atoms = protein->getAtoms();
  vector<Atom*> l_atoms = ligand->getAtoms();
  
  generate_topology(p_atoms, l_atoms);

}

void Topology::generate_topology(vector<Atom*> protein, vector<Atom*> ligand)
{
  protein_atoms = protein;
  ligand_atoms = ligand;
  
  all_atoms = protein;
  for(unsigned int i=0; i<ligand.size();i++)
    all_atoms.push_back(ligand[i]);
  
  generate_topology(all_atoms);

}

void Topology::generate_topology(vector<Atom*> all_atoms)
{
  printf("top for %d atoms\n",all_atoms.size());

  /*** Find coordinate box of protein-ligand complex ***/
  xmin = xmax = all_atoms[0]->x;
  ymin = ymax = all_atoms[0]->y;
  zmin = zmax = all_atoms[0]->z;

  for(unsigned int i=1; i<all_atoms.size();i++)
    {
      if(all_atoms[i]->x < xmin)
	xmin = all_atoms[i]->x;
      if(all_atoms[i]->y < ymin)
	ymin = all_atoms[i]->y;
      if(all_atoms[i]->z < zmin)
	zmin = all_atoms[i]->z;

      if(all_atoms[i]->x > xmax)
	xmax = all_atoms[i]->x;
      if(all_atoms[i]->y > ymax)
	ymax = all_atoms[i]->y;
      if(all_atoms[i]->z > zmax)
	zmax = all_atoms[i]->z;

    }
  
  printf("Complex box: [%f:%f][%f:%f][%f:%f]\n",
	 xmin,xmax,ymin,ymax,zmin,zmax);

  printf("Distances: (%f,%f,%f)\n",xmax-xmin,ymax-ymin,zmax-zmin);
	 


  /*** add buffer to box limits ***/
  xmin = xmin - buffer;
  ymin = ymin - buffer;
  zmin = zmin - buffer;
  
  xmax = xmax + buffer;
  ymax = ymax + buffer;
  zmax = zmax + buffer;

  printf("Buffered complex box: [%f:%f][%f:%f][%f:%f]\n",
	 xmin,xmax,ymin,ymax,zmin,zmax);

  /*** find number of grid points in each dimension ***/
  no_points[0] = (int) ceil((xmax-xmin)/spacing);
  no_points[1] = (int) ceil((ymax-ymin)/spacing);
  no_points[2] = (int) ceil((zmax-zmin)/spacing);
  
  printf("number of points (%d,%d,%d)\n",no_points[0],no_points[1],no_points[2]);

  /*** generate grid matrix ***/
  vector<int> temp_z(no_points[2],0);
  vector<vector<int> > temp_y(no_points[1],temp_z);
  grid.assign(no_points[0],temp_y);


  /*** set vdw radii ***/
  Generalparameters GP(rep);
  GP.set_vdw_radii(all_atoms);

  /*** set up topology ***/
  float threshold,sqd,sqThreshold;
  vector<float> coor;
  // bool go;

 
  
  printf("Setting up topology\n");

  vector<int> atom_indexes;
  vector<int> top_start(3,0);
  vector<int> top_end(3,0);
  
  int no_points_search = (int) ( (3.0*alpha_factor/spacing)+1 );
  
  for(unsigned int a=0; a<all_atoms.size(); a++)
    {  
      
      atom_indexes = get_indexes(all_atoms[a]->x,all_atoms[a]->y,all_atoms[a]->z);
      threshold = all_atoms[a]->vdw_dist*alpha_factor;
      sqThreshold = threshold*threshold;
      
      for(unsigned int c=0;c<3;c++)
	{
	  top_start[c] = atom_indexes[c] - no_points_search;
	  top_end[c]   = atom_indexes[c] + no_points_search;
	  
	  if(top_start[c] < 0)
	    top_start[c] = 0;
	  
	  if(top_end[c] > no_points[c])
	    top_end[c] = no_points[c];
	}
      
      for(int x=top_start[0]; x<top_end[0]; x++)
	for(int y=top_start[1]; y<top_end[1]; y++)
	  for(int z=top_start[2]; z<top_end[2]; z++)
	    {
	      coor = get_coordinates(x,y,z);
	      sqd = 
		(all_atoms[a]->x - coor[0])*(all_atoms[a]->x - coor[0])+
		(all_atoms[a]->y - coor[1])*(all_atoms[a]->y - coor[1])+
		(all_atoms[a]->z - coor[2])*(all_atoms[a]->z - coor[2]);
	      
	      if(sqd < sqThreshold)
		grid[x][y][z] = 1;
	    }
    }



  /*** old unoptimised algorithm ***/  
//   int a;
//       for(unsigned int x=0; x<no_points[0]; x++)
// 	for(unsigned int y=0; y<no_points[1]; y++)
// 	  for(unsigned int z=0; z<no_points[2]; z++)
// 	{
// 	  go = true;
// 	  a = 0;
//   	  coor = get_coordinates(x,y,z);
	  
// 	  while(go)
// 	    {
// 	      threshold = all_atoms[a]->vdw_dist*alpha_factor;
// 	      sqThreshold = threshold*threshold;
      

// 	      if(all_atoms[a]->x - coor[0] < threshold)
// 		if(all_atoms[a]->y - coor[1] < threshold)
// 		  if(all_atoms[a]->z - coor[2] < threshold)
// 		    {
// 		      sqd = 
// 			(all_atoms[a]->x - coor[0])*(all_atoms[a]->x - coor[0])+
// 			(all_atoms[a]->y - coor[1])*(all_atoms[a]->y - coor[1])+
// 			(all_atoms[a]->z - coor[2])*(all_atoms[a]->z - coor[2]);
		      
// 		      if(sqd < sqThreshold)
// 			{
// 			  grid[x][y][z] = 1;
// 			  go = false;
// 			}

// 		    }
// 	      a++;
// 	      if(a == all_atoms.size())
// 		go =false;

//   	    }
// 	}

//     
  printf("Setting up topology - done\n");

 //  for(unsigned int z=0; z<no_points[2]; z++)
//     {
//       printf("Cut at z=%d\n",z);
//       for(unsigned int x=0; x<no_points[0]; x++)
// 	{
// 	  for(unsigned int y=0; y<no_points[1]; y++)
// 	    {
// 	      if(grid[x][y][z] == 1)
// 		printf("# ");
// 	      else
// 		printf("| ");
// 	    }
// 	  printf("\n");
// 	}
//       printf("\n");
//     }


}


float Topology::find_surface_distace(Atom * atom)
{
  /*** find surface distances ***/

  vector<int> indexes = get_indexes(atom->x, atom->y, atom->z);

  int ib = index_buffer;

  float res = -1.0;

  while(res < -0.001)
    {
      for(unsigned int i=0;i<3;i++)
	{
	  start[i] = indexes[i]-ib;
	  end[i]   = indexes[i]+ib;
	}
      
      res = search_box(atom);
      
      //      printf("box size %d gave %f\n",ib,res);

      ib += 10;

      if(ib>50)
	{
	  printf("WARNING: Uanble to determine surface distance of atom %s\n",atom->element.c_str());
	  res = 0.0;
	}
    }

 
  return res;

}


vector<float> Topology::get_protein_surface_distances()
{
  vector<float> res;
  for(unsigned int i=0; i<protein_atoms.size(); i++)
    res.push_back(find_surface_distace(protein_atoms[i]));

  return res;

}
vector<float> Topology::get_ligand_surface_distances()
{
  vector<float> res;
  for(unsigned int i=0; i<ligand_atoms.size(); i++)
    res.push_back(find_surface_distace(ligand_atoms[i]));

  return res;


}



inline float Topology::search_box(Atom * atom)
{
  float min_dist,cur_dist;
  vector<float> coor;

  for(unsigned int i=0;i<3;i++)
    {
      if(start[i] < 0)
	start[i] = 0;
      
      if(end[i] > no_points[i])
	end[i] = no_points[i];
    }



  min_dist = 100000;
  
  for(int x=start[0]; x<end[0]; x++)
    for(int y=start[1]; y<end[1]; y++)
      for(int z=start[2]; z<end[2]; z++)
	{
	  if(grid[x][y][z] == 0)
	    {
	      coor = get_coordinates(x,y,z);

	      if(atom->x - coor[0] < min_dist)
		if(atom->y - coor[1] < min_dist)
		  if(atom->z - coor[2] < min_dist)
		    {
		      cur_dist = sqrt( 
				      (atom->x - coor[0])*(atom->x - coor[0])+
				      (atom->y - coor[1])*(atom->y - coor[1])+
				      (atom->z - coor[2])*(atom->z - coor[2]));
		      
		      if(cur_dist < min_dist)
			{
			  min_dist = cur_dist;
		
			}
		      
		    }
	    }
	}

  if(min_dist>1000)
    min_dist = -1.0;

  return min_dist;

}




inline vector<float> Topology::get_coordinates(int index_x, int index_y, int index_z)
{
  vector<float> res(3,0);
  res[0] = xmin + index_x*spacing;
  res[1] = ymin + index_y*spacing;
  res[2] = zmin + index_z*spacing;

  return res;
}


inline vector<int> Topology::get_indexes(float x, float y, float z)
{
  vector<int> res(3,0);
  res[0] = (int) floor((x - xmin)/spacing);
  res[1] = (int) floor((y - ymin)/spacing);
  res[2] = (int) floor((z - zmin)/spacing);

  //  printf("Coordinates (%f,%f,%f) translates into (%d,%d,%d)\n",x,y,z,res[0],res[1],res[2]);


  return res;
}



void Topology::set_parameters()
{

  spacing = 1.0;
  buffer = 5.0;
  alpha_factor = 2.0;

  index_buffer = 10;

  start.assign(3,0);
  end.assign(3,0);
  no_points.assign(3,0);
}

void Topology::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: topology.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());

  rep = reporter;
}

