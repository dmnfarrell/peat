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

/*
$Id: fragmental_volume.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "fragmental_volume.h"


Fragmental_Volume::Fragmental_Volume():Method()
{
  set_parameters();
}


Fragmental_Volume::Fragmental_Volume(FILE * reporter, string mode):Method(reporter)
{
  /*** write out to out file **/
  operation_mode = mode;
  writeDescription(reporter);
  set_parameters();
}

Fragmental_Volume::Fragmental_Volume(FILE * reporter, vector<Soup*> A, string mode):Method(reporter)
{
  /*** write out to out file **/
  operation_mode = mode;
  writeDescription(reporter);
  set_parameters();

  
  

  calculate(A);

  
  //Sampling test
  /*
  points_pr_cubic_A = 1;
  vector<float> acc_vols,tot_vols;
  vector<int> ppc;

  while (points_pr_cubic_A<1000001)
    {
      cout<<"********** points_pr_cubic_A "<<points_pr_cubic_A <<" ***********"<<endl;
      calculate(A);
      
      ppc.push_back(points_pr_cubic_A);
      acc_vols.push_back(total_volume);
      tot_vols.push_back(integrated_total_volume);

      points_pr_cubic_A=points_pr_cubic_A*10;

     }
  
  printf("Points pr cubic A | Accumulated volume |      Total volume  |            Error\n");
  printf("-------------------------------------------------------------------------------\n");
  for(unsigned int i=0;i<acc_vols.size();i++)
    printf("%17d   %18.5f   %17.5f   %17.5f\n",ppc[i],acc_vols[i],tot_vols[i],tot_vols[i]-acc_vols[i]);

  */
 }

Fragmental_Volume::Fragmental_Volume(FILE * reporter, SoupObject* A, string mode):Method(reporter)
 {
   /*** write out to out file **/
   operation_mode = mode;
   writeDescription(reporter);
   set_parameters();
   calculate(A);

 }



Fragmental_Volume::~Fragmental_Volume()
{fclose(out);}



void Fragmental_Volume::calculate(vector<Soup*>  A)
 {

   int a = 0;
   for(unsigned int i=0; i<A.size(); i++)
     {
       if (operation_mode == "ellipse")
	 calculate_protein_fraction_in_ellipse(A[i],focal_point_atoms1[a],focal_point_atoms2[a]);
       else
	 calculate(A[i]);
 
       a++;
     }
}

void Fragmental_Volume::calculate(Soup * A)
 {

   // sybyl atom types
   A->set_sybyl_atom_types();

   vector<vector<SoupObject *> > objects = A->getSoupObjects();

   for(unsigned int i=0; i<objects.size(); i++)
     for(unsigned int j=0; j<objects[i].size(); j++)
       calculate(objects[i][j]);
	
 }


void Fragmental_Volume::set_precalculated_fragmental_volumes(Soup * A)
 {

   // sybyl atom types
   A->set_sybyl_atom_types();

   vector<vector<SoupObject *> > objects = A->getSoupObjects();

   //apply precalculated volumes to protein atoms
   for(unsigned int j=0; j<objects[0].size(); j++)
     apply_precalculated_volumes(objects[0][j]);

   //calculate ligands
   // for(unsigned int i=1; i<objects.size(); i++)
   for(unsigned int j=0; j<objects[1].size(); j++)
     calculate(objects[1][j]);
   
 }


void Fragmental_Volume::apply_precalculated_volumes(SoupObject *A)
{
  vector<Atom*> atoms = A->getAtoms();
  
  vector<Atom*>::iterator atom = atoms.begin();

  while(atom != atoms.end())
    {
      if ((*atom)->element != "H")
	{
	  string residue = (*atom)->residue;
	  if(residue == "HIE" || residue == "HID" )
	    residue = "HIS";
	  if(residue == "CYX")
	    residue = "CYS";
	  
	  string name = (*atom)->name;
	  if (name=="O''")
	    name="OXT";
	  if (name=="O'")
	    name="O";
	  
	  string key = residue+"_"+name;


	  (*atom)->fragmental_volume=precalculated_volumes[key];
	  
	  if ((*atom)->fragmental_volume < 1.0)
	    {
	      (*atom)->fragmental_volume = 10.0;
	      //cout<<"Warning: Fragmental volume set to the standard value:"<< (*atom)->fragmental_volume<<"A^3 for ";
	      //(*atom)->display();
	    }
	}
      atom++;
    }

  return;
}


void Fragmental_Volume::calculate(SoupObject* A)
 {

   total_volume = 0.0;
   accumulated_sphere_volume = 0.0;
   accumulated_group_fragmental_volume = 0.0;

   // get atoms and set radii
   atoms = A->getAtoms();
   atom_radii.clear();
   for(unsigned int i=0; i<atoms.size(); i++)
     atom_radii[atoms[i]] = 
       atomic_radius(atoms[i]->element)*
       atomic_radius(atoms[i]->element);

   // Do bonds
   BM.calculate_list(atoms);
   bonds = BM.getResult();

   // calculate fragmental volumes
   for(unsigned int i=0; i<atoms.size(); i++)
     {
       // make list of spacially close (ra+rb<dist) atoms 
       find_touching_atoms(i);
       
       atoms[i]->fragmental_volume = atom_fragmental_volume(i);
       total_volume += atoms[i]->fragmental_volume;
     }
   
   // find group volumes
   for(unsigned int i=0; i<atoms.size(); i++)
     if(atoms[i]->element != "H")
       {
	 set_bonding_atoms(i);
	 float temp = 0;
	 for(unsigned int b=0; b<bonding_atoms.size(); b++)
	   temp+=bonding_atoms[b]->fragmental_volume;

	 atoms[i]->group_fragmental_volume=temp;
	 accumulated_group_fragmental_volume+=temp;
       }
     else
       atoms[i]->group_fragmental_volume=0.0;

   /*
   
   printf("   Atom  | Group type | Volume   | Sphere Volume | Group Volume\n");
   printf("---------------------------------------------------------------\n");
   for(unsigned int i=0; i<atoms.size(); i++)
     {
       sphere_volume = sphere_volume = 4.0/3.0*3.14159265*pow(atomic_radius(atoms[i]->element),3);
       accumulated_sphere_volume+=sphere_volume;
       set_bonding_atoms(i);

       printf("%3s-%4s | %-5s      | %7.4f  |  %7.4f      | %7.4f\n",
	      atoms[i]->residue.c_str(),
	      atoms[i]->name.c_str(),
	      atoms[i]->sybyl_type.c_str(),
	      //	      atoms[i]->element.c_str(),bonding_atoms.size()-1,
	      atoms[i]->fragmental_volume,
	      sphere_volume,
	      atoms[i]->group_fragmental_volume);
     }

   // volume of entire molcule
   integrated_total_volume = calculate_molecular_volume();

   // print
   cout<<A->identify()<<":"<<endl;
   cout<<"------------------------------------------------"<<endl;
   cout<<"Accumulated fragmental volume:      "<<total_volume<<endl;
   cout<<"Accumulated sphere volume:          "<<accumulated_sphere_volume<<endl;
   cout<<"Integrated total volue:             "<<integrated_total_volume<<endl;
   cout<<"Integration error:                  "<<integrated_total_volume - total_volume<<endl;
   cout<<"Accumulated group fragmental volume "<<accumulated_group_fragmental_volume<<endl<<endl;


   string name = A->identify();
   int loc2 = name.rfind(".");
   int loc1 = name.rfind("/")+1;
   if (loc1 == string::npos)
     loc1=0;
    

   fprintf(out,"# Fragmental volumes for %3s\n",(A->identify()).substr(0,3).c_str());
   fprintf(out,"# Point density is %d points / Aangstroem^3\n\n",points_pr_cubic_A);
   
   for(unsigned int i=0; i<atoms.size(); i++)
     if(atoms[i]->element != "H")
       {
	 set_bonding_atoms(i);

	 fprintf(out,"%3s-%4s  %-5s   %7.4f\n",
		 atoms[i]->residue.c_str(),
		 atoms[i]->name.c_str(),
		 atoms[i]->sybyl_type.c_str(),
		 atoms[i]->group_fragmental_volume);
       }
   */
   return;
}


float Fragmental_Volume::calculate_molecular_volume(Soup * A)
{

   vector<vector<SoupObject *> > objects = A->getSoupObjects();
   accumulated_group_fragmental_volume = 0.0;

   //calculate ligand
   if(objects[1].size()>0)
     calculate(objects[1][0]);

   
   return accumulated_group_fragmental_volume;


}


float Fragmental_Volume::calculate_molecular_volume_by_integration()
{
   // calculate volume of entire molecule

  // box centre
  float xcen=0,ycen=0,zcen=0;
  for(vector<Atom*>::iterator it=atoms.begin(); it != atoms.end(); it++)
    {
      xcen += (*it)->x;
      ycen += (*it)->y;
      zcen += (*it)->z;
    }
  xcen = xcen/atoms.size();
  ycen = ycen/atoms.size();
  zcen = zcen/atoms.size();

  // box size
  float size = 0.0;
  for(vector<Atom*>::iterator it=atoms.begin(); it != atoms.end(); it++)
    {
      float dist = 
	sqrt((xcen-(*it)->x)*(xcen-(*it)->x)+
	     (ycen-(*it)->y)*(ycen-(*it)->y)+
	     (zcen-(*it)->z)*(zcen-(*it)->z))
	+atomic_radius((*it)->element);
      
      if(dist>size)
	size=dist;
    }
  size = size*2.1;

  
  // find out how many points to use
  int N = ((int) (size*size*size*points_pr_cubic_A));

  //cout <<"MV: Size is:   "<<size<<" using "<<N<<" points"<<endl;

  // box boundaries
  float 
    xmin=xcen-size/2,
    ymin=ycen-size/2,
    zmin=zcen-size/2;
  /*
    xmax=xcen+size/2,
    ymax=ycen+size/2,
    ,zmax=zcen+size/2;
  */
  int n=0;

  float xran,yran,zran,dx,dy,dz,dist,ra;;
  float rmax = (float) RAND_MAX;
  vector<Atom*>::iterator it;
  bool inside;

  // Do monte carlo
  for(int j=0; j<N; j++)
    {
      // make a random point
      inside = false;
      xran = xmin + size*((float) rand())/rmax;
      yran = ymin + size*((float) rand())/rmax;
      zran = zmin + size*((float) rand())/rmax;

      // check if point is inside any atom 
      it=atoms.begin();
      while(it!=atoms.end() && !inside)
	{
	  ra = atom_radii[(*it)];
	  dx = xran-(*it)->x;
	  dy = yran-(*it)->y;
	  dz = zran-(*it)->z;
	  
	  dist = dx*dx+dy*dy+dz*dz;

	  if(dist<ra)
	    inside = true;
	  
	  it++;
	}
      
      if (inside)
	n++;
    }
    
  float res = size*size*size*n/N;

  return res;
}



void Fragmental_Volume::find_touching_atoms(int i)
{
  touching_atoms.clear();
  for(unsigned int j=0; j<atoms.size(); j++)
    {
      float dist = sqrt(
			(atoms[i]->x-atoms[j]->x)*(atoms[i]->x-atoms[j]->x)+
			(atoms[i]->y-atoms[j]->y)*(atoms[i]->y-atoms[j]->y)+
			(atoms[i]->z-atoms[j]->z)*(atoms[i]->z-atoms[j]->z));



      if (atoms[i] != atoms[j])
	if (dist<(atomic_radius(atoms[i]->element)+atomic_radius(atoms[j]->element)))
	  touching_atoms.push_back(atoms[j]);
      
    }
}


float Fragmental_Volume::atom_fragmental_volume(int i)
{

  // box centre
  float xcen=atoms[i]->x,ycen=atoms[i]->y,zcen=atoms[i]->z;

  // box size
  float size = sqrt(atom_radii[atoms[i]])*2.1;
  
  // find out how many points to use
  int N = ((int) (size*size*size*points_pr_cubic_A));

  //  cout <<" Size is:   "<<size<<" using "<<N<<" points"<<endl;

  // box boundaries
  float 
    xmin=xcen-size/2,
    ymin=ycen-size/2,
    zmin=zcen-size/2;
  /*,zmax=zcen+size/2
     xmax=xcen+size/2,
     ymax=ycen+size/2;
  */
  int n=0;
  float xran,yran,zran;
  float rmax = (float) RAND_MAX;

  // Do monte carlo
  for(int j=0; j<N; j++)
    {
      xran = xmin + size*((float) rand())/rmax;
      yran = ymin + size*((float) rand())/rmax;
      zran = zmin + size*((float) rand())/rmax;
      
      if (inside(xran,yran,zran,i))
	  n+=1;
      
     }
    
  // print result
  float res = size*size*size*n/N;
  
  return res;
}

inline bool Fragmental_Volume::inside(float x, float y, float z, int i)
{
  bool res = false;
  
  float dx,dy,dz,dist,ra,rb,da,D,costheta;

  // check if point is inside the atom 
  ra = atom_radii[atoms[i]];
  dx = x-atoms[i]->x;
  dy = y-atoms[i]->y;
  dz = z-atoms[i]->z;
  
  dist = dx*dx+dy*dy+dz*dz;

  if(dist<ra)
      res= true;

   // check that the point is not closer to any bonded heavy atom
   //
   // Using plane of "circle intersections"
   //

   if (res)
     {
       vector<float> norm,p,q;
       vector<Atom*>::iterator it = touching_atoms.begin();
       
       while (it != touching_atoms.end() && res)
	 {
	   norm.clear();
	   p.clear();
	   q.clear();

	   // find plane normal vector
	   norm.push_back((*it)->x - atoms[i]->x);
	   norm.push_back((*it)->y - atoms[i]->y);
	   norm.push_back((*it)->z - atoms[i]->z);

	   // lenght of normal vector
	   D = sqrt(norm[0]*norm[0]+
		    norm[1]*norm[1]+
		    norm[2]*norm[2]);


	   // normalise normal vector
	   for(unsigned int c=0;c<3;c++)
	     norm[c]=norm[c]/D;

	   // find distance from center of current atom to plane
	   rb = atom_radii[(*it)];
	   da = 0.5 * (D + (ra-rb)/D);

	   // find point on plane
	   p.push_back(atoms[i]->x + da*norm[0]);
	   p.push_back(atoms[i]->y + da*norm[1]);
	   p.push_back(atoms[i]->z + da*norm[2]);

	   // find vector from p to the random point
	   q.push_back(x - p[0]);
	   q.push_back(y - p[1]);
	   q.push_back(z - p[2]);

	   // find cos(angle btw p and q)
	   costheta = (norm[0]*q[0] + norm[1]*q[1] + norm[2]*q[2])/(D*sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]));

	   if(costheta > 0)
	       res= false;

	   it++;
	 }

     }

   return res;

 }


void Fragmental_Volume::set_bonding_atoms(int i)
 {
   // clear list
   bonding_atoms.clear();

   // insert the heavy atom 
   bonding_atoms.push_back(atoms[i]);

   // find bonding atoms
   for(unsigned int j=0; j<bonds.size(); j++)
     if (bonds[j][0]==i)
       if (atoms[bonds[j][1]]->element == "H")
	   bonding_atoms.push_back(atoms[bonds[j][1]]);

   return;
 }




void Fragmental_Volume::calculate_protein_fraction_in_ellipse(Soup * A, int a, int b)
{
  /*

    Estimates the fraction of space in an ellipse that is occucpied by protein
    Atoms with indices a and b are located in the two focal points of the ellipse

    Ellipse radius: sum of distance from border to the two focal points

  */

  // get atoms from soup object
  atoms = A->getAtoms();
  cout<<"Protein fraction in ellipse on "<<A->name<<" - focal point atoms "<<atoms[a]->name<<" and "<<atoms[b]->name<<endl;

  // ellipse setup
  float distance = sqrt((atoms[a]->x - atoms[b]->x)*(atoms[a]->x - atoms[b]->x)+
			(atoms[a]->y - atoms[b]->y)*(atoms[a]->y - atoms[b]->y)+
			(atoms[a]->z - atoms[b]->z)*(atoms[a]->z - atoms[b]->z));
 
  float ellipse_radius = ellipse_radius_factor*distance;
  float ellipse_volume = 1.0/6.0 * 3.14159265 * ellipse_radius * (ellipse_radius * ellipse_radius - distance * distance);
  
  // set radii
  atom_radii.clear();
  for(unsigned int i=0; i<atoms.size(); i++)
    atom_radii[atoms[i]] = 
      atomic_radius(atoms[i]->element)*
      atomic_radius(atoms[i]->element);

  // box centre
  float 
    xcen=(atoms[a]->x + atoms[b]->x)/2.0,
    ycen=(atoms[a]->y + atoms[b]->y)/2.0,
    zcen=(atoms[a]->z + atoms[b]->z)/2.0;

  // box size
  float size = ellipse_radius * 1.1;
 
  // find out how many points to use
  int N = ((int) (size*size*size*points_pr_cubic_A));

  // box boundaries
  float 
    xmin=xcen-size/2,
    ymin=ycen-size/2,
    zmin=zcen-size/2;
  /*
    xmax=xcen+size/2,
    ymax=ycen+size/2,
    zmax=zcen+size/2;
  */

  // make a list of atoms in and around the ellipse for speed optimization
  vector<Atom * >::iterator it;
  vector<Atom*> close_atoms;
  it=atoms.begin();

  int n=0,n_ellipse=0;
  float xran,yran,zran,dist_a,dist_b,ra,dx,dy,dz,dist;
  float rmax = (float) RAND_MAX;
  bool inside = false;
  
  while(it!=atoms.end())
    {
      
      dist_a = sqrt((atoms[a]->x - (*it)->x)*(atoms[a]->x - (*it)->x)+
		    (atoms[a]->y - (*it)->y)*(atoms[a]->y - (*it)->y)+
		    (atoms[a]->z - (*it)->z)*(atoms[a]->z - (*it)->z));

      dist_b = sqrt((atoms[b]->x - (*it)->x)*(atoms[b]->x - (*it)->x)+
		    (atoms[b]->y - (*it)->y)*(atoms[b]->y - (*it)->y)+
		    (atoms[b]->z - (*it)->z)*(atoms[b]->z - (*it)->z));
    
      if (dist_a+dist_b<=ellipse_radius+5.0) //buffer of 5.0 Aangstroem
	close_atoms.push_back((*it));

      it++;
    }


  // Do monte carlo
  for(int j=0; j<N; j++)
    {
      inside = false;
      xran = xmin + size*((float) rand())/rmax;
      yran = ymin + size*((float) rand())/rmax;
      zran = zmin + size*((float) rand())/rmax;
        
      // check if point is inside ellipse
      dist_a = sqrt((atoms[a]->x - xran)*(atoms[a]->x - xran)+
		    (atoms[a]->y - yran)*(atoms[a]->y - yran)+
		    (atoms[a]->z - zran)*(atoms[a]->z - zran));

      dist_b = sqrt((atoms[b]->x - xran)*(atoms[b]->x - xran)+
		    (atoms[b]->y - yran)*(atoms[b]->y - yran)+
		    (atoms[b]->z - zran)*(atoms[b]->z - zran));
    

      // check if point is inside any atom 
      if (dist_a+dist_b<=ellipse_radius)
	{
	  n_ellipse++;
	  it=close_atoms.begin();
	  while(it!=close_atoms.end() && !inside)
	    {
	      ra = atom_radii[(*it)];
	      dx = xran-(*it)->x;
	      dy = yran-(*it)->y;
	      dz = zran-(*it)->z;
	      
	      dist = dx*dx+dy*dy+dz*dz;
	      
	      if(dist<ra)
		inside = true;
	      
	      it++;
	    }
	}

      if (inside)
	  n++;
     }
    
  // print result
  float protein_volume = size*size*size*n/N;
  float ellipse_volume_integrated = size*size*size*n_ellipse/N;
  float fraction_occupied_by_protein = protein_volume/ellipse_volume;

  cout<<"Integrated eliipse volume is: "<<ellipse_volume_integrated<<endl;
  cout<<"Protein volume inside ellipse is: "<<protein_volume<<endl;
  cout<<"Fraction of the volume inside the ellipse occupied by protein: "<<fraction_occupied_by_protein<<endl;
    
  // make file 
  FILE * out;

  string filename;
  stringstream stream;
  stream<<A->name<<"_"<<a<<"_"<<b<<".txt";
  filename = stream.str();

  out = fopen(filename.c_str(),"w");

  fprintf(out,"# Ellipse fragmental volumes for %3s\n",A->name.c_str());
  fprintf(out,"# Point density is %d points / Aangstroem^3\n\n",points_pr_cubic_A);

  fprintf(out,"Focal atom one: %3s-%4s  (%-5d) coordinate  %7.4f, %7.4f, %7.4f\n",
	  atoms[a]->residue.c_str(),
	  atoms[a]->name.c_str(),
	  atoms[a]->residue_no,
	  atoms[a]->x,atoms[a]->y,atoms[a]->z);

  fprintf(out,"Focal atom two: %3s-%4s  (%-5d) coordinate  %7.4f, %7.4f, %7.4f\n\n",
	  atoms[b]->residue.c_str(),
	  atoms[b]->name.c_str(),
	  atoms[b]->residue_no,
	  atoms[b]->x,atoms[b]->y,atoms[b]->z);
   
  fprintf(out,"Focal point distance: %7.4f\n",distance);
  fprintf(out,"Ellipse radius: %7.4f\n",ellipse_radius);
  fprintf(out,"Ellipse volume: %7.4f\n",ellipse_volume);
  fprintf(out,"Integrated ellipse volume is: %7.4f\n",ellipse_volume_integrated);
  fprintf(out,"Protein volume inside ellipse is: %7.4f\n",protein_volume);
  fprintf(out,"Fraction of the volume inside the ellipse occupied by protein: %7.4f\n",fraction_occupied_by_protein);

  fclose(out);
  printf("Done making file\n");

  return;
}




void Fragmental_Volume::set_mode(string mode){operation_mode = mode;}


float Fragmental_Volume::atomic_radius(string element)
{
  float res = 2.0;
  if(radii.count(element)>0)
    res = radii[element];

  return res;
}

void Fragmental_Volume::set_parameters()
{
  
  // set atomic radii
  radii["H"] = 1.08;
  radii["C"] = 1.53;
  radii["O"] = 1.36;
  radii["N"] = 1.45;
  radii["F"] = 1.30;
  radii["CL"] = 1.65;
  radii["BR"] = 1.80;
  radii["I"] = 2.05;
  radii["S"] = 1.70;
  radii["P"] = 1.75;
  
  points_pr_cubic_A = 1000; //10000
  
  //
  // get atomic parameters for ellipse mode
  //
  //   cout<<"operation mode: "<<operation_mode<<endl;

   if (operation_mode == "ellipse")
     {
       
       ifstream parameters;

       parameters.open("ellipse.cfg");
       
       if(! parameters.is_open())
	 {
	   cout << "Error:File 'ellipse.cfg' could not be found\n";
	   exit(0);
	 }
       
       //printf("Setting parameters\n");
  
       string dummy;
       
       while (!parameters.eof())
	 {
	   
	   parameters >> dummy;
	   
	   if(dummy == "FOCAL_POINT_ATOM1:")
	     {
	       parameters >> dummy;
	       while(dummy != "end")
		 {
		   focal_point_atoms1.push_back(atoi(dummy.c_str()));
		   parameters >> dummy;
		 }
	     }
	   
	   if(dummy == "FOCAL_POINT_ATOM2:")
	     {
	       parameters >> dummy;
	       while(dummy != "end")
		 {
		   focal_point_atoms2.push_back(atoi(dummy.c_str()));
		   parameters >> dummy;
		 }
	     }
	   
	   
	   
	   if(dummy == "ELLIPSE_RADIUS_FACTOR:")
	     parameters >> ellipse_radius_factor;
	   
	 } //end read file

       //printing all configurations

       printf("Configurations for ellipse:\n");
       printf("\n");
       printf("\n");
       
       printf("FOCAL_POINT_ATOM1:        ");
       for(unsigned int i=0; i<focal_point_atoms1.size(); i++)
	 printf("%d ",focal_point_atoms1[i]);
       printf("end\n");
       
       printf("FOCAL_POINT_ATOM2:        ");
       for(unsigned int i=0; i<focal_point_atoms2.size(); i++)
	 printf("%d ",focal_point_atoms2[i]);
       printf("end\n");
       
       printf("ELLIPSE_RADIUS_FACTOR:    %f\n", ellipse_radius_factor);
       printf("\n");
       
       printf("Version: $Id: fragmental_volume.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
       printf("\n");

     
     }

   // read in precalculated volumes
   ifstream in;
   char line[256];
   string s;
   in.open("parameters/precalculated_volumes.txt");

   if(! in.is_open())
     {
       cout << "Error:File 'parameters/precalculated_volumes.txt' could not be found\n";
       exit(0);
     }

   while(!in.eof())
     {
       in.getline(line,256);
       //cout<<line<<endl;
       s = string(line);
      
       if (s.size()>3 && s.substr(0,1) != "#")
	 {
	   string key = s.substr(0,8);
	   for(unsigned int i=0; i<key.size(); )
	     if(string(key,i,1) ==" ")
	       key = key.erase(i, 1);
	     else
	       i++;

	   float volume = atof(s.substr(9,21).c_str());



	   precalculated_volumes[key]= volume;
	  
	   //cout<<"Atom "<<key<<" has the volume "<< precalculated_volumes[key]<<endl;
	}
      

   }






   // open file for printing atomic fragmental volumes
   out = fopen("atomic_fragmental_volumes.txt","w");

 }

 void Fragmental_Volume::writeDescription(FILE * reporter)
 {
   /*** write out to out file **/
   version =string("$Id: fragmental_volume.cpp 6326 2010-08-26 16:08:21Z nielsen $");

   description = 
      string("");

   fprintf(reporter,"%s\n",version.c_str());	
   fprintf(reporter,"%s\n",description.c_str());
 }
