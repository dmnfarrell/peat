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


#include "vdw_line.h"

Vdw_Line::Vdw_Line():Method() {
  set_parameters();
} 


Vdw_Line::Vdw_Line(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  
  set_parameters();


}


//Make new constructor so that vdw_line can be started directly without using energy

Vdw_Line::Vdw_Line(FILE * resultsFile, FILE * reporter, vector<Soup *> A, vector<Soup *> B):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();

  /*** read parameters ***/
  if(parameter_type == "AMBER")
    {
      PrepinReader prepA(reporter, A);
      PrepinReader prepB(reporter, B);
      printf("Using AMBER parameters\n");
    }
  else if(parameter_type == "GP")
    {
      Generalparameters prepA(reporter);
      for(unsigned int i=0; i<A.size(); i++)
	prepA.set_vdw_radii(A[i]->getAtoms());
      Generalparameters prepB(reporter);
      for(unsigned int i=0; i<B.size(); i++)
	prepB.set_vdw_radii(B[i]->getAtoms());
      printf("Using GP parameters\n");
    }

  double result = 0;

  if(A.size() != B.size())
    {
      cout << "WARNING: An equal amount of protein and ligand structures must be given to the van der Waals function!\n";
      return;
    }
  
  for(unsigned int i=0; i<A.size(); i++)
    {
      result = calculate(A[i],B[i],0);

      string decoy_number;
      if(B[i]->name.find_last_of("_")== string::npos)
      	decoy_number = "-1";
      else
	decoy_number = B[i]->name.substr(B[i]->name.find_last_of("_")+1);
      
      fprintf(resultsFile,"%-4s %-4s %8f\n",
	      A[i]->name.substr(0,4).c_str(),
	      decoy_number.c_str(),
	      result);
    }  
}



Vdw_Line::Vdw_Line(FILE * resultsFile, FILE * reporter, vector<Soup *> A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();

  /*** read parameters ***/
  if(parameter_type == "AMBER")
    {
      PrepinReader prepA(reporter, A);
      printf("Using AMBER parameters\n");
    }
  else if(parameter_type == "GP")
    {
      Generalparameters prepA(reporter);
      for(unsigned int i=0; i<A.size(); i++)
	prepA.set_vdw_radii(A[i]->getAtoms());
      printf("Using GP parameters\n");
    }

  double result = 0;

  
  for(unsigned int i=0; i<A.size(); i++)
    {
      result = calculate(A[i],0);

      string decoy_number;
      decoy_number = "-1";
            
      fprintf(resultsFile,"%-4s %-4s %8f\n",
	      A[i]->name.substr(0,4).c_str(),
	      decoy_number.c_str(),
	      result);
    }  
}






Vdw_Line::~Vdw_Line(){}


double Vdw_Line::calculate(Soup * A, Soup * B, double cutoff) {
  /*** calculate inter-soup vdw energies ***/
  vector<Atom *> atomsA = A->get_dry_atoms();  
  vector<Atom *> atomsB = B->get_dry_atoms();  
  
  result = 0;
  float temp_result,distance;

  /*** No cutoff - set cutoff to a silly and high value ***/
  if(cutoff == 0)
    cutoff = 100000000;

  /*** calculate bonds ***/
  if(!A->bonds_assigned)
    BM.assign_bonds(A);
  if(!B->bonds_assigned)
    BM.assign_bonds(B);

  for(unsigned int ai=0; ai< atomsA.size(); ai++)
    if(atomsA[ai]->element != "H") // make sure not to include hydrogen atoms
      for(unsigned int bi=0; bi< atomsB.size(); bi++) 
	if(atomsB[bi]->element != "H") // make sure not to include hydrogen atoms
	  {
	    distance = dist.calculate(atomsA[ai],atomsB[bi],true);
	    if(distance < cutoff)
	      if(find(atomsA[ai]->hydrogen_bonds.begin(),
		      atomsA[ai]->hydrogen_bonds.end(), 
		      atomsB[bi]) == atomsA[ai]->hydrogen_bonds.end()) //make sure atoms are not hydrogen bonded
		{
		  temp_result = calculate_vdw_line(atomsA[ai], atomsB[bi], distance);
		  if(include_surf_factor == 1)
		    temp_result = temp_result * (atomsA[ai]->desolvation_factor + atomsB[bi]->desolvation_factor)/2; //using average desolvation_factor factor of the two atoms
	    
		  result += temp_result;
	
	  }
      }
   return result;

}


double Vdw_Line::calculate(Soup * A, double cutoff)
{
  /*** calculate intra-soup vdw energies ***/
  vector<Atom *> atoms = A->get_dry_atoms();  
  //
  // Set up boxes
  //
  Boxes BOX(atoms,8.0);
  vector<int> close_atoms;
  //
  //
  //
  result = 0;
  float temp_result,distance;

  /*** No cutoff - set cutoff to a silly and high value ***/
  if(cutoff == 0)
    cutoff = 100000000;


  /*** calculate bonds ***/
  if(!A->bonds_assigned) {
    BM.assign_bonds(A);
  }

  for(unsigned int ai=0; ai< atoms.size(); ai++) {
    // 
    // Don't calculate for hydrogen atoms
    //
    if(atoms[ai]->element != "H") {
      close_atoms=BOX.get_close_atoms(atoms[ai]);
      for (unsigned int jcount=0;jcount<close_atoms.size();jcount++) {
	int bi=close_atoms[jcount];
	//for(unsigned int bi=0; bi< ai; bi++) {
	// make sure not to include hydrogen atoms
	if(atoms[bi]->element != "H") {
	  distance = dist.calculate(atoms[ai],atoms[bi],true);
	  if(distance < cutoff) {
	    //make sure atoms are not bonded
	    if(find(atoms[ai]->bonded_atoms.begin(),
		    atoms[ai]->bonded_atoms.end(), 
		    atoms[bi]) == atoms[ai]->bonded_atoms.end()) {
	      //
	      // make sure atoms are not hydrogen bonded
	      //
	      if(find(atoms[ai]->hydrogen_bonds.begin(),
		      atoms[ai]->hydrogen_bonds.end(), 
		      atoms[bi]) == atoms[ai]->hydrogen_bonds.end()) {
		temp_result = calculate_vdw_line(atoms[ai], atoms[bi], distance);
		if(include_surf_factor == 1) {
		  temp_result = temp_result * 0.5*(atoms[ai]->desolvation_factor + atoms[bi]->desolvation_factor)/2; //using average desolvation_factor factor of the two atoms
		}
		result += temp_result;
	      }
	    }
	  }
	}
      }
    }
  }
  return result;
}


inline double Vdw_Line::calculate_vdw_line(Atom * atom1, Atom * atom2, double dist) {
  //
  // This is the function used for calculating vdw energies
  //
  R = atom1->vdw_dist + atom2->vdw_dist;
  E = sqrt(atom1->vdw_energy * atom2->vdw_energy);

  if (LJ ==1)  {
    printf ("Doing LJ\n");
    fraction = R/dist;
    fraction2 = fraction*fraction;
    fraction3 = fraction2*fraction;
    fraction6 = fraction3*fraction3;
    fraction12 = fraction6*fraction6;
    res = E*(fraction12-2*fraction6);
    if (res > top) {
      res = top;
    }
    return res;
  }
  //
  // Now doing something else?
  //
  res = 0;
  //printf ("R: %6.3f, close: %6.3f\n",R,R-0.5*width);

  if(dist<R-0.5*width) {
    res = (-E - top)/(R-0.5 * width) * dist + top;
    //res=top*fabs(dist-(R-0.5*width));
    //printf ("Close: %6.3f\n",res);
  }
  else if( R-0.5*width <= dist && dist < R+0.5*width) {
    res = -E;
  }
  else if(R+0.5*width <= dist && dist < R+0.5*width+deltaR ) {
    res = E/deltaR * dist - E/deltaR *(R+0.5*width+deltaR);
  }


  residue1 = atom1->residue_no;
  residue2 = atom2->residue_no;

  if (residue1 == residue2) {
    res = res*same_residue_factor;
  }
  else if( (residue1 - residue2) == 1 || (residue1 - residue2) == -1) {
    res = res*neighbour_residue_factor;
  }
  return res;
}





void Vdw_Line::set_parameters()
{

  ifstream parameters;
  parameters.open("vdw_line.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'vdw_line.cfg' could not be found\n";
      exit(0);
    }


  string dummy;

  while (!parameters.eof())
    {

      parameters >> dummy;

      if(dummy == "MAX_ENERGY:")
	parameters >> top;
      if(dummy == "DELTA_R_FLAT:")
	parameters >> width;
      if(dummy == "DELTA_R_CUTOFF:")
	parameters >> deltaR;
      if(dummy == "PARAMETER_TYPE:")
	parameters >> parameter_type;
      if(dummy == "INCLUDE_SURF_FACTOR:")
	parameters >> include_surf_factor;
      if(dummy == "SAME_RESIDUE_FACTOR:")
	parameters >> same_residue_factor;
      if(dummy == "NEIGHBOUR_RESIDUE_FACTOR:")
	parameters >> neighbour_residue_factor;
      if(dummy == "USE_LJ:")
	parameters >> LJ;

    } //end read file

  //printing all configurations
  /*
  printf("Configurations for vdW Line:\n");
  printf("\n");
  printf("\n");

  printf("MAX_ENERGY:                  %f\n", top);
  printf("\n");

  printf("DELTA_R_FLAT:                %f\n", width);
  printf("\n");

  printf("DELTA_R_CUTOFF:              %f\n", deltaR);
  printf("\n");

  printf("PARAMETER_TYPE:              %s\n", parameter_type.c_str());
  printf("\n");

  printf("INCLUDE_SURF_FACTOR:         %d\n", include_surf_factor);
  printf("\n");

  printf("SAME_RESIDUE_FACTOR:         %f\n", same_residue_factor);
  printf("\n");

  printf("NEIGHBOUR_RESIDUE_FACTOR:    %f\n", neighbour_residue_factor);
  printf("\n");

  printf("Version: $Id: vdw_line.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");
  */
}





void Vdw_Line::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: vdw_line.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("Calculates van der Waals interactions");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());

  rep = reporter;
}
