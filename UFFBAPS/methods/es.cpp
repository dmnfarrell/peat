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
$Id: es.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "es.h"

Es::Es():Method()
{
  set_parameters();
}

Es::Es(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  set_parameters();
}

//Extra constructor so that es can be started directly without using energy

Es::Es(FILE * resultsFile, FILE * reporter, vector<Soup *> A, vector<Soup *> B):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  set_parameters();

  /*** read parameters ***/
  if(parameter_type == "AMBER")
    {
      PrepinReader prepA(reporter, A);
      PrepinReader prepB(reporter, B);
    }
  else if(parameter_type == "GP")
    {
      Generalparameters prepA(reporter, A);
      Generalparameters prepB(reporter, B);
    }
  else if(parameter_type == "SIMPLE")
    {
      Simple_Parameters prepA(reporter, A);
      Simple_Parameters prepB(reporter, B);
    }



  double result = 0;

  if(A.size() != B.size())
    {
      cout << "WARNING: An equal amount of protein and ligand structures must be given to the ES function!\n";
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


Es::Es(FILE * resultsFile, FILE * reporter, vector<Soup *> A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  set_parameters();

  /*** read parameters ***/
  if(parameter_type == "AMBER")
    {
      PrepinReader prepA(reporter, A);
    }
  else if(parameter_type == "GP")
    {
      Generalparameters prepA(reporter, A);
    }
  else if(parameter_type == "SIMPLE")
    {
      Simple_Parameters prepA(reporter, A);
    }



  double result = 0;

  
  for(unsigned int i=0; i<A.size(); i++)
    {

      result = calculate(A[i],0);

      fprintf(resultsFile,"%-4s %8f\n",
	      A[i]->name.substr(0,4).c_str(),
	      result);

    }  

}



Es::~Es(){}


double Es::calculate(Soup * A, Soup * B, double cutoff){

  vector<Atom *> atomsA = A->get_dry_atoms();  
  vector<Atom *> atomsB = B->get_dry_atoms();  

  result = 0;

  /*** No cutoff - set cutoff to a silly and high value ***/
  if(cutoff == 0)
    cutoff = 100000000;

  /*** calculate bonds ***/
  if(!A->bonds_assigned)
    BM.assign_bonds(A);
  if(!B->bonds_assigned)
    BM.assign_bonds(B);

  float distance;

  for(unsigned int ai=0; ai< atomsA.size(); ai++)
    if(atomsA[ai]->charge != 0)
      for(unsigned int bi=0; bi< atomsB.size(); bi++) 
	if(atomsB[bi]->charge != 0)
	  {
	    distance = dist.calculate(atomsA[ai],atomsB[bi],true);
	    if(distance < cutoff)
	      if(find(atomsA[ai]->hydrogen_bonds.begin(),atomsA[ai]->hydrogen_bonds.end(), atomsB[bi]) == atomsA[ai]->hydrogen_bonds.end()) //make sure atoms are not hydrogen bonded
		{

		  result += calculate_es(atomsA[ai], atomsB[bi], distance);
		  //		  if(include_surf_factor == 1)
		  //		    result = result * (atomsA[ai]->desolvation_factor + atomsB[bi]->desolvation_factor)/2;
		    		  
		}



      }
   return result;

}


double Es::calculate(Soup * A, double cutoff)
{
  vector<Atom *> atoms = A->get_dry_atoms();  
  result = 0;

  /*** No cutoff - set cutoff to a silly and high value ***/
  if(cutoff == 0)
    cutoff = 100000000;

  /*** calculate bonds ***/
  if(!A->bonds_assigned)
    BM.assign_bonds(A);

  float distance;

  for(unsigned int ai=0; ai< atoms.size(); ai++)
    if(atoms[ai]->charge != 0)
      for(unsigned int bi=0; bi< ai; bi++)
	if(atoms[bi]->charge != 0)
	  {
	    distance = dist.calculate(atoms[ai],atoms[bi],true);
	    if(distance < cutoff)
	      if(find(atoms[ai]->bonded_atoms.begin(),atoms[ai]->bonded_atoms.end(), atoms[bi]) == atoms[ai]->bonded_atoms.end()) //make sure atoms are not bonded
		if(find(atoms[ai]->hydrogen_bonds.begin(),atoms[ai]->hydrogen_bonds.end(), atoms[bi]) == atoms[ai]->hydrogen_bonds.end()) //make sure atoms are not hydrogen bonded
		  {
		    result += calculate_es(atoms[ai], atoms[bi], distance);
		    
		    
		  }
	  }

   return result;

}



double Es::calculate_es(Atom * atom1, Atom * atom2, double dist)
{
  /// 2.3066178170179668e-28 m J
  /// 2.3066178170179666e-18 A J
  /// 1389077.5420576576 A J/mol
  /// 1389.0775420576576 A kJ/mol


  D = 1389.0775; // e*e/(4*pi*eps_0) [A kJ/mol]
  res =  D * atom1->charge * atom2->charge / (permittivity(atom1, atom2, dist) *dist);

  residue1 = atom1->residue_no;
  residue2 = atom2->residue_no;

  if (residue1 == residue2)
    res = res*same_residue_factor;
  else if( (residue1 - residue2) == 1 || (residue1 - residue2) == -1)
    res = res*neighbour_residue_factor;

  if(include_surf_factor == 1)
    res = res * (atom1->desolvation_factor + atom2->desolvation_factor)/2;

  if(include_debye_huckel == 1)
    res = res * exp(-0.073315 * dist); //5.66*sqrt(I/T) I = 0.05, T = 298 FOLDX article

  //  cout<<"Found "<<res<<" between\n       ";
  //  atom1->display();
  //  cout<<"   and ";
  //  atom2->display();
    

  return res;
}


double Es::permittivity(Atom * atom1, Atom * atom2, double r)
{

  double eps;
  
  if(permittivity_type == "CONST")
    eps = constant;
  else if(permittivity_type == "CONST_DIST")
    eps = constant*r;
  else if(permittivity_type == "DESOLVATION")
    eps = 80 - 72 *(atom1->desolvation_factor + atom2->desolvation_factor)/2;
  else
    {
      eps = 0;
      printf("WARNING: Unknown permittivity type, %s\n",permittivity_type.c_str());
    }
    

  //  cout<<"eps "<<eps<<endl;
  return eps;
}


void Es::set_parameters()
{

  ifstream parameters;
  parameters.open("es.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'es.cfg' could not be found\n";
      exit(0);
    }


  string dummy;

  while (!parameters.eof())
    {

      parameters >> dummy;

      if(dummy == "PERMITTIVITY_TYPE:")
	parameters >> permittivity_type;
      if(dummy == "PARAMETER_TYPE:")
	parameters >> parameter_type;
      if(dummy == "CONSTANT:")
	parameters >> constant;
      if(dummy == "SAME_RESIDUE_FACTOR:")
	parameters >> same_residue_factor;
      if(dummy == "NEIGHBOUR_RESIDUE_FACTOR:")
	parameters >> neighbour_residue_factor;
      if(dummy == "INCLUDE_SURF_FACTOR:")
	parameters >> include_surf_factor;
      if(dummy == "INCLUDE_DEBYE_HUCKEL:")
	parameters >> include_debye_huckel;


    } //end read file

  //printing all configurations
  /*
  printf("Configurations for es:\n");
  printf("\n");
  printf("\n");

  printf("PERMITTIVITY_TYPE:             %s\n", permittivity_type.c_str());
  printf("\n");

  printf("PARAMETER_TYPE:                %s\n", parameter_type.c_str());
  printf("\n");

  printf("CONSTANT:                      %f\n", constant);
  printf("\n");

  printf("SAME_RESIDUE_FACTOR:           %f\n", same_residue_factor);
  printf("\n");

  printf("NEIGHBOUR_RESIDUE_FACTOR:      %f\n", neighbour_residue_factor);
  printf("\n");

  printf("INCLUDE_SURF_FACTOR:           %d\n", include_surf_factor);
  printf("\n");

  printf("INCLUDE_DEBYE_HUCKEL:          %d\n", include_debye_huckel);
  printf("\n");

  printf("Version: $Id: es.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");
  */
}



void Es::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: es.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("Calculates electrostatic interactions");

  rep =reporter;
  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
