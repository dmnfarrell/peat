
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
$Id: ligand_entropy.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "ligand_entropy.h"


LigandEntropy::LigandEntropy():Method()
{
  /*** set parameters ***/
  set_parameters();

}

LigandEntropy::LigandEntropy(FILE * reporter):Method(reporter)
{
  /*** set parameters ***/
  set_parameters();

 
  /*** write out to out file **/
  writeDescription(reporter);


  
}

LigandEntropy::LigandEntropy(FILE * resultsFile, FILE * reporter, vector<Soup*> A):Method(reporter)
{
  /*** set parameters ***/
  set_parameters();


  /*** write out to out file **/
  writeDescription(reporter);

  for(unsigned int i=0; i<A.size(); i++)
    {
      float result = calculate(A[i]);
      printf("result:'%f'\n",result);

      string decoy_number;
      if(A[i]->name.find_last_of("_")== string::npos)
      	decoy_number = "-1";
      else
	decoy_number = A[i]->name.substr(A[i]->name.find_last_of("_")+1);
      
      fprintf(resultsFile,"%-4s %-4s %8f\n",
	      A[i]->name.substr(0,4).c_str(),
	      decoy_number.c_str(),
	      result);


    }  


}

LigandEntropy::LigandEntropy(FILE * resultsFile, FILE * reporter, vector<Soup*> A, vector<Soup*> B):Method(reporter)
{
  /*** set parameters ***/
  set_parameters();


  /*** write out to out file **/
  writeDescription(reporter);

  for(unsigned int i=0; i<A.size(); i++)
    {
      float result = calculate_using_desolvation(A[i],B[i]);
      printf("result:'%f'\n",result);

      string decoy_number;
      if(A[i]->name.find_last_of("_")== string::npos)
      	decoy_number = "-1";
      else
	decoy_number = A[i]->name.substr(A[i]->name.find_last_of("_")+1);
      
      fprintf(resultsFile,"%-4s %-4s %8f\n",
	      A[i]->name.substr(0,4).c_str(),
	      decoy_number.c_str(),
	      result);


    }  


}



LigandEntropy::~LigandEntropy(){}

float LigandEntropy::calculate(Soup * A)
{
  return calculate(A->get_dry_atoms());
}

float LigandEntropy::calculate_using_desolvation(Soup * protein, Soup * ligand)
{

  float total_desolvation = DE.calculate_contact_between_A_and_B(ligand, protein); // calculate desolvation of A when also taking B into account
  vector<Atom*> ligand_atoms = ligand->getAtoms();
  float number_atoms = (float) ligand_atoms.size();

  float avr_desolvation = total_desolvation/number_atoms;


  // release of water atoms
  // molar concentration [mol/L] * Avogadro's constant [1/mol] * 1000 L/M^3 * 1e30 A^3/m^3
  float water_number_density = 0.033367; //55.4086*6.022e23*1000*1e-30;



  float ligand_volume = FV.calculate_molecular_volume(ligand);
  float displaced_waters = floor(ligand_volume*water_number_density);
  //cout<<"Volume is "<<ligand_volume<<" giving "<<displaced_waters<<" displaced waters"<< endl;

  
  float entropy_change = avr_desolvation - displaced_waters;


  

  return entropy_change;

}


float LigandEntropy::calculate(vector<Atom*> A)
{

  atoms = A;
  float entropy = 0;

  if(method == "LOG_ATOMS")
    entropy = log((double) atoms.size());
  else if(method == "NO_ATOMS")
    entropy = (double) atoms.size();
  else if(method == "NO_ROT_BONDS")
    entropy = no_rot_bonds();
  else if(method == "NO_ROT_BONDS_COMBINATORIAL")
    entropy = pow(2,no_rot_bonds());
  else if(method == "LOG_ROT_BONDS")
    entropy = log(no_rot_bonds());
  else
    printf("WARNING: Unknown ligand entropy method, '%s'\n", method.c_str());

  //Boltzmann's constant
  //float kb = 1.3807e-23; // J/K
  float entropy_unit = 0.003; //J/(mol K)

  entropy = temperature*entropy_unit*entropy;


  return entropy;
}


float LigandEntropy::no_rot_bonds()
{
  /*** calculate bonds ***/
  BondMatrix BM;
  BM.calculate(atoms);
  bonds = BM.getResult();

  /*** locate rings ***/
  int no_atoms = atoms.size();
  determine_ring_memberships(no_atoms);

  /*** print out info ***/
  /*  for(int i=0;i<bonds.size();i++)
    for(int j=0;j<bonds[i].size();j++)
      if(bonds[i][j] > 0)
	{
	  printf("Atoms '%4s-%3d' and '%4s-%3d' are ",
		 atoms[i]->name.c_str(),i+1,atoms[j]->name.c_str(),j+1);

	  if(is_ring_bond(i,j))
	    printf("Ring  ");
	  else
	    {

	      if(bonds[i][j] ==1 )
		printf("Single");
	      if(bonds[i][j] ==2 )
		printf("Double");
	    }
	  printf(" bound\n");
	}

  for(int i=0;i<bonds.size();i++)
    {
      printf("\n%2d",i+1);
      for(int j=0;j<bonds.size();j++)
	{
	  if(bonds[j][i] ==1 )
	    printf("- ");
	  else if(bonds[j][i] ==2 )
	    printf("= ");
	  else
	    printf(". ");

	}
    }
  printf("\n ");
  for(int i=0;i<bonds.size();i++)
    printf("%2d",i+1);
  printf("\n\n\n");

  */
  /*** find number of rotatable bonds ***/

  float res =0.0;
  
  for(unsigned int i=0;i<bonds.size();i++)
    for(unsigned int j=0;j<i+1;j++)
      if(bonds[j][i] ==1 )  // must be single bond
	{
	  if(!is_ring_bond(i,j))  // must not be part of a ring
	    {
	      if(!is_terminal(i))  // atom i not terminal
		{
		  if(!is_terminal(j))  //atom j not terminal
		    {
		      res = res +1.0;
		      //printf("Atoms '%4s-%3d' and '%4s-%3d' are bound by a ROTATBLE bond\n",
		      //	     atoms[i]->name.c_str(),i+1,atoms[j]->name.c_str(),j+1);
		      
		    }
		  //  else
		    //printf("Atoms '%4s-%3d' and '%4s-%3d' second atom is terminal\n",
		    //	   atoms[i]->name.c_str(),i+1,atoms[j]->name.c_str(),j+1);
		}
	      //      else
		//	printf("Atoms '%4s-%3d' and '%4s-%3d' first atom is terminal\n",
		//      atoms[i]->name.c_str(),i+1,atoms[j]->name.c_str(),j+1);
	    }
	  //  else
	    //printf("Atoms '%4s-%3d' and '%4s-%3d' ring members\n",
	    //	   atoms[i]->name.c_str(),i+1,atoms[j]->name.c_str(),j+1);
	}
  //      else
  //if(bonds[j][i] !=0 )
  //  printf("Atoms '%4s-%3d' and '%4s-%3d' are not single bonded\n",
  //	 atoms[i]->name.c_str(),i+1,atoms[j]->name.c_str(),j+1);
	

  printf("%f rotatable bonds found\n",res);

  
  return res;

}





bool LigandEntropy::is_ring_bond(int i, int j)
{
  bool res = false; 

  for(unsigned int k=0;k<ring_bonds.size();k++)
    if(ring_bonds[k][0] == i && ring_bonds[k][1] == j)
      res = true;
    

  return res;
}

bool LigandEntropy::is_terminal(int i)
{
  /*** is this atom terminal - ie only has one bond to a non-hydrogen atom***/
  int no = 0;
  bool res = true;

  for(unsigned int j=0;j<bonds[i].size();j++)
    if(bonds[i][j] > 0)
      if(atoms[j]->element != "H")
	no++;
  
  if(no>1)
    res = false;

  return res;
}

void LigandEntropy::determine_ring_memberships(int no_atoms)
{
  
  for(int i=0;i<no_atoms;i++)
    !is_ring_member(i,-1,i,1);
  
}


bool LigandEntropy::is_ring_member(int current, int last, int target, int no)
{

  /*** check if we have done this too many times ***/
  if(no > 6)
    return false;
    
  /*** no terminal atoms allowed in ring ***/
  if(is_terminal(current))
    return false;

  /*** find bonded atoms ***/
  vector<int> bonded;
  for(unsigned int i=0;i<bonds.size();i++)
    if(bonds[i][current] > 0)
      if((int)i != last)
	bonded.push_back(i);

  /*** see if target has been found ***/
  bool res = false;
  vector<int>::iterator result;
  result = find( bonded.begin(), bonded.end(), target);            

  if(result == bonded.end()) 
    {
      for(unsigned int i=0; i<bonded.size();i++)
	if(is_ring_member(bonded[i], current, target,no+1))
	  res = true;
    }
  else
    {
      res = true;
      vector<int> temp;
      temp.push_back(current);
      temp.push_back(last);
      ring_bonds.push_back(temp);
    }

  return res;
}              



void LigandEntropy::set_parameters()
{

  ifstream parameters;
  parameters.open("ligand_entropy.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'ligand_entropy.cfg' could not be found\n";
      exit(0);
    }


  string dummy;

  while (!parameters.eof())
    {

      parameters >> dummy;

      if(dummy == "METHOD:")
	parameters >> method;

      if(dummy == "TEMPERATURE:")
	parameters >> temperature;

    } //end read file

  //printing all configurations
  /*
  printf("Configurations for Ligand Entropy:\n");
  printf("\n");
  printf("\n");

  printf("METHOD:                 %s\n", method.c_str());
  printf("\n");

  printf("TEMPERATURE:            %f\n", temperature);
  printf("\n");

  printf("Version: $Id: ligand_entropy.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");
  */
}




void LigandEntropy::writeDescription(FILE * reporter)
{

  rep = reporter;
  /*** write out to out file **/
  version =string("$Id: ligand_entropy.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
