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
$Id: crystallized_water_bridges.cpp 18 2005-11-15 09:56:08Z chresten $
*/

#include "crystallized_water_bridges.h"


Crystallized_Water_Bridges::Crystallized_Water_Bridges():Method()
{
  set_parameters();  
}

Crystallized_Water_Bridges::Crystallized_Water_Bridges(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();  
}

Crystallized_Water_Bridges::Crystallized_Water_Bridges(FILE * reporter, FILE * resultsFile, vector<Soup *> A, vector<Soup *> B):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();  

  double result = 0;

  if(A.size() != B.size())
    {
      cout << "WARNING: An equal amount of protein and ligand structures must be given to the water bridge function!\n";
      return;
    }
  
  for(unsigned int i=0; i<A.size(); i++)
    {
      result = calculate(A[i],B[i]);
      
      
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

Crystallized_Water_Bridges::Crystallized_Water_Bridges(FILE * resultsFile, FILE * reporter, vector<Soup *> A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();  
    

  double result = 0;

  for(unsigned int i=0; i<A.size(); i++)
    {
      result = calculate(A[i]);
      printf("Water bridge for %s gave %8.2f\n",A[i]->name.c_str(), result);
    }  

}


Crystallized_Water_Bridges::~Crystallized_Water_Bridges(){}




float Crystallized_Water_Bridges::calculate(Soup * A)
{
  internal_dry_atoms = A->get_dry_atoms();

  float res = 0.0;

  vector<Water*> waters = A->getWaters();
  vector<Water*>::iterator water;

  for(water = waters.begin(); water != waters.end(); water++)
    {
      // find water oxygen atom
      Atom o;
      Atom * po = &o;
      (*water)->get_oxygen(po);

      if(po->name != "Dummy")
	  res += check_for_interal_water_bridge(po);
    }

  return water_bridge_energy*res;
}


float Crystallized_Water_Bridges::calculate(Soup * A, Soup * B)
{
  // remeber to get waters from both soups
  dry_atoms_A = A->get_dry_atoms();
  dry_atoms_B = B->get_dry_atoms();

  float res = 0.0;


  vector<Water*> waters = A->getWaters();
  vector<Water*> temp = B->getWaters();
  waters.insert(waters.end(), temp.begin(), temp.end());

  vector<Water*>::iterator water;

  for(water = waters.begin(); water != waters.end(); water++)
    {
      // find water oxygen atom
      Atom o;
      Atom * po = &o;
      (*water)->get_oxygen(po);

      if(po->name != "Dummy")
	  res += check_for_water_bridge(po);
    }

  return water_bridge_energy*res;

}



float Crystallized_Water_Bridges::check_for_interal_water_bridge(Atom * oxygen)
{
  float res = 0,d;

  vector<Atom*>::iterator atom;
  vector<int> residues;


  // find all atoms that h-bond to the water oxygen
  for(atom = internal_dry_atoms.begin(); atom != internal_dry_atoms.end(); atom++)
    {
      d = dist.calculate(*atom, oxygen, false);
      if(HB.water_bridge_hbond(*atom, oxygen, d) || HB.water_bridge_hbond(oxygen, *atom, d) )
	if (find(residues.begin(), residues.end(), (*atom)->residue_no) == residues.end()) // make sure we don't count multiple bonds to a residue
	  residues.push_back((*atom)->residue_no);
    }	  
  

  vector<int>::iterator res1, res2;

  //  printf("\nWater bridges for:  ");
  //  oxygen->display();
	
  for(res1 = residues.begin(); res1 != residues.end(); res1++)
    for(res2 = res1+1; res2 != residues.end(); res2++)
      {
	//printf("  %3d-%3d: ",*res1,*res2);
	if (*res1-*res2 == 0)
	  {
	    printf("  This should not happen!!!\n");

	  }

	if (abs(*res1-*res2) == 1)
	  {
	    res+= neighbour_residue_factor;
	    //printf("%8.2f\n",neighbour_residue_factor);

	  }
	else
	  {
	    res+= 1.0;
	    //printf("%8.2f\n",1.0);

	  }
	    
      }


  return res;
}








float Crystallized_Water_Bridges::check_for_water_bridge(Atom * oxygen)
{

  float res = 0,d;

  vector<Atom*>::iterator atom;
  vector<int> residues_A, residues_B;

  // find all atoms in A that h-bond to the water oxygen 
  for(atom = dry_atoms_A.begin(); atom != dry_atoms_A.end(); atom++)
    {
      d = dist.calculate(*atom, oxygen, false);
      if(HB.water_bridge_hbond(*atom, oxygen, d) || HB.water_bridge_hbond(oxygen, *atom, d) )
	if (find(residues_A.begin(), residues_A.end(), (*atom)->residue_no) == residues_A.end()) // make sure we don't count multiple bonds to a residue
	  residues_A.push_back((*atom)->residue_no);
    }	  

  // find all atoms in B that h-bond to the water oxygen 
  for(atom = dry_atoms_B.begin(); atom != dry_atoms_B.end(); atom++)
    {
      d = dist.calculate(*atom, oxygen, false);
      if(HB.water_bridge_hbond(*atom, oxygen, d) || HB.water_bridge_hbond(oxygen, *atom, d) )
	if (find(residues_B.begin(), residues_B.end(), (*atom)->residue_no) == residues_B.end()) // make sure we don't count multiple bonds to a residue
	  residues_B.push_back((*atom)->residue_no);
    }	  
  

  vector<int>::iterator res1, res2;

  // printf("\nWater bridges for:  ");
  // oxygen->display();
	
  for(res1 = residues_A.begin(); res1 != residues_A.end(); res1++)
    for(res2 = residues_B.begin(); res2 != residues_B.end(); res2++)
      {
	//printf("  %3d-%3d: ",*res1,*res2);
	res+= 1.0;
	//printf("%8.2f\n",1.0);

	
	    
      }


  return res;




}


void Crystallized_Water_Bridges::set_parameters()
{
  same_residue_factor = 0.0;
  neighbour_residue_factor = 0.5;
  water_bridge_energy = -1.5; //kJ/mol

}




void Crystallized_Water_Bridges::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: crystallized_water_bridges.cpp 18 2005-11-15 09:56:08Z chresten $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
