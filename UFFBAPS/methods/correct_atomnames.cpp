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
$Id: correct_atomnames.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "correct_atomnames.h"


Correct_Atomnames::Correct_Atomnames(FILE * reporter, vector<Soup*> A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  set_parameters();


  go(A);
}

Correct_Atomnames::Correct_Atomnames():Method()
{
  set_parameters();
}

Correct_Atomnames::~Correct_Atomnames(){}


void Correct_Atomnames::go(vector<Soup*> A)
{
  vector<ProteinChain*> pcs;

  for(unsigned int s=0; s<A.size(); s++)
    {
      printf("Correcting atom names for %s\n",A[s]->name.c_str());
      //extract protein chains
      pcs = A[s]->getProteinChains();
      
      //correct each protein chain
      for(unsigned int i=0; i<pcs.size(); i++)
	{
	  //correct names
	  atoms = pcs[i]->getAtoms();
	  correct_heavy_atoms();
	  correct_hydrogen_atoms();
	}
    }
}

void Correct_Atomnames::correct_heavy_atoms()
{
  /*
    Corrects heavy atom names
  */

  //corrects WhatIf names


  for(unsigned int i=0; i< atoms.size(); i++)
    {
      if(atoms[i]->name == "O''")
	atoms[i]->name = "OXT";
      if(atoms[i]->name == "O'")
	atoms[i]->name = "O";
    }


}


void Correct_Atomnames::correct_hydrogen_atoms()
{
  /*
    Corrects hydrogen atom names
    Assumes that heavy atom names are correct

    Hydrogen name: H{heavy atom tag}{branch id}

    only 1 hydrogen bounded to base heavy atom: no branch id

    more than 1 hydrogen to base heavy atom: branch id for hydrogen i: number of heavy atoms bounded - 1 + i

  */


  vector<int> bound_hydrogens;
  string heavy_atom_tag, name_base;

  for(unsigned int i=0; i< atoms.size(); i++)
    {
      if(atoms[i]->element != "H") //only look at heavy atoms
	{
	  //find bounded hydrogens
	  bound_hydrogens = find_bounded_hydrogens(i);
	  
	  //generate name base
	  heavy_atom_tag = atoms[i]->name.substr(1);
	  name_base = "H";
	  name_base.append(heavy_atom_tag);

	  //only 1 hydrogen bounded
	  if(bound_hydrogens.size()==1)
	    {  
	      atoms[bound_hydrogens[0]]->name = name_base;
	    }
	  //more than 1 hydrogen bounded
	  else if(bound_hydrogens.size() > 1)
	    {
	      //set counter start value
	      int counter = get_number_of_bounded_heavy_atoms(i) - 1;
	      
	      for(unsigned int h=0; h<bound_hydrogens.size(); h++)
		{
		  char no[6];
		  sprintf(no, "%d", counter);
		  atoms[bound_hydrogens[h]]->name = name_base+no;
		  counter++;
		}
	    }
	}//if H
    }//for loop

}


vector<int> Correct_Atomnames::find_bounded_hydrogens(int a)
{
  vector<int> res;

  for(unsigned int i=0; i< atoms.size(); i++)
    if(atoms[i]->element == "H")
      if(get_distance(i,a) < squareHdist)
	res.push_back(i);
  
  return res;
}

int Correct_Atomnames::get_number_of_bounded_heavy_atoms(int a)
{
  int res=0;
  for(unsigned int i=0; i< atoms.size(); i++)
    if(atoms[i]->element != "H")
      {
	if(atoms[i]->element == "S" && atoms[a]->element == "S")
	  {
	    if(get_distance(i,a) < squareSSdist)
	      res ++;
	  }
	else
	  if(get_distance(i,a) < squaredefaultDist)
	    res ++;
      }

  return res;
}

float Correct_Atomnames::get_distance(int i, int j) //to save memory instead of doing a distance matrix
{
  float res = 
    (atoms[i]->x - atoms[j]->x)*(atoms[i]->x - atoms[j]->x)+
    (atoms[i]->y - atoms[j]->y)*(atoms[i]->y - atoms[j]->y)+
    (atoms[i]->z - atoms[j]->z)*(atoms[i]->z - atoms[j]->z);
  
  return res;

}


void Correct_Atomnames::set_parameters()
{
  SSdist = 3.5;
  Hdist = 1.5; 
  defaultDist = 2.0;

  squareSSdist = SSdist*SSdist;
  squareHdist = Hdist*Hdist;
  squaredefaultDist= defaultDist*defaultDist;

}

void Correct_Atomnames::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: correct_atomnames.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());


}
