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
#include "proteinchain.h"

ProteinChain::~ProteinChain(){}

ProteinChain::ProteinChain(){}

ProteinChain::ProteinChain(vector<string> lines)
{
  int currentResidue = -100000; //silly value to ensure its not the first residue, eventhough it might be negative
  string residue_insertion_id = " ";
  vector<string> residueLines;
  
  /*** Splitting chain into residues ***/

  for(unsigned int i = 0; i< lines.size(); i++)
    {
      int tempResidue= atoi(string(lines[i],22,4).c_str());
      string temp_residue_insertion_id = string(lines[i],26,1);
      // printf("Found residue: %s\n", tempResidue.c_str());
      if(tempResidue != currentResidue || (temp_residue_insertion_id != residue_insertion_id && temp_residue_insertion_id != " ")) 
	// Make a new residue if the residue sequence number changes, or if a new residue insertion id code is found (but only if the new
	// residue insertion id code is not " " - this is to compensate for WHATIF sometimes forgetting residue insertion ids when protonating)
	{
	  if(currentResidue == -100000)
	    {
	      /*
	      printf("First split - id before: %s and after: %s\n",
		     residue_insertion_id.c_str(),
		     temp_residue_insertion_id.c_str());
	      */
	      currentResidue = tempResidue;
	      residue_insertion_id = temp_residue_insertion_id;
	    }
	  else
	    {
	      if(residueLines.size() > 0)
		{
		  Residue X(residueLines);
		  //		  for(int k=0;k<X.atoms.size();k++)
		  //		    atoms.push_back(X.atoms[k]);
		  residues.push_back(X);
		}
	      residueLines.clear();
	      
	      /*      if (temp_residue_insertion_id != residue_insertion_id)
		printf("Residue split based on insertion code for residue %d-%s and %d-%s - Id before: %s and after: %s\n",
		       currentResidue,
		       string(lines[i-1], 17, 3).c_str(),
		       tempResidue,
		       string(lines[i], 17, 3).c_str(),
		       residue_insertion_id.c_str(),
		       temp_residue_insertion_id.c_str());
	      */
	      currentResidue = tempResidue;
	      residue_insertion_id = temp_residue_insertion_id;

	      //printf("split line -> %s\n",lines[i].c_str());
	    }
	}
      
      if(string(lines[i], 16, 1) == " " || string(lines[i], 16, 1) == "A") //check for "alternative location tag"
	residueLines.push_back(lines[i]);
    }

  /*** last residue ***/
  if(residueLines.size() != 0)
    {
      Residue X(residueLines);
      //      for(int k=0;k<X.atoms.size();k++)
      //	atoms.push_back(X.atoms[k]);
      residues.push_back(X);
      residueLines.clear();
    }

}

void ProteinChain::listAtoms()
{
  for(unsigned int i=0; i<residues.size(); i++)
    residues[i].listAtoms();
}


vector<Atom*> ProteinChain::getAtoms()
{
  vector<Atom*> temp;
  vector<Atom*> temp2;
  for(unsigned int i=0; i<residues.size(); i++)
    {
      temp2 = residues[i].getAtoms();
      for(unsigned int j=0; j<temp2.size(); j++)
	temp.push_back(temp2[j]);
    }
  return temp;
}

int ProteinChain::numberOfAtoms()
{
  return getAtoms().size();
}



vector<Residue*> ProteinChain::getResidues()
{
  vector<Residue*> temp;
  for(unsigned int i = 0; i<residues.size(); i++)
    temp.push_back(&residues[i]);

  return temp;
}


vector<Residue> ProteinChain::copyResidues()
{
  return residues;
}

vector<Residue>* ProteinChain::get_pointer_to_vector_of_residues()
{return &residues;}



vector<Residue>::iterator ProteinChain::delete_residue(int i)
{
  vector<Residue>::iterator it = residues.begin()+i;
  return residues.erase(it);

}

vector<Residue>::iterator ProteinChain::delete_residue(vector<Residue>::iterator it)
{
  return residues.erase(it);
}



void ProteinChain::add_residue(Residue res)
{
  residues.push_back(res);
  return;
}
