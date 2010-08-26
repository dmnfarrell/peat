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
#include "soupmanip.h"

// Toolbox for handeling soups and soup objects

using namespace std;


vector<SoupObject *> convert_soup_to_objects(Soup * soup)
{
  vector<SoupObject *> res;

  for(unsigned int i=0; i<soup->getProteinChains().size();i++)
    res.push_back(soup->getProteinChains()[i]);
  for(unsigned int i=0; i<soup->getWaters().size();i++)
    res.push_back(soup->getWaters()[i]);
  for(unsigned int i=0; i<soup->getMolecules().size();i++)
    res.push_back(soup->getMolecules()[i]);

  return res;
}

vector<SoupObject *> convert_soup_to_objects(Soup * soup, bool incl_protein, bool incl_molecules, bool incl_water)
{
  vector<SoupObject *> res;

  if(incl_protein)
    for(unsigned int i=0; i<soup->getProteinChains().size();i++)
      res.push_back(soup->getProteinChains()[i]);
  
  if(incl_water)
    for(unsigned int i=0; i<soup->getWaters().size();i++)
      res.push_back(soup->getWaters()[i]);

  if(incl_molecules)
    for(unsigned int i=0; i<soup->getMolecules().size();i++)
      res.push_back(soup->getMolecules()[i]);

  return res;
}



vector<SoupObject *> convert_soups_to_objects(vector<Soup *> soups)
{

  vector<SoupObject *> res;
  for(unsigned int s=0; s<soups.size(); s++)
    {
      for(unsigned int i=0; i<soups[s]->getProteinChains().size();i++)
	res.push_back(soups[s]->getProteinChains()[i]);
      for(unsigned int i=0; i<soups[s]->getWaters().size();i++)
	res.push_back(soups[s]->getWaters()[i]);
      for(unsigned int i=0; i<soups[s]->getMolecules().size();i++)
	res.push_back(soups[s]->getMolecules()[i]);
    }

  return res;
}

vector<Atom *> convert_objects_to_atoms(vector<SoupObject *> objects)
{

  vector<Atom *> res;
  vector<Atom *> temp;
  for(unsigned int s=0; s<objects.size(); s++)
    {
      temp = objects[s]->getAtoms();
      for(unsigned int i=0; i<temp.size();i++)
	res.push_back(temp[i]);
    }

  return res;
}



vector<Residue> cut_residues_from_protein_chains_using_residue_number(vector<ProteinChain*> pc, int number)
{
  vector<Residue> res;

  for(vector<ProteinChain*>::iterator ipc = pc.begin(); ipc != pc.end(); ipc++)
    {
      vector<Residue>* residues = (*ipc)->get_pointer_to_vector_of_residues();
      
      for(vector<Residue>::iterator ir = residues->begin(); ir != residues->end(); )
	{
	  if((*ir).residueNumber == number)
	    {
	      res.push_back(*ir);
	      ir = (*ipc)->delete_residue(ir);
	    }
	  else
	    ir++;

	}
    }

  return res;
}


vector<Atom*> get_atoms_from_protein_chains(vector<ProteinChain*> pcs)
{
  vector<Atom*> res;

  for(vector<ProteinChain*>::iterator pc = pcs.begin(); pc != pcs.end(); pc++)
    {
      vector<Atom*> temp = (*pc)->getAtoms();
      for(vector<Atom*>::iterator atom = temp.begin(); atom != temp.end(); atom++)
	res.push_back(*atom);
    }

  return res;
}


vector<float> setup_boxes (vector<Atom*> atoms, float box_size)
{

  float 
    minx= 10000.0,miny= 10000.0,minz= 10000.0,
    maxx=-10000.0,maxy=-10000.0,maxz=-10000.0;

  vector<Atom*>::iterator atom;

  // find min and max values
  for (atom=atoms.begin(); atom!=atoms.end(); atom++)
    {
      //printf("%8.2f,%8.2f,%8.2f\n",
      //     (*atom)->x, (*atom)->y, (*atom)->z);


      if(minx > (*atom)->x)
	minx = (*atom)->x;

      if(miny > (*atom)->y)
	miny = (*atom)->y;

      if(minz > (*atom)->z)
	minz = (*atom)->z;

      if(maxx < (*atom)->x)
	maxx = (*atom)->x;

      if(maxy < (*atom)->y)
	maxy = (*atom)->y;

      if(maxz < (*atom)->z)
	maxz = (*atom)->z;

    }

  printf("x:[%8.2f:%8.2f]\n",minx,maxx);
  printf("y:[%8.2f:%8.2f]\n",miny,maxy);
  printf("z:[%8.2f:%8.2f]\n",minz,maxz);

  // put atoms into boxes
  for (atom=atoms.begin(); atom!=atoms.end(); atom++)
    {
      (*atom)->bx = (int) ( ((*atom)->x - minx)/box_size );
      (*atom)->by = (int) ( ((*atom)->y - miny)/box_size );
      (*atom)->bz = (int) ( ((*atom)->z - minz)/box_size );
    }

  
  vector<float> res;
  res.push_back(minx);
  res.push_back(miny);
  res.push_back(minz);
  res.push_back(maxx);
  res.push_back(maxy);
  res.push_back(maxz);

  return res;

}
