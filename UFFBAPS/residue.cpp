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
#include "residue.h"


Residue::Residue(vector<string> lines)
{
  residueType = string(lines[0], 17, 3);
  residueNumber = atoi((lines[0].substr(22, 4)).c_str());


  for(unsigned int i =0; i < lines.size(); i++)
    {
      //string name = string(lines[i], 17, 3);
      //int number = atoi((lines[i].substr(22, 4)).c_str());

      /*
      if(name != residueType)
	printf("Warning: residue is corrupted - found '%d-%s' and '%d-%s' in same residue!\n",number,name.c_str(),residueNumber,residueType.c_str());
      
      if(number != residueNumber)
	printf("Warning: residue is corrupted - found '%d-%d' and '%d-%d' in same residue!\n",number,name.c_str(),residueNumber,residueType.c_str());
      */

      Atom X(lines[i]);
      atoms.push_back(X);
      
    }



}
Residue::~Residue(){}

void Residue::listAtoms()
{
  for(unsigned int i=0; i<atoms.size(); i++)
    atoms[i].display();

}

vector<Atom*> Residue::getAtoms()
{
  vector<Atom*> temp;
  for(unsigned int i = 0; i<atoms.size(); i++)
    temp.push_back(&atoms[i]);

  return temp;
}

void Residue::delete_backbone_atoms()
{
  string bb[] = {"C", "CA", "N", "O"};
  vector<string> backbone_atoms(bb,bb+4);



  for(vector<Atom>::iterator a = atoms.begin(); a != atoms.end();)
    if (find(backbone_atoms.begin(), backbone_atoms.end(), (*a).name ) != backbone_atoms.end())
      {
	cout<<"deleting atom "<<a->name<<endl;
	a = atoms.erase(a);
      }
    else
      a++;
  

}


void Residue::find_atom_named(string name, Atom * a)
{
  // This method assumes that all atoms in a residue has unique names!

  vector<Atom*> atoms = getAtoms();

  for(vector<Atom*>::iterator atom = atoms.begin(); atom != atoms.end(); atom++)
    if((*atom)->name==name)
      {
	*a = **atom;
	return;
      }

  return;
}
			     
