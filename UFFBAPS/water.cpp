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

#include "water.h"

Water::Water(string line)
{
  Atom X(line);
  atoms.push_back(X);

}

Water::Water(vector<string> lines)
{
  for(unsigned int i=0; i<lines.size(); i++)
    {
      Atom X(lines[i]);
      atoms.push_back(X);
    }

}

int Water::numberOfAtoms()
{
  return atoms.size();
}

void Water::listAtoms()
{
  for(unsigned int i=0; i<atoms.size(); i++)
    atoms[i].display();
}

vector<Atom *> Water::getAtoms()
{
  vector<Atom*> temp;
  for(unsigned int i = 0; i<atoms.size(); i++)
    temp.push_back(&atoms[i]);

  return temp;

}

void Water::get_oxygen(Atom * a)
{
  for(unsigned int i=0; i<atoms.size(); i++)
    if (atoms[i].element == "O")
      {
	*a = atoms[i];
	return;
      }

  return;
}


Water::~Water(){}


