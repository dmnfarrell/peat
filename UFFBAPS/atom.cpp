#include "atom.h"
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

Atom::Atom()
{
  name = "Dummy";
  hydrogen_bond_could_be_replaced_by_water=false;
  return;
}


Atom::Atom(Vector v, int residue_number)
{
  name = "V";
  x = v.x;
  y = v.y;
  z = v.z;

  
  residue = "HOH";
  residue_no = residue_number;
  element = "O";
  vdw_dist = 0;
  vdw_energy = 0;
  charge = 0;
  hydrogen_bond_could_be_replaced_by_water = false;
  occupancy = 0.0;
  tempFactor = 0.0;
  return;
}


Atom::Atom(string line)
{

  is_hydrogen = false;
  hydrogen_bond_could_be_replaced_by_water = false;

  string identifier(line,0,6);

  if(identifier == "ATOM  " || identifier == "HETATM") //PDB file
    {
      /*** Parsing name and element ***/
      name =removeWhite(string(line,12,4));
      residue = string(line, 17, 3);
      residue_no = atoi(string(line, 22, 4).c_str());
      element = removeNumbers(removeWhite(string(line,12,2)));

      if (element == "H")
	is_hydrogen = true;

      sybyl_type = "";

      /*** Parsing coordinates ***/
      x = atof(string(line,30,8).c_str());
      y = atof(string(line,38,8).c_str());
      z = atof(string(line,46,8).c_str());
      
      /*** Parsing occupacy and tempFactor ***/
      occupancy = atof(string(line,54,6).c_str());
      tempFactor = atof(string(line,60,6).c_str());
      charge = 0;

    }
  else //MOL2 file
    {
      
      /*** break line into words ***/
      vector<string> words;
      words.clear();
      char temp[300];
      istringstream stream(line.c_str());
      while(stream>>temp)
	words.push_back(removeWhite(string(temp)));
    
      /*** setting name and element ***/
      name = words.at(1);
      string tempElem = words.at(5);
      sybyl_type = tempElem;
      int loc = tempElem.find( ".", 0 );
      element = string(tempElem, 0, loc);
      
      if (element == "H")
	is_hydrogen = true;

      charge = atof(words.at(8).c_str());
      residue_no = 0;

      /*** setting coordinates ***/
      x = atof(words.at(2).c_str());
      y = atof(words.at(3).c_str());
      z = atof(words.at(4).c_str());

      /*** correcting GAFF types ***/
      char letter = element.substr(0,1)[0];
      if(isupper(letter) == 0)
      	{
	  if(element == "br")
	    element = "BR";
	  else if(element =="Cl")
	    element = "CL";
	  if(element == "ca")
	    element = "CA";
	  else if(element =="na")
	    element = "NA";

	  else
	    {
	      //Element is capitalized version of first letter
	      element.clear();
	      element.push_back(toupper(letter));
	    }
	}

    }


  /*** testing element ***/
  string validElements[] ={"C","H","O","N","CL","F","S","CA","NA","P","BR","I","ZN","MG","U","MN","CD","FE","K","NI","CO","CU"};
  bool valid = false;
  for(int i=0;i<22;i++)
    if(element == validElements[i])
      valid = true;
  //  printf("element found:%s!\n",element.c_str());
  if(!valid)
    {
      if(element.substr(0,1) == "H")
	{
	  //	  printf("Warning: Guessing that '%s' is hydrogen!\n",element.c_str());
	  element = "H";
	}
      else if(element == "Cl")
	  element = "CL";
      else if(element == "Br")
	  element = "BR";
      else
	printf("Warning: Invalid element found:%s!\n",element.c_str());

    }
}


void Atom::display()
{
  try{
    printf("Atom: %4s (%3s %4d) %7.3f %7.3f %7.3f %2s %3.5f %3.5f %3.5f\n", 
	   name.c_str(),
	   residue.c_str(),
	   residue_no,
	   x,y,z,element.c_str(),
	   vdw_dist, vdw_energy, charge);
    /*
  printf("%4s %4s %4s", 
	 name.c_str(),adv_type.c_str(),element.c_str());
    */
  }
  catch(...){
    printf("WARNING: Can't display atom\n");
  }

}

Vector Atom::get_pos()
{
  return Vector(x,y,z);
}
