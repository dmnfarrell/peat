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

#include "soup.h"


Soup::Soup()
{
  isActive=true;
  bonds_assigned=false;
  sybyl_types_assigned=false;
}

Soup::~Soup(){}

void Soup::addToSoup(string filename)
{
  //  printf("Adding contents of %s to the soup.\n",filename.c_str());

  vector<string> lines;
  string line;

  ifstream file(filename.c_str());

  if(!file.is_open())
    {
      printf("File %s not found!\n",filename.c_str());
      exit(0);
    }

  while(!file.eof())
    {
      getline(file, line);
      lines.push_back(line);
    }

  printf("File %s has been read:\n",filename.c_str());
  
  int loc = filename.rfind( ".");
  string fileExtension(filename, loc, filename.size()-loc);
 
  if(fileExtension == ".pdb")
    interpretatePDBLines(lines,filename);
  else if(fileExtension == ".mol2")
    interpretateMOL2Lines(lines,filename);
  else
    {
      printf("Unknown file extension: '%s'\n",fileExtension.c_str());
      exit(0);
    }
}

void Soup::addToSoup(vector<string> lines, string filename)
{
  //printf("Adding some lines from %s to the soup.\n",filename.c_str());

  int loc = filename.rfind( ".");
  string fileExtension(filename, loc, filename.size()-loc);

  if(fileExtension == ".pdb")
    interpretatePDBLines(lines,filename);
  else if(fileExtension == ".mol2")
    interpretateMOL2Lines(lines,filename);
  else
    {
      printf("Unknown file extension: '%s'\n",fileExtension.c_str());
      exit(0);
    }
 
}



void Soup::soupStat()
{

  /*** counting atoms ***/
  int noProteinAtoms = 0;
  int noWaterAtoms = 0;
  int noOtherAtoms = 0;
  int noAtoms = 0;
  int noObjects = 0;

  for(unsigned int i=0; i<proteinChains.size();i++)
    noProteinAtoms += proteinChains[i].numberOfAtoms();
  
  for(unsigned int i=0; i<waters.size();i++)
    noWaterAtoms += waters[i].numberOfAtoms();
  
  for(unsigned int i=0; i<molecules.size();i++)
    noOtherAtoms += molecules[i].numberOfAtoms();
  
  noAtoms = noProteinAtoms + noWaterAtoms + noOtherAtoms;
  noObjects = proteinChains.size() + waters.size() + molecules.size();

  printf("---------------------------------------------------\n");
  printf("The soup '%s' contains:\n",name.c_str());
  printf("Object                Number                Atoms\n");
  printf("---------------------------------------------------\n");
  printf("Protein chain         %6d               %6d\n",proteinChains.size(), noProteinAtoms);
  printf("Water molecule        %6d               %6d\n",waters.size(), noWaterAtoms);
  printf("Other molecules/ions  %6d               %6d\n",molecules.size(), noOtherAtoms);  
  printf("---------------------------------------------------\n");
  printf("In total              %6d               %6d\n",noObjects, noAtoms);  
  printf("---------------------------------------------------\n");
}

void Soup::listAtoms()
{
  printf("Listing all atoms in soup:\n");

  printf("Proteins:\n");
  for(unsigned int i=0; i<proteinChains.size();i++)
    proteinChains[i].listAtoms();
  printf("Water:\n");
  for(unsigned int i=0; i<waters.size();i++)
    waters[i].listAtoms();
    printf("Other:\n");
  for(unsigned int i=0; i<molecules.size();i++)
    molecules[i].listAtoms();
}


vector<Atom *> Soup::getAtoms()
{
  vector<Atom*> temp;
  vector<Atom*> temp2;
  for(unsigned int i=0; i<proteinChains.size();i++)
    {
      temp2 = proteinChains[i].getAtoms();
      for(unsigned int j=0; j<temp2.size(); j++)
	temp.push_back(temp2[j]);
    }
  for(unsigned int i=0; i<waters.size();i++)
    {
      temp2 = waters[i].getAtoms();
      for(unsigned int j=0; j<temp2.size(); j++)
	temp.push_back(temp2[j]);
    }
  for(unsigned int i=0; i<molecules.size();i++)
    {
      temp2 = molecules[i].getAtoms();
      for(unsigned int j=0; j<temp2.size(); j++)
	temp.push_back(temp2[j]);
    }

  return temp;
}


vector<Atom *> Soup::get_dry_atoms() //return non-water atoms
{
  vector<Atom*> temp;
  vector<Atom*> temp2;
  for(unsigned int i=0; i<proteinChains.size();i++)
    {
      temp2 = proteinChains[i].getAtoms();
      for(unsigned int j=0; j<temp2.size(); j++)
	temp.push_back(temp2[j]);
    }
  for(unsigned int i=0; i<molecules.size();i++)
    {
      temp2 = molecules[i].getAtoms();
      for(unsigned int j=0; j<temp2.size(); j++)
	temp.push_back(temp2[j]);
    }

  return temp;
}



vector<vector<SoupObject *> > Soup::getSoupObjects()
{
  vector<vector<SoupObject *> > res;
  vector<SoupObject *> temp;
  
  for(unsigned int i = 0; i<proteinChains.size(); i++)
    temp.push_back(&proteinChains[i]);
  res.push_back(temp); temp.clear();

  for(unsigned int i = 0; i<molecules.size(); i++)
    temp.push_back(&molecules[i]);
  res.push_back(temp); temp.clear();
  
  for(unsigned int i = 0; i<waters.size(); i++)
    temp.push_back(&waters[i]);
  res.push_back(temp); temp.clear();
  
  return res;
}


vector<ProteinChain *> Soup::getProteinChains()
{
  vector<ProteinChain *> temp;
  for(unsigned int i = 0; i<proteinChains.size(); i++)
    temp.push_back(&proteinChains[i]);

  return temp;
}

vector<Water *> Soup::getWaters()
{
  vector<Water*> temp;
  for(unsigned int i = 0; i<waters.size(); i++)
    temp.push_back(&waters[i]);

  return temp;

}
vector<Molecule *> Soup::getMolecules()
{
  vector<Molecule*> temp;
  for(unsigned int i = 0; i<molecules.size(); i++)
    temp.push_back(&molecules[i]);

  return temp;

}

vector<ProteinChain> Soup::copyProteinChains()
{
  return proteinChains;
}




void Soup::add_protein_chain(ProteinChain object)
{
  //cout<<"Adding "<<object.identify()<<endl;
  proteinChains.push_back(object);
  return;
}

void Soup::add_molecule(Molecule object)
{
  //cout<<"Adding "<<object.identify()<<endl;
  molecules.push_back(object);
  return;
}

void Soup::add_water(Water object)
{
  //cout<<"Adding "<<object.identify()<<endl;
  waters.push_back(object);
  return;
}

void Soup::delete_protein_chain(int i)
{
  vector<ProteinChain>::iterator it = proteinChains.begin()+i;
  cout<<"deleting "<<(*it).identify()<<endl;
  proteinChains.erase(it);
  return;
}


void Soup::clearSoup()
{
  proteinChains.clear();
  waters.clear();
  molecules.clear();

  return;
}

/*** assumes that mol2 files only contain molecules - not protein or water ***/
void Soup::interpretateMOL2Lines(vector<string> lines, string filename)
{
  /*** number of molecules already in soup ***/
  int oldMolecules = molecules.size();
  
  vector<string> tempLines;
  vector<vector<string> > moleculeLines;
  bool reading = false;

  /*** divide lines into molecules ***/
  for(unsigned int i=0; i<lines.size(); i++)
    {
      string atomTag(lines[i],0,13);
      string anyTag(lines[i],0,1);
      //      printf("Reading anyTag:'%s' and atomTag:'%s'\n", anyTag.c_str(), atomTag.c_str());
 
      if(anyTag=="@" && tempLines.size() != 0)
	{
	  moleculeLines.push_back(tempLines);
	  reading = false;
	  tempLines.clear();
	}

     if(reading)
	tempLines.push_back(lines[i]);

      if(atomTag == "@<TRIPOS>ATOM")
	 reading = true;
    }  

  /*** last molecule ***/
  if(tempLines.size() != 0)
    moleculeLines.push_back(tempLines);

  //printf("Making %d molecules\n",moleculeLines.size());
  /*** make molecules from lines ***/
  for(unsigned int i=0; i<moleculeLines.size(); i++)
    {
      Molecule X(moleculeLines[i]);
      molecules.push_back(X);
    }

 /*** soup object ids ***/
 for(unsigned int i=oldMolecules;i<molecules.size();i++)
   {
     char no[6];
     sprintf(no, "%d", i);
     molecules.at(i).id=filename + string(": Molecule no ")+string(no);
    }

}

void Soup::interpretatePDBLines(vector<string> lines, string filename)
{

  /*** determining number of objects already in soup ***/
  int oldProteinChains = proteinChains.size();
  int oldMolecules = molecules.size();
  int oldWaters = waters.size();

  vector<string> tempLines;
  vector<vector<string> > chains;

  /*** devide lines into chains ***/
  string currentType = "-";
  for(unsigned int i=0;i<lines.size();i++)
    {
      string identifier(lines[i],0,6);
      if(identifier ==currentType) //adding to current chain
	tempLines.push_back(lines[i]);
      else if(tempLines.size() != 0)
	{
	  chains.push_back(tempLines);
	  tempLines.clear();
	  currentType ="-";
	}
      if(identifier == "ATOM  " && currentType == "-")//starting atom chain
	{
	  tempLines.push_back(lines[i]);
	  currentType ="ATOM  ";
	}

      if(identifier == "HETATM" && currentType == "-") //starting hetatm chain
	{
	  tempLines.push_back(lines[i]);
	  currentType = "HETATM";
	}


    }

  /*** doing the last chain ***/
  if(tempLines.size() != 0) 
    {
      chains.push_back(tempLines);
      tempLines.clear();
    }

  /*** test chain division ***/
  /*printf("Found following chains:\n");
  for(unsigned int i=0;i<chains.size();i++)
      printf("  Chain %d: has %d elements \n", i, chains[i].size());
  */
  /*** make soup objects from chains ***/
  for(unsigned int i=0; i<chains.size(); i++)
    {
      vector<string> chain =chains[i];
      string identifier(chain[0],0,6);

      /*** ATOM  - Protein chain ***/
      if(identifier == "ATOM  ")
	{
	  //check that all ids are ATOM
	  bool allAtom = true;
	  for(unsigned int j=0; j<chain.size(); j++)
	    {
	      string identifier2(chain[j],0,6);
	      if(identifier2 != "ATOM  ")
		allAtom = false;  
	    }

	  if(allAtom)
	    {
	      //printf("All lines in chain %d have ATOM identifiers \n",i);
	      //printf(" - will be treated as protein chain\n");
	      ProteinChain X(chain);
	      proteinChains.push_back(X);
	    }
	  else
	    {
	      printf("Warning: Not all lines in chain %d have ATOM ids\n",i);
	      ProteinChain X(chain);
	      proteinChains.push_back(X);

	    }


	}
      /*** HETATM  - water, ions, other molecules ***/
      else 
	{
	  //printf("Chain %d contains HETATMs \n",i);
	  
	  string currentMolecule = "dummy";
	  vector<string> moleculeLines;
	  bool isWater = true;

	  for(unsigned int j=0; j<chain.size(); j++)
	    {
	      if(currentMolecule != string(chain[j],22,4))
		{
		  if(currentMolecule == "dummy")
		    currentMolecule = string(chain[j],22,4);
		  else
		    {
		      /*** create water molecule ***/
		      if(isWater)
			{
			  Water X(moleculeLines);
			  waters.push_back(X);
			  moleculeLines.clear();
			  currentMolecule = string(chain[j],22,4);
			}
		      /*** Create ion/other molecule ***/
		      else
			{
			  //printf("Found other molecule: %s \n",string(moleculeLines[0],17,3).c_str());
			  Molecule X(moleculeLines);
			  molecules.push_back(X);
			  moleculeLines.clear();
			  currentMolecule = string(chain[j],22,4);
			  isWater = true;
			}
		    }
		}
	      moleculeLines.push_back(chain[j]);
	      string thisName(chain[j],17,3);
	      if(thisName != "HOH")
		isWater = false;  
	      

	    } //loop over chain 


	  /*** last molecule ***/

	  if(moleculeLines.size() != 0)
	    {
	      if(isWater)
		{
		  Water X(moleculeLines);
		  waters.push_back(X);
		  moleculeLines.clear();
		}
	      else
		{
		  printf("Found other molecule: %s \n",string(chain[chain.size()-1],17,3).c_str());
		  Molecule X(moleculeLines);
		  molecules.push_back(X);
		  moleculeLines.clear();
		}
	    }	


	} // if HETATM      
    } // loop over chains  


  /*** soupObject ids ***/


  for(unsigned int i=oldProteinChains;i<proteinChains.size();i++)
    {
      char no[6];
      sprintf(no, "%d", i);
      proteinChains.at(i).id=filename + string(": Protein chain no ")+string(no);
    }
  for(unsigned int i=oldMolecules;i<molecules.size();i++)
    {
      char no[6];
      sprintf(no, "%d", i);
      molecules.at(i).id=filename + string(": Molecule no ")+string(no);
    }
  for(unsigned int i=oldWaters;i<waters.size();i++)
    {
      char no[6];
      sprintf(no, "%d", i);
      waters.at(i).id=filename + string(": Water no ")+string(no);
    }


}


void Soup::set_sybyl_atom_types()
{

  // Setup sybyl atom types
  sybyl_types.clear();
  ifstream in;
  char line[256];
  string s;
  in.open("parameters/sybyl_types.txt");

  if(! in.is_open())
    {
      cout << "Error:File 'parameters/sybyl_types.txt' could not be found\n";
      exit(0);
    }

  while(!in.eof())
    {
      in.getline(line,256);
      //cout<<line<<endl;
      s = string(line);
      
      if (s.size()>3)
	{
	  string key = s.substr(0,8);
	  for(unsigned int i=0; i<key.size(); )
	    if(string(key,i,1) ==" ")
	      key = key.erase(i, 1);
	    else
	      i++;

	  string atom_name = s.substr(8,14);

	  for(unsigned int i=0; i<atom_name.size(); )
	    if(string(atom_name,i,1) ==" ")
	      atom_name = atom_name.erase(i, 1);
	    else
	      i++;


	  sybyl_types[key]= atom_name.c_str();
	  
	  //	  cout<<"Atom "<<key<<" has the sybyl type "<< sybyl_types[key]<<endl;
	}
      

   }

  vector<Atom*> atoms = getAtoms();
  vector<Atom*>::iterator it = atoms.begin();
  
  while(it != atoms.end())
    {
      if((*it)->residue == "HIE" || (*it)->residue == "HID" )
	(*it)->residue = "HIS";
      if((*it)->residue == "CYX")
	(*it)->residue = "CYS";
      
      if ((*it)->name=="O''")
	(*it)->name="OXT";
      if ((*it)->name=="O'")
	(*it)->name="O";
  
      string k = (*it)->residue+"-"+(*it)->name;
      (*it)->generic_key = k;

      if ((*it)->sybyl_type=="" && (*it)->element != "H")
	(*it)->sybyl_type = sybyl_types[k];
      //cout<<"KEY:"<<k<<" gives: "<<(*it)->sybyl_type << "   in"<< name<<endl;
      
      
      it++;
    }

  sybyl_types_assigned=true;  
  return;

}


// operators

Soup Soup::operator+(Soup &other)
{
  Soup res;

  // add names
  if (name != "")
    res.name = name+" + "+other.name;
  else
    res.name = other.name;

  // add protein chains
  vector<ProteinChain *> these_protein_chains = getProteinChains();
  for(unsigned int i=0; i<these_protein_chains.size(); i++)
    res.add_protein_chain(*these_protein_chains[i]);
  
  vector<ProteinChain *> other_protein_chains = other.getProteinChains();
  for(unsigned int i=0; i<other_protein_chains.size(); i++)
    res.add_protein_chain(*other_protein_chains[i]);

  // add molecules
  vector<Molecule *> these_molecules = getMolecules();
  for(unsigned int i=0; i<these_molecules.size(); i++)
    res.add_molecule(*these_molecules[i]);
  
  vector<Molecule *> other_molecules = other.getMolecules();
  for(unsigned int i=0; i<other_molecules.size(); i++)
    res.add_molecule(*other_molecules[i]);

  // add waters
  vector<Water *> these_waters = getWaters();
  for(unsigned int i=0; i<these_waters.size(); i++)
    res.add_water(*these_waters[i]);
  
  vector<Water *> other_waters = other.getWaters();
  for(unsigned int i=0; i<other_waters.size(); i++)
    res.add_water(*other_waters[i]);


  return res;


}



