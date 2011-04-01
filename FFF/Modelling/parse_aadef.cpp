/*
 #
 # FFF - Flexible Force Field
 # Copyright (C) 2010 Jens Erik Nielsen
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
#include "parse_aadef.h"

parse_aadefs::parse_aadefs(vector<string> aadeflines) {
  //
  // mode==0: Amino Acid definition file
  // mode==1: Rotamer library
  //
  //
  // Split the amino acid definition files up into sections separated by "*"
  // and send each block of lines to parse_aadef (parse_aadef parses 1 amino acid definition)
  //      
  printf ("Parsing Topology files containing %d lines\n",static_cast<int>(aadeflines.size()));
  bool done=false;
  unsigned int count=0;
  string line,name;
  vector<string> lines;
  //
  // Init array for amino acids
  //
  amino_acids.resize(0);
  //
  // Skip comments
  //
  bool aapart=true;
  while (!done) {
    if (!aadeflines[count].substr(0,2).compare("//") || !aadeflines[count].substr(0,1).compare("#"))  {
      //
      if (aadeflines[count].size()>25) {
	if (aadeflines[count].substr(0,25).compare("#$$$Entering OTHER.DAT$$$")==0) {
	  aapart=false;
	}
      }
      count++;
      continue;
    }
    done=true;
  }
  done=false;
  printf ("reading first record\n");
  //
  // First record starts
  //
  if (strip(aadeflines[count]).substr(0,1).compare("*")==0) {
    count++;
    name=strip(aadeflines[count]);
    count++;
    while (!done) {
      line=aadeflines[count];
      //
      // Deal with comments
      //
      if (!line.substr(0,2).compare("//") || !line.substr(0,1).compare("#"))  {
	count++;
	printf ("Line:%s %d\n",line.c_str(),line.size());
	if (line.size()>25) {
	  if (line.substr(0,25).compare("#$$$Entering OTHER.DAT$$$")==0) {
	    aapart=false;
	  }
	}
	if (count>=aadeflines.size()) {
	  break;
	}
	continue;
      }
      //
      if (line.substr(0,1).compare("*")==0 || count==aadeflines.size()-1) {
	if (count==aadeflines.size()-1) {
	  lines.push_back(line);
	  printf ("Parsing last record: %s\n",name.c_str());
	  if (aapart) {
	    amino_acids.push_back(name);
	  }
	  aadefs.push_back(parse_aadef(name,lines));
	  return;
	}
	//
	//printf ("Parsing: %s\n",name.c_str());
	if (aapart) {
	  amino_acids.push_back(name);
	}
	aadefs.push_back(parse_aadef(name,lines));
	lines.resize(0);
	count++;
	name=strip(aadeflines[count]);
      }
      else {
	lines.push_back(line);
      }
      count++;
    }
  }
  else if (strip(aadeflines[count]).substr(0,1).compare("#")==0) {
    count++;
  }
  else {
    cout << "ERROR: Wrong format for Amino Acid definition file\n";
    exit(0);
  }
  //
  // Parse the last record
  //
  printf ("Parsing last rescord: %s\n",name.c_str());
  if (aapart) {
    amino_acids.push_back(name);
  }
  aadefs.push_back(parse_aadef(name,lines));
  lines.resize(0);
  cout << "Parsing done\n";
  return;
}

//
// -------------------------------------------------------------
//

parse_aadef::parse_aadef() {
  //
  // Dummy constructor
  //
}

//
// ------------------------------------------------------------
//

parse_aadef::parse_aadef (string exname,vector<string> lines) {
  //
  // Parses a single amino acid definition
  // and instantiates atom_class for each atom
  //
  atoms.resize(0);
  c_alpha=-1;
  n_mainchain=-1;
  c_mainchain=-1;
  string line;
  double x,y,z;
  int atomcounter=0;
  this->name=exname;
  for (int i=0;i<static_cast<int>(lines.size()-2);i++) {
    //
    // Get the X, Y and Z coordinates and the atom name
    //
    line=lines[i];
    char buffer[120];
    string X(line,6,7);
    string Y(line,14,7);
    string Z(line,22,7);
    X.copy(buffer,X.length());
    x=atof(buffer);
    Y.copy(buffer,Y.length());
    y=atof(buffer);
    Z.copy(buffer,Z.length());
    z=atof(buffer);
    string atname(line,1,6);
    atname=strip(atname);
    //
    // Use the general atom_class to store coordinates
    //
    atom_class INS(atname,x,y,z);
    atoms.push_back(INS);
    //
    // See if we should set one of the special pointers
    //
    if (strip(atname)=="CA") {
      this->c_alpha=atomcounter;
    }
    else {
      if (strip(atname)=="N") {
	this->n_mainchain=atomcounter;
      }
      else {
	if (strip(atname)=="C") {
	  this->c_mainchain=atomcounter;
	}
      }
    }
    atomcounter++;
  }
  //
  // Parse the bonds information - this info is in a single line
  // <atom1> <atom2> <atom3> <atom4>
  // Bonds are between 1 and 2, and between 3 and 4
  //
  line=lines[lines.size()-2];
  //printf ("line: %s\n",line.c_str());
  vector<int> numbers;
  //
  // Get the numbers out of the line
  //
  numbers=getints(line);
  for (int i=0;i<static_cast<int>(numbers.size());i++) {
    int first=numbers[i];
    int second=numbers[i+1];
    i++;
    atoms[first-1].bonds.push_back(second-1);
    atoms[second-1].bonds.push_back(first-1);
  }
  //
  // Parse the definition of side chain dihedral angles
  //
  // This info is in the last line of the entry
  // <atom1> <atom2> <atom3> <atom4> <atom5> ...
  // where the dihedral is defined by atoms 1-4
  //  1-2-Dihedral-3-4
  //
  line=lines[lines.size()-1];
  numbers.resize(0);
  numbers=getints(line);
  numchi=numbers[0];
  if (numchi>0) {
    for (int i=1;i<static_cast<int>(numbers.size());i++) {
      dihedral_atoms.push_back(atoms[numbers[i]-1]);
    }
  }
  //
  // An extra safety check
  //
  if (static_cast<unsigned int>(numchi*4)!=dihedral_atoms.size()) {
    cout << "Corrupt entry for torsion angles: " << line << endl;
    exit(0);
  }
  //
  // Construct the information on the # of bonds from CA
  //
  if (c_alpha!=-1) {
    for (int atom=0;atom<static_cast<int>(atoms.size());atom++) {
      //
      // Start at CA, and count number of bonds until we reach the atom
      //
      if (atoms[atom].is_backbone()) {
	//
	// For backbone atoms we set the value to -1
	//
	bonds_from_CA[strip(atoms[atom].name)]=-1;
      } 
      else {
	//
	// Get the real value for all other atoms
	//
	vector<int> unchecked_atoms;
	vector<int> bonds_away;
	unchecked_atoms.push_back(c_alpha);
	bonds_away.push_back(0);
	int exam_atom=0;
	while (true) {
	  if (unchecked_atoms[exam_atom]==atom) {
	    bonds_from_CA[strip(atoms[atom].name)]=bonds_away[exam_atom];
	    //printf ("Bonds from CA:%4s  %2d\n-------------\n",atoms[atom].name.c_str(),
	    //		  bonds_from_CA[strip(atoms[atom].name)]);
	    break;
	  }
	  //
	  // Fill in the next atoms to be checked
	  //
	  for (int bond=0;
	       bond<static_cast<int>(atoms[unchecked_atoms[exam_atom]].bonds.size())
		 ;bond++) {
	    //
	    // Check that we didn't examine this atom already
	    //
	    int atom_to_check=atoms[unchecked_atoms[exam_atom]].bonds[bond];
	    bool already_done=false;
	    for (int test=0;test<static_cast<int>(unchecked_atoms.size());test++) {
	      if (atom_to_check==unchecked_atoms[test]) {
		already_done=true;
		break;
	      }
	    }
	    if (!already_done) {
	      unchecked_atoms.push_back(atoms[unchecked_atoms[exam_atom]].bonds[bond]);
	      bonds_away.push_back(bonds_away[exam_atom]+1);
	    }
	  }
	  exam_atom=exam_atom+1;
	}
      }
    }
  }
  return;
}


