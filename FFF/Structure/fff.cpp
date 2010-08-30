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
#include "fff.h"
using namespace std;

void FFF::read_pdb (const std::string pdbfilename) {
  //
  // Initialise class. Read a pdb file
  //
  string line;
  vector<string> lines;
  //
  //  Read the file
  //
  printf("Starting to read pdb file; %s\n",pdbfilename.c_str());
  ifstream file(pdbfilename.c_str());
  if (!file) {
    printf("PDB file not found!: %s\n",pdbfilename.c_str());
    exit(0);
  }
  pdbname=pdbfilename; //store pdbfile as protein-name
  // Put the entire file in lines
  while(getline(file,line)) {
     lines.push_back(line);
  }
  file.close();
  printf("File has been read.\n");
  //
  // Call the parse_lines routine
  //
  parse_lines(lines);
    return;
}
//
// -----------------------------------------
//

void FFF::remove_atoms_with_tag(const std::string tag) {
  update_all_atoms();
  vector<atom_class> delete_atoms;
  for (unsigned int atom=0;atom<all_atoms.size();atom++) {
    //printf ("looking for tag: '%s' with: '%s'\n",tag.c_str(),((*(all_atoms[atom])).tag.c_str()));
    //printf ("%d - %d\n",tag.size(),((*(all_atoms[atom])).tag).size());
    //(*(all_atoms[atom])).print_detail();
    if (strip(tag)==(strip(all_atoms[atom]->tag))) {
      delete_atoms.push_back((*(all_atoms[atom])));
      //remove_atom((*(all_atoms[atom])));
    }
  }
  for (unsigned int atom=0;atom<delete_atoms.size();atom++) {
    remove_atom(delete_atoms[atom]);
  }
  check_soup();
  update_all_atoms();
  return;
}

//
// ----
//

void FFF::remove_atom(atom_class& atom) {
  // Remove this atom from the residue
  //printf ("calling remove_atom with");
  //atom.print_detail();
  int inchain=atom.inchain;
  int inresidue=atom.inresidue;
  for (unsigned int latom=0;latom<chains[inchain].residues[inresidue].atoms.size();latom++) {
    if (chains[inchain].residues[inresidue].atoms[latom]==atom) {
      //printf ("Erasing: ");
      //(*(chains[inchain].residues[inresidue].atoms.begin()+latom)).print_detail();
      chains[inchain].residues[inresidue].atoms.erase(chains[inchain].residues[inresidue].atoms.begin()+latom);
      //chains[inchain].residues[inresidue].update();
      break;
    }
  }
  return;
}

//
// -----------------------------------------
//

void FFF::parse_lines(const std::vector<std::string> lines) {
  //
  // Parse all the lines we got
  //
  //printf("Parsing lines\n");
  //
  // Count the number of atoms and the number of residues
  nresidues=0;
  natoms=0;
  string lastnum="   ";
  string line;
  for (int i=0;i<static_cast<int>(lines.size());i++) {
    line=lines[i];
      if (line.length()>=6) {
    string ID(line,0,6);
    if (ID=="ATOM  " || ID=="HETATM") {
      natoms++;
      string thisnum(line,22,4);
      if (thisnum!=lastnum) {
	     nresidues++;
	     lastnum=thisnum;
      }
    }
      }
  }
  printf("Number of residues=%d, Number of atoms=%d\n",nresidues,natoms);
  //
  // Look for 'TER' records or changes in the chain ID
  //
  vector<string> chain_lines;
    string oldchainID("None");
    for(int i=0;i<static_cast<int>(lines.size());i++) {
        line=lines[i];
        if (line.length()<1) {
                continue;
       }
        //printf ("Line: '%s'\n",line.c_str());
        vector<string> sp=split(strip(line));
        if (sp[0]!="ATOM" && sp[0]!="HETATM" && sp[0]!="TER") {
            continue;
        }
        string SS(line,0,3);
        string chainID;
        if (sp[0]=="ATOM" || sp[0]=="HETATM") {
            string XX(line,21,1);
            chainID=strip(XX);
        } else {
            chainID=oldchainID;
        }
        //
        if (SS=="TER" || chainID.compare(oldchainID)!=0) {
            oldchainID=chainID;
            //
            if (chain_lines.size()!=0) {
                chain_class X(this,chain_lines);
                chains.push_back(X);
            }
            chain_lines.resize(0);
            if (SS!="TER") {
                string ID(line,0,6);
                if (ID=="ATOM  " || ID=="HETATM") {
                    chain_lines.push_back(line);
                }
            }
        }
        else {
            string ID(line,0,6);
            if (ID=="ATOM  " || ID=="HETATM") {
                chain_lines.push_back(line);
            }
        }
    }
    //
    // In case there are no TER records in the file
    //
    printf ("Pushing back last chain\n");
    if (chain_lines.size()!=0) {
        chain_class X(this,chain_lines);
        chains.push_back(X);
    }
    printf ("PDB file has been parsed\n");
    check_soup();
    return;
}

//
// --------------------------------------
//

void FFF::write_pqr(const std::string filename) {
  //
  // Write PQR file using the current assigned charges and radii 
  //
  char pdbfile[200];
  strcpy(pdbfile,filename.c_str());
  FILE* outfile;
  //
  // Does the dir exist?
  //
  vector<string> spath=splitpath(filename);
  if (not isdirectory(spath[0]) and spath[0]!="") {
    printf ("Directory for writing does not exist: '%s'\n",spath[0].c_str());
    exit(0);
  }
  outfile=fopen(pdbfile, "w");
  vector<string> pdblines=make_pdblines("PQR");
  for (unsigned int count=0;count<pdblines.size();count++) {
    fprintf(outfile,"%s",pdblines[count].c_str());
  }
  fclose(outfile);
  return;
}

void FFF::writepdb(const std::string filename) {
  write_pdb(filename);
  return;
}

void FFF::write_pdb(const std::string filename) {
  //
  // Write PDB file
  //
  char pdbfile[200];
  strcpy(pdbfile,filename.c_str());
  FILE* outfile;
  //
  // Does the dir exist?
  //
  vector<string> spath=splitpath(filename);
  if (not isdirectory(spath[0]) and spath[0]!="") {
    printf ("Directory for writing does not exist: '%s'\n",spath[0].c_str());
    exit(0);
  }
  outfile=fopen(pdbfile, "w");
  vector<string> pdblines=make_pdblines("PDB");
  for (unsigned int count=0;count<pdblines.size();count++) {
    fprintf(outfile,"%s",pdblines[count].c_str());
  }
  fclose(outfile);
  return;
}

vector<string> FFF::make_pdblines(string type) {
  //
  // Create a vector with all the lines
  //
  vector<string> pdblines;
  //
  // Loop over all atoms and residues and write
  //
  int atomcounter=0;
  for (unsigned int chain=0;chain<chains.size();chain++) {
    for (unsigned int resnum=0;resnum<chains[chain].residues.size();resnum++) {
      for (unsigned int atom=0;atom<chains[chain].residues[resnum].atoms.size();atom++) {
	char buffer [100];
	atomcounter++;
	sprintf(buffer,"ATOM  %4d",atomcounter);
	string s1(buffer);
	//
	// Special checks for the atomname
	//
	string atname(strip(chains[chain].residues[resnum].atoms[atom].name),0,1);
	if (atname=="1" || atname=="2" || atname=="3") {
	  sprintf(buffer,"  %-4s",strip(chains[chain].residues[resnum].atoms[atom].name).c_str());
	} else {
	  sprintf(buffer,"   %-3s",strip(chains[chain].residues[resnum].atoms[atom].name).c_str());
	}
	s1=s1+string(buffer);
	// Chain identifier is missing 
	// Residue name
	sprintf(buffer," %3s %1s",strip(chains[chain].residues[resnum].name).c_str(),
		chains[chain].name.c_str());
	s1=s1+string(buffer);
	//
	// Handle residue numbers > 999
	//
	if (strip(chains[chain].residues[resnum].pdbnum).size()==4) {
	  sprintf(buffer,"%4s ",strip(chains[chain].residues[resnum].pdbnum).c_str());
	} else {
	  sprintf(buffer," %3s ",strip(chains[chain].residues[resnum].pdbnum).c_str());
	}
	s1=s1+string(buffer);
	// Coordinates
	sprintf(buffer,"    %7.3f %7.3f %7.3f",
		chains[chain].residues[resnum].atoms[atom].x,
		chains[chain].residues[resnum].atoms[atom].y,
		chains[chain].residues[resnum].atoms[atom].z);
	s1=s1+string(buffer);
	//
	// Occupancy
	//
	if (type=="PDB") {
	  sprintf(buffer,"  %4.2f %5.2f\n",1.00,chains[chain].residues[resnum].atoms[atom].bfactor);
	} else if (type=="PQR") {
	  sprintf(buffer,"  %4.2f %5.2f\n",chains[chain].residues[resnum].atoms[atom].charge,
		  chains[chain].residues[resnum].atoms[atom].radius);
	} else {
	  printf ("Incorrect type for make_pdblines: %s\n",type.c_str());
	  exit(0);
	}
	s1=s1+string(buffer);
	pdblines.push_back(s1);
      }
    }
  }
  return pdblines;
}

//
// --------------------------------------
//

void FFF::check_soup() {
  //
  // Update some general arrays and check the soup
  //
    for (unsigned chain=0;chain<chains.size();chain++) {
        chains[chain].number=chain;
        chains[chain].update();
    }
  return;
}


void FFF::update_BOXLJ() {
    // Update the boxes for Lennard-Jones and Steric_Clash calculations
    BOXLJ = new Boxes(all_atoms,6.0);
    return;
}

void FFF::update_all_atoms() {
    //
    // Loop over all atoms and do lots of stuff
    //
    //printf ("in update_all atoms\n");
    all_atoms.resize(0);
    int count=0;
    unsigned int numchains=chains.size();
    for (unsigned int chain=0;chain<numchains;chain++) {
        //
        // Set the number of each chain
        //
        chains[chain].number=chain;
        //
        // Loop over all residues
        //
        unsigned int numres=chains[chain].residues.size();
        for (unsigned residue=0;residue<numres;residue++) {
	  // set the is_aa flag
	  chains[chain].residues[residue].is_aa=false;
	  for (unsigned int i=0;i<amino_acids.size();i++) {
	    if (chains[chain].residues[residue].name==amino_acids[i]) {
	      chains[chain].residues[residue].is_aa=true;
	    }
	  }
	  //
	  unsigned int numatoms=chains[chain].residues[residue].atoms.size();
            for (unsigned int atom=0;atom<numatoms;atom++) {
                //
                // Fill all_atoms
                //
                all_atoms.push_back(&chains[chain].residues[residue].atoms[atom]);
                //
                // Update the index
                //
                chains[chain].residues[residue].atoms[atom].index=count;
                count=count+1;
            }
        }
    }
    //printf ("all_atoms.size(): %d\n",static_cast<int>(all_atoms.size()));
    return ;
}

bool FFF::could_be_Hbonded(const atom_class& atom1, const atom_class& atom2) {
    //
    // Simple function for identifying a potential Hbond
    // No counting of valences or any angles
    //
    bool Hbonded=false;
    bool doac=false;
    if (atom1.is_donor && atom2.is_acceptor) {
        doac=true;
    } else if (atom2.is_donor && atom1.is_acceptor) {
        doac=true;
    }
    if (doac) {
        double distance=Dist(atom1,atom2);
        if (distance<3.5) {
            Hbonded=true;
        }
    }
    return Hbonded;
}

bool FFF::are_bonded(const atom_class& atom1, const atom_class& atom2) {
    // Are these two atoms bonded?
    bool bonded=false;
    //atom1.print_detail();
    //atom2.print_detail();
    //printf ("---");
    //printf ("Bond1\n");
    //printf ("cov bonds: %d %d\n",atom1.covalent_bonds.size(),
    //                            atom2.covalent_bonds.size());
    for (unsigned int j=0;j<atom1.covalent_bonds.size();j++) {
        //printf ("j; %d \n",j);
        //(*(atom1.covalent_bonds[j])).print();
         //(*(atom1.covalent_bonds[j])).print_detail();
        if ((*(atom1.covalent_bonds[j]))==atom2) {
            bonded=true;
            break;
        }
    }
    //printf ("Bond2\n");
    bool bonded2=false;
    for (unsigned int j=0;j<atom2.covalent_bonds.size();j++) {
        //printf ("j: %d size: %d\n",j,atom2.covalent_bonds.size());
        if ((*(atom2.covalent_bonds[j]))==atom1) {
            bonded2=true;
            break;
        }
    }
    //printf ("here\n");
    if (bonded!=bonded2) {
        printf ("Bonded one way but not the other\n");
        atom1.print_detail();
        atom2.print_detail();
        printf ("distance; %7.3f \n",Dist(atom1,atom2));
        exit(0);
    }
    return bonded;
}

bool FFF::are_1_3_bonded(const atom_class& atom1,const atom_class atom2) {
    // This function is wrong and should be rewritten
    bool bonded=false;
    for (unsigned int j=0;j<atom1.covalent_bonds.size();j++) {
        if (are_bonded(*(atom1.covalent_bonds[j]),atom2)) {
            bonded=true;
            break;
            }
    }
    return bonded;
}

bool FFF::are_1_4_bonded(const atom_class& atom1,const atom_class atom2) {
     // This function is wrong and should be rewritten
    bool bonded=false;
    for (unsigned int j=0;j<atom1.covalent_bonds.size();j++) {
        if (are_1_3_bonded(*(atom1.covalent_bonds[j]),atom2)) {
            bonded=true;
            break;
        }
    }
    return bonded;
}

//
// --------------------------------------
//

void FFF::soup_stat() {
  unsigned int numchains=chains.size();
  printf ("Number of chains: %4d \n",static_cast<int>(numchains));
  for (unsigned chain=0;chain<numchains;chain++) {
    unsigned int numres=chains[chain].residues.size();
    printf ("Chain: %4d, Number of residues: %d \n",chain,numres);
    for (unsigned residue=0;residue<numres;residue++) {
      unsigned int numatoms=chains[chain].residues[residue].atoms.size();
      printf ("Chain: %d, Residue: %4d, pdbnum: %4s, name: %4s, number of atoms: %d\n",
	      chains[chain].residues[residue].inchain,residue,chains[chain].residues[residue].pdbnum.c_str(),
	      chains[chain].residues[residue].name.c_str(),numatoms);
      for (unsigned int atom=0;atom<numatoms;atom++) {
	chains[chain].residues[residue].atoms[atom].print_detail();
      }
    }
  }
  return;
}
	
  
vector<int> FFF::find_residue(const std::string chainname,const std::string pdbresnumber) {
  vector<int> result;
    for (unsigned int chain=0;chain<chains.size();chain++) {
      if (strip(chainname)==strip(chains[chain].name)) {
            for (unsigned int residue=0;residue<chains[chain].residues.size();residue++) {
                if (strip(pdbresnumber)==strip(chains[chain].residues[residue].pdbnum)) {
		  result.push_back(chain);
		  result.push_back(residue);
		  return result;
                }
            }
        }
    }
    printf ("Could not find residue: %s %s\n",chainname.c_str(),pdbresnumber.c_str());
    return result;
}
