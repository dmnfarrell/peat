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


#include "prepinreader.h"


PrepinReader::PrepinReader(FILE * reporter, vector<Soup *> A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  go(A);
}

PrepinReader::PrepinReader():Method(){}

PrepinReader::PrepinReader(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

}

PrepinReader::~PrepinReader(){}


void PrepinReader::go(vector<Soup *> A)
{
  for(unsigned int s=0; s<A.size(); s++)
    {
      string name = A[s]->name;
      int no_pc =  A[s]->getProteinChains().size();
      int no_mo =  A[s]->getMolecules().size();
      int no_wa =  A[s]->getWaters().size();
      
      printf("Preparing parameters for soup '%s'\n", name.c_str());
      /*      
      printf("Soup '%s' has:\n", name.c_str());
      printf("  '%d' protein chains\n", no_pc);
      printf("  '%d' molecules\n", no_mo);
      printf("  '%d' waters\n", no_wa);
      */

      //making parameters for all objects in soup
      for(int pc=0; pc<no_pc; pc++)
	do_proteinChain(A[s]->getProteinChains()[pc]);
      for(int mo=0; mo<no_mo; mo++)
	do_molecule(A[s]->getMolecules()[mo], A[s]->name);
      for(int wa=0; wa<no_wa; wa++)
	do_water(A[s]->getWaters()[wa]);
      
      A[s]->assigned_parameters = "AMBER";
      /*
      printf("Soup '%s' has the atoms:\n", name.c_str());
      
      vector<Atom*> atoms = A[s]->getAtoms();
      
      for(unsigned int a=0; a<atoms.size(); a++)
	printf("%4s %2s %3s %2s %4.3f %4.3f %4.3f\n",
	       atoms[a]->name.c_str(),
	       atoms[a]->element.c_str(),
	       atoms[a]->residue.c_str(),
	       atoms[a]->adv_type.c_str(),
	       atoms[a]->vdw_dist,
	       atoms[a]->vdw_energy,
	       atoms[a]->charge);
      */
    }
}




void PrepinReader::find_protein_chain_charges(vector<ProteinChain*> P)
{
  //Used by generalparameters to set protein charges

  ProteinChain * temp;

  for(unsigned int pi=0; pi<P.size(); pi++)
    {
      temp = P[pi];
      /*** N-terminus ***/
      read_prep("all_aminont94.in");
      compare(temp->getResidues()[0]->getAtoms());

      /*** everything between ***/
      read_prep("all_amino94.in");
      for(unsigned int i=1; i<temp->getResidues().size()-1; i++)
	compare(temp->getResidues()[i]->getAtoms());

      /*** C-terminus ***/
      read_prep("all_aminoct94.in");
      compare(temp->getResidues()[temp->getResidues().size()-1]->getAtoms());
    }
}


void PrepinReader::do_proteinChain(ProteinChain * p)
{

  /*** N-terminus ***/
  read_prep("all_aminont94.in");
  compare(p->getResidues()[0]->getAtoms());


  /*** everything between ***/
  read_prep("all_amino94.in");
  for(unsigned int i=1; i<p->getResidues().size()-1; i++)
    compare(p->getResidues()[i]->getAtoms());

  /*** C-terminus ***/
  read_prep("all_aminoct94.in");
  compare(p->getResidues()[p->getResidues().size()-1]->getAtoms());

  /*** setting parameters ***/
  read_parm("parm94.dat");
  set_parameters(p->getAtoms());
} 

void PrepinReader::do_molecule(Molecule * m, string soup_name)
{
  string guess = soup_name+"and-bcc.prepi";
  printf("Guessing that the name of the prepin file is: %s\n",guess.c_str());

  read_prep(guess);
  compare(m->getAtoms());

  read_parm("gaff.dat");
  set_parameters(m->getAtoms());
}

void PrepinReader::do_water(Water * w)
{
  read_prep("opls_uni.in");
  compare(w->getAtoms());
  
  read_parm("parm94.dat");
  set_parameters(w->getAtoms());
}


void PrepinReader::read_prep(string fn)
{
  /*** clear maps ***/
  advtype.clear();
  pc.clear();
  keys.clear();

  /*** read file ***/
  //printf("Reading file '%s'\n", fn.c_str());
  string filename="parameters/"+fn;
  string line;
  file.open(filename.c_str());

  if(!file) 
    {
      printf("WARNING: could not read parameter file '%s'!",fn.c_str());
      return;
    }  
  bool go = false;

  while(!file.eof())
    {

      getline(file, line);

      /*** break line into words ***/
      words.clear();
      char temp[300];
      istringstream stream(line.c_str());
      while(stream>>temp)
	words.push_back(removeWhite(string(temp)));

      /*** new residue ***/
      if(words[1] == "INT")
	{
	  go = true;
	  amino = words[0];
	}
      if(words.size() == 0)
	go = false;
      
      if(go == true)
	 if(words.size() > 4 && words[1] != "DUMM")
	   {
	     string key1;
	     key1 = words[1];
	     
	     /*** to avoid confusion when reading molecule prepins***/
	     if(amino=="SUB")
	       amino="";

	     string collectKey = amino + " " + key1;
	  
	     keys.push_back(collectKey);
	     
	     string tmpAdvType;
	     tmpAdvType = words[2];
	     advtype[collectKey]=tmpAdvType;
	     
	     double tmpPartialCharge;
	     istringstream instr(words[10]);
	     instr >>  tmpPartialCharge;
	     pc[collectKey]=tmpPartialCharge;
	     
	     //	     cout << "Found key:" << collectKey << " " << tmpAdvType << " " << tmpPartialCharge << endl;

	   }
	  
	
    }
  file.close();
}

void PrepinReader::read_parm(string fn)
{
  nonbondKeys.clear();
  nonbondEnergy.clear();
  nonbondR.clear();

  printf("Reading parm file '%s'\n", fn.c_str());
  string filename="parameters/"+fn;
  string line,tag;

  file.open(filename.c_str());

  while(!file.eof())
    {
      getline(file, line);
      tag = string(line, 0, 4);
      if(tag =="MOD4")
	{
	  while(tag !="END")
	    {
	      getline(file, line);
	      tag = string(line, 0, 3);
	      if(line.size() > 5)
		{
		  string name;
		  double energy,r;

		  istringstream stream(line);

		  stream >> name;

		  nonbondKeys.push_back(name);
		  stream >> r;
		  stream >> energy;

		  nonbondEnergy[name]=energy;
		  nonbondR[name]=r;

		  //		  cout <<"Found parms :"<<name<<" r:"<<nonbondR[name]<<" e:"<<nonbondEnergy[name]<<endl;
		}
	    }
	}
    }
  file.close();
}


void PrepinReader::compare(vector<Atom *> atoms)
{
  /*  if(atoms.size()> 0)
    cout << "Comparing residue: "<<atoms[0]->residue<<endl;
  else
    cout << "WARNING: Comparing empty residue"<<endl;
  */
  if(atoms.size() == 0)
    cout << "WARNING: Comparing empty residue"<<endl;

  //Setting protonations states
  set_protonation_state_HIS(atoms);
  set_protonation_state_LYS(atoms);
  set_protonation_state_GLU(atoms);
  set_protonation_state_ASP(atoms);
  set_protonation_state_CYS(atoms);


  for(unsigned int i = 0; i< atoms.size();i++)
    {
      //correct H names
      //      atoms[i]->name = correct_hydrogen_name(atoms[i]->name);

      //correct WHATIF names
      //      atoms[i]->name = correct_WI_names(atoms[i]->name, atoms[i]->residue);

      //correct 1-3 errors
      //      atoms[i]->name = correct_1_3_error(atoms[i]->name, atoms[i]->residue);

      string searchKey = make_search_key(atoms[i]->residue, atoms[i]->name);

      if(look_up_key(searchKey))
	{
	  atoms[i]->adv_type = advtype[searchKey];
	  atoms[i]->charge = pc[searchKey];
	  
	  //	  cout <<"adv, pc set: "<<atoms[i]->adv_type<<" "<<atoms[i]->charge<<endl;
	}
      else
	{
	  /*** gap ? ***/
	  //  if(atoms[i]->name != "OXT")
	  //  if(atoms[i]->name != "H2")
	  //    if(atoms[i]->name != "H3")
	  
	  cout << "WARNING: atom NOT FOUND in parameter file. searchKey:  "<< searchKey << endl;
	  
	  atoms[i]->adv_type = "TEST";
	  atoms[i]->charge = 0;
	  
	}
    }
  
}

void PrepinReader::set_parameters(vector<Atom *> atoms)
{
  vector<string> C_sp2_redundant, N_sp2_redundant;

  C_sp2_redundant.push_back("C*");
  C_sp2_redundant.push_back("CB");
  C_sp2_redundant.push_back("CN");
  C_sp2_redundant.push_back("CW");
  C_sp2_redundant.push_back("CC");
  C_sp2_redundant.push_back("CR");
  C_sp2_redundant.push_back("CV");
  C_sp2_redundant.push_back("CK");
  C_sp2_redundant.push_back("CQ");

  N_sp2_redundant.push_back("N*");
  N_sp2_redundant.push_back("N2");
  N_sp2_redundant.push_back("NB");
  N_sp2_redundant.push_back("NA");
  N_sp2_redundant.push_back("NC");

  
  for(unsigned int i = 0; i< atoms.size();i++)
    {
      /*** checking for redundant C_sp2 type***/
      if(atoms[i]->element == "C")
	for(unsigned int c=0; c<C_sp2_redundant.size(); c++)
	  if(C_sp2_redundant[c] == atoms[i]->adv_type)
	    atoms[i]->adv_type = "C";

      /*** checking for redundant N_sp2 type***/
      if(atoms[i]->element == "N")
	for(unsigned int c=0; c<N_sp2_redundant.size(); c++)
	  if(N_sp2_redundant[c] == atoms[i]->adv_type)
	    atoms[i]->adv_type = "N";


      bool found = look_up_non_bound_key(atoms[i]->adv_type);
      
      if(found)
	{
	  atoms[i]->vdw_dist = nonbondR[atoms[i]->adv_type];
	  atoms[i]->vdw_energy = nonbondEnergy[atoms[i]->adv_type];
	  //cout <<"nonR, nonE set: "<<atoms[i]->vdw_dist<<" "<<atoms[i]->vdw_energy<<endl;
	}
      else if(atoms[i]->adv_type == "TEST") //ignored atom
	{
	  atoms[i]->vdw_dist = 0;
	  atoms[i]->vdw_energy = 0;
	}
      else
	{
	  atoms[i]->vdw_dist = 0;
	  atoms[i]->vdw_energy = 0;
	  cout<<"WARNING: No vdw parameters found for advanced type '"<<atoms[i]->adv_type<<"'"<<endl;
	}

    }

}


string PrepinReader::make_search_key(string residue, string name)
{

  string searchKey = residue + " " + name;
  return searchKey;
}

bool PrepinReader::look_up_key(string key)
{
  bool found = false;
  for(unsigned int j=0; j< keys.size() ; j++)
    {
      if(keys[j] == key)
	found = true;
    }

  return found;
}

bool PrepinReader::look_up_non_bound_key(string key)
{
  bool found = false;
  for(unsigned int j=0; j< nonbondKeys.size() ; j++)
    {
      if(nonbondKeys[j] == key)
	found = true;
    }

  return found;
}

void PrepinReader::set_protonation_state_HIS(vector<Atom *> atoms)
{
  //
  //Sorting HIS residues into the HIP, HIE and HID catagories
  //
  if(atoms.size() > 0)
    if(atoms[0]->residue == "HIS")
      {
	int no_H_delta = 0;
	int no_H_epsilon = 0;

	for(unsigned int i = 0; i< atoms.size();i++)
	  { 
	    //correcting atom name before doing comparison
	    //	    atoms[i]->name = correct_hydrogen_name(atoms[i]->name);

	    string comp = string(atoms[i]->name,0,2);
	    if(comp == "HD")
	      no_H_delta++;
	    else if(comp == "HE")
	      no_H_epsilon++;
	  }
	string new_res = "HIS";
	if(no_H_delta == 2 && no_H_epsilon < 2)
	  new_res = "HID";
	else if(no_H_delta < 2 && no_H_epsilon == 2)
	  new_res = "HIE";
	else if(no_H_delta == 2 && no_H_epsilon == 2)
	  new_res = "HIP";
	else
	  {
	    cout <<"WARNING: Failed to define protonation state for HIS residue, using HID"<<endl;
	    new_res = "HID";
	  }
	//cout<<"HIS residue set to be "<<new_res<<endl;
	for(unsigned int i = 0; i< atoms.size();i++)
	  atoms[i]->residue = new_res;
      }

}

void PrepinReader::set_protonation_state_LYS(vector<Atom *> atoms)
{
  //
  //Sorting LYS residues into the LYS and LYN catagories
  //
  if(atoms.size() > 0)
    if(atoms[0]->residue == "LYS")
      {
	int no_H_zeta = 0;

	for(unsigned int i = 0; i< atoms.size();i++)
	  { 
	    //correcting atom name before doing comparison
	    //	    atoms[i]->name = correct_hydrogen_name(atoms[i]->name);

	    string comp = string(atoms[i]->name,0,2);
	    if(comp == "HZ")
	      no_H_zeta++;
	  }
	string new_res = "LYS";
	if(no_H_zeta == 2)
	  new_res = "LYN";
	else if(no_H_zeta == 3)
	  new_res = "LYS";
	else
	  cout <<"WARNING: Failed to define protonation state for LYS residue, using LYS"<<endl;
	for(unsigned int i = 0; i< atoms.size();i++)
	  atoms[i]->residue = new_res;
      }


}
void PrepinReader::set_protonation_state_GLU(vector<Atom *> atoms)
{
  //
  //Sorting GLU residues into the GLU and GLH catagories
  //
  if(atoms.size() > 0)
    if(atoms[0]->residue == "GLU")
      {
	int no_H_epsilon = 0;

	for(unsigned int i = 0; i< atoms.size();i++)
	  { 
	    //correcting atom name before doing comparison
	    //	    atoms[i]->name = correct_hydrogen_name(atoms[i]->name);

	    string comp = string(atoms[i]->name,0,2);
	    if(comp == "HE")
	      no_H_epsilon++;
	  }
	string new_res = "GLU";
	if(no_H_epsilon == 1)
	  new_res = "GLH";
	else if(no_H_epsilon == 0)
	  new_res = "GLU";
	else
	  cout <<"WARNING: Failed to define protonation state for GLU residue, using GLU"<<endl;
	for(unsigned int i = 0; i< atoms.size();i++)
	  atoms[i]->residue = new_res;
      }


}

void PrepinReader::set_protonation_state_ASP(vector<Atom *> atoms)
{
  //
  //Sorting GLU residues into the ASP and ASH catagories
  //
  if(atoms.size() > 0)
    if(atoms[0]->residue == "ASP")
      {
	int no_H_delta = 0;

	for(unsigned int i = 0; i< atoms.size();i++)
	  { 
	    //correcting atom name before doing comparison
	    //	    atoms[i]->name = correct_hydrogen_name(atoms[i]->name);

	    string comp = string(atoms[i]->name,0,2);
	    if(comp == "HD")
	      no_H_delta++;
	  }
	string new_res = "ASP";
	if(no_H_delta == 1)
	  new_res = "ASH";
	else if(no_H_delta == 0)
	  new_res = "ASP";
	else
	  cout <<"WARNING: Failed to define protonation state for ASP residue, using ASP"<<endl;
	for(unsigned int i = 0; i< atoms.size();i++)
	  atoms[i]->residue = new_res;
      }


}

void PrepinReader::set_protonation_state_CYS(vector<Atom *> atoms)
{
  //
  //Sorting CYS residues into the CYS, CYM and CYX catagories
  //
  if(atoms.size() > 0)
    if(atoms[0]->residue == "CYS")
      {
	int no_H_gamma = 0;

	for(unsigned int i = 0; i< atoms.size();i++)
	  { 
	    //correcting atom name before doing comparison
	    //	    atoms[i]->name = correct_hydrogen_name(atoms[i]->name);

	    string comp = string(atoms[i]->name,0,2);
	    if(comp == "HG")
	      no_H_gamma++;
	  }
	string new_res = "CYS";
	if(no_H_gamma == 1)
	  new_res = "CYS";
	else if(no_H_gamma == 0)
	  new_res = "CYX";
	//	  new_res = "CYM"; currently skipped due to missing charges in all_amino94.in
	else
	  cout <<"WARNING: Failed to define protonation state for GLU residue, using CYS"<<endl;
	for(unsigned int i = 0; i< atoms.size();i++)
	  atoms[i]->residue = new_res;
      }


}


string PrepinReader::correct_hydrogen_name(string name)
{
  //
  // Fixing hydrogen naming
  //
  string new_name = name;
  if(isdigit(name[0]) && name.c_str()[1] == 'H')
    {
      new_name = string(name,1,name.size()-1);
      new_name.append(name,0,1);
    }

  return new_name;
}

string PrepinReader::correct_1_3_error(string name, string residue)
{

  vector<string> has_HD1;
  has_HD1.push_back("HID");
  has_HD1.push_back("HIP");
  has_HD1.push_back("HIS");
  has_HD1.push_back("TRP");
  has_HD1.push_back("PHE");
  has_HD1.push_back("TYR");

  vector<string> has_HE1;
  has_HE1.push_back("HID");
  has_HE1.push_back("HIE");
  has_HE1.push_back("HIP");
  has_HE1.push_back("HIS");
  has_HE1.push_back("TRP");
  has_HE1.push_back("PHE");
  has_HE1.push_back("TYR");
  has_HE1.push_back("MET");

  if(name == "HA1")
    name = "HA3";
  if(name == "HB1" && residue != "ALA")
    name = "HB3";
  if(name == "HG1" && residue != "THR")
    name = "HG3";
  if(name == "HD1" && find(has_HD1.begin(), has_HD1.end(), residue) == has_HD1.end())
    name = "HD3";
  if(name == "HE1" && find(has_HE1.begin(), has_HE1.end(), residue) == has_HE1.end())
    name = "HE3";

  if(name == "HG11" && residue != "VAL")
    name = "HG13";

  return name;
}

string PrepinReader::correct_WI_names(string name, string residue)
{
  if(name == "O''")
    name = "OXT";
  if(name == "O'")
    name = "O";
  if(name == "H1")
    if(!look_up_key(make_search_key(residue, name)))
      {
	cout << "Assuming gap in file - Changing H1 to H of residue "<<residue<<", will ignore OXT, H2 and H3 "<<endl;
	name = "H";
      }


  return name;
}



void PrepinReader::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: prepinreader.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("Reads AMBER atom names and charges from AMBER prepin files");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
