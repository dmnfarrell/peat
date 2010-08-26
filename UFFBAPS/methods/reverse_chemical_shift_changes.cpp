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


#include "reverse_chemical_shift_changes.h"


Reverse_Chemical_Shift_Changes::Reverse_Chemical_Shift_Changes(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();  


}

Reverse_Chemical_Shift_Changes::Reverse_Chemical_Shift_Changes(FILE * reporter,vector<Soup*> A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  cout<<"checkpoint 1"<<endl;
  set_parameters();
  cout<<"checkpoint 2"<<endl;
  calculate_dielectrics(A);
}


Reverse_Chemical_Shift_Changes::~Reverse_Chemical_Shift_Changes(){}



void Reverse_Chemical_Shift_Changes::calculate_dielectrics(vector<Soup*> A)
{
  for(vector<Soup*>::iterator soup = A.begin(); soup != A.end(); soup++)
    calculate_dielectrics(*soup);
}


void Reverse_Chemical_Shift_Changes::calculate_dielectrics(Soup * A)
{

  A->set_sybyl_atom_types();
  calculate_dielectrics(A->getAtoms());

}

void Reverse_Chemical_Shift_Changes::calculate_dielectrics(vector<Atom*> atoms_in)
{

  atoms = atoms_in;
  read_in_chemical_shifts(data_file);  
  
  FILE * mfile = fopen("dielectrics.txt","w");
  fprintf(mfile,"Residue, Residue number, Dielectric, Distance, electrostatic field, electrostatic potential\n");

  for(vector<amide>::iterator cur_amide = amides.begin(); cur_amide != amides.end(); cur_amide++ )
    {
 
     // charge, distance, vectors 
     (*cur_amide).dielectric = CSC.induced_shift(affector_atom, (*cur_amide).N,(*cur_amide).H, (*cur_amide).C); //units of e/(A*A)

     // convert to SI units
     (*cur_amide).dielectric = (*cur_amide).dielectric * e / (1e-10*1e-10); //units of C/(m*m)

     // constants
     (*cur_amide).dielectric = (*cur_amide).dielectric * alpha/(4*pi*eps_0*(*cur_amide).experimental_chemical_shift); //units: ppm/au /(C*C/(J*m)*ppm)*C/(m*m) = J /(C*m) * 1/au

     // a u
     (*cur_amide).dielectric = (*cur_amide).dielectric * e*a0/Eh; //units: ppm


     cout<<"Distance: "<<(*cur_amide).N->residue_no<<(*cur_amide).N->residue<< "  "<<D.calculate(affector_atom, (*cur_amide).N, true) <<" to ";
     (*cur_amide).N->display();
     cout<<"Distance: "<< (*cur_amide).H->residue_no<< (*cur_amide).H->residue<< "  "<<D.calculate(affector_atom, (*cur_amide).H, true) <<" to ";
     (*cur_amide).H->display();
     cout<<"Affector atom: ";
     affector_atom->display();

     float distance = sqrt((((*cur_amide).N->x+(*cur_amide).H->x)/2 - affector_atom->x)*(((*cur_amide).N->x+(*cur_amide).H->x)/2 - affector_atom->x)+
			   (((*cur_amide).N->y+(*cur_amide).H->y)/2 - affector_atom->y)*(((*cur_amide).N->y+(*cur_amide).H->y)/2 - affector_atom->y)+
			   (((*cur_amide).N->z+(*cur_amide).H->z)/2 - affector_atom->z)*(((*cur_amide).N->z+(*cur_amide).H->z)/2 - affector_atom->z));

     // electric field
     (*cur_amide).electrostatic_field = 1/(4*pi*eps_0*(*cur_amide).dielectric)*affector_atom->charge*e/(distance*1e-10*distance*1e-10);

     // electrostatic potential
     (*cur_amide).electrostatic_potential = 1/(4*pi*eps_0*(*cur_amide).dielectric)*affector_atom->charge*e/(distance*1e-10);
     
     cout<<(*cur_amide).dielectric<<"   dist:"<< distance<<endl;
     fprintf(mfile,"%3s %4d %10.2e %10.2e %10.2e %10.2e\n",(*cur_amide).H->residue.c_str(),(*cur_amide).H->residue_no,(*cur_amide).dielectric, distance, (*cur_amide).electrostatic_field, (*cur_amide).electrostatic_potential);
    }
 
  fclose(mfile);
}


void Reverse_Chemical_Shift_Changes::read_in_chemical_shifts(string filename)
{
  ifstream list(filename.c_str());
  cout<<"hallllos"<<endl;
  if(! list.is_open())
    {
      cout << "Error:File '"<<filename<<"' could not be found\n";
      exit(0);
    }

  string line,affector_name;
  vector<string> words;

  //set affector atom
  getline(list, line);
  set_affector_atom(line);

  cout<<"spasser"<<endl;
  while (!list.eof())
    {
      // get line and split it into words
      getline(list, line);
      /*** replaces ',' with ' ' ***/
      for(int i=0; i<line.size(); i++)
	if(string(line,i,1) ==",")
	  {
	    line.replace(i, 1, string(" "));
	    i++;
	  }
      words.clear();
      char temp[300];
      istringstream stream(line.c_str());
      while(stream>>temp)
	{
	  words.push_back(string(temp));
	}

      // insert value in amide chem shift list
      if(words.size()>1)
	{
	  
	  // check if last word can be converted to float
	  istringstream i(words[words.size()-1]);
	  float chem_shift;
	  if (i >> chem_shift)
	    if(words.size()==2) //Only backbone amides!!!
	    {
	      int res_no = atoi(words[0].substr(1,words[0].size()-1).c_str());
	      insert_amide(res_no, chem_shift);

	    }
	}
    

    }
}


void Reverse_Chemical_Shift_Changes::insert_amide(int residue, float chem_shift)
{
  amide temp;

  cout<<chem_shift<<" for residue "<<residue<<endl;

  // set atoms
  for(vector<Atom*>::iterator atom = atoms.begin();atom!=atoms.end();atom++)
    {
      if((*atom)->residue_no == residue)
	{
	  if((*atom)->name == "N")
	    temp.N = *atom;
	  if((*atom)->name == "H")
	    temp.H = *atom;
   	}
      if((*atom)->residue_no == residue-1)
	if((*atom)->name == "C")
	  temp.C = *atom;

    }
  
  // set chem shift
  temp.experimental_chemical_shift = chem_shift; 

  cout<<"Inserted residue: "<<residue<<" with the chemical shift "<<chem_shift<<endl;
  temp.N->display();
  temp.H->display();
  temp.C->display();

  // store
  amides.push_back(temp);

}


void Reverse_Chemical_Shift_Changes::set_affector_atom(string residue)
{

  cout<<"looking for "<<residue<<endl;
  cout<<"split "<<residue.substr(10,4)<<endl;
  int residue_no = atoi(residue.substr(10,4).c_str());
  cout<<"residue no "<<residue_no<<endl;
  // find affector atom
  for(vector<Atom*>::iterator atom = atoms.begin(); atom != atoms.end(); atom++ )
    if((*atom)->residue_no == residue_no && CSC.titratable_residues.count((*atom)->generic_key)==1)
      {
	affector_atom = *atom;
	affector_atom_charge = CSC.titratable_residues[(*atom)->generic_key];
      }

  // set the charge
  affector_atom->charge = CSC.titratable_residues[affector_atom->generic_key];

 

  cout<<"Affector "<<residue_no<<" ";
  affector_atom->display();
}

void Reverse_Chemical_Shift_Changes::set_parameters()
{
  //  alpha = 70.0; //ppm/au
  alpha = 977.0; //A|| for C-N 
  pi = 3.14159265;
  eps_0 = 8.85419e-12;//C*C/(J*m)
  e = 1.602177e-19; //C
  a0 = 5.291772e-11; //m
  Eh = 4.359744e-18; //J



  

  ifstream parameters;
  parameters.open("reverse_chemical_shift_changes.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'reverse_chemical_shift_changes.cfg' could not be found\n";
      exit(0);
    }

  printf("Setting parameters\n");

  string dummy;

  while (!parameters.eof())
    {

      parameters >> dummy;

      if(dummy == "EXPERIMENTAL_DATA_FILE:")
	parameters >> data_file;

    } //end read file

  //printing all configurations

  printf("Configurations for reverse chemical shift changes:\n");
  printf("\n");
  printf("\n");

  printf("EXPERIMENTAL_DATA_FILE:             %s\n", data_file.c_str());
  printf("\n");

  printf("Version: $Id: es.cpp 160 2007-11-26 19:44:53Z chresten $\n");
  printf("\n");

}

void Reverse_Chemical_Shift_Changes::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: reverse_chemical_shift_changes.cpp 18 2005-11-15 09:56:08Z chresten $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
