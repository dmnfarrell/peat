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

/*
$Id: desolvation.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "desolvation.h"


Desolvation::Desolvation(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();
}

Desolvation::Desolvation():Method()
{
  set_parameters();
}


Desolvation::Desolvation(FILE * resultsFile,FILE * reporter, vector<Soup*> A, vector<Soup*> B):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();

  float result = 0;
  for(unsigned int i=0; i<A.size(); i++)
    {
      
      result = calculate(A[i], B[i]);
      
      fprintf(resultsFile,"%-4s %8f\n",
 	      A[i]->name.substr(0,4).c_str(),
	      result);
    }  
}
 
Desolvation::Desolvation(FILE * resultsFile,FILE * reporter, vector<Soup*> A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();

  float result = 0;
  for(unsigned int i=0; i<A.size(); i++)
    {
      result = calculate(A[i]);
      
      fprintf(resultsFile,"%-4s %8f\n",
 	      A[i]->name.substr(0,4).c_str(),
	      result);
    }  
}

Desolvation::~Desolvation(){
//fclose(out);
}

float Desolvation::calculate(Soup * A)
{
  cout <<"preparing a"<<endl<<flush;
  prepare_soup(A);
  cout <<"prepare end a"<<endl<<flush;

  // set desolvation_factors
  set_desolvation_factors(A);

  return sum_atoms(A->get_dry_atoms());
}

float Desolvation::calculate(Soup * A, Soup * B) // calculate desolvation for protein-ligand complexes 
{

  prepare_soup(A);
  prepare_soup(B);

  // set desolvation_factors
  set_desolvation_factors(A,B);


  vector<Atom*> all = A->get_dry_atoms();
  vector<Atom*> atomsB = B->get_dry_atoms();

  all.insert (all.end(), atomsB.begin(), atomsB.end()); // append atomsB

  return sum_atoms(all);
}


float Desolvation::calculate_contact_between_A_and_B(Soup * A, Soup * B) // calculate desolvation of A when also taking B into account
// A must be a ligand 
{

  prepare_soup(A);
  prepare_soup(B);

  /*** calculate bonds ***/
  if(!A->bonds_assigned)
    BM.assign_bonds(A);
  if(!B->bonds_assigned)
    BM.assign_bonds(B);

  vector<Atom*> atoms = A->get_dry_atoms();
  vector<Atom*> atomsB = B->get_dry_atoms();
  atoms.insert(atoms.end(), atomsB.begin(), atomsB.end()); // append atomsB

  set_desolvation_factors(atoms, true);

  vector<Atom*> all = A->get_dry_atoms();
  float res = 0;
  for(vector<Atom*>::iterator it = all.begin(); it < all.end(); it++)
    if((*it)->element != "H")
      res += (*it)->desolvation_factor;

  return res;
}



void Desolvation::prepare_soup(Soup * A)
{

  // find fragmental volumes - used to calculate occupancies
  FV.set_precalculated_fragmental_volumes(A);

  // set sybyl types
  A->set_sybyl_atom_types();


  return;
}

float Desolvation::sum_atoms(vector<Atom*> all)
{
  // Just sums up enerG for all atoms in A

  float res = 0;
  
  for(vector<Atom*>::iterator it = all.begin(); it < all.end(); it++)
    if((*it)->element != "H")
      {
	if(desolvation_parameters.find((*it)->sybyl_type) != desolvation_parameters.end())
	  res += atomic_desolvation((*it));
	//	else
	//	  if ((*it)->residue != "HOH")
	//	  1  printf("WARNING: Could not find desolvation parameter for %s-%s\n",(*it)->residue.c_str(),(*it)->name.c_str());
      }

  //  printf("Total desolvation for %s: %f\n",A->name.c_str(),res);
  return res;


}



float Desolvation::atomic_desolvation(Atom * atom)
{

  // desolvation pr atom
  float res = desolvation_parameters[atom->sybyl_type];
  // volume of group
  // res = res * atom->fragmental_volume;

  // desolvation_factor of atom
  res = res * atom->desolvation_factor;
  //    (atom->desolvation_factor - get_min_desolvation_factor(atom->generic_key.replace(3,1,"_")))/
  //  (get_max_desolvation_factor(atom->generic_key.replace(3,1,"_")) - get_min_desolvation_factor(atom->generic_key.replace(3,1,"_")));
  /*
  printf("Desolvation energy:%8.2f occFactor%8.2f (min:%8.2f, max:%8.2f, this:%8.2f) res:%8.2f for ",
	 desolvation_parameters[atom->sybyl_type],
	 (atom->desolvation_factor - get_min_desolvation_factor(atom->generic_key.replace(3,1,"_")))/
	 (get_max_desolvation_factor(atom->generic_key.replace(3,1,"_")) - get_min_desolvation_factor(atom->generic_key.replace(3,1,"_"))),
	 get_min_desolvation_factor(atom->generic_key.replace(3,1,"_")),
	 get_max_desolvation_factor(atom->generic_key.replace(3,1,"_")),
	 atom->desolvation_factor,
	 res);
  atom->display();
  */
  return res;
}


float Desolvation::get_min_desolvation_factor(string key)
{
  string::size_type loc = key.find( "-", 0 );
  if( loc != string::npos ) 
    key.replace(loc,1,"_");
  

  float res = avr_min_desolvation_factor;
  if (min_occupancies.find(key) != min_occupancies.end())
    if (min_occupancies[key] != max_occupancies[key]) // to avoid cases where only one value was found (CYS OXT)
      res = min_occupancies[key];

  return res;
}

float Desolvation::get_max_desolvation_factor(string key)
{
  string::size_type loc = key.find( "-", 0 );
  if( loc != string::npos ) 
    key.replace(loc,1,"_");
  
  float res = avr_max_desolvation_factor;
  if (max_occupancies.find(key) != max_occupancies.end())
    if (min_occupancies[key] != max_occupancies[key]) // to avoid cases where only one value was found (CYS OXT)
      res = max_occupancies[key];

  return res;
}


void Desolvation::set_desolvation_factors(Soup * A)
{
  /*** calculate bonds ***/
  if(!A->bonds_assigned)
    BM.assign_bonds(A);

  vector<Atom*> atoms = A->get_dry_atoms();
  set_desolvation_factors(atoms, false);

  return;
}

void Desolvation::set_desolvation_factors(Soup * A, Soup * B)
{

  /*** calculate bonds ***/
  if(!A->bonds_assigned)
    BM.assign_bonds(A);
  if(!B->bonds_assigned)
    BM.assign_bonds(B);

  vector<Atom*> atoms = A->get_dry_atoms();
  vector<Atom*> atomsB = B->get_dry_atoms();
  atoms.insert(atoms.end(), atomsB.begin(), atomsB.end()); // append atomsB

  set_desolvation_factors(atoms, false);

  return;
}



void Desolvation::set_desolvation_factors(vector<Atom*> atoms, bool include_hydrophilic) {
  // find hydrophobic atoms
  HP.find_hydrophobic_atoms(atoms);
  desolvation_distance = 6.0*6.0;
  float distance, desolvation_factor, sigma = 3.5*3.5, factor;

  //vector<Atom*>::iterator atom1 = atoms.begin(), atom2;
  Atom* jatom2;
  cout << "DE begin"<<flush;
  Boxes BOX(atoms,sqrt(desolvation_distance));
  vector<int> close_atoms;
  //
  for (unsigned int atom1=0;atom1<atoms.size();atom1++) {
    atoms[atom1]->desolvation_factor=0.0; // Init
    if (atoms[atom1]->element!="H") {
      close_atoms=BOX.get_close_atoms(atoms[atom1]);
      for (unsigned int jcount=0;jcount<close_atoms.size();jcount++) {
	int atom2=close_atoms[jcount];
	if (atoms[atom2]->element!="H" && atom1!=atom2) {
	  if (include_hydrophilic or atoms[atom2]->is_hydrophobic) {
	    if (atoms[atom1]->residue_no != atoms[atom2]->residue_no) {
	      //
	      // Make sure the atoms are not bonded
	      //
	      bool bonded=false;
	      for (unsigned int kk=0;kk<atoms[atom1]->bonded_atoms.size();kk++) {
		if (atoms[atom2]==atoms[atom1]->bonded_atoms[kk]) {
		  bonded=true;
		}
	      }
	      if (not bonded) {
		//
		// Calculate the distance
		//
		distance=
		  ((atoms[atom1])->x - (atoms[atom2])->x)*((atoms[atom1])->x - (atoms[atom2])->x)+
		  ((atoms[atom1])->y - (atoms[atom2])->y)*((atoms[atom1])->y - (atoms[atom2])->y)+
		  ((atoms[atom1])->z - (atoms[atom2])->z)*((atoms[atom1])->z - (atoms[atom2])->z);
		if (distance < desolvation_distance) {
		  factor=1.0;
		  if (abs(atoms[atom1]->residue_no - atoms[atom2]->residue_no) ==1) {
		    factor=0.5;
		  }
		  atoms[atom1]->desolvation_factor += factor*(atoms[atom2]->fragmental_volume)*exp(-distance/(2*sigma));
		}
	      }
	    }
	  }
	}
      }
      atoms[atom1]->desolvation_factor = (atoms[atom1]->desolvation_factor - get_min_desolvation_factor(atoms[atom1]->generic_key))/
	(get_max_desolvation_factor((atoms[atom1])->generic_key) - get_min_desolvation_factor((atoms[atom1])->generic_key));
      
      atoms[atom1]->tempFactor = atoms[atom1]->desolvation_factor; //to make nice pictures in pymol	  
    }
    else {
      atoms[atom1]->tempFactor = 0.0; //to make nice pictures in pymol	
    }
  }

 //  while(atom1 != atoms.end()) {
//       if((*atom1)->element != "H")
// 	{
// 	  (*atom1)->desolvation_factor = 0.0;
// 	  // Get the close atoms
// 	  close_atoms=BOX.get_close_atoms((*atom1));
// 	  for (unsigned int jcount=0;jcount<close_atoms.size();jcount++) {
// 	    //atom2 = atoms.begin();
// 	    //while(atom2 != atoms.end()) {
// 	    int atom2_number=close_atoms[jcount];
// 	    atom2=&(atoms[atom2_number]);
// 	      if(include_hydrophilic or (*atom2)->is_hydrophobic) // Only hydrophobic atoms contribute to desolvation
// 		if((*atom2)->element != "H" && atom1 != atom2) //disregard hydrogens and self
// 		  if((*atom1)->residue_no != (*atom2)->residue_no) //disregard same residue atoms
// 		    if(find((*atom1)->bonded_atoms.begin(),
// 			    (*atom1)->bonded_atoms.end(),
// 			    (*atom2)) == (*atom1)->bonded_atoms.end()) //make sure atoms are not bonded
// 		      {
// 			distance = 
// 			  ((*atom1)->x - (*atom2)->x)*((*atom1)->x - (*atom2)->x)+
// 			  ((*atom1)->y - (*atom2)->y)*((*atom1)->y - (*atom2)->y)+
// 			  ((*atom1)->z - (*atom2)->z)*((*atom1)->z - (*atom2)->z); 
		      
// 			if(distance < desolvation_distance)
// 			  {
// 			    factor = 1.0;
// 			    if( abs((*atom1)->residue_no - (*atom2)->residue_no) == 1 )
// 			      factor = 0.5;

// 			    (*atom1)->desolvation_factor += factor*(*atom2)->fragmental_volume*exp(-distance/(2*sigma));
			    
// 			  }
			
// 		      }
// 	      //atom2++;
// 	  }
	  

// 	  //////////////////////
// 	  /*

// 	  printf("Desolvation_Factor: %8.2f for ",(*atom1)->desolvation_factor);
// 	  (*atom1)->display();

// 	  // write out UNSCALED occupancies for development statistics
// 	  // 
// 	  fprintf(out,"%3s:%3d-%4s  %-5s   %7.2f\n",
// 		  (*atom1)->residue.c_str(),
// 		  (*atom1)->residue_no,
// 		  (*atom1)->name.c_str(),
// 		  (*atom1)->sybyl_type.c_str(),
// 		  (*atom1)->desolvation_factor);
// 	  */	  
// 	  ////////////////////
	  
// 	  (*atom1)->desolvation_factor = ((*atom1)->desolvation_factor - get_min_desolvation_factor((*atom1)->generic_key))/
// 	    (get_max_desolvation_factor((*atom1)->generic_key) - get_min_desolvation_factor((*atom1)->generic_key));

// 	  (*atom1)->tempFactor = (*atom1)->desolvation_factor; //to make nice pictures in pymol	  
	  
// 	  // write out SCALED stats for surface/core analysis
// 	  /*
// 	  fprintf(out,"%3s:%3d-%4s  %-5s   %7.2f\n",
// 		  (*atom1)->residue.c_str(),
// 		  (*atom1)->residue_no,
// 		  (*atom1)->name.c_str(),
// 		  (*atom1)->sybyl_type.c_str(),
// 		  (*atom1)->desolvation_factor);
// 	  */
// 	}
//       else
// 	(*atom1)->tempFactor = 0.0; //to make nice pictures in pymol	  

//       atom1++;
//     }
  cout  <<"DE end"<<endl<<flush;

  return;

}



void Desolvation::set_parameters()
{
  ifstream in;
  in.open("parameters/desolvation_parameters.txt");

  if(! in.is_open())
    {
      cout << "Error:File 'parameters/desolvation_parameters.txt' could not be found\n";
      exit(0);
    }

  char line[256];
  string s;
  while(!in.eof())
    {
      in.getline(line,256);
      s = string(line);
      if (s.substr(0,4) == "ADP:")
	{
	  string atom_name = s.substr(5,6);

	  for(unsigned int i=0; i<atom_name.size(); )
	    if(string(atom_name,i,1) ==" ")
	      atom_name = atom_name.erase(i, 1);
	    else
	      i++;

	

	  desolvation_parameters[atom_name] = atof(s.substr(13,8).c_str());
	  
	  // cout<<"Atom "<<atom_name<<" has desolvation energy "<<desolvation_parameters[atom_name]<<endl;
	  
	}
    }
  in.close();
  
  
  // read in min and max occupancies
  in.open("parameters/occupancies.txt");

  if(! in.is_open())
    {
      cout << "Error:File 'parameters/occupancies.txt' could not be found\n";
      exit(0);
    }

  while(!in.eof())
    {
      in.getline(line,256);
      //cout<<line<<endl;
      s = string(line);
      
      if (s.size()>3 && s.substr(0,1) != "#")
	{
	  string key = s.substr(0,8);
	  for(unsigned int i=0; i<key.size(); )
	    if(string(key,i,1) ==" ")
	      key = key.erase(i, 1);
	    else
	      i++;

	  float min = atof(s.substr(32,40).c_str());
	  float max = atof(s.substr(41,51).c_str());

	  min_occupancies[key]= min;
	  max_occupancies[key]= max;
	  
	  //cout<<"Atom "<<key<<" has the min and max occupancies "<< min_occupancies[key]<<"   "<<max_occupancies[key]<<endl;
	}
      

   }

  // set average and default min and max occupancies
  avr_min_desolvation_factor = 0.0;
  for(map<string,float>::iterator it=min_occupancies.begin(); it != min_occupancies.end(); it++)
    avr_min_desolvation_factor += it->second;
  avr_min_desolvation_factor = avr_min_desolvation_factor/min_occupancies.size();

  avr_max_desolvation_factor = 0.0;
  for(map<string,float>::iterator it=max_occupancies.begin(); it != max_occupancies.end(); it++)
    avr_max_desolvation_factor += it->second;
  avr_max_desolvation_factor = avr_max_desolvation_factor/max_occupancies.size();

  //printf("Min and max occupancies are %f and %f\n",avr_min_desolvation_factor, avr_max_desolvation_factor);




  // open file for printing atomic occupancies
  //out = fopen("atomic_occupancies.txt","w");

  
}


void Desolvation::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: desolvation.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
