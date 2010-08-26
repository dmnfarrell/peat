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
$Id: Chemical_Shift_Changes.cpp 18 2005-11-15 09:56:08Z chresten $
*/

#include "chemical_shift_changes.h"


Chemical_Shift_Changes::Chemical_Shift_Changes():Method()
{
  set_parameters();
}


Chemical_Shift_Changes::Chemical_Shift_Changes(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();
}

Chemical_Shift_Changes::Chemical_Shift_Changes(FILE * reporter, vector<Soup*> A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();

  for(vector<Soup*>::iterator soup = A.begin(); soup!=A.end(); soup++)
    amide_shifts(*soup);
  
}


Chemical_Shift_Changes::~Chemical_Shift_Changes(){}


void Chemical_Shift_Changes::amide_shifts(Soup * A)
{
  cout<<"Calculating amide shifts for "<<A->name<<endl;
  A->set_sybyl_atom_types();
  BM.assign_bonds(A);


  // for each titratable group
  //    for each N-H
  //       set atom bfactor
  // write out pdb with bfactors

  //set all tempFactors = 0
  vector<Atom*> atoms = get_atoms_from_protein_chains(A->getProteinChains());
  vector<Atom*>::iterator atom = atoms.begin();
  vector<vector<Atom*> >::iterator amide;
  while(atom != atoms.end())
    {
      (*atom)->tempFactor=0;
      atom++;
    }

  // find amides
  generate_list_of_amides(atoms);  
  atom = atoms.begin();
  
  float chem_shift, max;
  
  while(atom != atoms.end())
    {

      if (titratable_residues.count((*atom)->generic_key )>0)
	{
	  (*atom)->display();
	   max=-1000;
	  // set all amide chem shifts
	  amide = amides.begin();
	  while(amide != amides.end())
	    {
	    if(amide->at(0) != (*atom))
	      {
		/*
		cout<<"   ";
		amide->at(0)->display();
		cout<<"   ";
		amide->at(1)->display();
		cout<<"   ";
		amide->at(2)->display();
		*/

		chem_shift = induced_shift((*atom),amide->at(0), amide->at(1), amide->at(2));//units of e/(A*A)
		chem_shift = 1/(4*pi*eps_0*eps_rel)*alpha * e * chem_shift; //J m /C*C ppm/au C /(A*A) = J m /C ppm/au  /(A*A) = (J m ppm )/(C au A*A)
		chem_shift = chem_shift * 1e10*1e10;// (J m ppm )/(C au A*A) * A*A/(m*m) = (J ppm )/(C au m) = J / (C * m) * ppm/au
		chem_shift = chem_shift * (e*a0)/Eh;// J / (C * m) * ppm *C*m/J)= ppm 
		cout<<chem_shift<<endl;
		//amide->at(0)->tempFactor = -chem_shift;
		amide->at(1)->tempFactor = chem_shift;
		/*
		if(fabs(chem_shift) > max)
		  max = fabs(chem_shift);
		*/
	      }
	    else
	      {
		//amide->at(0)->tempFactor = 0;
		amide->at(1)->tempFactor = 0;
	      }
	    amide++;
	    }
	  // normalise chem shifts
	  //vector<Atom*>::iterator atom1 = atoms.begin();
	  // cout<<"max:"<<max<<endl;
	  /*	  while(atom1 != atoms.end())
	    {
	      (*atom1)->tempFactor = (*atom1)->tempFactor / max;
	      atom1++;
	    }
	  */


	  // print out pdb
	  
	  ostringstream o;
	  o << (*atom)->residue_no;

	  string fn = A->name + "_"+(*atom)->residue+"_"+o.str()+".pdb";
	  PDB_writer.write(fn,A);
	  
	  write_out_chemical_shifts(fn+".m",atoms);
	  
	}
     


      atom++;
    }

}

void Chemical_Shift_Changes::write_out_chemical_shifts(string fn, vector<Atom*> atoms)
{
  //writing NH bfactors to file
  
  FILE * mfile = fopen(fn.c_str(),"w");
  for(vector<Atom*>::iterator atom = atoms.begin(); atom != atoms.end(); atom++)
    if((*atom)->element == "H")
      {
	vector<Atom*>::iterator N = has_element_bonded(*atom,"N");
	if(N != (*atom)->bonded_atoms.end())
	  if(no_hydrogen_bounded(*N) == 1)
	    {
	      fprintf(mfile, "%1s%-4d %4s %7.2f\n", names3to1[(*atom)->residue].c_str(), (*atom)->residue_no, (*atom)->name.c_str(), (*atom)->tempFactor);
	    }
      }
  fclose(mfile);
  
}



void Chemical_Shift_Changes::generate_list_of_amides(vector<Atom*> atoms)
{
  amides.clear();
  vector<Atom*> temp;  
  vector<Atom*>::iterator atom = atoms.begin();
  
  while(atom != atoms.end())
    {
      if((*atom)->name=="N")
	if(no_hydrogen_bounded(*atom) == 1)
	  {
	    vector<Atom*>::iterator bonded_H =  has_element_bonded(*atom,"H");
	    vector<Atom*>::iterator bonded_C =  has_atom_name_bonded(*atom,"C");
	    
	    temp.clear();
	    temp.push_back((*atom));
	    temp.push_back((*bonded_H));
	    temp.push_back((*bonded_C));
	    amides.push_back(temp);
	    
	    cout<<"Amide: ";
	    (*atom)->display();
	    cout<<"       ";
	    (*bonded_H)->display();
	    cout<<"       ";
	    (*bonded_C)->display();
	    
	  }
      
      atom++;
    }
  
}

vector<Atom*>::iterator Chemical_Shift_Changes::has_element_bonded(Atom * atom, string element)
{

  vector<Atom*>::iterator res = (atom->bonded_atoms).end();
  vector<Atom*>::iterator batom = (atom->bonded_atoms).begin();

  while(batom != (atom->bonded_atoms).end())
    {
      if((*batom)->element == element)
	{
	  res = batom;
	}
      batom++;
    }

  return res;
}

vector<Atom*>::iterator Chemical_Shift_Changes::has_atom_name_bonded(Atom * atom, string name)
{

  vector<Atom*>::iterator res = (atom->bonded_atoms).end();
  vector<Atom*>::iterator batom = (atom->bonded_atoms).begin();

  while(batom != (atom->bonded_atoms).end())
    {
      if((*batom)->name == name)
	{
	  res = batom;
	}
      batom++;
    }

  return res;
}



int Chemical_Shift_Changes::no_hydrogen_bounded(Atom * atom)
{
  int res = 0;
  vector<Atom*>::iterator batom = atom->bonded_atoms.begin();
  while(batom != atom->bonded_atoms.end())
    {
      if((*batom)->element == "H")
	res++;
      batom++;
    }
  
  return res;
}



float Chemical_Shift_Changes::induced_shift(Atom * charged_atom, Atom * N, Atom * H, Atom * C)
{

  // v1: N -> charged_atom
  // v2: N -> C
  //
  //


  
  vector<float> v1,v2;
  v1.push_back(charged_atom->x - N->x);
  v1.push_back(charged_atom->y - N->y);
  v1.push_back(charged_atom->z - N->z);

  v2.push_back(C->x - N->x);
  v2.push_back(C->y - N->y);
  v2.push_back(C->z - N->z);

  // normalise vectors v1 and v2
  float l1,l2;
  l1 = sqrt(v1[0]*v1[0]+
	    v1[1]*v1[1]+
	    v1[2]*v1[2]);

  l2 = sqrt(v2[0]*v2[0]+
	    v2[1]*v2[1]+
	    v2[2]*v2[2]);

  v1[0] = v1[0]/l1;
  v1[1] = v1[1]/l1;
  v1[2] = v1[2]/l1;

  v2[0] = v2[0]/l2;
  v2[1] = v2[1]/l2;
  v2[2] = v2[2]/l2;
  

  float res = 0;
  res += v1[0]*v2[0];
  res += v1[1]*v2[1];
  res += v1[2]*v2[2];

  res = res/(l1*l1); //distance squared
  res = res * titratable_residues[charged_atom->generic_key];//sign of charge

  return res; //units of e/(A*A)

 



  /*  
  vector<float> v1,v2;
  v1.push_back((N->x+H->x)/2 - charged_atom->x);///avr btw N and H
  v1.push_back((N->y+H->y)/2 - charged_atom->y);
  v1.push_back((N->z+H->z)/2 - charged_atom->z);

  v2.push_back(H->x - N->x);
  v2.push_back(H->y - N->y);
  v2.push_back(H->z - N->z);

  float l1,l2;
  l1 = sqrt(v1[0]*v1[0]+
	    v1[1]*v1[1]+
	    v1[2]*v1[2]);

  l2 = sqrt(v2[0]*v2[0]+
	    v2[1]*v2[1]+
	    v2[2]*v2[2]);

  v1[0] = v1[0]/l1;
  v1[1] = v1[1]/l1;
  v1[2] = v1[2]/l1;

  v2[0] = v2[0]/l2;
  v2[1] = v2[1]/l2;
  v2[2] = v2[2]/l2;
  

  float res = 0;
  res += v1[0]*v2[0];
  res += v1[1]*v2[1];
  res += v1[2]*v2[2];

  res = res/(l1*l1); //distance squared
  res = res * titratable_residues[charged_atom->generic_key];//sign of charge

  return res; //units of e/(A*A)

  */

}

float Chemical_Shift_Changes::angle(vector<float> v1, vector<float> v2){return 0.0;}


void Chemical_Shift_Changes::set_parameters()
{
  titratable_residues["ASP-CG"] = -1;
  titratable_residues["GLU-CD"] = -1;
  titratable_residues["HIS-CG"] =  1;
  titratable_residues["LYS-NZ"] =  1;
  titratable_residues["ARG-CZ"] =  1;


  
  names3to1["ALA"] = "A";
  names3to1["ARG"] = "R";
  names3to1["ASN"] = "N";
  names3to1["ASP"] = "D";
  names3to1["CYS"] = "C";
  names3to1["GLU"] = "E";
  names3to1["GLN"] = "Q";
  names3to1["GLY"] = "G";
  names3to1["HIS"] = "H";
  names3to1["ILE"] = "I";
  names3to1["LEU"] = "L";
  names3to1["LYS"] = "K";
  names3to1["MET"] = "M";
  names3to1["PHE"] = "F";
  names3to1["PRO"] = "P";
  names3to1["SER"] = "S";
  names3to1["THR"] = "T";
  names3to1["TRP"] = "W";
  names3to1["TYR"] = "Y";
  names3to1["VAL"] = "V";

  for(map<string, string>::iterator i = names3to1.begin(); i != names3to1.end(); i++)
    names1to3[(*i).second] = (*i).first;

  //alpha = 70.0; //C-H bonds!
  alpha = 977.0; //A|| for C-N 
  pi = 3.14159265;
  eps_0 = 8.85419e-12;//C*C/(J*m)
  eps_rel = 8; // 1
  e = 1.602177e-19; //C
  a0 = 5.291772e-11; //m
  Eh = 4.359744e-18; //J

}


void Chemical_Shift_Changes::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: Chemical_Shift_Changes.cpp 18 2005-11-15 09:56:08Z chresten $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
