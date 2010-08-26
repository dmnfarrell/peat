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
$Id: generalparameters.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "generalparameters.h"


Generalparameters::Generalparameters(FILE * reporter, vector<Soup *> A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  Correct_Atomnames correct_names(reporter,A);

  set_parameters();
  
  apply_parameters(A);

}

Generalparameters::Generalparameters(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  set_parameters();

}

Generalparameters::Generalparameters():Method()
{
  set_parameters();
}


Generalparameters::~Generalparameters(){}


void Generalparameters::apply_parameters(vector<Soup *> A)
{




  //make a list of all atoms in A
  vector<Atom*> all_atoms;
  vector<ProteinChain*> all_pcs;
  vector<Molecule*> all_mols;
  for(unsigned int s=0; s< A.size(); s++)
    if(A[s]->assigned_parameters != "GP") //make sure that parameters are not assigned more than once
      {
	vector<Atom*> temp_atoms = A[s]->getAtoms();
	for(unsigned int i=0; i< temp_atoms.size(); i++)
	  all_atoms.push_back(temp_atoms[i]);
	
	vector<ProteinChain*> temp_pcs = A[s]->getProteinChains();
	for(unsigned int i=0; i< temp_pcs.size(); i++)
	  all_pcs.push_back(temp_pcs[i]);
	
	vector<Molecule*> temp_mols = A[s]->getMolecules();
	for(unsigned int i=0; i< temp_mols.size(); i++)
	  if(temp_mols[i]->getAtoms().size()==1) //only include molecules with one atom => monoatomic ions
	    all_mols.push_back(temp_mols[i]);

	
	A[s]->assigned_parameters = "GP";
	
      }

  //apply parameters
  set_vdw_radii(all_atoms);
  set_protein_charges(all_pcs);
  set_ion_charges(all_mols);
  check_parameters(all_atoms);

  //print all atoms
  //print_all(A);

}


void Generalparameters::set_vdw_radii(vector<Atom *> all_atoms)
{
  for(unsigned int i=0; i<all_atoms.size(); i++)
    {
      if(vdw_radius.count(all_atoms[i]->element) == 1)
	all_atoms[i]->vdw_dist = vdw_radius[all_atoms[i]->element];
      else
	{
	  printf("WARNING: van der Waals radius of '%s' atom is set standard value of 2.0 A\n",all_atoms[i]->element.c_str());
	  all_atoms[i]->vdw_dist = 2.0;
	}
      
      //set the vdw energy to 1.00 kJ/mol
      all_atoms[i]->vdw_energy = 1.00; //kJ/mol //////////////////////////
      //all_atoms[i]->vdw_dist *= 1.5; // No hydrogens -> multiply all radii with 1.5

    }
}

void Generalparameters::set_protein_charges(vector<ProteinChain *> P)
{
  PrepinReader proteinCharges(rep);
  proteinCharges.find_protein_chain_charges(P);
}



void Generalparameters::set_ion_charges(vector<Molecule *> all_mols)
{
  for(unsigned int i=0; i<all_mols.size(); i++)
    {
      vector<Atom*> temp = all_mols[i]->getAtoms();
      for(unsigned int j=0; j<temp.size(); j++)
	{
	  if(ion_charge.count(temp[j]->element) == 1)
	    temp[j]->charge = ion_charge[temp[j]->element];
	  else
	    {
	      printf("WARNING: charge of ion '%s' could not be set\n",temp[j]->element.c_str());
	    }
	}
    }

}


void Generalparameters::check_parameters(vector<Atom *> all_atoms)
{}

void Generalparameters::print_all(vector<Soup *> A)
{
  for(unsigned int s=0; s< A.size(); s++)
    {
      printf("Soup '%s' has the atoms:\n",A[s]->name.c_str());
      printf("--------------------------------------------\n");
      printf("| NAME|  E|    radius|    charge|vdw_energy|\n");
      printf("--------------------------------------------\n");
      vector<Atom*> atoms = A[s]->getAtoms();

      for(unsigned int i=0; i< atoms.size(); i++)
	printf("|%5s|%3s|%10f|%10f|%10f|\n", 
	       atoms[i]->name.c_str(),
	       atoms[i]->element.c_str(),
	       atoms[i]->vdw_dist,
	       atoms[i]->charge,
	       atoms[i]->vdw_energy);
      printf("-------------------------------------------\n");
      
    }



}

void Generalparameters::set_parameters()
{

  //Ref: A Bondi: van der Waals Volumes and Radii
  //J. Phys. chem. 1964, 68(3), pp 441-451

  vdw_radius["AG"] =  1.72;
  vdw_radius["AR"] =  1.88;
  vdw_radius["AS"] =  1.85;
  vdw_radius["AU"] =  1.66;
  vdw_radius["BR"] =  1.85;      
  vdw_radius["C"]  =  1.70;    
  vdw_radius["CD"] =  1.58;     
  vdw_radius["CL"] =  1.75;
  vdw_radius["CU"] =  1.40;      
  vdw_radius["F"]  =  1.47;     
  vdw_radius["GA"] =  1.87;     
  vdw_radius["H"]  =  1.20;
  vdw_radius["HE"] =  1.40;      
  vdw_radius["HG"] =  1.55;     
  vdw_radius["I"]  =  1.98;     
  vdw_radius["IN"] =  1.93;
  vdw_radius["K"]  =  2.75;      
  vdw_radius["KR"] =  2.02;     
  vdw_radius["LI"] =  1.82;     
  vdw_radius["MG"] =  1.73;
  vdw_radius["N"]  =  1.55;      
  vdw_radius["NA"] =  2.27;     
  vdw_radius["NE"] =  1.54;     
  vdw_radius["NI"] =  1.63;
  vdw_radius["O"]  =  1.52;      
  vdw_radius["P"]  =  1.80;     
  vdw_radius["PB"] =  2.02;     
  vdw_radius["PD"] =  1.63;
  vdw_radius["PT"] =  1.72;
  vdw_radius["S"]  =  1.80;
  vdw_radius["SE"] =  1.90;
  vdw_radius["SI"] =  2.10;
  vdw_radius["SN"] =  2.17;
  vdw_radius["TE"] =  2.06;     
  vdw_radius["TL"] =  1.96;     
  vdw_radius["U"]  =  1.86;
  vdw_radius["XE"] =  2.16;      
  vdw_radius["ZN"] =  1.39;

  //halogens
  ion_charge["CL"] = -1;
  ion_charge["F"]  = -1;
  ion_charge["BR"] = -1;
  ion_charge["I"]  = -1;

  //Alkali metals
  ion_charge["NA"] = 1;
  ion_charge["K"]  = 1;

  //Alkaline earth metals
  ion_charge["CA"] = 2.0;
  ion_charge["MG"] = 2.0;

  //metals
  ion_charge["ZN"] = 2.0;
  ion_charge["MN"] = 2.0; //3.0
  ion_charge["FE"] = 2.0; //3.0
  ion_charge["NI"] = 2.0; //3.0
  ion_charge["CO"] = 2.0; //3.0
  ion_charge["CU"] = 2.0; //3.0
  ion_charge["CD"] = 2.0; //3.0

  //others
  ion_charge["P"] = -3.0;
  ion_charge["U"] = 4.0; //????????????


}

void Generalparameters::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: generalparameters.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());

  rep = reporter;
}
