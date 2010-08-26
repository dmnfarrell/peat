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


#include "simple_parameters.h"


Simple_Parameters::Simple_Parameters()
{
  set_parameters();
}


Simple_Parameters::Simple_Parameters(FILE * reporter):Method(reporter)
{

  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();
}

Simple_Parameters::Simple_Parameters(FILE * reporter, vector<Soup*> A):Method(reporter)
{

  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();

  prepare(A);
}

Simple_Parameters::~Simple_Parameters(){}

void Simple_Parameters::prepare(vector<Soup*> A)
{
  for(vector<Soup*>::iterator it = A.begin(); it != A.end(); it++)
    prepare(*it);
}

void Simple_Parameters::prepare(Soup * A)
{

  if (A->assigned_parameters != "SP") //make sure parameters are assigned only once
    {  
      if (!(A->sybyl_types_assigned))
	A->set_sybyl_atom_types(); //this also sets generic_keys

      set_electrostatic_parameters(A);
      set_dispersive_parameters(A);

      A->assigned_parameters = "SP";
    }

}

void Simple_Parameters::set_electrostatic_parameters(Soup * A)
{
  // protein chains, molecules, waters
  vector<vector<SoupObject *> > objects = A->getSoupObjects();  
  vector<SoupObject *> protein_chains = objects[0];
  vector<SoupObject *> molecules = objects[1];
  vector<SoupObject *> waters = objects[2];

  //protein chains
  for(vector<SoupObject *>::iterator protein_chain = protein_chains.begin(); protein_chain != protein_chains.end(); protein_chain++)
    set_electrostatic_parameters_for_protein_chain( (static_cast<ProteinChain *> (*protein_chain) ));

  //molecules
  for(vector<SoupObject *>::iterator molecule = molecules.begin(); molecule != molecules.end(); molecule++)
    set_electrostatic_parameters_for_molecule( (static_cast<Molecule *> (*molecule) ));

  //waters
  for(vector<SoupObject *>::iterator water = waters.begin(); water != waters.end(); water++)
    set_electrostatic_parameters_for_water( (static_cast<Water *> (*water) ));

}

void Simple_Parameters::set_electrostatic_parameters_for_protein_chain(ProteinChain * p)
{
  vector<Residue*> residues = p->getResidues();
  vector<Atom*> atoms;

  for(vector<Residue*>::iterator residue = residues.begin(); residue != residues.end(); residue++)
    {
      atoms = (*residue)->getAtoms();
      
      for(vector<Atom*>::iterator atom = atoms.begin(); atom != atoms.end(); atom++)
	if((*atom)->element != "H")
	  {
	    // set n-term
	    if (residue == residues.begin() && (*atom)->name=="N")
	      {
		(*atom)->charge = atomic_charges["NTERM"];
		//printf("N-TERM charge: %4.1f for ",atomic_charges["NTERM"]);
		//(*atom)->display();
	      }
	    
	    // set c-term
	    if (residue == residues.end()-1 && ((*atom)->name=="OXT" || (*atom)->name=="O''"))
	      {
		(*atom)->charge = atomic_charges["CTERM"];
		//printf("C-TERM charge: %4.1f for ",atomic_charges["CTERM"]);
		//(*atom)->display();
	      }
	    
	    // side chains
	    if(atomic_charges.find((*atom)->generic_key) != atomic_charges.end())
	      {
		(*atom)->charge = atomic_charges[(*atom)->generic_key];
		//printf("SChain charge: %4.1f for ", atomic_charges[(*atom)->generic_key]);
		//(*atom)->display();
	      }
	  }

    }
  
  return;
}

void Simple_Parameters::set_electrostatic_parameters_for_molecule(Molecule * m)
{
  //bzzzz - for now we are just assuming that molecules are read in from mol2 files with partial charges

  return;
}


void Simple_Parameters::set_electrostatic_parameters_for_water(Water * w)
{
  //The partial charge of water atoms is ignored

  vector<Atom*> atoms = w->getAtoms();
 
  for(vector<Atom*>::iterator atom = atoms.begin(); atom != atoms.end(); atom++)
    (*atom)->charge = 0.0;


  return;
}


void Simple_Parameters::set_dispersive_parameters(Soup * A)
{
  vector<Atom*> atoms = A->getAtoms();
  /*
    for(vector<Atom*>::iterator atom = atoms.begin(); atom != atoms.end(); atom++)
    {

      //(*atom)->vdw_dist = 1.5;
      //(*atom)->vdw_energy = 0.1; //Arbitrary values for testing
     }
  */
  // use the vdw paramters in general_parameters method
  
  GP.set_vdw_radii(atoms);



  return;

}




void Simple_Parameters::set_parameters()
{
  atomic_charges["ARG-CZ"]  =  1.0;
  atomic_charges["ASP-CG"]  = -1.0;
  atomic_charges["GLU-CD"]  = -1.0;
  atomic_charges["HIS-NE2"] =  0.5;
  atomic_charges["LYS-NZ"]  =  1.0;
  atomic_charges["NTERM"]   =  1.0;
  atomic_charges["CTERM"]   = -1.0;
  


}

void Simple_Parameters::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: simple_parameters.cpp 18 2005-11-15 09:56:08Z chresten $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
