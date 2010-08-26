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


#include "pdbwriter.h"

Pdbwriter::Pdbwriter():Method(){}

Pdbwriter::Pdbwriter(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  
}

Pdbwriter::Pdbwriter(FILE * reporter, string fn, Soup * A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  write(fn, A);

}

Pdbwriter::~Pdbwriter(){}


void Pdbwriter::write(string fn, Soup * A)
{

  pdb = fopen( fn.c_str(), "w" );
  atom_no = 1;
  resi_no = 1;

  fprintf(pdb, "Generated using genialNavn\n");

  vector<ProteinChain *> pc = A->getProteinChains();
  vector<Water *> wa = A->getWaters();
  vector<Molecule *> mo = A->getMolecules();
  vector<Residue *> residues;
  vector<Atom *> a;



  /*** write protein chains ***/
  for(unsigned int i=0; i<pc.size(); i++)
    {
      residues = pc[i]->getResidues();
      for(unsigned int r=0; r<residues.size(); r++)
	{
	  make_atoms("ATOM  ", residues[r]->getAtoms());
	  resi_no++;
	}
    }

  /*** write waters ***/


  /*** write molecules ***/


  fclose(pdb);

}





void Pdbwriter::make_atoms(string tag, vector<Atom *> a)
{



  for(unsigned int k=0; k< a.size(); k++)
    {
     
      fprintf(pdb,"%s %4d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
	      tag.c_str(), atom_no, a[k]->name.c_str(), a[k]->residue.c_str(), resi_no, a[k]->x, a[k]->y, a[k]->z,a[k]->occupancy,a[k]->tempFactor);
     
      atom_no++;

    }

}













void Pdbwriter::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: pdbwriter.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
