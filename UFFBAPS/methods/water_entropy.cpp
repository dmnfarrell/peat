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


#include "water_entropy.h"


Water_Entropy::Water_Entropy():Method()
{
  set_parameters();  


}

Water_Entropy::Water_Entropy(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  set_parameters();  



}

Water_Entropy::Water_Entropy(FILE * resultsFile, FILE * reporter, vector<Soup *> A):Method(reporter)
{

  /*** write out to out file **/
  writeDescription(reporter);
 
  set_parameters();


  float result;

  vector<Atom*> ligand_atoms;
  vector<vector<float> > distances;

  Distances dist(reporter);
  BondMatrix BM(reporter);


  for(unsigned int i=0; i<A.size(); i++)
    {

      /*** make atom vectors ***/
      ligand_atoms = A[i]->getAtoms();

      /*** find bonds ***/
      BM.assign_bonds(A[i]);

      result = calculate(ligand_atoms);

      string decoy_number;
      if(A[i]->name.find_last_of("_")== string::npos)
      	decoy_number = "-1";
      else
	decoy_number = A[i]->name.substr(A[i]->name.find_last_of("_")+1);

      fprintf(resultsFile,"%-4s %-4s %8f\n",
	      A[i]->name.substr(0,4).c_str(),
	      decoy_number.c_str(),
	      result);


    }  

}


Water_Entropy::~Water_Entropy(){}


float Water_Entropy::calculate(vector<Atom*> ligand_atoms)
{


  /*** find hydrophobic protein atoms ***/
  HP.find_hydrophobic_atoms(ligand_atoms);
  
  float res = 0.0;
  
  for(unsigned int l=0; l<ligand_atoms.size(); l++)
    if(ligand_atoms[l]->is_hydrophobic)
      res = res +1.0;

  res = -res *4.184*0.1; //each hydrophobic atom counts 0.1 kcal/mol
  return res;
}




void Water_Entropy::set_parameters()
{
}






void Water_Entropy::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: water_entropy.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());


  rep = reporter;
}
