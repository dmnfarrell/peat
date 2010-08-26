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


#include "protein_entropy.h"


ProteinEntropy::ProteinEntropy():Method()
{
  /*** set parameters ***/
  set_parameters();
}


ProteinEntropy::ProteinEntropy(FILE * reporter):Method(reporter)
{
  /*** set parameters ***/
  set_parameters();

 
  /*** write out to out file **/
  writeDescription(reporter);


  
}

ProteinEntropy::ProteinEntropy(FILE * resultsFile, FILE * reporter, vector<Soup*> A, vector<Soup*> B):Method(reporter)
{
  /*** set parameters ***/
  set_parameters();


  /*** write out to out file **/
  writeDescription(reporter);

  if(A.size() != B.size())
    {
      cout << "WARNING: An equal amount of protein and ligand structures must be given to the ProteinEntropy function!\n";
      return;
    }

  vector<vector<float> > distances;  

  for(unsigned int i=0; i<A.size(); i++)
    {
      float result = calculate(A[i],B[i]);

      string decoy_number;
      if(B[i]->name.find_last_of("_")== string::npos)
      	decoy_number = "-1";
      else
	decoy_number = B[i]->name.substr(B[i]->name.find_last_of("_")+1);
      
      fprintf(resultsFile,"%-4s %-4s %8f\n",
	      A[i]->name.substr(0,4).c_str(),
	      decoy_number.c_str(),
	      result);
      
      
      //      fprintf(resultsFile,"%-10s %f\n",convert_soup_to_objects(A[i])[0]->identify().substr(0,4).c_str(),result);
    }  


}

ProteinEntropy::~ProteinEntropy(){}


float ProteinEntropy::calculate(Soup * A, Soup * B)
{

  vector<ProteinChain*> pcs = A->getProteinChains();
  vector<Residue*> res;
  vector<Atom*> ligand_atoms = B->get_dry_atoms();
  float entropy = 0;

  //for each protein chain 
  for(unsigned int i=0; i<pcs.size(); i++)
    {
      res = pcs[i]->getResidues();
      for(unsigned int r=0; r<res.size(); r++)
	entropy += check_residue(res[r], ligand_atoms);
    }

  return temperature*entropy*4.184;
}



float ProteinEntropy::check_residue(Residue * residue, vector<Atom*> ligand_atoms)
{
  vector<Atom*> res_atoms = residue->getAtoms();

  Distances dist;
  dist.calculate(res_atoms,ligand_atoms,true);
  vector<vector<float> >  dists = dist.getResult();
  
   //is any ligand atom closer than threshold to any residue atom?
  bool close = false;
  for(unsigned int r=0; r<res_atoms.size(); r++)
    if (res_atoms[r]->element != "H") // disregard hydrogen
      if (find(backbone_atoms.begin(), backbone_atoms.end(), res_atoms[r]->name ) == backbone_atoms.end()) //disregard back bone atoms
	for(unsigned int l=0; l<ligand_atoms.size(); l++)
	  if(dists[r][l] < threshold)
	    close = true;
  
  //what is the entropy loss of this residue?
  float res = 0;
  if(close)
    {  
      if(method == "count")
	res = 1;
      else if(method == "abagyan")
	res = abagyan[residue->residueType];
      else
	{
	  res = 1;
	  printf("Warning: Invalid method '%s' in Protein Entropy\n",method.c_str());
	}
    }  

  return res;
}


float ProteinEntropy::calculate(Soup * A)
{
  // Calculate protein entropy intra-soup
  //
  // Entropic penalties for all residues that are closer than threshold to another residue 
  // 
  // Neighbour residues should be scaled ~ 0.5
  //
  //
  /*  
  vector<ProteinChain*> pcs = A->getProteinChains();
  vector<Residue*> residues;
  vector<Atom*> ligand_atoms = B->get_dry_atoms();
  float entropy = 0;

  //for each protein chain 
  for(unsigned int i=0; i<pcs.size(); i++)
    {
      res = pcs[i]->getResidues();
      for(unsigned int r=0; r<res.size(); r++)
	entropy += check_residue(res[r], ligand_atoms);
    }

  return temperature*entropy;
  */

  return 0.0;
}


void ProteinEntropy::set_parameters()
{


  /*** Set up side-chain entropi library ***/
  /*** ref: Abagyan, J Mol Biol 235:983-1002 ***/
  /*** units: kcal/(mol*K) ***/


  abagyan["ALA"] = 0.0000000;
  abagyan["ARG"] = 0.0071000;
  abagyan["ASN"] = 0.0027000;
  abagyan["ASP"] = 0.0020333;
  abagyan["CYS"] = 0.0038000;
  abagyan["GLN"] = 0.0067333;
  abagyan["GLU"] = 0.0055000;
  abagyan["GLY"] = 0.0000000;
  abagyan["HIS"] = 0.0033000;
  abagyan["ILE"] = 0.0025000;
  abagyan["LEU"] = 0.0025000;
  abagyan["LYS"] = 0.0073667;
  abagyan["MET"] = 0.0051000;
  abagyan["PHE"] = 0.0029333;
  abagyan["PRO"] = 0.0000000;
  abagyan["SER"] = 0.0020667;
  abagyan["THR"] = 0.0020333;
  abagyan["TRP"] = 0.0032333;
  abagyan["TYR"] = 0.0033000;
  abagyan["VAL"] = 0.0016667;


  // make list of backbone atom names
  string bb[] = {"C", "CA", "N", "O"};
  backbone_atoms = vector<string>(bb,bb+4);

  ifstream parameters;
  parameters.open("protein_entropy.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'protein_entropy.cfg' could not be found\n";
      exit(0);
    }


  string dummy;

  while (!parameters.eof())
    {

      parameters >> dummy;


      if(dummy == "THRESHOLD:")
	parameters >> threshold;
      if(dummy == "TEMPERATURE:")
	parameters >> temperature;
      if(dummy == "METHOD:")
	parameters >> method;


    } //end read file

  //printing all configurations
  /*
  printf("Configurations for Protein Entropy:\n");
  printf("\n");
  printf("\n");

  printf("THRESHOLD:              %f\n", threshold);
  printf("\n");

  printf("TEMPERATURE:            %f\n", temperature);
  printf("\n");

  printf("METHOD:                 %s\n", method.c_str());
  printf("\n");

  printf("Version: $Id: protein_entropy.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");
  */
}




void ProteinEntropy::writeDescription(FILE * reporter)
{

  rep = reporter;
  /*** write out to out file **/
  version =string("$Id: protein_entropy.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
