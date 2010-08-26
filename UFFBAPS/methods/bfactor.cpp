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
$Id: bfactor.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "bfactor.h"


Bfactor::Bfactor(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  
}

Bfactor::Bfactor(FILE * reporter, vector<Soup*> proteins):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  
  calculate(proteins);

}


Bfactor::Bfactor(FILE * reporter, vector<Soup*> proteins, vector<Soup*> ligands):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  
  calculate(proteins,ligands);

}



Bfactor::~Bfactor(){}


float Bfactor::surface_dist_factor(float protein_distance, float ligand_distance)
{
  float min_dist = min(protein_distance, ligand_distance);
  
  float alpha = 0.5;
  
  if(4.0 < min_dist)
      alpha = 1-1/(min_dist - 2);

  return alpha;
}


void Bfactor::calculate(vector<Soup*> proteins, vector<Soup*> ligands)
{

  if(proteins.size() != ligands.size())
    {
      printf("WARNING: Bfactors must equal numbers of proteins and ligands as input\n");
      exit(0);
    }

  
  FILE * dataBindingSite = fopen("bfactorsBindingSite.gnu","w");
  FILE * dataElsewhere = fopen("bfactorsNotBindingSite.gnu","w");
  FILE * data = fopen("bfactorsAll.gnu","w");
  FILE * dataRandom = fopen("bfactorsRand4pct.gnu","w");
  FILE * dataHistoBinding = fopen("histoBindingSite.gnu","w");
  FILE * dataHistoAll = fopen("histoAll.gnu","w");

  //histogram
  int no_bins = 20;
  //  int min = 0;
  //  int max = 20;

  vector<float> avr_b_b;
  avr_b_b.assign(no_bins,0.0);
  vector<float> avr_b_all;
  avr_b_all.assign(no_bins,0.0);


  vector<int> no_vals_b;
  no_vals_b.assign(no_bins,0);
  vector<int> no_vals_all;
  no_vals_all.assign(no_bins,0);

  
  Topology top(rep);
  Distances dist(rep);

  
  vector<Atom*> protein_atoms;
  vector<Atom*> ligand_atoms;


  for(unsigned int i=0; i<proteins.size(); i++)
    {
      protein_atoms = proteins[i]->getAtoms();
      ligand_atoms  = ligands[i]->getAtoms();

      dist.calculate(protein_atoms,ligand_atoms,false);
      dists = dist.getResult();

      top.generate_topology(protein_atoms,ligand_atoms);
      
      for(unsigned int p=0; p<protein_atoms.size(); p++)
	if(protein_atoms[p]->element != "H")
	  if(protein_atoms[p]->tempFactor > 0.0001)
	    { 
	      float surf_dist = top.find_surface_distace(protein_atoms[p]);
	      //all
	      fprintf(data, "%8.5f  %8.5f\n",surf_dist,protein_atoms[p]->tempFactor);
	      
	      //do histogram stuff
	      int bin_no = (int) surf_dist;
	      avr_b_all.at(bin_no) += protein_atoms[p]->tempFactor;
	      no_vals_all.at(bin_no) ++;
	      
	      
	      //4 pct
	      if(p%25 == 0)
		fprintf(dataRandom, "%8.5f  %8.5f\n",surf_dist,protein_atoms[p]->tempFactor);
	      
	      
	      //binding site
	      if(is_close(p))
		{
		  printf("Protein atom %4d has surf. dist. %8f and bfactor %8f - at binding site\n",
			 p+1,
			 surf_dist,
			 protein_atoms[p]->tempFactor);
		  
		  fprintf(dataBindingSite, "%8.5f  %8.5f\n",surf_dist,protein_atoms[p]->tempFactor);
		  
		  //do histogram stuff
		  int bin_no = (int) surf_dist;
		  avr_b_b.at(bin_no) += protein_atoms[p]->tempFactor;
		no_vals_b.at(bin_no) ++;
		
		}
	      
	      //not binding site
	      else
		{
		  
		  printf("Protein atom %4d has surf. dist. %8f and bfactor %8f - somewhere else\n",
			 p+1,
			 surf_dist,
			 protein_atoms[p]->tempFactor);
		  
		  fprintf(dataElsewhere, "%8.5f  %8.5f\n",surf_dist,protein_atoms[p]->tempFactor);
		}
	  }
      
    }
  

  //histogram
  for(unsigned int i=0; i<avr_b_b.size(); i++)
    {
      avr_b_b[i] = avr_b_b[i]/no_vals_b[i];
      fprintf(dataHistoBinding,"%8f  %8f\n", (float) i,avr_b_b[i]);

      avr_b_all[i] = avr_b_all[i]/no_vals_all[i];
      fprintf(dataHistoAll,"%8f  %8f\n", (float) i,avr_b_all[i]);

      
    }
}

void Bfactor::calculate(vector<Soup*> proteins)
{
  
  FILE * data = fopen("bfactorsAll.gnu","w");
  // FILE * dataRandom = fopen("bfactorsRand4pct.gnu","w");

  Topology top(rep);
  
  vector<Atom*> protein_atoms;
  
  for(unsigned int i=0; i<proteins.size(); i++)
    {
      protein_atoms = proteins[i]->getAtoms();

      top.generate_topology(protein_atoms);
      
      for(unsigned int p=0; p<protein_atoms.size(); p++)
	if(protein_atoms[p]->element != "H")
	  {
	    float surf_dist = top.find_surface_distace(protein_atoms[p]);
	    printf("Protein atom %4d has surf. dist. %8f and bfactor %8f\n",
		   p+1,
		   surf_dist,
		 protein_atoms[p]->tempFactor);
	    
	    fprintf(data, "%8.5f  %8.5f\n",surf_dist,protein_atoms[p]->tempFactor);
	  } 
    }

}


inline bool Bfactor::is_close(int p)
{
  bool res = false;

  for(unsigned int l=0; l<dists[p].size(); l++)
    if(dists[p][l]<25)  //5 A squared 
      res = true;

  return res;
}


void Bfactor::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: bfactor.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  rep = reporter;
  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
