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
$Id: hp.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "hp.h"


Hp::Hp(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  set_parameters();  

}

Hp::Hp():Method()
{
  set_parameters();  
}


Hp::Hp(FILE * resultsFile, FILE * reporter, vector<Soup *> A, vector<Soup *> B):Method(reporter)
{

  /*** write out to out file **/
  writeDescription(reporter);
 
  set_parameters();

  if(A.size() != B.size())
    {
      cout << "WARNING: An equal amount of protein and ligand structures must be given to the HP function!\n";
      return;
    }

  float result;

  for(unsigned int i=0; i<A.size(); i++)
    {
      result = calculate(A[i],B[i]);

      string decoy_number;
      if(B[i]->name.find_last_of("_")== string::npos)
      	decoy_number = "-1";
      else
	decoy_number = B[i]->name.substr(B[i]->name.find_last_of("_")+1);

      fprintf(resultsFile,"%-4s %-4s %8f\n",
	      A[i]->name.substr(0,4).c_str(),
	      decoy_number.c_str(),
	      result);


    }  

}


Hp::Hp(FILE * resultsFile, FILE * reporter, vector<Soup *> A):Method(reporter)
{

  /*** write out to out file **/
  writeDescription(reporter);
 
  set_parameters();

  float result;

  for(unsigned int i=0; i<A.size(); i++)
    {
      result = calculate(A[i]);

      fprintf(resultsFile,"%-4s %8f\n",
	      A[i]->name.substr(0,4).c_str(),
	      result);
    }  
}


Hp::~Hp(){}



float Hp::calculate(Soup * A, Soup * B)
{
  float result;

  vector<Atom*> protein_atoms, ligand_atoms;

  vector<vector<float> > distances;

  /*** make atom vectors ***/
  protein_atoms = A->getAtoms();
  ligand_atoms = B->getAtoms();

  /*** find bonds ***/
  if(!A->bonds_assigned)
    BM.assign_bonds(A);
  if(!B->bonds_assigned)
    BM.assign_bonds(B);

  /*** find inter-molecule distances ***/
  dist.calculate(protein_atoms, ligand_atoms, false);
  distances = dist.getResult();

  result = calculate(protein_atoms, ligand_atoms, distances, A->assigned_parameters, B->assigned_parameters);

  return result;
}


float Hp::calculate(vector<Atom*> protein_atoms, 
		    vector<Atom*> ligand_atoms,
		    vector<vector<float> > distances, 
		    string assigned_parametersA,
		    string assigned_parametersB)
{

//   if(include_surf_factor == 1)
//     {
//       Topology top(rep);
//       top.generate_topology(protein_atoms, ligand_atoms);
//       ligand_surf_dists = top.get_ligand_surface_distances();
//       protein_surf_dists = top.get_protein_surface_distances();
//     }


  /*** set vdw radii ***/
  if (assigned_parametersA=="")
    {
      printf("No parameters have been set for protein atoms, using GP vdw distances\n");
      Generalparameters GP(rep);
      GP.set_vdw_radii(protein_atoms);
    }
  if(assigned_parametersB=="")
    {
      printf("No parameters have been set for ligand atoms, using GP vdw distances\n");
      Generalparameters GP(rep);
      GP.set_vdw_radii(ligand_atoms);
    }

  float res = 0;
  float temp_result;

  /*** find hydrophobic protein atoms ***/
  find_hydrophobic_atoms(protein_atoms);
  find_hydrophobic_atoms(ligand_atoms);
  
  for(unsigned int p=0; p<protein_atoms.size(); p++)
    if(protein_atoms[p]->is_hydrophobic)
      for(unsigned int l=0; l<ligand_atoms.size(); l++)
	if(ligand_atoms[l]->is_hydrophobic)
	  {
	    temp_result = interaction(sqrt(distances[p][l]), protein_atoms[p]->vdw_dist + ligand_atoms[l]->vdw_dist);
	    if(include_surf_factor == 1)
	      temp_result = temp_result*(protein_atoms[p]->desolvation_factor + ligand_atoms[l]->desolvation_factor)/2; //using average desolvation_factor factor of the two atoms
	    
	    res += temp_result;

	   
// 	    printf("Hydrophobic protein atom %s-%d and Ligand atom %s-%d at distance %f has interaction %f\n",
// 		   protein_atoms[p]->element.c_str(),
// 		   p+1,
// 		   ligand_atoms[l]->element.c_str(),
// 		   l+1,
// 		   sqrt(distances[p][l]),
//		   interaction(sqrt(distances[p][l]), protein_atoms[p]->vdw_dist + ligand_atoms[l]->vdw_dist));
	   }
  


//    printf("Testing interaction:\n");
//    float d =0.0;
//    float r1 = 1.70;
//    float r2 = 1.20;
//    while(d<10)
//      {
//        printf("At %3f interaction is %f\n",d,interaction(d, r1+r2));
//        d =d +0.1;
//      }
    
  return res;
}


float Hp::calculate(Soup * A)
{
 
  float result;

  vector<Atom*> atoms;
  vector<vector<int> > bonds;

  /*** make atom vector ***/
  atoms = A->getAtoms();
  
  /*** find bonds ***/
  if(!A->bonds_assigned)
    BM.assign_bonds(A);

  result = calculate(atoms, A->assigned_parameters);

  return result;
}


float Hp::calculate(vector<Atom*> atoms,
		    string assigned_parameters)
{

  /*** set vdw radii ***/
  if (assigned_parameters=="")
    {
      printf("No parameters have been set for atoms, using GP vdw distances\n");
      Generalparameters GP(rep);
      GP.set_vdw_radii(atoms);
    }

  float res = 0;
  float temp_result;

  Distances dist(rep);
  
  /*** find hydrophobic protein atoms ***/
  find_hydrophobic_atoms(atoms);
  
  for(unsigned int p=0; p<atoms.size(); p++)
    if(atoms[p]->is_hydrophobic)
      for(unsigned int l=0; l<p; l++)
	if(atoms[l]->is_hydrophobic)
	  {
	    temp_result = interaction(dist.calculate(atoms[p],atoms[l],true), atoms[p]->vdw_dist + atoms[l]->vdw_dist);
	    res += temp_result;
	            
	    printf("Hydrophobic atoms %s-%d and %s-%d at distance %f has interaction %f\n",
		   atoms[p]->element.c_str(),
		   p+1,
		   atoms[l]->element.c_str(),
		   l+1,
		   dist.calculate(atoms[p],atoms[l],true),
		   temp_result);
	  }
  

   
  return res;
}





inline float Hp::interaction(float dist, float sum_vdw_radii)
{

  float R1 = sum_vdw_radii + 0.5;
  float R2 = R1 + 3.0;
  float res = -4.184 * 0.01; //each hydrophobic interaction counts 0.01 kcal/mol

  if(R1 < dist && dist < R2 )
    res = res * (1 - (dist - R1)/(3.0));
  else if ( dist >= R2)
    res = 0;

  return res;

}


void Hp::find_hydrophobic_atoms(vector<Atom*> atoms)
{

  bool hydrophobic;

  for(vector<Atom*>::iterator atom = atoms.begin(); atom != atoms.end(); atom++)
    {
      hydrophobic = false;

      if(include_non_hetero_carbon == 1 && is_non_hetero_carbon(*atom))
	hydrophobic = true;

      if(include_non_hetero_hydrogen == 1 && is_non_hetero_hydrogen(*atom))
	hydrophobic = true;

      if(include_halogens == 1 && is_non_ionic_halogen(*atom))
	hydrophobic = true;

      (*atom)->is_hydrophobic = hydrophobic;


      //      cout<<"Hydrophobic atom: "<<(*atom)->is_hydrophobic<<" ";
      //      (*atom)->display();
      
    }


  return;
}


bool Hp::is_non_hetero_hydrogen(Atom * atom)
{
  bool res = false;
 
  if(atom->element != "H")
    return false;

  if(atom->bonded_atoms.size()==1)
    {
      if(is_non_hetero_carbon(*(atom->bonded_atoms.begin())))
	 res =true;
    }
  else
    {
      cout<<"WARNING: This hydrogen does not have exactly one bond: ";
      atom->display();
      for(unsigned int j =0; j<atom->bonded_atoms.size(); j++)
	{
	  cout<<"   ";
	  atom->bonded_atoms[j]->display();
	}
    }

  /*
    for(unsigned int i=0; i<local_bonds.size(); i++)
    if(local_bonds[i][0] == a)
    if(is_non_hetero_carbon(local_bonds[i][1]))
    res = true;
  */
  return res;


}


bool Hp::is_non_hetero_carbon(Atom * atom)
{
  bool res = true;
 

  if(atom->element != "C")
    return false;

  for(vector<Atom*>::iterator bonded_atom = atom->bonded_atoms.begin(); 
      bonded_atom != atom->bonded_atoms.end(); bonded_atom++)
    if((*bonded_atom)->element != "C" && (*bonded_atom)->element != "H")
      res = false;
  
  /*

    for(unsigned int i=0; i<local_bonds.size(); i++)
    if(local_bonds[i][0] == a) 
      if(local_atoms[local_bonds[i][1]]->element != "C" && local_atoms[local_bonds[i][1]]->element != "H")
	res = false;
  */
  return res;

}


bool Hp::is_non_ionic_halogen(Atom * atom)
{
  bool res = false;

  /*** check if halogen ***/
  vector<string>::iterator result;
  result = find(halogens.begin(), halogens.end(), atom->element);            

  if(result == halogens.end()) 
    return false;

  /*** check that atom is not an ion ***/
  if(atom->bonded_atoms.size()>0)
    res = true;

  /*
  for(unsigned int i=0; i<local_bonds[a].size(); i++)
    if(local_bonds[i][0] == a) 
      res = true;
  */
  return res;
}


void Hp::set_parameters()
{


  halogens.push_back("F");
  halogens.push_back("CL");
  halogens.push_back("BR");
  halogens.push_back("I");

  ifstream parameters;
  parameters.open("hp.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'hp.cfg' could not be found\n";
      exit(0);
    }


  string dummy;

  while (!parameters.eof())
    {

      parameters >> dummy;

      if(dummy == "INCLUDE_NON_HETERO_CARBON:")
	parameters >> include_non_hetero_carbon;
      if(dummy == "INCLUDE_NON_HETERO_HYDROGEN:")
	parameters >> include_non_hetero_hydrogen;
      if(dummy == "INCLUDE_HALOGENS:")
	parameters >> include_halogens;
      if(dummy == "INCLUDE_SURF_FACTOR:")
	parameters >> include_surf_factor;

    } //end read file

  //printing all configurations
  /*
  printf("Configurations for hydrophobicity:\n");
  printf("\n");
  printf("\n");

  printf("INCLUDE_NON_HETERO_CARBON:     %d\n", include_non_hetero_carbon);
  printf("\n");

  printf("INCLUDE_NON_HETERO_HYDROGEN:   %d\n", include_non_hetero_hydrogen);
  printf("\n");

  printf("INCLUDE_HALOGENS:              %d\n", include_halogens);
  printf("\n");

  printf("INCLUDE_SURF_FACTOR:           %d\n", include_surf_factor);
  printf("\n");

  printf("Version: $Id: hp.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");
  */


  

}






void Hp::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: hp.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());


  rep = reporter;
}
