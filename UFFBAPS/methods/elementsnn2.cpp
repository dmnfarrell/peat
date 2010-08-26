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
$Id: elementsnn2.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "elementsnn2.h"


Elementsnn2::Elementsnn2(FILE * reporter)//:Nnmodel()
{

  /*** write out to out file **/
  writeDescription(reporter);
  
  printf("Using the elements 2 network!\n");

  set_parameters();
  if(include_hb_info == 1)
    HB.set_parameters();
 

}

Elementsnn2::~Elementsnn2(){}


void Elementsnn2::create_network()
{

  network.create_network(topology);
  network.set_learning_rate(learning_rate);
 
  network.set_transfer_functions(tf);
  network.set_linear_coef(linear_coef);
  printf("%s\n",network.print_network().c_str());
}


void Elementsnn2::build_inputs(pose * A)
{


  protein_atoms = convert_objects_to_atoms(A->protein);
  ligand_atoms = convert_objects_to_atoms(A->ligand);
  
  vector<float> res;


  if(include_surf_dist == 1)
    {
      Topology top(rep);
      top.generate_topology(protein_atoms, ligand_atoms);
      ligand_surf_dists = top.get_ligand_surface_distances();
      protein_surf_dists = top.get_protein_surface_distances();

    }

  if(A->pairs.size() == 0)
    {
      A->pairs = first_build(protein_atoms, ligand_atoms);

      //bonds 
      if(no_protein_bonds_incl > 0)
	{
	  //find protein bonds
	  BondMatrix proteinBonds(rep);
	  proteinBonds.calculate(protein_atoms);
	  A->protein_bonds = proteinBonds.getResult();
	}
      
      if(no_ligand_bonds_incl > 0)
	{
	  //find ligand bonds
	  BondMatrix ligandBonds(rep);
	  ligandBonds.calculate(ligand_atoms);
	  A->ligand_bonds = ligandBonds.getResult();
	}
    }

  if(no_protein_bonds_incl > 0)
    protein_bonds = A->protein_bonds;
	
  if(no_ligand_bonds_incl > 0)
    ligand_bonds = A->ligand_bonds;
	
    

  int p,l;

  for(unsigned int i=0; i<A->pairs.size(); i++)
    {
      p = A->pairs[i][0];
      l = A->pairs[i][1];

      //input for distance
      float dist = sqrt((protein_atoms[p]->x - ligand_atoms[l]->x)*(protein_atoms[p]->x - ligand_atoms[l]->x)+
			(protein_atoms[p]->y - ligand_atoms[l]->y)*(protein_atoms[p]->y - ligand_atoms[l]->y)+
			(protein_atoms[p]->z - ligand_atoms[l]->z)*(protein_atoms[p]->z - ligand_atoms[l]->z));
   

      res.clear();
      
      //input for protein atom
      res = protein_types[protein_atoms[p]->element.c_str()];
      
      if(no_protein_bonds_incl == 1)
	{
	  vector<float> protein_freqs = get_protein_freqs(p);
	  res.insert(res.end(), protein_freqs.begin(), protein_freqs.end()); 
	}
      
      //input for ligand atom
      res.insert(res.end(), ligand_types[ligand_atoms[l]->element.c_str()].begin(), ligand_types[ligand_atoms[l]->element.c_str()].end()); 

      
      if(no_ligand_bonds_incl == 1)
	{
	  vector<float> ligand_freqs = get_ligand_freqs(l);
	  res.insert(res.end(), ligand_freqs.begin(), ligand_freqs.end()); 
	}
      
      res.insert(res.end(),dist);
      
      //input for context
      if(include_surf_dist == 1)
	{
	  float pair_surf_dist = (protein_surf_dists[p] + ligand_surf_dists[l])/2;
	  res.insert(res.end(),pair_surf_dist);
	  //		    printf("Surf: %f\n",pair_surf_dist);
	}
      
      if(include_hb_info == 1)
	{
	  float hb = 0.0;
	  if(HB.is_hbond(protein_atoms[p],ligand_atoms[l],dist)||HB.is_hbond(ligand_atoms[l],protein_atoms[p],dist))
	    hb=1.0;
	  res.insert(res.end(),hb);
	  //printf("Atoms %4s and %4s hbond:%f\n",protein_atoms[p]->name.c_str()
	  //	 ,ligand_atoms[l]->name.c_str(), hb);
	  
	}
      
      if(res.size() == (unsigned int) no_inputs)
	{
	  A->inputs.push_back(res);
	  /*printf("inputs: ");
	  for (vector<float>::iterator itr = res.begin(); itr!=res.end(); itr++)
	    printf("%1.1f ",*itr );
	  printf("\n");
	  */
	}
      else
	{
	  printf("Warning: could not create input for atom pair '%s'-'%s'\n",
		 protein_atoms[p]->element.c_str(),ligand_atoms[l]->element.c_str());
	  
	  /*		    printf("%d \n", no_inputs);
	    for(unsigned int r=0; r<res.size(); r++)
	    printf("%f ", res[r]);
	    printf("\n");
	  */
	  
	}
      
    }
 


}


vector<vector<int> >  Elementsnn2::first_build(vector<Atom*> protein_atoms, vector<Atom*> ligand_atoms)
{
  vector<int> temp;
  vector<vector<int> > res;


  for(unsigned int p=0; p<protein_atoms.size(); p++)
    for(unsigned int l=0; l<ligand_atoms.size(); l++)
      {
	temp.clear();
	//input for distance
	float dist = sqrt((protein_atoms[p]->x - ligand_atoms[l]->x)*(protein_atoms[p]->x - ligand_atoms[l]->x)+
			  (protein_atoms[p]->y - ligand_atoms[l]->y)*(protein_atoms[p]->y - ligand_atoms[l]->y)+
			  (protein_atoms[p]->z - ligand_atoms[l]->z)*(protein_atoms[p]->z - ligand_atoms[l]->z));
   
	if(dist < threshold)
	  if(count(ignore_elements_p.begin(), ignore_elements_p.end(),protein_atoms[p]->element) < 1)
	    if(count(ignore_elements_l.begin(), ignore_elements_l.end(),ligand_atoms[l]->element) < 1)
	      {
		temp.push_back(p);
		temp.push_back(l);
		res.push_back(temp);
	      }
      }

  return res;
}

void Elementsnn2::define_input()
{

  int one_index = 0;
  
  vector<float> temp;
  //map protein type vectors
  for(int i=0; i<no_protein_types; i++)
    {
      for(int e=0; e<no_protein_types; e++)
	{
	  if(e == one_index)
	    temp.push_back(1);
	  else
	    temp.push_back(0);
	}
      protein_types[include_elements_p[i]] =  temp;
      temp.clear();
      one_index++;
    }

  //map ligand type vectors  
  one_index = 0;
  for(int i=0; i<no_ligand_types; i++)
    {
      for(int e=0; e<no_ligand_types; e++)
	{
	  if(e == one_index)
	    temp.push_back(1);
	  else
	    temp.push_back(0);
	}
      ligand_types[include_elements_l[i]] =  temp;
      temp.clear();
      one_index++;
    }


  //print out mapping
  for(int i=0; i<no_protein_types; i++)
    {
      printf("Protein element '%2s' has been mapped to input ",include_elements_p[i].c_str());
      for(int e=0; e<no_protein_types; e++)
	printf("'%.0f'",protein_types[include_elements_p[i]][e]);
      cout<<"\n";
    }

  printf("Ignored protein elements: ");
  for(unsigned int i=0; i<ignore_elements_p.size(); i++)
    printf("'%s' ", ignore_elements_p[i].c_str());
  printf("\n");

  for(int i=0; i<no_ligand_types; i++)
    {
      printf("Ligand element '%2s' has been mapped to input ",include_elements_l[i].c_str());
      for(int e=0; e<no_ligand_types; e++)
	printf("'%.0f'",ligand_types[include_elements_l[i]][e]);
      cout<<"\n";
    }

  printf("Ignored ligand elements: ");
  for(unsigned int i=0; i<ignore_elements_l.size(); i++)
    printf("'%s' ", ignore_elements_l[i].c_str());
  printf("\n");



}




void Elementsnn2::test_pose(pose * A)
{
  float energy = 0;


  /*** sum up energy contributions from all atom pairs ***/
  for(unsigned int i=0; i < A->inputs.size(); i++)
    {
      vector<float> out = network.test(A->inputs[i]);
      energy += out[0];
    }

  /*** set prediction and the error ***/
  A->prediction = energy;
  A->error = A->target - A->prediction;
  A->rel_error = A->error / A->target;

}
void Elementsnn2::train_pose(pose * A)
{
  vector<float> error(1,A->error);
  for(unsigned int i=0; i<A->inputs.size();i++)
    network.train_ext_error(A->inputs[i],error);

  /*** print current network ***/
  //printf("%s\n",network.print_network().c_str());

}

void Elementsnn2::update_weights()
{
  network.update_weights();
}

vector<float> Elementsnn2::get_protein_freqs(int p)
{

  vector<float> temp_types, freqs;
  freqs.assign(include_elements_p.size(),0);

  //find bonded atoms
  for(unsigned int i=0; i<protein_bonds[p].size(); i++)
    {
      if(protein_bonds[p][i]>0)
	if(count(ignore_elements_p.begin(), ignore_elements_p.end(),protein_atoms[i]->element) < 1)
	  {
	    temp_types = protein_types[protein_atoms[i]->element.c_str()];
	    for(unsigned int j=0; j<freqs.size();j++)
		freqs[j]+=temp_types[j];

	  }
    }

  //normalise frequencies
  float total =0;
  for(unsigned int j=0; j<freqs.size();j++)
    total += freqs[j];


  if(total > 0.1)
    {
      for(unsigned int j=0; j<freqs.size();j++)
	freqs[j] = freqs[j]/total;
    }
  else
    printf("WARNING: Protein atom %d-%s has been found to have no bonds in elementsNN2\n",p+1,protein_atoms[p]->element.c_str());

  return freqs;
}


vector<float> Elementsnn2::get_ligand_freqs(int p)
{

  vector<float> temp_types, freqs;
  freqs.assign(include_elements_l.size(),0);
 
  //find bonded atoms
  for(unsigned int i=0; i<ligand_bonds[p].size(); i++)
    {
      if(ligand_bonds[p][i]>0)
	if(count(ignore_elements_l.begin(), ignore_elements_l.end(),ligand_atoms[i]->element) < 1)
	  {
	    temp_types = ligand_types[ligand_atoms[i]->element.c_str()];
	    for(unsigned int j=0; j<freqs.size();j++)
	      freqs[j]+=temp_types[j];
	  }
    }

  //normalise frequencies
  float total =0;
  for(unsigned int j=0; j<freqs.size();j++)
    total += freqs[j];

  if(total > 0.1)
    {
      for(unsigned int j=0; j<freqs.size();j++)
	freqs[j] = freqs[j]/total;
    }
  else
    printf("WARNING: Ligand atom %d-%s has been found to have no bonds in elementsNN2\n",p+1,ligand_atoms[p]->element.c_str());
  
  return freqs;
}



void Elementsnn2::write_network(string name)
{
  /*** save network to file ***/
  network.write_file(name);
}


void Elementsnn2::read_network(string name)
{
  /*** save network to file ***/
  network.read_file(name);
  network.set_learning_rate(learning_rate);
  network.set_transfer_functions(tf);
  network.set_linear_coef(linear_coef);
  
}


void Elementsnn2::set_parameters()
{
  ifstream parameters;
  parameters.open("elementsnn2.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'elementsnn2.cfg' could not be found\n";
      exit(0);
    }

  printf("Setting parameters\n");

  string dummy;

  while (!parameters.eof())
    {

      parameters >> dummy;

      if(dummy == "DATASET:")
	parameters >> dataset;
      if(dummy == "NO_EPOCHS:")
	parameters >> no_epochs;
      if(dummy == "TARGET_FACTOR:")
	parameters >> target_factor;
      if(dummy == "MAX_EPOCHS_S_RECORD:")
	parameters >> max_epochs_since_record;

      if(dummy == "TRAIN_PCT:")
	parameters >> train_pct;
      if(dummy == "VALIDATE_PCT:")
	parameters >> validate_pct;
      if(dummy == "TEST_PCT:")
	parameters >> test_pct;

      if(dummy == "LEARNING_RATE:")
	parameters >> learning_rate;

      if(dummy == "SMOOTHING:")
	parameters >> smoothing;


      if(dummy == "TOPOLOGY:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      topology.push_back(atoi(dummy.c_str()));
	      parameters >> dummy;
	    }
	}

      if(dummy == "TRANSER_FUNCTIONS:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      tf.push_back(atoi(dummy.c_str()));
	      parameters >> dummy;
	    }
	}

      if(dummy == "LINEAR_COEF:")
	parameters >> linear_coef;


      if(dummy == "THRESHOLD:")
	parameters >> threshold;

      if(dummy == "INCLUDE_ELEMENTS_P:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      include_elements_p.push_back(dummy);
	      parameters >> dummy;
	    }
	}

      if(dummy == "IGNORE_ELEMENTS_P:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      ignore_elements_p.push_back(dummy);
	      parameters >> dummy;
	    }
	}
      if(dummy == "INCLUDE_ELEMENTS_L:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      include_elements_l.push_back(dummy);
	      parameters >> dummy;
	    }
	}

      if(dummy == "IGNORE_ELEMENTS_L:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      ignore_elements_l.push_back(dummy);
	      parameters >> dummy;
	    }
	}
      if(dummy == "INCLUDE_SURF_DIST:")
	parameters >> include_surf_dist;

      if(dummy == "INCLUDE_HB_INFO:")
	parameters >> include_hb_info;


      if(dummy == "NO_INCLUDE_BONDS_P:")
	parameters >> no_protein_bonds_incl;
      if(dummy == "NO_INCLUDE_BONDS_L:")
	parameters >> no_ligand_bonds_incl;


      if(dummy == "INCLUDE_CO_FACTORS:")
	parameters >> include_co_factors;
      if(dummy == "INCLUDE_WATER:")
	parameters >> include_water;

    } //end read file

  //printing all configurations

  printf("Configurations for Elements NN:\n");
  printf("\n");
  printf("\n");

  printf("DATASET:              %s\n",dataset.c_str());
  printf("NO_EPOCHS:            %d\n",no_epochs);
  printf("TARGET_FACTOR:        %f\n",target_factor);
  printf("MAX_EPOCHS_S_RECORD:  %d\n",max_epochs_since_record);
  printf("\n");

  printf("TRAIN_PCT:            %.2f\n",train_pct);
  printf("VALIDATE_PCT:         %.2f\n",validate_pct);
  printf("TEST_PCT:             %.2f\n",test_pct);
  printf("\n");

  printf("LEARNING_RATE:        %.10f\n",learning_rate);
  printf("SMOOTHING:            %d\n",smoothing);

  printf("LINEAR_COEF:          %f\n",linear_coef);
  printf("\n");

  printf("TOPOLOGY:             ");
  for(unsigned int i=0; i<topology.size(); i++)
    printf("%d ",topology[i]);
  printf("end\n");

  printf("TRANSER_FUNCTIONS:    ");
  for(unsigned int i=0; i<tf.size(); i++)
    printf("%d ",tf[i]);
  printf("end\n");


  printf("THRESHOLD:            %.2f\n",threshold);
  printf("\n");
  printf("INCLUDE_ELEMENTS_P:   ");
  for(unsigned int i=0; i<include_elements_p.size(); i++)
    printf("%s ",include_elements_p[i].c_str());
  printf("end\n");

  printf("IGNORE_ELEMENTS_P:    ");
  for(unsigned int i=0; i<ignore_elements_p.size(); i++)
    printf("%s ",ignore_elements_p[i].c_str());
  printf("end\n");
  printf("\n");

  printf("INCLUDE_ELEMENTS_L:   ");
  for(unsigned int i=0; i<include_elements_l.size(); i++)
    printf("%s ",include_elements_l[i].c_str());
  printf("end\n");

  printf("IGNORE_ELEMENTS_L:    ");
  for(unsigned int i=0; i<ignore_elements_l.size(); i++)
    printf("%s ",ignore_elements_l[i].c_str());
  printf("end\n");
  printf("\n");

  printf("INCLUDE_SURF_DIST:    %d\n",include_surf_dist);
  printf("INCLUDE_HB_INFO:      %d\n",include_hb_info);

  printf("NO_INCLUDE_BONDS_P:   %d\n",no_protein_bonds_incl);
  printf("NO_INCLUDE_BONDS_L:   %d\n",no_ligand_bonds_incl);
  printf("\n");

  printf("INCLUDE_CO_FACTORS:   %d\n",include_co_factors);
  printf("INCLUDE_WATER:        %d\n",include_water);
  printf("\n");

  printf("Version: $Id: elementsnn2.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");


  //finding number of inputs
  no_protein_types = include_elements_p.size();
  no_ligand_types = include_elements_l.size();
  no_inputs = 
    no_protein_types * (1+no_protein_bonds_incl) + 
    no_ligand_types  * (1+no_ligand_bonds_incl) +
    1;

  if(include_surf_dist == 1)
    no_inputs++;

  if(include_hb_info == 1)
    no_inputs++;


  //correcting topology 0 -> no_inputs
  for(unsigned int i=0; i<topology.size(); i++)
    if(topology[i] == 0)
      topology[i] = no_inputs;


}

int Elementsnn2::get_no_epochs(){return no_epochs;}
int Elementsnn2::get_smoothing(){return smoothing;}
float Elementsnn2::get_target_factor(){return target_factor;}
int Elementsnn2::get_max_epochs_since_record(){return max_epochs_since_record;}
string Elementsnn2::get_dataset(){return dataset;}

float Elementsnn2::get_train_pct(){return train_pct;}
float Elementsnn2::get_validate_pct(){return validate_pct;}
float Elementsnn2::get_test_pct(){return test_pct;}

bool Elementsnn2::get_include_co_factors()
{
  if(include_co_factors == 1)
    return true;
  else
    return false;
}
bool Elementsnn2::get_include_water()
{
  if(include_water == 1)
    return true;
  else
    return false;
}


void Elementsnn2::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: elementsnn2.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  rep = reporter;

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
