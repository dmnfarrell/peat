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
$Id: distancenn.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "distancenn.h"


DistanceNN::DistanceNN(FILE * reporter):Nnmodel()
{

  /*** write out to out file **/
  writeDescription(reporter);
 
  printf("Using the distance network!\n");

  set_parameters();
}

DistanceNN::~DistanceNN(){}


void DistanceNN::create_network()
{

  //instantiate networks
  NN temp1;
  vector<NN> temp2;
  for(int i=0; i<no_types; i++)
    temp2.push_back(temp1);
  for(int i=0; i<no_types; i++)
    networks.push_back(temp2);

  //set up networks
  for(int i=0; i<no_types; i++)  
    for(int j=0; j<no_types; j++)
      { 
	networks[i][j].create_network(topology);
	networks[i][j].set_learning_rate(learning_rate);

	networks[i][j].set_transfer_functions(tf);
	networks[i][j].set_linear_coef(linear_coef);
	printf("[%d][%d]:\n %s\n",i,j,networks[i][j].print_network().c_str());
      }
}


void DistanceNN::build_inputs(pose * A)
{

  vector<Atom *> protein_atoms = convert_objects_to_atoms(A->protein);
  vector<Atom *> ligand_atoms = convert_objects_to_atoms(A->ligand);
  vector<float> res;
  vector<int> index;
  for(unsigned int p=0; p<protein_atoms.size(); p++)
    for(unsigned int l=0; l<ligand_atoms.size(); l++)
      {

	//input for distance
	float dist = sqrt((protein_atoms[p]->x - ligand_atoms[l]->x)*(protein_atoms[p]->x - ligand_atoms[l]->x)+
			  (protein_atoms[p]->y - ligand_atoms[l]->y)*(protein_atoms[p]->y - ligand_atoms[l]->y)+
			  (protein_atoms[p]->z - ligand_atoms[l]->z)*(protein_atoms[p]->z - ligand_atoms[l]->z));
   
	if(dist < threshold)
	  //check if elements should be ignored 
	  if(count(ignore_elements.begin(), ignore_elements.end(),protein_atoms[p]->element) < 1)
	    if(count(ignore_elements.begin(), ignore_elements.end(),ligand_atoms[l]->element) < 1)
	      {
		res.clear();
		index.clear();
		//temp		
		if(dist < 1)
		  printf("Warning: bumb registered: '%s' and '%s'",protein_atoms[p]->name.c_str(),ligand_atoms[l]->name.c_str() );

		//temp end
		float invD1 = 1/dist;

		for(unsigned int i=0; i<powers.size(); i++)
		  res.push_back(pow(invD1, powers[i]));		  
				
		if(element_indexes.count(protein_atoms[p]->element.c_str())>0)
		  index.push_back(element_indexes[protein_atoms[p]->element.c_str()]);
		if(element_indexes.count(ligand_atoms[l]->element.c_str())>0)
		  index.push_back(element_indexes[ligand_atoms[l]->element.c_str()]);
		

		if(res.size() == (unsigned int) no_inputs && index.size() == 2)
		  {	    
		    A->inputs.push_back(res);
		    A->indexes.push_back(index);
		  }
		else
		  printf("Warning: could not create input for atom pair '%s'-'%s'\n",
			 protein_atoms[p]->element.c_str(),ligand_atoms[l]->element.c_str());
		
	      }
      } 
}

void DistanceNN::define_input()
{

  int index=0;
  no_types=include_elements.size();

  //map vectors
  for(unsigned int i=0; i<include_elements.size(); i++)
    {
      element_indexes[include_elements[i]] = index;
      index++;
    }

  //print out mapping
  for(unsigned int i=0; i<include_elements.size(); i++)
      cout<<"Element '"<<include_elements[i]<<"' has been mapped to index '"<<element_indexes[include_elements[i]]<<"'\n";


}




void DistanceNN::test_pose(pose * A)
{
  float energy = 0;
  /*** sum up energy contributions from all atom pairs ***/
  for(unsigned int i=0; i < A->inputs.size(); i++)
    {
      vector<float> out = networks[A->indexes[i][0]][A->indexes[i][1]].test(A->inputs[i]);
      energy += out[0];
    }

  /*** set prediction and the error ***/
  A->prediction = energy;
  A->error = A->target - A->prediction;
  A->rel_error = A->error / A->target;

}
void DistanceNN::train_pose(pose * A)
{
  vector<float> error(1,A->error);
  for(unsigned int i=0; i<A->inputs.size();i++)
    {

      networks[A->indexes[i][0]][A->indexes[i][1]].train_ext_error(A->inputs[i],error);

    }

  /*** print current network ***/
  // printf("%s\n",network.print_network().c_str());

}



void DistanceNN::write_network(string name)
{
  /*** save network to file ***/
 for(int i=0; i<no_types; i++)  
    for(int j=0; j<no_types; j++)
      {  
	ostringstream si,sj;
	si << i;
	sj << j;
	networks[i][j].write_file(name+si.str()+"_"+sj.str());
      }
}


void DistanceNN::read_network(string name)
{
  /*** save network to file ***/
 for(int i=0; i<no_types; i++)  
    for(int j=0; j<no_types; j++)
      {
	ostringstream si,sj;
	si << i;
	sj << j;
	networks[i][j].read_file(name+si.str()+"_"+sj.str());
	networks[i][j].set_learning_rate(learning_rate);
	networks[i][j].set_transfer_functions(tf);
	networks[i][j].set_linear_coef(linear_coef);


      }

}


void DistanceNN::set_parameters()
{
  ifstream parameters;
  parameters.open("distancenn.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'distancenn.cfg' could not be found\n";
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

      if(dummy == "LINEAR_COEF:")
	parameters >> linear_coef;

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

      if(dummy == "INPUT_POWERS:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      powers.push_back(atoi(dummy.c_str()));
	      parameters >> dummy;	  
	    }
	  no_inputs = powers.size();
	}
      
      if(dummy == "THRESHOLD:")
	parameters >> threshold;

      if(dummy == "INCLUDE_ELEMENTS:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      include_elements.push_back(dummy);
	      parameters >> dummy;
	    }
	}

      if(dummy == "IGNORE_ELEMENTS:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      ignore_elements.push_back(dummy);
	      parameters >> dummy;
	    }
	}

      if(dummy == "INCLUDE_CO_FACTORS:")
	parameters >> include_co_factors;
      if(dummy == "INCLUDE_WATER:")
	parameters >> include_water;

    } //end read file

  //printing all configurations

  printf("Configurations for Distance NN:\n");
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

  printf("INPUT_POWERS:         ");
  for(unsigned int i=0; i<powers.size(); i++)
    printf("%d ",powers[i]);
  printf("end\n");

  printf("\n");
  printf("THRESHOLD:            %.2f\n",threshold);

  printf("INCLUDE_ELEMENTS:     ");
  for(unsigned int i=0; i<include_elements.size(); i++)
    printf("%s ",include_elements[i].c_str());
  printf("end\n");

  printf("IGNORE_ELEMENTS:      ");
  for(unsigned int i=0; i<ignore_elements.size(); i++)
    printf("%s ",ignore_elements[i].c_str());
  printf("end\n");
  printf("\n");

  printf("INCLUDE_CO_FACTORS:   %d\n",include_co_factors);
  printf("INCLUDE_WATER:        %d\n",include_water);
  printf("\n");

  printf("Version: $Id: distancenn.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");

  //correcting topology 0 -> no_inputs
  for(unsigned int i=0; i<topology.size(); i++)
    if(topology[i] == 0)
      topology[i] = no_inputs;

}

int DistanceNN::get_no_epochs(){return no_epochs;}
float DistanceNN::get_target_factor(){return target_factor;}
int DistanceNN::get_max_epochs_since_record(){return max_epochs_since_record;}
string DistanceNN::get_dataset(){return dataset;}

float DistanceNN::get_train_pct(){return train_pct;}
float DistanceNN::get_validate_pct(){return validate_pct;}
float DistanceNN::get_test_pct(){return test_pct;}

bool DistanceNN::get_include_co_factors()
{
  if(include_co_factors == 1)
    return true;
  else
    return false;
}
bool DistanceNN::get_include_water()
{
  if(include_water == 1)
    return true;
  else
    return false;
}



void DistanceNN::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: distancenn.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
