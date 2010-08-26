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
$Id: elementsnn.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "elementsnn.h"


ElementsNN::ElementsNN(FILE * reporter):Nnmodel()
{

  /*** write out to out file **/
  writeDescription(reporter);
 
  printf("Using the elements network!\n");

  set_parameters();
 

}

ElementsNN::~ElementsNN(){}


void ElementsNN::create_network()
{
  vector<int> topology;
  topology.push_back(no_inputs);
  topology.push_back(1);


  network.create_network(topology);
  network.set_learning_rate(learning_rate);
  printf("%s\n",network.print_network().c_str());
}


void ElementsNN::build_inputs(pose * A)
{

  vector<Atom *> protein_atoms = convert_objects_to_atoms(A->protein);
  vector<Atom *> ligand_atoms = convert_objects_to_atoms(A->ligand);
  vector<float> res;
  for(unsigned int p=0; p<protein_atoms.size(); p++)
    for(unsigned int l=0; l<ligand_atoms.size(); l++)
      {

	//input for distance
	float dist = sqrt((protein_atoms[p]->x - ligand_atoms[l]->x)*(protein_atoms[p]->x - ligand_atoms[l]->x)+
			  (protein_atoms[p]->y - ligand_atoms[l]->y)*(protein_atoms[p]->y - ligand_atoms[l]->y)+
			  (protein_atoms[p]->z - ligand_atoms[l]->z)*(protein_atoms[p]->z - ligand_atoms[l]->z));
   
	if(dist < threshold)
	  if(count(ignore_elements.begin(), ignore_elements.end(),protein_atoms[p]->element) < 1)
	    if(count(ignore_elements.begin(), ignore_elements.end(),ligand_atoms[l]->element) < 1)
	      {
		res.clear();
	    
		//input for protein atom
		res = types[protein_atoms[p]->element.c_str()];
	    
		//input for ligand atom
		res.insert(res.end(), types[ligand_atoms[l]->element.c_str()].begin(), types[ligand_atoms[l]->element.c_str()].end()); 
	      
		res.insert(res.end(),dist);
	      
		//input for context
		//res.insert(res.end(),); 
	      
		if(res.size() == (unsigned int) no_inputs)
		  A->inputs.push_back(res);
		else
		  printf("Warning: could not create input for atom pair '%s'-'%s'\n",
			 protein_atoms[p]->element.c_str(),ligand_atoms[l]->element.c_str());

	    }
      } 


}

void ElementsNN::define_input()
{

  int one_index = 0;
  no_types = include_elements.size();
  no_inputs = 2*no_types + 1 + 0;  

  vector<float> temp;
  //map vectors
  for(int i=0; i<no_types; i++)
    {
      for(int e=0; e<no_types; e++)
	{
	  if(e == one_index)
	    temp.push_back(1);
	  else
	    temp.push_back(0);
	}
      types[include_elements[i]] =  temp;
      temp.clear();
      one_index++;
    }

  //print out mapping
  for(int i=0; i<no_types; i++)
    {
      printf("Element '%2s' has been mapped to input ",include_elements[i].c_str());
      for(int e=0; e<no_types; e++)
	printf("'%.0f'",types[include_elements[i]][e]);
      cout<<"\n";
    }

  printf("IGNORE_ELEMENTS SIZE: %d\n",ignore_elements.size());
 
  printf("Ignore elements: ");
  for(unsigned int i=0; i<ignore_elements.size(); i++)
    printf("'%s' ", ignore_elements[i].c_str());
  printf("\n");
}




void ElementsNN::test_pose(pose * A)
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
void ElementsNN::train_pose(pose * A)
{
  vector<float> error(1,A->error);
  for(unsigned int i=0; i<A->inputs.size();i++)
    network.train_ext_error(A->inputs[i],error);

  /*** print current network ***/
  // printf("%s\n",network.print_network().c_str());

}



void ElementsNN::write_network(string name)
{
  /*** save network to file ***/
  network.write_file(name);
}


void ElementsNN::read_network(string name)
{
  /*** save network to file ***/
  network.read_file(name);
  network.set_learning_rate(learning_rate);

}


void ElementsNN::set_parameters()
{
  ifstream parameters;
  parameters.open("elementsnn.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'elementsnn.cfg' could not be found\n";
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


      if(dummy == "TOPOLOGY:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      topology.push_back(atoi(dummy.c_str()));
	      parameters >> dummy;
	    }
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
  printf("\n");

  printf("TOPOLOGY:             ");
  for(unsigned int i=0; i<topology.size(); i++)
    printf("%d ",topology[i]);
  printf("end\n");

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

  printf("Version: $Id: elementsnn.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");

  //correcting topology 0 -> no_inputs
  for(unsigned int i=0; i<topology.size(); i++)
    if(topology[i] == 0)
      topology[i] = no_inputs;


}

int ElementsNN::get_no_epochs(){return no_epochs;}
float ElementsNN::get_target_factor(){return target_factor;}
int ElementsNN::get_max_epochs_since_record(){return max_epochs_since_record;}
string ElementsNN::get_dataset(){return dataset;}

float ElementsNN::get_train_pct(){return train_pct;}
float ElementsNN::get_validate_pct(){return validate_pct;}
float ElementsNN::get_test_pct(){return test_pct;}

bool ElementsNN::get_include_co_factors()
{
  if(include_co_factors == 1)
    return true;
  else
    return false;
}
bool ElementsNN::get_include_water()
{
  if(include_water == 1)
    return true;
  else
    return false;
}


void ElementsNN::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: elementsnn.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
