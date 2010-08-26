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
$Id: nn_ls.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "nn_ls.h"


NN_LS::NN_LS(FILE * reporter):Nnmodel()
{

  /*** write out to out file **/
  writeDescription(reporter);
 
  printf("Using the NN_LS network!\n");

  set_parameters();
}

NN_LS::~NN_LS(){}


void NN_LS::create_network()
{

  //set up network
  network.create_network(topology);
  network.set_learning_rate(learning_rate);
  network.set_transfer_functions(tf);
  network.set_linear_coef(linear_coef);
  printf("Network:\n %s\n",network.print_network().c_str());
}


void NN_LS::build_inputs(pose * A)
{

  vector<Atom *> protein_atoms = convert_objects_to_atoms(A->protein);
  vector<Atom *> ligand_atoms = convert_objects_to_atoms(A->ligand);
  vector<float> res;



  LigandEntropy LS(rep);
  float no_rot_bonds = 0.0;/////////////////////////////////////////// NOT USED ANYMORE !!!!!!!!!!!! LS.no_rot_bonds(ligand_atoms);

  res.push_back((float) protein_atoms.size());
  res.push_back((float) ligand_atoms.size());
  res.push_back(no_rot_bonds);


  if(res.size() == (unsigned int)no_inputs)
    {	    
      A->inputs.push_back(res);
      
      printf("Inputs: ");
      for(unsigned int i=0; i< A->inputs[A->inputs.size()-1].size(); i++)
	printf("'%f' ", A->inputs[A->inputs.size()-1][i]);
		
		    			
    }
  else
    printf("Warning: could not create input for '%s'\n",A->protein[0]->identify().c_str());
  
}

void NN_LS::define_input(){}


void NN_LS::test_pose(pose * A)
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

void NN_LS::train_pose(pose * A)
{
  vector<float> error(1,A->error);
  for(unsigned int i=0; i<A->inputs.size();i++)
    {
      network.train_ext_error(A->inputs[i],error);
    }

  /*** print current network ***/
  // printf("%s\n",network.print_network().c_str());

}



void NN_LS::write_network(string name)
{
  /*** save network to file ***/
  network.write_file(name);
}


void NN_LS::read_network(string name)
{
  /*** read network from file ***/
  network.read_file(name);
  network.set_learning_rate(learning_rate);
  network.set_transfer_functions(tf);
  network.set_linear_coef(linear_coef);
}


void NN_LS::set_parameters()
{
  no_inputs =3;

  ifstream parameters;
  parameters.open("nn_ls.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'nn_ls.cfg' could not be found\n";
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


    } //end read file

  //printing all configurations

  printf("Configurations for es NN:\n");
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
  printf("\n");

  printf("Version: $Id: nn_ls.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");

  //correcting topology 0 -> no_inputs
  for(unsigned int i=0; i<topology.size(); i++)
    if(topology[i] == 0)
      topology[i] = no_inputs;

}

int NN_LS::get_no_epochs(){return no_epochs;}
float NN_LS::get_target_factor(){return target_factor;}
int NN_LS::get_max_epochs_since_record(){return max_epochs_since_record;}
string NN_LS::get_dataset(){return dataset;}

float NN_LS::get_train_pct(){return train_pct;}
float NN_LS::get_validate_pct(){return validate_pct;}
float NN_LS::get_test_pct(){return test_pct;}


void NN_LS::writeDescription(FILE * reporter)
{

  rep=reporter;
  /*** write out to out file **/
  version =string("$Id: nn_ls.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
