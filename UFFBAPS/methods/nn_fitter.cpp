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
$Id: nn_fitter.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "nn_fitter.h"


NN_FITTER::NN_FITTER(FILE * resultsFile,FILE * reporter)
{

  /*** write out to out file **/
  writeDescription(reporter,resultsFile);
 
  printf("Using the NN_FITTER network!\n");

  set_parameters();
  
  build_inputs();
  set_targets(dataset);
  organise();
}

NN_FITTER::~NN_FITTER(){}


void NN_FITTER::create_network()
{

  //set up network
  network.create_network(topology);
  network.set_learning_rate(learning_rate);
  network.set_transfer_functions(tf);
  network.set_linear_coef(linear_coef);
  printf("Network:\n %s\n",network.print_network().c_str());
}


void NN_FITTER::build_inputs()
{
  
  ifstream input_lines;
  input_lines.open("input_list.txt");
  vector<string> words;
  vector<float> values;
  string dummy;

  while (!input_lines.eof())
    {
    
      getline(input_lines,dummy);
      string checkComment(dummy,0,1);
      if(checkComment != "#")
	{
	  /*** break line into words ***/
	  words.clear();
	  char temp[300];
	  istringstream stream(dummy.c_str());
	  while(stream>>temp)
	    words.push_back(removeWhite(string(temp)));
		    
	  /*** check for lines with only white space ***/
	  if(words.size() > 0)
	    {
	  
	      values.clear();
	      vector<string>::iterator iterator = words.begin();
	      iterator++;
	      while(iterator != words.end())
		{
		  values.push_back(atof((*iterator).c_str()));
		  iterator++;
		}
	      
	      inputs[words[0]] = values;
	    
	      vector<string>::iterator its;
	      for(its=include_contacts.begin();its<include_contacts.end();its++)
		inputs[words[0]].push_back(get_contacts(words[0],*its));
	    }
	}
    }
  
  cout << "printing inputs\n";
  for( map<string, vector<float> >::iterator iter = inputs.begin(); iter != inputs.end(); iter++ ) 
    {
      cout <<(*iter).first << " " ;
      for(unsigned int i=0;i<(*iter).second.size();i++)
 	cout<<setw(10)<< (*iter).second[i] << " ";
      cout << endl;   
    }
  
  
}

float NN_FITTER::get_contacts(string pdb, string dist)
{
  ifstream input_lines;
  
  string temp = "contacts";
  temp.append(dist.c_str());
  temp.append("A.txt");
  //  cout<<"name: "<<temp<<endl;
  input_lines.open(temp.c_str());
  vector<string> words;
  string dummy;
  float res=0.0;

  while (!input_lines.eof())
    {
      getline(input_lines,dummy);
      string cur_pdb(dummy,0,4);
      if(cur_pdb == pdb)
	{
	  /*** break line into words ***/
	  words.clear();
	  char temp[300];
	  istringstream stream(dummy.c_str());
	  while(stream>>temp)
	    words.push_back(removeWhite(string(temp)));
	
	  res = atof(words[2].c_str());
	  
	}
    }

  return res;
}
void NN_FITTER::set_targets(string dataset)
{
  Targets target(rep, dataset);
  map<string,vector<float> >::iterator it;

  for(it = inputs.begin();it!=inputs.end();it++)
    {
      targets[(*it).first] = target_factor * target.get_G((*it).first);
      //     printf("Target of '%s' is %f\n",
      //	     (*it).first.c_str(),
      //	     targets[(*it).first]);
    }
}


void NN_FITTER::organise()
{

  int no_test_cases = (int) ceil(inputs.size()*(test_pct)/100.0);
  int no_train_cases = (int) (inputs.size()*(train_pct)/100.0);
  int no_cross_validations = 1;
  if(test_pct != 0)
    no_cross_validations = (int) ceil(100.0/test_pct);

  average_test_error = 0;

  map<string,vector<float> >::iterator it;
  /*** for each cross validation ***/
  for(int i=0;i<no_cross_validations;i++)
    {
      /*** organise lists ***/
      int no_test=0,no_train=0;
      test_names.clear();
      test_set.clear();
      test_set_targets.clear();
      training_set.clear();
      training_set_targets.clear();
      validation_set.clear();
      validation_set_targets.clear();

      for(it = inputs.begin();it!=inputs.end();it++)
	{
	  vector<float> target(1,targets[(*it).first]);
	  
	  //test 
	  if (no_test < no_test_cases && count(tested_names.begin(),tested_names.end(),(*it).first)==0)
	    {
	      //	      cout<<"Adding to test set: "<<(*it).first<<endl;
	      test_set.push_back((*it).second);
	      test_set_targets.push_back(target);
	      test_names.push_back((*it).first);
	      tested_names.push_back((*it).first);
	      no_test++;
	    }
	  //training
	  else if (no_train < no_train_cases)
	    {
	      //	      cout<<"Adding to training set: "<<(*it).first<<endl;
	      training_set.push_back((*it).second);
	      training_set_targets.push_back(target);
	      no_train++;
	    }
	  //validation
	  else
	    {
	      //      cout<<"Adding to validation set: "<<(*it).first<<endl;
	      validation_set.push_back((*it).second);
	      validation_set_targets.push_back(target);
	    }
	  
	}
 
      cout<<"no_test_cases:  "<<no_test_cases<<endl;
      cout<<"no_train_cases: "<<no_train_cases<<endl;
      cout<<"In test set:    "<<test_set.size()<<endl;
      cout<<"In training set:"<<training_set.size()<<endl;
      cout<<"In validati set:"<<validation_set.size()<<endl;
      cout<<"In all sets:    "<<validation_set.size()+training_set.size()+test_set.size()<<endl;
      
      go();
    }
  average_test_error = average_test_error/((float) no_cross_validations);
  cout<<"Overall AUE on test: "<<average_test_error<<endl;

}


void NN_FITTER::go()
{


  bool converged = false;
  float record = 100000;
  float aue_v,aue_tr,aue_te;
  int epochs_since_record = 0,no_epoch=1;

  create_network();

  while (!converged)
    {
      aue_tr = train();
      aue_v = validate();
      if (aue_v<record)
	{
	  record=aue_v;
	  epochs_since_record = 0;
	  // cout<<"Record broken: "<<aue_v<<endl;
	  //printf("%s\n",network.print_network().c_str());
	  write_network("record.net");

	}

      /*** check for conversion ***/
      if(max_epochs_since_record<=epochs_since_record)
	converged=true;
      if (no_epochs<=no_epoch)
	converged=true;
      
      if(no_epoch%10==0)
	printf("At epoch %d, the error is %f (validation) and %f (training)\n",no_epoch,aue_v,aue_tr);
      no_epoch++;
      epochs_since_record++;
    }
  
  read_network("record.net");
  printf("%s\n",network.print_network().c_str());
  

  aue_te = test();
  average_test_error += aue_te;
  printf("Last training error is:     %f\n",aue_tr);
  printf("Record validation error is: %f\n",record);
  printf("Avr testing error is:       %f\n",aue_te);

}


float NN_FITTER::train()
{
 
  float aue = network.train(training_set,training_set_targets);
 
  /*** print current network ***/
  //  printf("%s\n",network.print_network().c_str());

  return aue/fabs(target_factor);
}

float NN_FITTER::validate()
{
  float aue = 0.0;
  int no = 0;
  vector<vector<float> >::iterator v_input  = validation_set.begin();
  vector<vector<float> >::iterator v_target = validation_set_targets.begin();
  
  while(v_input<validation_set.end() && v_target<validation_set_targets.end())
    {
      aue += fabs(network.test(*v_input).at(0)-((*v_target).at(0)));
      v_input++;
      v_target++;
      no++;
    }

  return aue/fabs(((float) no)*target_factor);
}


float NN_FITTER::test()
{
  float aue = 0.0,prediction;
  int no = 0;
  vector<vector<float> >::iterator t_input  = test_set.begin();
  vector<vector<float> >::iterator t_target = test_set_targets.begin();
  vector<string>::iterator t_name = test_names.begin();

  //  printf("Testing:\n");
  //  printf("-------------------------------------------------\n");
  //  printf("Prediction   Target        Error\n");
  //  printf("-------------------------------------------------\n");

  while(t_input<test_set.end() && t_target<test_set_targets.end())
    {
      prediction = network.test(*t_input).at(0);
      aue += fabs(prediction-(*t_target).at(0));
      //  printf("%5.2f        %5.2f        %5.2f\n",
      //	     prediction/target_factor,
      //	     (*t_target).at(0)/target_factor,
      //	     (prediction-(*t_target).at(0))/target_factor);

      fprintf(results,"%-4s %-4s %8f\n",
	      (*t_name).c_str(),
	      "-1",
	      prediction/target_factor);


      t_input++;
      t_target++;
      t_name++;
      no++;
    }
  printf("-------------------------------------------------\n");
  return aue/fabs(((float) no)*target_factor);
}



void NN_FITTER::write_network(string name)
{
  /*** save network to file ***/
  network.write_file(name);
}


void NN_FITTER::read_network(string name)
{
  printf("TFsize(read): %d\n",tf.size());
  /*** read network from file ***/
  network.read_file(name);
  printf("Network:\n %s\n",network.print_network().c_str());

  network.set_learning_rate(learning_rate);
  network.set_transfer_functions(tf);
  network.set_linear_coef(linear_coef);
}


void NN_FITTER::set_parameters()
{


  ifstream parameters;
  parameters.open("nn_fitter.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'nn_fitter.cfg' could not be found\n";
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
      
      if(dummy == "NO_INPUTS:")
	parameters >> no_inputs;
      
      if(dummy == "INCLUDE_CONTACTS:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      include_contacts.push_back(dummy);
	      parameters >> dummy;
	    }
	}


    } //end read file

  //printing all configurations

  printf("Configurations for fitter NN:\n");
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

  printf("INCLUDE_CONTACTS:    ");
  for(unsigned int i=0; i<include_contacts.size(); i++)
    printf("%s ",include_contacts[i].c_str());
  printf("end\n");
  printf("\n");

  printf("NO_INPUTS:            %d\n",no_inputs);
  printf("\n");
  printf("\n");
  printf("Version: $Id: nn_fitter.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");

  //correcting topology 0 -> no_inputs
  for(unsigned int i=0; i<topology.size(); i++)
    if(topology[i] == 0)
      topology[i] = no_inputs;


  /*** reading input file ***/
  string line;

  ifstream file("input_list.txt");

  if(!file.is_open())
    {
      printf("File input_list.txt not found!\n");
      exit(0);
    }

  while(!file.eof())
    {
      getline(file, line);
      input_lines.push_back(line);
    }


}


void NN_FITTER::writeDescription(FILE * reporter, FILE * resultsFile)
{

  rep=reporter;
  results=resultsFile;
  /*** write out to out file **/
  version =string("$Id: nn_fitter.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
