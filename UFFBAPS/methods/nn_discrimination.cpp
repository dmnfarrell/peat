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
$Id: nn_discrimination.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "nn_discrimination.h"


NN_DISCRIMINATION::NN_DISCRIMINATION(FILE * resultsFile,FILE * reporter)
{

  /*** write out to out file **/
  writeDescription(reporter,resultsFile);
 
  printf("Using the NN_DISCRIMINATION network!\n");
  
  target_factor=1.0;

  build_inputs();
  set_parameters();
  
  error_file.open("error.gnu");

  organise();
}

NN_DISCRIMINATION::~NN_DISCRIMINATION(){}


void NN_DISCRIMINATION::create_network()
{
  printf("TFsize(create): %d\n",tf.size());
  //set up network
  network.create_network(topology);
  network.set_learning_rate(learning_rate);
  network.set_transfer_functions(tf);
  // network.set_linear_coef(linear_coef);
  printf("Network:\n %s\n",network.print_network().c_str());
}


void NN_DISCRIMINATION::build_inputs()
{
  
  ifstream input_lines;
  input_lines.open("input_list.txt");
  vector<string> words;
  vector<float> values;
  string dummy;

  // read inputs from inputs.txt file

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
	      iterator++;
	      iterator++;
	      
	      while(iterator != words.end())
		{
		  values.push_back(atof((*iterator).c_str()));
		  iterator++;
		}
	      // save inputs
	      if (inputs.count(words[0])==0)
		{
		  map<string,vector<float>,DecoyCmp > temp;
		  inputs[words[0]] = temp;
		}

	      inputs[words[0]][words[1]] = values;
	      
	      // save rmsd
	      if (rmsds.count(words[0]) == 0)
		{
		  map<string,float> temp;
		  rmsds[words[0]] = temp;
		}
	      rmsds[words[0]][words[1]] = atof(words[2].c_str());

	    }
	}
    }


  map<string, map<string, vector<float>,DecoyCmp > >::iterator iter = inputs.begin();
  map<string, vector<float>,DecoyCmp >::iterator iter2 = (*iter).second.begin();
  int no_terms = (*iter2).second.size();
  no_inputs = 2*no_terms;

  // Normalize inputs
  


  cout << "printing inputs\n";


  for(iter=inputs.begin(); iter != inputs.end(); iter++ ) //elements of first map
    {
     
      /*
	(*iter).first:   pdb (string)
	(*iter).second:  map<decoy, inputs>
	(*iter2).first:  decoy (string)
	(*iter2).second: inputs (vector<float>)

      */

      vector<float> mins(no_terms, 1e7);
      vector<float> maxs(no_terms, -1e7);

      cout<<"Reading: "<<(*iter).first<<endl;

      for(iter2=(*iter).second.begin(); iter2 != (*iter).second.end(); iter2++ ) //elements of second map
	{
	  
	  for(unsigned int i=0;i<(*iter2).second.size();i++)
	    { 
	      
	      if (mins[i]>(*iter2).second[i])
		{
		  mins[i] = (*iter2).second[i];
		}
	      if (maxs[i]<(*iter2).second[i])
		{
		  maxs[i] = (*iter2).second[i];
		}
	      
	    }
	}


      /*      cout<<"mins: ";
      for(unsigned int i=0;i<mins.size();i++)
	cout<<mins[i]<<" ";
      cout<<endl;
      
      cout<<"maxs: ";
      for(unsigned int i=0;i<maxs.size();i++)
	cout<<maxs[i]<<" ";
      cout<<endl;
      
      */
      for(unsigned int i=0;i<mins.size();i++)
	if(maxs[i]-mins[i]<1e-7)
	  maxs[i]+=1.0;


      // normalising
      for(iter2=(*iter).second.begin(); iter2 != (*iter).second.end(); iter2++ ) //elements of second map
	for(unsigned int i=0;i<(*iter2).second.size();i++)
	  (*iter2).second[i] = ((*iter2).second[i]-mins[i])/(maxs[i]-mins[i]);    
  
      
    }
  

  // checking normalization
  for(iter=inputs.begin(); iter != inputs.end(); iter++ ) //elements of first map
    for(iter2=(*iter).second.begin(); iter2 != (*iter).second.end(); iter2++ ) //elements of second map
      {

	for(unsigned int i=0;i<(*iter2).second.size();i++)
	  if((*iter2).second[i]<0 || (*iter2).second[i]>1)
	    cout <<"WARNING: "<<(*iter).first<<" "<<(*iter2).first<<"("<<(*iter2).second.size()<<")"<<setw(10)<< (*iter2).second[i] << " is not normalised properly" << endl;   
      }
  
}


vector<float> NN_DISCRIMINATION::get_target(string pdb, string decoy1, string decoy2)
{
  vector<float> res(2,0);



  /****************/
  //res[0] = 1.0;

  
  if(rmsds[pdb][decoy1] < rmsds[pdb][decoy2])
    res[1] = 1.0;
  else
    res[0] = 1.0;
  
  return res;
}


void NN_DISCRIMINATION::organise()
{

  int no_test_cases = (int) ceil(inputs.size()*(test_pct)/100.0);
  int no_train_cases = (int) (inputs.size()*(train_pct)/100.0);
  int no_cross_validations = 1;
  if(test_pct != 0)
    no_cross_validations = (int) ceil(100.0/test_pct);

  average_test_error = 0;
  map<string,map<string, vector<float>,DecoyCmp > >::iterator it;

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
	  //test 
	  if (no_test < no_test_cases && count(tested_names.begin(),tested_names.end(),(*it).first)==0)
	    {
	      cout<<"Adding to test set: "<<(*it).first<<endl;
	      prepare_set(0, (*it).first);

	      tested_names.push_back((*it).first);
	      no_test++;
	    }
	  //training
	  else if (no_train < no_train_cases)
	    {
	      cout<<"Adding to training set: "<<(*it).first<<endl;
	      prepare_set(1, (*it).first);

	      no_train++;
	    }
	  //validation
	  else
	    {
	      cout<<"Adding to validation set: "<<(*it).first<<endl;
	      prepare_set(2, (*it).first);
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


void NN_DISCRIMINATION::prepare_set(int set, string pdb)
{

  map<string, vector<float> >::iterator it1, it2;
  vector<string> temp;

  for(it1 = inputs[pdb].begin();it1!=inputs[pdb].end();it1++)
    for(it2 = inputs[pdb].begin();it2!=inputs[pdb].end();it2++)
      if ((*it1).first != (*it2).first)
	{
	  //	  cout << "Preparing "<<(*it1).first<<", "<<(*it2).first<<" for "<<pdb<<endl;
	  
	  
	  vector<float> res = (*it1).second;
	  res.insert(res.end(), (*it2).second.begin(), (*it2).second.end());
	  
	  /*	  
	    cout<<"Inputs:";
	  for (int i=0;i<res.size();i++)
	    cout<<res[i]<<" ";
	  cout<<endl;
	  */

	  switch(set)
	    {
	    case 0:
	      test_set.push_back(res);
	      test_set_targets.push_back(get_target(pdb,(*it1).first,(*it2).first));
	      temp.clear();
	      temp.push_back(pdb);
	      temp.push_back((*it1).first);
	      temp.push_back((*it2).first);
	      test_names.push_back(temp);
	      break;
	    case 1:
	      training_set.push_back(res);
	      training_set_targets.push_back(get_target(pdb,(*it1).first,(*it2).first));
	      break;
	    case 2:
	      validation_set.push_back(res);
	      validation_set_targets.push_back(get_target(pdb,(*it1).first,(*it2).first));
	      break;
	    default:
	      cout<<"Unknown set: "<<set<<endl;
	      break;
	    }
	  
	  /*
	  vector<float> temp = get_target(pdb,(*it1).first,(*it2).first);
	  cout<<"targets:";
	  for (int i=0;i<temp.size();i++)
	    cout<<temp[i]<<" ";
	  cout<<endl;
	  */


	}


  return;

}





void NN_DISCRIMINATION::go()
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
	  cout<<"Record broken: "<<aue_v<<endl;
	  //	  printf("New record (check if match with read!)\n%s\n",network.print_network().c_str());
	  write_network("record.net");

	}

      /*** check for conversion ***/
      if(max_epochs_since_record<=epochs_since_record)
	converged=true;
      if (no_epochs<=no_epoch)
	converged=true;
      
      if(no_epoch%1==0)
	printf("At epoch %d, the error is %f (validation) and %f (training)\n",no_epoch,aue_v,aue_tr);


      error_file <<no_epoch<<"  "<<aue_tr<<"  "<<aue_v<<endl;
      error_file.flush();
      // printf("%s\n",network.print_network().c_str());
      no_epoch++;
      epochs_since_record++;
    }
  
  //  printf("Before retrieving:\n%s\n",network.print_network().c_str());
  read_network("record.net");
  //  printf("Record network retrieved:\n%s\n",network.print_network().c_str());
  
  
  aue_te = test();
  cout<<"After test"<<endl;
  average_test_error += aue_te;
  printf("Last training error is:     %f\n",aue_tr);
  printf("Record validation error is: %f\n",record);
  printf("Avr testing error is:       %f\n",aue_te);

}


float NN_DISCRIMINATION::train()
{
  vector<vector<float> >::iterator input_iter,target_iter;
  input_iter = training_set.begin();
  target_iter = training_set_targets.begin();


  float aue = 0.0;
  int smoothing = 1;
  int no = 0;

  while (input_iter != training_set.end() && target_iter != training_set_targets.end())
    {
      vector<vector<float> > input(input_iter, input_iter+smoothing), target(target_iter, target_iter+smoothing);
 
      aue+=network.train(input, target);
      network.update_weights();

      input_iter+=smoothing;
      target_iter+=smoothing;
      no++;

    }
  /*** print current network ***/
  //  printf("%s\n",network.print_network().c_str());
  return aue/(((float) no)*target_factor);
}

float NN_DISCRIMINATION::validate()
{
  float aue = 0.0;
  int no = 0;
  vector<float> res;


  vector<vector<float> >::iterator v_input  = validation_set.begin();
  vector<vector<float> >::iterator v_target = validation_set_targets.begin();
  
  while(v_input<validation_set.end() && v_target<validation_set_targets.end())
    {
      res = network.test(*v_input);
      for (unsigned int i=0;i<res.size();i++)
	{
	  aue += fabs(res[i]-(*v_target)[i]);
	  no++;
	}
      v_input++;
      v_target++;

    }

  return aue/(((float) no)*target_factor);
}


float NN_DISCRIMINATION::test()
{
  float aue = 0.0;
  vector<float> prediction;
  int no = 0;
  vector<vector<float> >::iterator t_input  = test_set.begin();
  vector<vector<float> >::iterator t_target = test_set_targets.begin();
  vector<vector<string> >::iterator t_name = test_names.begin();

  map<vector<string>, vector<float>,DecoyDecoyCmp > results;

  /*  printf("Testing:\n");
  printf("-------------------------------------------------\n");
  printf("Prediction            Target \n");
  printf("-------------------------------------------------\n");
  */
  while(t_input<test_set.end() && t_target<test_set_targets.end())
    {
      prediction = network.test(*t_input);

      for (unsigned int i=0;i<prediction.size();i++)
	{
	  aue += fabs(prediction[i]-(*t_target)[i]);
	  no++;
	}


      /*
      printf("Predictions: %s %4s %4s:  %5.2f %5.2f        %5.2f %5.2f\n",
      	     (*t_name)[0].c_str(),
	     (*t_name)[1].c_str(),
	     (*t_name)[2].c_str(),
	     prediction[0]/target_factor,
      	     prediction[1]/target_factor,
      	     (*t_target).at(0)/target_factor,
     	     (*t_target).at(1)/target_factor);

            fprintf(results,"%-4s %-4s %8f\n",
	      (*t_name).c_str(),
	      "-1",
	      prediction/target_factor);
      */

      results[*t_name] = prediction;

      t_input++;
      t_target++;
      t_name++;
      no++;
    }
  // printf("-------------------------------------------------\n");




  rank_results(results);



  return aue/fabs(((float) no)*target_factor);
}



void NN_DISCRIMINATION::rank_results(map<vector<string>, vector<float>,DecoyDecoyCmp> results)
{

  map<vector<string>, vector<float>, DecoyDecoyCmp>::iterator iter = results.begin();


  map<string, float, PDBDecoyCmp> simple_score;


  while(iter != results.end())
    {

      //  string temp = "less";
      string name = (*iter).first[0]+" "+(*iter).first[1];
      /* cout<<"Scoring "<<(*iter).first[0]
	  <<" "<<(*iter).first[1]
	  <<" "<<(*iter).first[2]
	  <<" -> "<<(*iter).second[0]
	  <<" and "<<(*iter).second[1]
	  <<endl;
      */
      if(simple_score.count(name)==0)
	simple_score[name]=-1.0;


      if((*iter).second[0]<(*iter).second[1])
	{
	  //temp = "more";
	  
	  //  cout<<"in loop with "<<(*iter).second[0]<<" and "<<(*iter).second[1] <<", now "<<simple_score.size()<<" keys"<<endl;

	  simple_score[name]-=1.0;
	  
	} 
      /*  cout<<"In "
	  <<(*iter).first[0]
	  <<" "
	  <<(*iter).first[1]
	  <<" is "
	  <<temp
	  <<" than "
	  <<(*iter).first[2]
	  <<endl;
      */
      iter++;
    }


  map<string, float, PDBDecoyCmp>::iterator score = simple_score.begin();

  while(score != simple_score.end())
    {

      //      cout<<"Now writting "<<(*score).first<<" of "<<simple_score.size()<<endl;
      fprintf(resultsF,"%-9s %f\n",(*score).first.c_str(), (*score).second);	
      score++;
    }


}




void NN_DISCRIMINATION::write_network(string name)
{
  /*** save network to file ***/
  network.write_file(name);
}


void NN_DISCRIMINATION::read_network(string name)
{
  /*** read network from file ***/
  network.read_file(name);

  network.set_learning_rate(learning_rate);
  network.set_transfer_functions(tf);
  //network.set_linear_coef(linear_coef);
}


void NN_DISCRIMINATION::set_parameters()
{


  ifstream parameters;
  parameters.open("nn_discrimination.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'nn_discrimination.cfg' could not be found\n";
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

  printf("Configurations for fitter NN:\n");
  printf("\n");
  printf("\n");

  printf("DATASET:              %s\n",dataset.c_str());
  printf("NO_EPOCHS:            %d\n",no_epochs);
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

  printf("TRANSER_FUNCTIONS:    ");
  for(unsigned int i=0; i<tf.size(); i++)
    printf("%d ",tf[i]);
  printf("end\n");
  printf("\n");

  printf("\n");
  printf("\n");
  printf("Version: $Id: nn_discrimination.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");

  //correcting topology 0 -> no_inputs
  for(unsigned int i=0; i<topology.size(); i++)
    if(topology[i] == 0)
      topology[i] = no_inputs;

}


void NN_DISCRIMINATION::writeDescription(FILE * reporter, FILE * resultsFile)
{

  rep=reporter;
  resultsF=resultsFile;
  /*** write out to out file **/
  version =string("$Id: nn_discrimination.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
