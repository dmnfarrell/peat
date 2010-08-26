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
$Id: nnpot.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "nnpot.h"


Nnpot::Nnpot(FILE * resultsFile, FILE * reporter, 
	     vector<Soup *>  protein_set, vector<Soup *>  ligand_set,  int model_id, 
	     bool do_training):Method(reporter)

{
  /*************************************** 
gets a vector of proteins and a vector of lignads

if train
   do n-fold cross-validation, where n = round_up(100/test_pct)
      divide into train, validation, and test set
      train on train set until error on validation set is minimised
      test on test set and write to results.txt
   repeat until all poses have been tested
      



if !train
   put everything in test set
   load network from file
   test


  ***************************************/


 /*** write out to out file **/
  writeDescription(reporter);
  results_f = resultsFile;

  /*** check that equal amount of proteins and ligands are given ***/
  if(protein_set.size() != ligand_set.size())
    {

	if(ligand_set.size() % protein_set.size() == 0)
	{
	  printf("Assuming that decoys are included - reorganising lists\n");
	  reorganise_lists(protein_set, ligand_set);
	  ligand_set = new_ligands;
	  
	}
	else
	  {
	    printf("Warning: NNpot needs equal amounts of protein and ligand structures in order to work!\n");
	    exit(0);
	  }
    }



  /*** select model to use ***/

  if(model_id == 1)
    {
      printf("Elements chosen!\n");
      ElementsNN temp_model(reporter);

      model = &temp_model;
      organise(do_training, protein_set, ligand_set);
    }
  else if(model_id == 2)
    {
      printf("Distances chosen!\n");
      DistanceNN temp_model(reporter);

      model = &temp_model;
      organise(do_training, protein_set, ligand_set);
    }

  else if(model_id == 3)
    {
      printf("NN ES chosen!\n");

      Generalparameters prepA(reporter, protein_set);
      Generalparameters prepB(reporter, ligand_set);

      NN_ES temp_model(reporter);

      model = &temp_model;
      organise(do_training, protein_set, ligand_set);
    }
  else if(model_id == 4)
    {
      printf("NN vdW chosen!\n");

      //      Generalparameters prepA(reporter, protein_set);
      //      Generalparameters prepB(reporter, ligand_set);

      NN_VDW temp_model(reporter);

      model = &temp_model;
      organise(do_training, protein_set, ligand_set);
    }

  else if(model_id == 5)
    {
      printf("NN vdW chosen!\n");
      
      
      NN_HBOND temp_model(reporter);

      model = &temp_model;
      organise(do_training, protein_set, ligand_set);
    }
  else if(model_id == 6)
    {
      printf("Elements 2 chosen!\n");
      Elementsnn2 temp_model(reporter);

      model = &temp_model;
      organise(do_training, protein_set, ligand_set);
    }
  else if(model_id == 7)
    {
      printf("NN ligand entropy chosen!\n");
      NN_LS temp_model(reporter);

      model = &temp_model;
      organise(do_training, protein_set, ligand_set);
    }
  else if(model_id == 8)
    {
      printf("Generic nn!\n");
      GenericNN temp_model(reporter);

      model = &temp_model;
      organise(do_training, protein_set, ligand_set);
    }
  
  
}

Nnpot::~Nnpot(){}


void Nnpot::organise(bool do_training, vector<Soup *>  protein_set, vector<Soup *>  ligand_set)
{

  //get parameters from selected model
  no_epochs               = model->get_no_epochs();
  target_factor           = model->get_target_factor();
  max_epochs_since_record = model->get_max_epochs_since_record();
  dataset                 = model->get_dataset();
  train_pct               = model->get_train_pct();
  validate_pct            = model->get_validate_pct();
  test_pct                = model->get_test_pct();

  model->define_input();

  incl_co_factors         = model->get_include_co_factors();
  incl_water              = model->get_include_water();

  inputs_on_the_fly = true;

  //generating poses
  for(unsigned int i=0; i<protein_set.size(); i++)
    prepare_pose(protein_set[i], ligand_set[i]);
  
  //generating decoy poses
  for(unsigned int i=0; i<decoys.size(); i++)
    prepare_decoy_pose(protein_set, decoys[i]);


  //set parameters

  stats_file.open("stats.gnu");
  error_file.open("errors.gnu");


  /*** find experimental binding affinities ***/
  set_targets(dataset);

  /*** prepare for bookkeeping ***/
  vector<string> pdb_names;
  for(unsigned int i=0; i<all_poses.size(); i++)
    if(count(pdb_names.begin(), pdb_names.end(),all_poses[i].name.substr(0,4)) == 0)
      pdb_names.push_back(all_poses[i].name.substr(0,4));
  
  vector<int> pdbs_tested(pdb_names.size(),0);

  vector<int> tested(all_poses.size(),0);
  int no_crossfold_validation;
  if(test_pct != 0)
    no_crossfold_validation = (int) ceil(100.0/test_pct);
  else
    no_crossfold_validation = 1;
  
  int no_test_pdbs = (int) ceil(pdb_names.size()*(test_pct)/100.0);
  int no_train_pdbs = (int) (pdb_names.size()*(train_pct)/100.0);
      
  /*** do n-fold cross-validation ***/
  for(int run=0; run<no_crossfold_validation; run++)
    {
      vector<int> testing(all_poses.size(),0);
      int no_testing = 0;
      int no_training = 0;

      train_poses.clear();
      validate_poses.clear();
      test_poses.clear();

      /*** divide poses into training, validation, and test poses ***/
      for(unsigned int i=0; i<pdb_names.size(); i++)
	{
	  if(pdbs_tested[i] == 0 && no_testing < no_test_pdbs) 
	    {
	      for(unsigned int j=0; j<all_poses.size(); j++)
		if(all_poses[j].name.substr(0,4) == pdb_names[i])
		  {
		    test_poses.push_back(&all_poses[j]);//use this structure for testing
		    tested[j] = 1;
		    testing[j] = 1;
		  }
	      pdbs_tested[i] = 1;
	      no_testing++;
	    }
	  else
	    {
	      if(no_training < no_train_pdbs)
		{
		  for(unsigned int j=0; j<all_poses.size(); j++)
		    if(all_poses[j].name.substr(0,4) == pdb_names[i])
		      train_poses.push_back(&all_poses[j]); //use this structure for training
		  no_training++;
		}
	      else
		{
		  for(unsigned int j=0; j<all_poses.size(); j++)
		    if(all_poses[j].name.substr(0,4) == pdb_names[i])
		      validate_poses.push_back(&all_poses[j]); //use this structure for validation
		

		}
	
	    }
	}

      /*** print out poses ***/
      printf("Cross-validation no. %d has been prepared\n",run+1);
      for(unsigned int i=0; i<all_poses.size(); i++)
	printf("Test status of %s is %d\n", all_poses[i].name.c_str(),tested[i]);
      write_poses();
      
      /*** do training ***/
      if(train_poses.size() > 0)
	train();
      
      /*** test ***/
      if(train_poses.size() == 0)
	model->read_network("best_net");
      test();

      /*** write status ***/
      write_status();
      
    }
  
}



void Nnpot::train()
{

  model->create_network();



  int epoch = 1;
  float min_validate_error = 100000000;
  int epochs_since_record = 0;
  int train_number = 0;

  bool go = true;

  /*** for epochs ***/
  while(go)
    {
      /*** train on all train poses ***/
      for(unsigned int p=0; p < train_poses.size(); p++)
	{
	  if(inputs_on_the_fly)
	    model->build_inputs(train_poses[p]);

	  model->test_pose(train_poses[p]);
	  
	  /*** update weights in nn ? ***/
	  train_number++;
	  if (train_number % model->get_smoothing()==0)
	    model->update_weights();
	    
	  /*** find error and write to stats file and stdout ***/
	  stats_file << train_poses[p]->rel_error*100<<"  ";

	  /*** do training ***/
	  model->train_pose(train_poses[p]);
	  
	  if(inputs_on_the_fly)
	    train_poses[p]->inputs.clear();
	}
      /*** test on all validation poses ***/
      for(unsigned int p=0; p < validate_poses.size(); p++)
	{
	  if(inputs_on_the_fly)
	    model->build_inputs(validate_poses[p]);

	  model->test_pose(validate_poses[p]);

	  /*** find error and write to stats file and stdout ***/
	  stats_file << validate_poses[p]->rel_error*100<<"  ";

	  if(inputs_on_the_fly)
	    validate_poses[p]->inputs.clear();

	}

      if(epoch % 10 == 0 || epoch ==1)
	write_train_status(epoch);

      /*** calc mean unsigned error and print to error file ***/
      float train_error = 0;
      float validate_error = 0;
      
      error_file << epoch<<"    ";

      for(unsigned int p=0; p < train_poses.size(); p++)
	train_error += fabs(train_poses[p]->error);
      error_file << train_error/(fabs(target_factor) * train_poses.size())<<"    ";

      for(unsigned int p=0; p < validate_poses.size(); p++)
	validate_error += fabs(validate_poses[p]->error);
      error_file << validate_error/(fabs(target_factor) * validate_poses.size())<<"\n";
      
      
      /*** check if validation error record has been broken ***/


      if(validate_error < min_validate_error)
	{
	  printf("Record broken - new record is %f\n", validate_error/validate_poses.size());
	  min_validate_error = validate_error;
	  epochs_since_record = 0;
	  model->write_network("best_net");
	}
      else
	{
	  epochs_since_record++;
	  if(epochs_since_record > max_epochs_since_record)
	    {
	      model->read_network("best_net");
	      go = false;
	    }
	}

      error_file.flush();


      /*** flush stats file ***/
      stats_file <<"\n";
      stats_file.flush();

      epoch++;

      if(epoch >= no_epochs)
	go = false;
    }
  stats_file <<"#cross-validation done\n\n";
  

}


void Nnpot::test()
{
  printf("Testing started!\n");
  vector<int> erasable;

  for(unsigned int p=0; p<test_poses.size(); p++)
    {
      printf("Now looking for %s in %d decoys\n",test_poses[p]->name.substr(0,4).c_str(),decoy_poses.size());
      /*** do decoys ***/
      
      for(unsigned int d=0; d<decoy_poses.size(); d++)
	if(decoy_poses[d].name.substr(0,4) == test_poses[p]->name.substr(0,4))
	  {
	    /*** test ***/
	    model->build_inputs(&decoy_poses[d]);
	    model->test_pose(&decoy_poses[d]);

	    string decoy_number;
	    if(decoy_poses[d].name.find_last_of("_")== string::npos)
	      decoy_number = "-1";
	    else
	      decoy_number = decoy_poses[d].name.substr(decoy_poses[d].name.find_last_of("_")+1);
	    
	    //printf("Found decoy %s\n", decoy_poses[d].name.c_str());
	    
	    fprintf(results_f,"%-4s %-4s %8f\n",
		    decoy_poses[d].protein[0]->identify().substr(0,4).c_str(),
		    decoy_number.c_str(),
		    decoy_poses[d].prediction/target_factor);

	    /*** clean up to save memory ***/
	    //	    decoy_poses.erase(d);
	    erasable.push_back(d);

   	  }
      int move_back =0;
      for(unsigned int i=0; i<erasable.size(); i++)
	{
	  //	  printf("Deleting %s\n",decoy_poses[erasable[i]- move_back].name.c_str());	
	  decoy_poses.erase(decoy_poses.begin() + erasable[i]-move_back);
	  move_back++;
	}
      erasable.clear();

      /*** test experimental poses ***/
      if(inputs_on_the_fly)
	model->build_inputs(test_poses[p]);

      model->test_pose(test_poses[p]);

      if(inputs_on_the_fly)
	test_poses[p]->inputs.clear();


      /*** write to file and stdout ***/
      test_poses[p]->tested = true;
      test_poses[p]->test_prediction = test_poses[p]->prediction/target_factor;
      test_poses[p]->test_error = test_poses[p]->target/target_factor - test_poses[p]->test_prediction;

      string decoy_number;
      if(test_poses[p]->name.find_last_of("_")== string::npos)
	decoy_number = "-1";
      else
	decoy_number = test_poses[p]->name.substr(test_poses[p]->name.find_last_of("_")+1);
      
      fprintf(results_f,"%-4s %-4s %8f\n",
	      test_poses[p]->protein[0]->identify().substr(0,4).c_str(),
	      decoy_number.c_str(),
	      test_poses[p]->prediction/target_factor);
     

      printf("Result for %s is %5.3f and target is %5.3f\n",test_poses[p]->protein[0]->identify().substr(0,4).c_str(),
	     test_poses[p]->prediction/target_factor,test_poses[p]->target/target_factor);

    }
}




void Nnpot:: reorganise_lists(vector<Soup *>  protein_set, vector<Soup *>  ligand_set)
{
  /*** reorganise soup list so that an equal amount of soups excists in protein and ligand list
       and that all decoys are moved to a seperate list ***/


  string lig_name;

  for(unsigned int p=0; p<protein_set.size(); p++)
    {
      lig_name = protein_set[p]->name.substr(0,4)+"-lig_101";

      for(unsigned int l=0; l<ligand_set.size(); l++)
	if(ligand_set[l]->name == lig_name)
	    new_ligands.push_back(ligand_set[l]);
    }

  for(unsigned int l=0; l<ligand_set.size(); l++)
    if(ligand_set[l]->name.substr(9) != "101")
	decoys.push_back(ligand_set[l]);

  printf("Reorganisation complete - got %d new ligands and %d decoys\n",new_ligands.size(),decoys.size());

}


void Nnpot::write_train_status(int epoch)
{
  
  printf("***Now starting epoch %d ***\n",epoch);
  printf("Pose                                 Prediction       Target        Error   Rel. error\n");
	
  /*** print all train poses ***/
  for(unsigned int p=0; p < train_poses.size(); p++)
    printf("TRAIN:    %s  %12.2f %12.2f %12.2f %12.2f\n",
	   train_poses[p]->name.c_str(), train_poses[p]->prediction, train_poses[p]->target, 
	   train_poses[p]->error,train_poses[p]->rel_error*100);

  /*** print all validation poses ***/
  for(unsigned int p=0; p < validate_poses.size(); p++)
    printf("VALIDATE: %s  %12.2f %12.2f %12.2f %12.2f\n",
	   validate_poses[p]->name.c_str(), validate_poses[p]->prediction, validate_poses[p]->target, 
	   validate_poses[p]->error,validate_poses[p]->rel_error*100);

}



void Nnpot::write_status()
{ 
  /*** write networks to files ***/
  //  model->write_network("networks");

  printf("--------------------------------------------------------------------------------------\n");
  printf(" Pose    Prediction [kJ/mol]      Target [kJ/mol]       Error [kJ/mol]      Rel.error%%\n");
  printf("--------------------------------------------------------------------------------------\n");

  for(unsigned int p=0; p<all_poses.size(); p++)
    {
      if(all_poses[p].tested == true)
	  printf(" %s         %16.4f     %16.4f     %16.4f %12.2f\n",
		 all_poses[p].protein[0]->identify().substr(0,4).c_str(),
		 all_poses[p].test_prediction, 
		 all_poses[p].target/target_factor, 
		 all_poses[p].test_error,
		 all_poses[p].test_error/(all_poses[p].target/target_factor)*100);
      else
	  printf(" %s                      N/A         %12.4f                  N/A          N/A\n",
		 all_poses[p].protein[0]->identify().substr(0,4).c_str(),
		 all_poses[p].target/target_factor);


    }

}

void Nnpot::write_poses()
{
  /*** write train poses ***/
  for(unsigned int i=0; i<train_poses.size(); i++)
    printf("Train pose nr %d is %s and has %d input vectors\n",i,train_poses[i]->name.c_str(),train_poses[i]->inputs.size());
  /*** write validate poses ***/
  for(unsigned int i=0; i<validate_poses.size(); i++)
    printf("Validate pose nr %d is %s and has %d input vectors\n",i,validate_poses[i]->name.c_str(),validate_poses[i]->inputs.size());
  /*** write test poses ***/
  for(unsigned int i=0; i<test_poses.size(); i++)
    printf("Test pose nr %d is %s and has %d input vectors\n",i,test_poses[i]->name.c_str(),test_poses[i]->inputs.size());
  /*** sum up number of poses ***/
  printf("Total number of poses is: %d (train) + %d (validate) + %d (test) = %d\n",
	 train_poses.size(),validate_poses.size(),test_poses.size(),
	 train_poses.size()+validate_poses.size()+test_poses.size());
}


void Nnpot::prepare_pose(Soup * protein, Soup * ligand) //set: 1->train 2->validate 3->test
{
  pose res;
  //set pose name, protein, and ligand
  res.name = protein->name + " vs " + ligand->name;
  res.protein = convert_soup_to_objects(protein, true, incl_co_factors, incl_water);
  res.ligand =  convert_soup_to_objects(ligand);
  res.tested = false;

  if(!inputs_on_the_fly)
    model->build_inputs(&res);

  printf("Pose %s has been prepared with %d inputs\n",res.name.c_str(),res.inputs.size());

  all_poses.push_back(res);
}


void Nnpot::prepare_decoy_pose(vector<Soup *> proteins, Soup * ligand) 
{
 
  string pdb = ligand->name.substr(0,4);
  Soup * protein;

  /*** find corresponding protein ***/
  for(unsigned int p=0; p<proteins.size(); p++)
    if(proteins[p]->name.substr(0,4) == pdb)
      {
	protein = proteins[p];
	pose res;
	//set pose name, protein, and ligand
	res.name = protein->name + " vs " + ligand->name;
	res.protein = convert_soup_to_objects(protein, true, incl_co_factors, incl_water);
	res.ligand =  convert_soup_to_objects(ligand);
	res.tested = false;
	
	//  model->build_inputs(&res);
	
	//  printf("Decoy pose %s has been prepared with %d inputs\n",res.name.c_str(),res.inputs.size());
	
	decoy_poses.push_back(res);

      }

}



void Nnpot::set_targets(string dataset)
{
  Targets target(report, dataset);

  string decoy_number;
      
  
  for(unsigned int i=0; i<all_poses.size(); i++)
    {
      if(all_poses[i].name.find_last_of("_")== string::npos)
	all_poses[i].target = target_factor * target.get_G(all_poses[i].name.substr(0,4));
      else
	{
	  decoy_number = all_poses[i].name.substr(all_poses[i].name.find_last_of("_")+1);
	  float rmsd = target.get_rmsd(all_poses[i].name.substr(0,4), atoi(decoy_number.c_str()));
	  all_poses[i].target = target_factor * target.get_decoy_G(all_poses[i].name.substr(0,4),rmsd);
	}

      printf("Target of '%s' is %f\n",all_poses[i].name.c_str(),all_poses[i].target);
    }

  

}




void Nnpot::writeDescription(FILE * reporter)
{
  //Saves the reporter file internally so that it can be passed on to sub methods
  report = reporter;
  /*** write out to out file **/
  version =string("$Id: nnpot.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
