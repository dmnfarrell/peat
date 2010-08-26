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
$Id: genericnn.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "genericnn.h"

/*
Generic NN

Number of contacts at a given threshold distance - e.g. 5 A




*/

GenericNN::GenericNN(FILE * reporter):Nnmodel()
{

  /*** write out to out file **/
  writeDescription(reporter);
 
  printf("Using the generic network!\n");
 
  threshold = 10;

}

GenericNN::~GenericNN(){}


void GenericNN::create_network()
{
  vector<int> topology;
  topology.push_back(no_inputs);
  topology.push_back(10);
  topology.push_back(1);

  vector<int> tf;
  tf.push_back(0);
  tf.push_back(0);
  //  topology.push_back(1);

  float linear_coef = 0.01;
  float learning_rate = 100; //0


  network.create_network(topology);
  network.set_learning_rate(learning_rate);
  network.set_transfer_functions(tf);
  network.set_linear_coef(linear_coef);

  printf("%s\n",network.print_network().c_str());
}


void GenericNN::build_inputs(pose * A)
{
  string lelements[] = {"C", "N", "O", "F", "NA", "P", "S", "CL", "CA", "BR", "I"};
  vector<string> elements(lelements, lelements+11);
  
  float lbins[] = {10,8,6,5,4};
  vector<float> bins(lbins, lbins+5);  


  vector<string> combos;
  vector<string>::iterator elements_iter1;
  vector<string>::iterator elements_iter2;
  for(elements_iter1=elements.begin(); elements_iter1!=elements.end(); elements_iter1++ )
    for(elements_iter2=elements.begin(); elements_iter2!=elements.end(); elements_iter2++ )
      {
	string s = *elements_iter1 + " " + *elements_iter2;
	combos.push_back(s);
      }
  no_inputs = combos.size()*bins.size();
  map<float,map<string, int> > contacts;
  vector<float>::iterator bin_iter = bins.begin();
  
  while(bin_iter != bins.end())
    {
      map<string, int> temp;
      
      for(vector<string>::iterator s=combos.begin(); s!=combos.end(); s++)
	temp[*s]=0;

      contacts[*bin_iter] = temp;
      bin_iter++;
    }

  vector<float> res;

  vector<Atom *> protein_atoms = convert_objects_to_atoms(A->protein);
  vector<Atom *> ligand_atoms = convert_objects_to_atoms(A->ligand);

  for(unsigned int p=0; p<protein_atoms.size(); p++)
    for(unsigned int l=0; l<ligand_atoms.size(); l++)
      {
	//input for distance
	float dist = sqrt((protein_atoms[p]->x - ligand_atoms[l]->x)*
			  (protein_atoms[p]->x - ligand_atoms[l]->x)+
			  (protein_atoms[p]->y - ligand_atoms[l]->y)*
			  (protein_atoms[p]->y - ligand_atoms[l]->y)+
			  (protein_atoms[p]->z - ligand_atoms[l]->z)*
			  (protein_atoms[p]->z - ligand_atoms[l]->z));
	
	for(unsigned int b=1; b<bins.size(); b++)
	  if(dist < bins[b-1] && dist > bins[b])
	    {
	      string s = protein_atoms[p]->element + " " + ligand_atoms[l]->element;
	      if (contacts[bins[b-1]].count(s)!=0)
		contacts[bins[b-1]][s]++;

	    }
   
      
      } 

  /*** print contacts ***/
  cout<< "   ";
  for(vector<string>::iterator s=combos.begin(); s!=combos.end(); s++)
    printf("%6s ",(*s).c_str());
  cout<<endl;
  
  for(map<float,map<string, int> >::iterator iter=contacts.begin(); 
      iter!=contacts.end(); iter++ )
    {
      printf("%3.0f",(*iter).first);
      for(vector<string>::iterator s=combos.begin(); s!=combos.end(); s++)
	{
	  printf("%6d ",(*iter).second[*s]);
	}
      cout <<  endl;
      
    }
  cout<< "--------------------------"<<endl;
  
  for(vector<string>::iterator s=combos.begin(); s!=combos.end(); s++)
    for(vector<float>::iterator b=bins.begin(); b!=bins.end(); b++)
      res.push_back((float) contacts[*b][*s]);
	

  /*** normalize inputs ***/
  float total = 0;
  for(unsigned int i=0; i<res.size();i++)
    total += res[i];
  for(unsigned int i=0; i<res.size();i++)
    res[i] =res[i]/total;
  
  if(res.size() == (unsigned int) no_inputs)
    A->inputs.push_back(res);
  else
    printf("Warning: could not create input - no_input:%d, res.size():%d\n",
	   no_inputs, res.size());


}

void GenericNN::define_input()
{

 

}




void GenericNN::test_pose(pose * A)
{
  float energy = 0;
  /*** sum up energy contributions from all atom pairs ***/
  for(unsigned int i=0; i < A->inputs.size(); i++)
    {
      /*      cout<<"inputs ";
      for(unsigned int a=0; a<A->inputs[i].size();a++) 
	cout<<A->inputs[i][a]<<" ";
      cout<<endl;
      */
      vector<float> out = network.test(A->inputs[i]);
      //cout<<"no outs:"<<out.size()<<"  out[0]:"<<out[0]<<endl;
      energy += out[0];
    }

  /*** set prediction and the error ***/
  //  cout<<"prediction: "<<energy<<endl;
  A->prediction = energy;
  A->error = A->target - A->prediction;
  A->rel_error = A->error / A->target;

}
void GenericNN::train_pose(pose * A)
{
  vector<float> error(1,A->error);
  for(unsigned int i=0; i<A->inputs.size();i++)
    network.train_ext_error(A->inputs[i],error);

  /*** print current network ***/
  //  printf("%s\n",network.print_network().c_str());

}



void GenericNN::write_network(string name)
{
  /*** save network to file ***/
  network.write_file(name);
}


void GenericNN::read_network(string name)
{
  /*** save network to file ***/
  network.write_file(name);
}


void GenericNN::update_weights()
{
  network.update_weights();
}



int GenericNN::get_no_epochs(){return 10000;}
int GenericNN::get_max_epochs_since_record(){return 10000;}
int GenericNN::get_smoothing(){return 1;}
float GenericNN::get_target_factor(){return -0.01;}
string GenericNN::get_dataset(){return "wang";}

float GenericNN::get_train_pct(){return 50;}
float GenericNN::get_validate_pct(){return 20;}
float GenericNN::get_test_pct(){return 30;}

bool GenericNN::get_include_co_factors(){return 0;}
bool GenericNN::get_include_water(){return 0;}
  





void GenericNN::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: genericnn.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
