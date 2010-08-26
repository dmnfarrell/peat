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
$Id: newnn.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "elementsnn.h"


ElementsNN::ElementsNN(FILE * reporter):Nnmodel()
{

  /*** write out to out file **/
  writeDescription(reporter);
 
  printf("Using the elements network!\n");
 
  threshold = 10;

}

ElementsNN::~ElementsNN(){}


void ElementsNN::create_network()
{
  vector<int> topology;
  topology.push_back(no_inputs);
  topology.push_back(1);

  float learning_rate = 0.000000005; //0


  network.create_network(topology);
  network.set_learning_rate(learning_rate);
  printf("%s\n",network.print_network().c_str());
}


void ElementsNN::build_inputs(pose * A)
{

  vector<Atom *> protein_atoms = convert_objects_to_atoms(A->protein);
  vector<Atom *> ligand_atoms = convert_objects_to_atoms(A->ligand);
  vector<float> res;
  for(int p=0; p<protein_atoms.size(); p++)
    for(int l=0; l<ligand_atoms.size(); l++)
      {

	//input for distance
	float dist = sqrt((protein_atoms[p]->x - ligand_atoms[l]->x)*(protein_atoms[p]->x - ligand_atoms[l]->x)+
			  (protein_atoms[p]->y - ligand_atoms[l]->y)*(protein_atoms[p]->y - ligand_atoms[l]->y)+
			  (protein_atoms[p]->z - ligand_atoms[l]->z)*(protein_atoms[p]->z - ligand_atoms[l]->z));
   
	if(dist < threshold)
	  {
	    res.clear();
	    
	    //input for protein atom
	    res = types[protein_atoms[p]->element.c_str()];
	    
	    //input for ligand atom
	    res.insert(res.end(), types[ligand_atoms[l]->element.c_str()].begin(), types[ligand_atoms[l]->element.c_str()].end()); 
	    
	    res.insert(res.end(),dist);
	    
	    //input for context
	    //res.insert(res.end(),); 
	    
	    if(res.size() == no_inputs)
	      A->inputs.push_back(res);
	    else
	      printf("Warning: could not create input for atom pair '%s'-'%s'\n",
		     protein_atoms[p]->element.c_str(),ligand_atoms[l]->element.c_str());

	  }
      } 


}

void ElementsNN::define_input()
{

  no_types = 12;
  no_inputs = 2*no_types + 1 + 0;

  //use elements
  float  H [] = {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  float  C [] = {0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  float  N [] = {0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  float  O [] = {0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  float  F [] = {0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  float  NA[] = {0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0};
  float  P [] = {0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0};
  float  S [] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0};
  float  CL[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0};
  float  CA[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0};
  float  BR[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0};
  float  I [] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0};
  
  //convert to vectors  
  vector<float>  Hvector ( H, H + 12 );
  vector<float>  Cvector ( C, C + 12 );
  vector<float>  Nvector ( N, N + 12 );
  vector<float>  Ovector ( O, O + 12 );
  vector<float>  Fvector ( F, F + 12 );
  vector<float> NAvector (NA,NA + 12 );
  vector<float>  Pvector ( P, P + 12 );
  vector<float>  Svector ( S, S + 12 );
  vector<float> CLvector (CL,CL + 12 );
  vector<float> CAvector (CA,CA + 12 );
  vector<float> BRvector (BR,BR + 12 );
  vector<float>  Ivector ( I, I + 12 );
  
  //map vectors
  types["H" ] =  Hvector;
  types["C" ] =  Cvector;
  types["N" ] =  Nvector;
  types["O" ] =  Ovector;
  types["F" ] =  Fvector;
  types["NA"] = NAvector;
  types["P" ] =  Pvector;
  types["S" ] =  Svector;
  types["CL"] = CLvector;
  types["CA"] = CAvector;
  types["BR"] = BRvector;
  types["I" ] =  Ivector;
 
  



}




void ElementsNN::test_pose(pose * A)
{
  float energy = 0;
  /*** sum up energy contributions from all atom pairs ***/
  for(int i=0; i < A->inputs.size(); i++)
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
  for(int i=0; i<A->inputs.size();i++)
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
  network.write_file(name);
}










/*
  for(int p=0; p<protein_atoms.size(); p++)
    for(int l=0; l<ligand_atoms.size(); l++)
      {
	res.clear();
	
	//input for protein atom
	res = types[protein_atoms[p]->element.c_str()];

	//input for ligand atom
	res.insert(res.end(), types[ligand_atoms[l]->element.c_str()].begin(), types[ligand_atoms[l]->element.c_str()].end()); 

	//input for distance
	float dist = sqrt((protein_atoms[p]->x - ligand_atoms[l]->x)*(protein_atoms[p]->x - ligand_atoms[l]->x)+
			  (protein_atoms[p]->y - ligand_atoms[l]->y)*(protein_atoms[p]->y - ligand_atoms[l]->y)+
			  (protein_atoms[p]->z - ligand_atoms[l]->z)*(protein_atoms[p]->z - ligand_atoms[l]->z));

	res.insert(res.end(),dist);

	//input for context
	//res.insert(res.end(),); 

	/*	if(l==0)
	  {	
	    printf("--------------------------------------------------------------------------------------------------------------\n");
	    printf(" E- E:   H   C   N   O   F  NA   P   S  CL  CA  BR   I |  H   C   N   O   F  NA   P   S  CL  CA  BR   I | dist\n");
	    printf("--------------------------------------------------------------------------------------------------------------\n");
	       
	  }
	
	if(res.size() == no_inputs)
	  {
	    /*  printf("%2s-%2s: ",protein_atoms[p]->element.c_str(),ligand_atoms[l]->element.c_str());
	    for(int i=0; i<res.size(); i++ )
	      {
		if(i == 24)
		  printf("%3.9f ",res[i]);
		else
		  printf("%3.0f ",res[i]);
		if(i==11 || i == 23)
		  printf("|");
	      }
	    
	    float dist_check = sqrt(
				    (protein_atoms[p]->x - ligand_atoms[l]->x)*(protein_atoms[p]->x - ligand_atoms[l]->x)+
				    (protein_atoms[p]->y - ligand_atoms[l]->y)*(protein_atoms[p]->y - ligand_atoms[l]->y)+
				    (protein_atoms[p]->z - ligand_atoms[l]->z)*(protein_atoms[p]->z - ligand_atoms[l]->z));

	    printf(" %3.9f",dist_check);
	    if(res[24]-dist_check > 0.000001)
	      printf("DISTERROR!");
	    printf("\n");
	    
	    if(res[24] < 10)
	      A->inputs.push_back(res);
	  }
	else
	  {
	    printf("Warning: could not create input for atom pair '%s'-'%s'\n",
		   protein_atoms[p]->element.c_str(),ligand_atoms[l]->element.c_str());
	  }


*/






void ElementsNN::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: newnn.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
