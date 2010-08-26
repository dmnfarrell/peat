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


#include "NN.h"

NN::NN(){}


NN::NN(vector<int> layers) 
{
  create_network(layers);
}

NN::NN(string filename) //read a network from file
{
  read_file(filename);
}

NN::~NN(){}


void NN::create_network(vector<int> layers) //create a new network
{
  network.clear();
  //  srand(time(NULL));
  for (unsigned int i=1; i<layers.size(); i++)
    {
      if(layers[i] == 0)
	{
	  printf("You can't have a network layer with zero neurons, aborting!\n");
	  exit(0);
	}
      network.push_back(generate_layer(layers[i],layers[i-1])); //(out, in)
    }
  
  learning_rate = 0.001; //standard value - can be changed using the set_learning_rate method
  linear_coef = 5;

}


float NN::train(vector<vector<float> > input, vector<vector<float> > expectation) //train on a set of inputs and targets
{

  float aue=0;
  float instant_error=0;

  //checking that input has the same number of input and expectation vectors 
  if(input.size() != expectation.size())
    {
      printf("Error: Input and expectation vectors not of same lenght!\n");
      return 0;
    }

  for(unsigned int m=0; m<input.size(); m++)
    {

      if(expectation[m].size() != network[network.size()-1].no_neuron_out)
	{
	  printf("Error: The network has %d output neurons and must be given %d expectation values!\n"
		 ,network[network.size()-1].no_neuron_out, network[network.size()-1].no_neuron_out);
	  return 0;
	}

      
      test(input[m]);
 
      //calculates the error on the output layer
      network[network.size()-1].error.clear();
      for(unsigned int i=0; i<expectation[m].size(); i++)
	network[network.size()-1].error.push_back(expectation[m][i] - network[network.size()-1].output[i]);

      instant_error = 0;
      for(unsigned int i=0; i<network[network.size()-1].error.size(); i++)
	instant_error += fabs(network[network.size()-1].error[i]);

      instant_error = instant_error/network[network.size()-1].error.size();
      aue += fabs(instant_error);

      //updating weights
      backpropagate_error();
      calculate_delta_weights();
      //      update_weights();        
    }   


  aue = aue/input.size();
  return aue;

}



void NN::train_ext_error(vector<float> input, vector<float> errors) //train on inputs using externally found errors
{

  
  if(errors.size() != network[network.size()-1].no_neuron_out)
    {
      printf("Error: The network has %d output neurons and must be given %d expectation values!\n"
	     ,network[network.size()-1].no_neuron_out, network[network.size()-1].no_neuron_out);
      return;
    }

  test(input);
  

  //sets the error on the output layer
  network[network.size()-1].error = errors;


  //updating weights
  backpropagate_error();
  calculate_delta_weights();
  //  update_weights();        

}



vector<float> NN::test(vector<float> input)
{

  if(input.size() != network[0].no_neuron_in)
    {
      printf("Error: The network has %d input neurons and must be given %d input values!\n"
	     ,network[0].no_neuron_in, network[0].no_neuron_in);

      vector<float> zero(1,0);
      return zero;
    }

  network[0].output = calc_layer(input, &network[0]);

  for (unsigned int i=1; i<network.size(); i++)     
    network[i].output = calc_layer(network[i-1].output, &network[i]);





  /* *****************************

  cout<<"NN got: ";
  for(unsigned int i=0; i< input.size();i++)
    cout<<input[i]<<" ";

  cout<<endl;

  cout<<"Output: ";
  for(unsigned int i=0; i< network[network.size()-1].output.size();i++)
    cout<<network[network.size()-1].output[i]<<" ";

  cout<<endl;


  ******************************/

  return network[network.size()-1].output;
}

void NN::set_learning_rate(float learn)
{
  learning_rate = learn;
}
         
void NN::set_linear_coef(float coef)
{
  linear_coef = coef;
}


void NN::set_transfer_functions(vector<int> tfs)
{

  /*
    0: Sigmoid
    1: linear
  */

  if(tfs.size() != network.size())
    {
      printf("ERROR: Set transfer function: size of input vector must be equal to number network layers!\n");
      exit(0);
    }
  
  for(unsigned int i=0; i<network.size(); i++)
      network[i].transfer_function = tfs[i];
}


void NN::calculate_delta_weights() //calculates deltaWeights - but does not update the weights 
{
  for(unsigned int n=network.size()-1; n>=0; n--)
    {
      
      for(unsigned int row=0; row<network.at(n).W.size(); row++)
	{
	  for(unsigned int col=0; col<network.at(n).W.at(0).size(); col++)
	    { 
	      float deltaW = 
		learning_rate   //learning rate
		*network.at(n).error.at(row)  //error on this layer
		//	*network.at(n).output.at(row)*(1-network.at(n).output.at(row)) //output on this layer
		*do_diff_transfer(network.at(n).output.at(row), network.at(n).transfer_function)
		*network.at(n).input.at(col);  //input from previus layer
	      
	      
	      network.at(n).deltaW.at(row).at(col) += deltaW;
	      
	    } 
	}
    }


}

void NN::update_weights()
{
  for(unsigned int n=network.size()-1; n>=0; n--)
    for(unsigned int row=0; row<network.at(n).W.size(); row++)
      for(unsigned int col=0; col<network.at(n).W.at(0).size(); col++)
	{
	  network.at(n).W.at(row).at(col) += network.at(n).deltaW.at(row).at(col);
	  network.at(n).deltaW.at(row).at(col) = 0.0;
	}
}


void NN::backpropagate_error()
{
  //W^T*error*est*(1-est)
  vector<float> new_error;

  for(unsigned int n=network.size()-1; n>0; n--)
    {
      new_error.clear();
      for(unsigned int col=0; col<network.at(n).W.at(0).size(); col++)
	{
	  float ne = 0; 
	  for(unsigned int row=0; row<network.at(n).W.size(); row++)
	    {
	      ne += 
		network.at(n).W.at(row).at(col)  //weight
		*network.at(n).error.at(row)   //error
		*network.at(n).input.at(col)*(1-network.at(n).input.at(col)); //input
	    }
	  new_error.push_back(ne);
	}     
      network[n-1].error = new_error;
    }
}



string NN::print_network()
{
  stringstream res;
  res << "N E T W O R K   L A Y O U T \n\n";

  for(unsigned int n=0; n<network[0].W[0].size(); n++)
    res << " (I)";
  res << "\n\n";
 
 
 for(unsigned int l=0; l<network.size(); l++)
    {
      res << "Weights of layer " << setw(3)<< l+1 << " are [ " 
	  << setw(3)<< network[l].W.size() 
	  << " X " << setw(3)<< network[l].W[0].size() <<  " ]:\n";

      for(unsigned int i=0; i<network[l].W.size(); i++)
	{
	  for(unsigned int j=0; j<network[l].W[0].size(); j++)
	    res << " " << setw(12) << network[l].W[i][j];
	  res << "\n";
	}
      res << "\n";
      for(unsigned int n=0; n<network[l].W.size(); n++)
	{
	  if(l == network.size()-1)
	    res << " (O)";
	  else
	    res << " (H)";
	}
      res << "\n\n";
 
    }
 return res.str();
}

vector<float> NN::get_all_weights()
{

  vector<float> weights;

   for(unsigned int l=0; l<network.size(); l++)
     {
       for(unsigned int i=0; i<network[l].W.size(); i++)
	 {
	   for(unsigned int j=0; j<network[l].W[0].size(); j++)
	     weights.push_back(network[l].W[i][j]);
	}
     }
   return weights;

}


vector<float> NN::calc_layer(vector<float> input, nnlayer * layer)
{
  int n_out = layer->no_neuron_out; 
  int n_in = layer->no_neuron_in;

  layer->input = input;

  vector<float> result;
  result.assign(n_out, 0);


  for(int j=0; j<n_in; j++)
    if(input.at(j) != 0)
      for(int i=0; i<n_out; i++)
	result.at(i) += layer->W.at(i).at(j) * input.at(j);
  
  for(int i=0; i<n_out; i++)
    result[i] = do_transfer(result[i],layer->transfer_function);//sigmoid(result[i]);

  return result;
}


float NN::do_transfer(float in, int transfer_function)
{
  float res= 0;
  if(transfer_function == 0)
    res = sigmoid(in);
  else if(transfer_function == 1)
    res = linear_coef * in;
  else
    {
      printf("Error could not find transfer function with index '%d'\n",transfer_function);
      exit(0);
    }

  return res;
}


float NN::do_diff_transfer(float in, int transfer_function)
{
  /*
    Returns f'(x) as function of f(x) for each of the implemented transfer function
  */

  float res= 0;
  if(transfer_function == 0)
    res = in*(1-in); //sigmoid
  else if(transfer_function == 1)
    res = linear_coef;
  else
    {
      printf("Error could not find transfer function with index '%d'",transfer_function);
      exit(0);
    }

  return res;


}



float NN::sigmoid(float in)//the sigmoid funtion
{
  float res = 1/(1+exp(-in));
  return res;
}

float NN::sigmoid_diff(float in)//the differentiated sigmoid - NOT USED?
{
  float sig = sigmoid(in);
  float res = sig*(1-sig);
  return res;
}




nnlayer NN::generate_layer(int row, int col)//randomizes weights
{

  nnlayer res;
 
 
  vector<float> temp;
  vector<float> temp2;
  vector<vector<float> > weights;
  vector<vector<float> > deltaWeights;

  for(int i=0; i<row; i++)
    {
      for(int j=0; j<col; j++)
	{
	  temp.push_back( ( (float)(rand()%10000) ) /100000 - 0.05);
	  temp2.push_back(0.0);
	}
      weights.push_back(temp);
      deltaWeights.push_back(temp2);
      temp.clear();
      temp2.clear();
    }

  res.W = weights;
  res.deltaW = deltaWeights;

  res.no_neuron_out = row;
  res.no_neuron_in  = col;
  res.transfer_function = 0;

  return res;
}

void NN::write_file(string filename) //write network to file
{
  ofstream file (filename.c_str());

  if(!file.is_open())
    {
      printf("Error: Could not open file '%s'!",filename.c_str());
      return;
    }

  file << "*********************************************************\n";
  file << "This is an autogenerated file containing a neural network\n";
  file << "*********************************************************\n\n";

  file << print_network();

  file.close();
  
}

void NN::read_file(string filename) //clears current network and reads a new one from file
{

  network.clear();

  ifstream file (filename.c_str());

  if(!file.is_open())
    {
      printf("Error: Could not open file '%s'!",filename.c_str());
      return;
    }

  
  string line;
  while(file)
    {
      getline(file, line);
      
      if(line.substr(0,16) == "Weights of layer")
	{
	  //int layer_no = atoi(line.substr(17,3).c_str());
	  int no_rows  = atoi(line.substr(27,3).c_str());
	  int no_cols  = atoi(line.substr(33,3).c_str());
	  
	  nnlayer res;
	  vector<float> temp(no_cols,0);
	  vector<vector<float> > weights(no_rows,temp);

	  for(int i=0; i<no_rows; i++)
	    {
	      getline(file, line);
	      int pos = 0;
	      for(int j=0; j<no_cols; j++)
		{
		  weights.at(i).at(j) = atof(line.substr(pos,13).c_str());
		  pos += 13;
		}	  

	      res.W = weights;
	      res.no_neuron_out = no_rows;
	      res.no_neuron_in  = no_cols;

	    }
	  network.push_back(res);
	}
    }//end while

  file.close();
}


