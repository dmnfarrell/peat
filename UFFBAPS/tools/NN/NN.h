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

#ifndef NN_H
#define NN_H
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

using namespace std;


struct nnlayer{
  vector<vector<float> > W; 
  vector<vector<float> > deltaW; //for storing deltaWwights when smoothing 
  //[no_neuron_out][no_neuron_in]
  //[row]          [column]

  unsigned int no_neuron_in; 
  unsigned int no_neuron_out;

  vector<float> output; //has length equal to no_neuron_out
  vector<float> error; //has length equal to no_neuron_out
  vector<float> input; //has length equal to no_neuron_in

  int transfer_function;
};

class NN{

 public:
  NN();
  NN(vector<int> layers);
  NN(string filename);
  ~NN();

  float train(vector<vector<float> > input, vector<vector<float> > expectation);
  void train_ext_error(vector<float> input, vector<float> errors);
  vector<float> test(vector<float> input);

  void create_network(vector<int> layers); //create a new network
  void write_file(string filename);
  void read_file(string filename);
  string print_network();
  void set_learning_rate(float learn);
  void set_transfer_functions(vector<int> tfs);
  void set_linear_coef(float coef);
  vector<float> get_all_weights();
  void update_weights();
 private:
  vector<nnlayer> network; // hidden_layer,output_layer;

  float do_transfer(float in, int transfer_function);
  float do_diff_transfer(float in, int transfer_function);

  float learning_rate;
  float sigmoid(float in);
  float sigmoid_diff(float in);
  nnlayer generate_layer(int row, int col);
  vector<float> calc_layer(vector<float> input, nnlayer * layer);
  vector<vector<float> > train_test(vector<float> input);
  void calculate_delta_weights();

  void backpropagate_error();
  float linear_coef;


};

#endif
