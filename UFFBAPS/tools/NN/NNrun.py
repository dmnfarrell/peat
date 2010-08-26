

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

import sys
import math

class NNpython:
    "A test-only python implementation of the NN c++ class"

    def __init__(self):
        self.network = []

    def test(self,input):
        
        if(len(input) !=  self.network[0].cols):
            print "Error: The network has %d input neurons and must be given %d input values!\n" %(self.network[0].cols, self.network[0].cols)
            
            zero = []
            return zero
        

        self.network[0].output = self.network[0].calc_layer(input)
        
        for i in range(len(self.network)-1):
            self.network[i+1].output = self.network[i+1].calc_layer(self.network[i].output)
            
        return self.network[len(self.network)-1].output


    def print_network(self):
        print "P Y T H O N   N E T W O R K   L A Y O U T \n"

        print "Size of nn is %d."%len(self.network)
        for n in range(len(self.network[0].W[0])):
            sys.stdout.write(" (I)")
        print "\n\n"
  
        for l in range(len(self.network)):
            print "Weights of layer ", l+1, " are [ ", len(self.network[l].W)," X ",len(self.network[l].W[0])," ]:"

            for i in range(len(self.network[l].W)):
                for j in range(len(self.network[l].W[0])):
                    sys.stdout.write("  %f" %self.network[l].W[i][j])
                print " ";
	
            print " ";
            for n in range(len(self.network[l].W)):
                if(l == len(self.network)-1):
                    sys.stdout.write(" (O)")
                else:
                    sys.stdout.write(" (H)")
	
        print "\n"
 

    def read_file(self,name):

        file = open(name)

        #read file
        lines = file.readlines()
        for i in range(len(lines)):
            #find weights in files
            if(lines[i][0:16] == "Weights of layer"):
                layer_no = int(lines[i][17:20])
                no_rows  = int(lines[i][27:30])
                no_cols  = int(lines[i][33:36])
                print "Generating layer number %d with [%d x %d] matrix" %(layer_no,no_rows,no_cols)
                res = layer()
                
                #make a matrix
                weights = []
		for row in range(no_rows):
                    weights.insert(row, [])
                    for col in range(no_cols):
                        weights[row].insert(col, 0)

                #read weights
                for a in range(no_rows):
                    pos = 0
                    for j in range(no_cols):
                        weights[a][j] = float(lines[i+a+1][pos:pos+13])
                        pos += 13
 
                res.W = weights
                res.rows = no_rows
                res.cols  = no_cols
                self.network.append(res)
               # del res

        file.close()


class layer:

    def __init__(self):
        "A NN network layer"
        
        self.W = []
        
        self.rows = 0
        self.cols = 0
        self.output = []
    
    def set_weights(self,matrix, irows, icols):
        W = matrix
        rows = irows
        cols = icols


    def calc_layer(self,input):

        result = []
        for i in range(self.rows):
            result.append(0)
        for i in range(self.rows):
            for j in range(self.cols):
                result[i] += self.W[i][j] * input[j]

        for i in range(self.rows):
            result[i] = self.sigmoid(result[i])

        return result



    def sigmoid(self,inp):
        res = 1/(1+math.exp(-inp));
        return res
