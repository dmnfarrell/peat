#!/bin/env python

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




import NN
import os
import NNrun

stats = open('stats.txt','w')

#net topology
nettop =(2,4,1)

#make a new network 
network = NN.NN(nettop)
network.set_learning_rate(0.1)

#change transfer function
tf = (0,0)
network.set_transfer_functions(tf)

#print network to screen
print network.print_network()


#train exclusive or 
input = ([0,0],
         [1,0],
         [0,1],
         [1,1])
expect = ([0],
          [1],
          [1],
          [0])


#input = ([1,1],)
#expect = ([1],)


#res = network.test(input)

#train a lot of times
for x in range(0,30000):
   res = network.train(input,expect)

   #print net for every 200 training rounds
   if(x%200 == 0 ):
       stats.write("%5d "%x)
       for inp in input:
           for i in inp:
               stats.write("%f "%i)
           est = network.test(inp)
           for i in est:
               stats.write("%f "%i)
           
     
       stats.write("%f "%res)
       stats.write("\n")
       print network.print_network()
stats.close()

#make plot
os.system('gnuplot xor.gnu')
#os.system('evince file.eps')

#try out NNrun
network.write_file('net-test.txt')

pythonnet = NNrun.NNpython()
pythonnet.read_file('net-test.txt')
pythonnet.print_network()
#pythonnet.test(expect[0])
print "c++ network  vs  python network  diff"
for inp in input:
   print "%6f    -    %6f      = %10f"%(pythonnet.test(inp)[0],network.test(inp)[0],pythonnet.test(inp)[0]-network.test(inp)[0])
