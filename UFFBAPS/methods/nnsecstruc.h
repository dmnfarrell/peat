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
#ifndef NNSECSTRUC_H
#define NNSECSTRUC_H
#include <map>
#include <math.h>
#include <fstream>
#include "method.h"


using namespace std;

class Nnsecstruc:  Method{

 public:
  Nnsecstruc(FILE * reporter, vector<SoupObject *> trainSet, vector<SoupObject *> testSet);
  Nnsecstruc(FILE * reporter, vector<SoupObject *> testSet);
  ~Nnsecstruc();

 private:
  void writeDescription(FILE * reporter);

  void train(vector<SoupObject *> trainSet);
  void test(vector<SoupObject *> testSet);


  void setup();
  void setWeights();
  double randomNo();
  double sgm(double in);
  void count(int size);
  void flushHistogram();
  vector<vector<double> > compute(vector<Residue> residues, int res);
  vector<double> getAnswer(int res);

  vector< vector<int> > input;
  vector<double> hiddenOutput;
  vector<double> output;
  vector<vector<double> > predictions;
  vector< vector<double> > weight11;
  vector< vector<double> > weight21;
  vector<double> weight12;
  vector<double> weight22;
  vector< vector<double> > errors;
  vector<int> histogram;
  vector<int> helixResidues;
  vector<int> sheetResidues;

  vector<vector<int> > helixLibrary;
  vector<vector<int> > sheetLibrary;

  int allright,allwrong;

  map<string, int> conversion;

  double learningRate;


};

#endif
