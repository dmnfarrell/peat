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
#ifndef TASKMASTER_H
#define TASKMASTER_H
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include "soupobject.h"
#include "soupmanip.h"
#include "soup.h"
#include "atom.h"
#include "methods/simplelj.h"
#include "methods/distances.h"
#include "methods/peoe.h"
#include "methods/nnsecstruc.h"
#include "methods/complementarity.h"
#include "methods/up_complementarity.h"
#include "methods/nnpot.h"
#include "methods/prepinreader.h"
#include "methods/generalparameters.h"
#include "methods/simple_parameters.h"
#include "methods/energy.h"
#include "methods/vdw.h"
#include "methods/vdw_flat.h"
#include "methods/vdw_line.h"
#include "methods/es.h"
#include "methods/hbonds.h"
#include "methods/pdbwriter.h"
#include "methods/protein_entropy.h"
#include "methods/ligand_entropy.h"
#include "methods/water_entropy.h"
#include "methods/correct_atomnames.h"
#include "methods/topology.h"
#include "methods/metal.h"
#include "methods/hp.h"
#include "methods/bfactor.h"
#include "methods/nn_fitter.h"
#include "methods/genericnn.h"
#include "methods/nn_discrimination.h"
#include "methods/mutant_stability_change.h"
#include "methods/write_file.h"
#include "methods/desolvation.h"
#include "methods/fragmental_volume.h"
#include "methods/chemical_shift_changes.h"
#include "methods/reverse_chemical_shift_changes.h"
#include "methods/backbone_entropy.h"
#include "methods/water_bridges.h"

using namespace std;

class TaskMaster{

 public:
  TaskMaster();
  TaskMaster(string outfile, string results);
  ~TaskMaster();

  int job(vector<Soup *> A, vector<Soup *> B, string method);

 private:

  FILE * reporter;
  FILE * resultsFile;
  void determineMethod(vector<Soup *> A, vector<Soup *> B, string method);
  void prepareFiles(string ofile, string res);
  void closeFile();

};

#endif
