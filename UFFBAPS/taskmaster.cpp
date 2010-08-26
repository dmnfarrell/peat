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

#include "taskmaster.h"

TaskMaster::TaskMaster(){}
TaskMaster::TaskMaster(string outfile, string results){prepareFiles(outfile, results);}
TaskMaster::~TaskMaster(){closeFile();}

int TaskMaster::job(vector<Soup *> A, vector<Soup *> B, string method)
{

  determineMethod(A, B, method);
 
  return 1;
}



void TaskMaster::determineMethod(vector<Soup *> A, vector<Soup *> B, string method)
{
  /*** determine method and do method specific calls ***/
  transform(method.begin(), method.end(),method.begin(),::tolower);

  if(method == "simplelj")
    SimpleLJ calc(resultsFile,reporter,convert_soups_to_objects(A),convert_soups_to_objects(B));
  else if(method == "complementarity")
    Complementarity calc(resultsFile,reporter,A,B);
  else if(method == "up_complementarity")
    Up_Complementarity calc(resultsFile,reporter,A,B);
  else if(method == "peoe")
    Peoe calc(reporter,convert_soups_to_objects(A));

  else if(method == "distancenn-train")
    Nnpot calc(resultsFile,reporter,A,B,2,true);
  else if(method == "distancenn-test")
    Nnpot calc(resultsFile,reporter,A,B,2,false);

  else if(method == "elementnn-train")
    Nnpot calc(resultsFile,reporter,A,B,1,true);
  else if(method == "elementnn-test")
    Nnpot calc(resultsFile,reporter,A,B,1,false);

  else if(method == "nn_es")
    Nnpot calc(resultsFile,reporter,A,B,3,true);
  else if(method == "nn_vdw")
    Nnpot calc(resultsFile,reporter,A,B,4,true);
  else if(method == "nn_hbond")
    Nnpot calc(resultsFile,reporter,A,B,5,true);

  else if(method == "elementnn2")
    Nnpot calc(resultsFile,reporter,A,B,6,true);

  


  else if(method == "nn_ls")
    Nnpot calc(resultsFile,reporter,A,B,7,true);
  else if(method == "genericnn")
    Nnpot calc(resultsFile,reporter,A,B,8,true);

  else if(method == "nn_fitter")
    NN_FITTER calc(resultsFile,reporter);


  else if(method == "prepin")
    PrepinReader calc(reporter,A);
  else if(method == "gp")
    Generalparameters calc(reporter,A);
  else if(method == "sp")
    Simple_Parameters calc(reporter,A);
  else if(method == "correct_names")
    Correct_Atomnames calc(reporter,A);


  else if(method == "energy")
    Energy calc(reporter,A,B);
  else if(method == "writepdb")
    Pdbwriter calc(reporter,"out.pdb", A[0]);

  else if(method == "topology")
    Topology calc(resultsFile,reporter,A,B);

  else if(method == "vdw")
    Vdw calc(resultsFile,reporter,A,B);
  else if(method == "vdw_flat")
    Vdw_Flat calc(resultsFile,reporter,A,B);
  else if(method == "vdw_line")
    Vdw_Line calc(resultsFile,reporter,A,B);
  else if(method == "vdw_line_internal")
    Vdw_Line calc(resultsFile,reporter,A);

  else if(method == "es")
    Es calc(resultsFile,reporter,A,B);
  else if(method == "es_internal")
    Es calc(resultsFile,reporter,A);


  else if(method == "hp")
    Hp calc(resultsFile,reporter,A,B);
  else if(method == "hp_internal")
    Hp calc(resultsFile,reporter,A);

  else if(method == "hbonds")
    Hbonds calc(resultsFile,reporter,A,B);
  else if(method == "hbonds_internal")
    Hbonds calc(resultsFile,reporter,A);


  else if(method == "protein_entropy")
    ProteinEntropy calc(resultsFile,reporter,A,B);

  else if(method == "ligand_entropy")
    LigandEntropy calc(resultsFile,reporter,A);
  else if(method == "ligand_entropy_desolvation")
    LigandEntropy calc(resultsFile,reporter,A,B);

  else if(method == "water_entropy")
    Water_Entropy calc(resultsFile,reporter,A);

  else if(method == "metals")
    Metal calc(resultsFile,reporter,A,B);

  else if(method == "bfactor")
    Bfactor calc(reporter,A);
  else if(method == "bfactor_binding_site")
    Bfactor calc(reporter,A,B);

  else if(method == "nn_discrimination")
    NN_DISCRIMINATION calc(resultsFile,reporter);

  else if(method == "msc")
    MSC calc(reporter, resultsFile, A, B);
  else if(method == "write_pdb")
    Write_File calc(reporter,A,string("pdb"));
  else if(method == "write_pqr")
    Write_File calc(reporter,A,string("pqr"));
  else if(method == "desolvation")
    Desolvation calc(resultsFile,reporter,A,B);
  else if(method == "desolvation_internal")
    Desolvation calc(resultsFile,reporter,A);
  else if(method == "fragmental_volume")
    Fragmental_Volume calc(reporter,A,"normal");
  else if(method == "fragmental_volume_ellipse")
    Fragmental_Volume calc(reporter,A,"ellipse");

  else if(method == "chemical_shift_changes")
    Chemical_Shift_Changes calc(reporter,A);
  else if(method == "reverse_chemical_shift_changes")
    Reverse_Chemical_Shift_Changes calc(reporter,A);

  else if(method == "backbone_entropy_train")
    BackboneEntropy calc(reporter,A,true);
  else if(method == "backbone_entropy")
    BackboneEntropy calc(reporter,A,false);

  else if(method == "water_bridges")
    Water_Bridges calc(reporter, resultsFile, A, B);
  else if(method == "water_bridges_internal")
    Water_Bridges calc(reporter, resultsFile, A);
 

  else
    printf("Can not do job. Unknown method: '%s'!\n", method.c_str());
 
}


void TaskMaster::prepareFiles(string ofile, string res)
{

  reporter = fopen(ofile.c_str(),"w");
  resultsFile = fopen(res.c_str(),"w");

  if(reporter == NULL || resultsFile == NULL)
    {
      printf("Warning: A output file could not bee opened!\n");
      exit(0);
    }

  fprintf(reporter,"This is a output file!\n");


	

}
void TaskMaster::closeFile()
{
  

}

