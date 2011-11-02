/*
 #
 # FFF - Flexible Force Field
 # Copyright (C) 2010 Jens Erik Nielsen
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
#include "main.h"

int main(int argc, char *argv[]) {
    //
    // Do some simple tests
    //
    FFF Protein;
    Protein.read_pdb("2lzt_A.pdb");
    Rotamer_class ROT("parameters/small_lib");
    //pKa_class PKA(Protein,ROT,"parameters");
    //PKA.resolve_clashes();
    //Protein.write_pdb("resolved.pdb");
    //exit(0);
    //M.repair_all();
    //M.build_hydrogens();
    //Protein.write_pdb("repaired.pdb");
    //exit(0);
//
    model_class M(Protein,ROT,"parameters");
    for (unsigned int res=1;res<10;res++) {
      M.Mutate("A:0001","PHE",3,0.5);
      printf ("modelled\n");
      Protein.write_pdb("alascan/testCpp.pdb");
      printf ("wrote pdb\n");
      M.undo_mutation();
      printf ("undid mutation\n");
      //exit(0);
  }
    exit(0);
  //
  // Test that everything is ok
  //
  //Protein.soup_stat();
  //
  // Instantiate the model class
  //
  //
  //
  // Get a chiangle
  //
  //int chain=0;
  //int residue=12;
  //int chiangle=1;
  //printf ("Chiangle for residue: %d: %f\n",residue,M.getchi(chain,residue,chiangle));
  
  //
  // Mutate a residue
  //
  //M.Mutate(chain,residue,"ALA",0);
  //M.Mutate(chain,residue,"LYS",0);
  //
  // Write PDB file
  //
  //Protein.write_pdb("mutated.pdb");
  //
  cout << argc << endl;
  if (argc==1) {
    // Interpreter mode
    cout << "Interpreter mode selected" << endl;
    Control X;
    X.interpreter();
  } else {
    // command file mode
    string line;
    vector<string> filelines;
    //printf("Starting to read command file\n");
    ifstream file(argv[1]);
    if (!file) {
      printf("Command file not found!: %s\n",argv[0]);
      exit(0);
    }
    // Put the entire file in lines
    while(getline(file,line)) {
      filelines.push_back(line);
    }
    file.close();
    //printf("Command file has been read.\n");
    Control X(filelines);
  }
  return 0;
}


