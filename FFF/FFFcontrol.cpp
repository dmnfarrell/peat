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
#include "FFFcontrol.h"

Control::Control(vector<string> commandlines) {
  //
  // Execute the commands in teh file
  //
  cout << "here instead" << endl;
  running=1;
  for(unsigned int i=0; i<commandlines.size(); i++) {
    run_command(commandlines[i]);
  }
  if (running==1) {
    interpreter();
  }
  return;
}

Control::Control() {
  // 
  // Pass control to the interpreter
  //
  return;
}

void Control::interpreter() {
  //
  // Take input from std. in
  //
  string command;
  running = 1;
  while (running==1) {
    cout << "FFF> ";
    getline (cin, command);
    cout << "echo command" << command << endl;
    run_command(command);
  }
  return;
}

void Control::run_command(const std::string command) {
  //
  // Run a command
  //
  cout << "running:" << command << endl;
  if (command=="exit" or command=="bye" or command=="quit") {
    cout << "GoodBye" << endl;
    running=0;
  }
  //
  // Split the string
  //
  string object, name;
  object="";
  name="";
  vector<string> args=split(command);
  printf ("Command: %s, ",command.c_str());
  printf ("I found %d words\n",static_cast<int>(args.size()));
  string action=args[0];
  //
  // Execute the command - I'm sure this could be made nicer
  //
  // Load: Load <filename> <protein name>
  if (action=="Load") {
    Proteins.push_back(new FFF);
      unsigned int protnumber=Proteins.size()-1;
      (*(Proteins[protnumber])).read_pdb(args[1]);
    Protein_names.push_back(args[2]);
    
    Rotamer_class ROT("../Protool/bbdep02.May.sortlib");
    Model_instance.push_back(new model_class(*Proteins[protnumber],ROT,"AA.DAT"));
  } 
  // Soup stat
  else if (action=="Soup_stat") {
    (*Proteins[(Proteins.size()-1)]).soup_stat();
  }
  // Mutate: Mutate <protein_name> <chain> <residue> <new residue>
  else if (action=="Mutate") {
    printf ("Mutating now\n");
    
  }
}
