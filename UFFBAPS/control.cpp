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

#include "control.h"

Control::Control()
{
  taskMaster = TaskMaster(string("out.txt"), string("results.txt"));
  script_initiated = false;
}
Control::~Control()
{
  // write out script results if any
  if(A.size()>0)
    {
      Matrix_IO M;
      M.open("script.m");
      M.write_matrix("A", A);
      M.write_vector("b",b);
      M.write_string_vector("Names",names);
      M.close();
      cout<<"Done writting to octave file"<<endl;
    }


}

void Control::command()
{
  string commandLine;
  while(commandLine != "exit")
    {
      cout<<"> ";
      getline(cin, commandLine); 
      interpretator(commandLine);
    }
  
}


void Control::command(string line)
{
  printf("not implemented\n");
}
void Control::commands(vector<string> lines)
{
  printf("not implemented\n");
}
void Control::commands(string infile)
{
 
  /*** reads command file ***/
  printf("Control reads commands from %s.\n",infile.c_str());
  string line;
  ifstream file(infile.c_str());

  if(!file.is_open())
    {
      printf("Command file %s not found!\n",infile.c_str());
      exit(0);
    }
  while(!file.eof())
    {
      getline(file, line);
      lines.push_back(line);
    }

  for(unsigned int i=0; i<lines.size(); i++)
    interpretator(lines[i]);

}
  

void Control::interpretator(string line)
{

  /*** Transforms to lower case ***/
  //transform(line.begin(), line.end(), line.begin(),::tolower);

  /*** check for commentation ***/
  string checkComment(line,0,1);
  if(checkComment == "#" || line.size()<2)
    return;

  /*** replaces '{' and '}' with ' { ' and ' } ' ***/
  for(unsigned int i=0; i<line.size(); i++)
    {
      if(string(line,i,1) =="{")
	{
	  line.replace(i, 1, string(" { "),0,3);
	  i++;
	}
      else if(string(line,i,1)=="}")
	{
	  line.replace(i, 1, string(" } "),0,3);
	  i++;
	}
    }  

  /*** break command into words ***/
  vector<string> words;
  words.clear();
  char temp[300];
  istringstream stream(line.c_str());
  while(stream>>temp)
    words.push_back(removeWhite(string(temp)));

  /*** check for lines with only white space ***/
  if(words.size() < 1)
    return;
  /*
  printf("Control has the words:\n");
  for(int i=0;i<words.size();i++)
    printf("'%s' ",words[i].c_str());
  printf("\n");
  */
  
  try
    {
      string com = words.at(0);

      if(com == "load")
	{
	if(words.size() > 2)
	  load(words.at(1), words.at(2));
	else
	  {
	    char no[6];
	    sprintf(no, "%d", soups.size());	 
	    load(words.at(1), string("AutoSoup") + string(no));
	  }
	}
      else if(com == "loadtosoup")
	{
	  if(words.size() > 3)
	    loadToSoup(atoi(words.at(1).c_str()), words.at(2), words.at(3));
	  else
	    {
	      char no[6];
	      sprintf(no, "%d", soups.size());	 
	      load(words.at(1), string("AutoSoup") + string(no));
	    }
	}
      else if(com == "stat")
	stat();
      else if(com == "list")
	list();
      else if(com == "clear")
	clear();
      else if(com == "task")
	doTask(words);
      else if(com == "script")
	do_script(words);
      else if(com == "exit")
	return;
      else
	printf("Don't understand the command '%s'!\n",words.at(0).c_str());
    }
  catch (...)
    {
      printf("The command %s caused a critical error!\n",line.c_str());
    }


}


void Control::doTask(vector<string> words)
{
 
  /*** making lists of soups to use with method ***/
  
  int listNo=0;
  
  for(unsigned int i=2; i<words.size();i++)
    {
      
      /*** Start next list ***/
      if(words.at(i) =="{")
	listNo++;
      else if(words.at(i) !="}")
	addSoupToList(words.at(i),listNo);
      
      if(listNo>2)
	{
	  printf("Too many list found!\n");
	}
    }
  /*** call taskMaster ***/
  
  printf("Created lists:\n");
  for(unsigned int i=0;i<list1.size();i++)
    printf("  Soup %d of list 1 is: %s\n",i,list1.at(i)->name.c_str());
  for(unsigned int i=0;i<list2.size();i++)
    printf("  Soup %d of list 2 is: %s\n",i,list2.at(i)->name.c_str());
  
  taskMaster.job(list1, list2, words.at(1));

  list1.clear();
  list2.clear();
}


void Control::do_script(vector<string> words)
{
  
  //check if scripting is initiated
  if (!script_initiated)
    initiate_scripting();
  
  printf("\n");
  printf("  Script mode - Note that absolute stability values do not make sense!\n");
 
  string method = "",sign="";
  Soup this_soup,this_soup2;
  float result = 0, this_result = 0;
  int no_soups = 1;
  bool name_is_set = 0;
  vector<double> components(Stability.get_no_components(),0); //this line could be problematic if binding is included in scripting


  for(unsigned int i=1; i<words.size();i++)
    {
      //cout<<">"<<words.at(i)<<endl;
      // find out what the sign of the term is
      if (words.at(i) =="+" or words.at(i) =="-")
	{
	  if (method == "") // if method is set it means that we are adding soups not energy terms
	    sign = words.at(i);
	}
      // check which method to use for this term
      else if (words.at(i) =="stability" or words.at(i) =="binding" or words.at(i) =="binding_just_PL")
	{
	  method = words.at(i);
	  //printf("Got method %s\n",method.c_str());
	}

      // calculate the term
      else if(words.at(i) =="}")
	{
	  vector<double> these_components(8,0);
	  if (method == "stability")
	    {
	      this_result = Stability.calculate_stability(&this_soup);
	      printf("  Stability: %-40s       %1s %8.2f kJ/mol\n",this_soup.name.c_str(),sign.c_str(), this_result);
	      these_components = Stability.get_components();
	    }
	  else if (method == "binding_just_PL")
	    {
	      if (no_soups == 1)
		{
		  no_soups = 2;
		  //cout<<"Starting second soup"<<endl;
		  continue;
		}
	      
	      this_result = Binding.calculate_binding(&this_soup, &this_soup2);
	      printf("  Binding:   %-20s and %-20s  %1s %8.2f kJ/mol\n",this_soup.name.c_str(),this_soup2.name.c_str(),sign.c_str(), this_result);
	      these_components = Binding.get_components();
	    }
	  else if (method == "binding")
	    {
	      if (no_soups == 1)
		{
		  no_soups = 2;
		  //cout<<"Starting second soup"<<endl;
		  continue;
		}
	      
	      Soup holo, apo;
	      apo = this_soup;
	      holo = this_soup + this_soup2;
	      
	      float stab_apo, stab_holo, ligand_entropy;
	      vector<double> components_apo, components_holo;
	      stab_apo = Stability.calculate_stability(&apo);
	      components_apo  = Stability.get_components();
	      stab_holo = Stability.calculate_stability(&holo);
	      components_holo = Stability.get_components();

	      ligand_entropy = Ligand_entropy.calculate_using_desolvation(&apo, &this_soup2);

	      this_result = stab_holo - stab_apo + ligand_entropy;

	      printf("  Binding_ds:%-20s and %-20s  %1s %8.2f kJ/mol\n",this_soup.name.c_str(),this_soup2.name.c_str(),sign.c_str(), this_result);
	      

	      for (int k =0; k<components_holo.size(); k++)
		these_components[k] = components_holo[k] - components_apo[k];
	      these_components[6] = ligand_entropy; //include ligand entropy
	    }

	  
	  

	  

	  if (sign=="+" or sign =="")
	    {
	      result += this_result;
	      for(unsigned int c=0;c<components.size();c++)
		components[c] += these_components[c];
	    }
	  else if (sign=="-")
	    {
	      result -= this_result;
	      for(unsigned int c=0;c<components.size();c++)
		components[c] -= these_components[c];

	    }

	  this_soup.clearSoup();
	  this_soup.name = "";
	  this_soup2.clearSoup();
	  this_soup2.name = "";
	  method = "";
	  no_soups = 1;

	}
      // set target value
      else if(words.at(i) =="target")
	{
	  // printf("got target: %s\n",words.at(i+1).c_str());
	  b.push_back(atof(words.at(i+1).c_str()));
	  i++;//jump
	}
      // set calculation name
      else if(words.at(i) =="name")
	{
	  // printf("got name: %s\n",words.at(i+1).c_str());
	  names.push_back(words.at(i+1));
	  i++;//jump
	  name_is_set = 1;
	}

      // find soups to use for this term
      else if(words.at(i) !="{" and words.at(i) !="}")
	{
	  //      	  cout<<"Looking for "<<words.at(i)<<endl;
	  bool found = false;
	  for(unsigned int s=0; s<soups.size(); s++)
	    if(soups[s].name == words.at(i))
	      {
		//cout<<"adding to soup "<<no_soups<<endl;
		if (no_soups == 1)
		  this_soup = this_soup + soups[s];
		else if(no_soups == 2)
		  this_soup2 = this_soup2 + soups[s];
		found = true;
	      }
	  if (!found)
	    {
	      printf("Script error: invalid soup name: %s\n",words.at(i).c_str());
	      exit(0);
	    }
	}



    }

  if (name_is_set != 1)
    names.push_back("");

  printf("  ---------------------------------------------------------------------------\n");
  printf("  Result for %-30s                   %8.2f kJ/mol\n",names.at(names.size()-1).c_str(),result);
  printf("  ===========================================================================\n");
  A.push_back(components);

  
}

void Control::initiate_scripting()
{
  Stability.initiate();
  Binding.initiate();

  script_initiated = true;
  
  return;
}


void Control::addSoupToList(string name, int list) //type: 1:prot, 2:mol, 3:water
{

  Soup * temp;
  bool found = false;

  /*** find soup in question ***/
  for(unsigned int s=0; s<soups.size(); s++)
    if(soups[s].name == name)
      {
	temp = &soups[s];
	found = true;
      }
  /*** add object to correct list ***/
  if(found)
    {
      if(list==1)
	list1.push_back(temp);
      else if(list==2)
	list2.push_back(temp);
    }
  else
    printf("Soup '%s' could not be found!\n",name.c_str());
}

void Control::load(string file, string name)
{
  //model control
  //if no models in file, make one soup of entire file and push back in all

  /*** read pdb file ***/
  vector<string> lines = readFile(file);
  vector<vector<string> > modelLines;
  string::size_type found_pdb = file.find(".pdb",0);
  string::size_type found_mol2 = file.find(".mol2",0);

  if(found_pdb != string::npos)
    modelLines = divide_pdb(lines);
  else if(found_mol2 != string::npos)
    modelLines = divide_mol2(lines);
  else
    printf("Unknown file extension of file %s !", file.c_str());

  if(modelLines.size() != 0)
    {
      printf("%d models found in %s - creating %d soups\n",modelLines.size(),file.c_str(),modelLines.size());
      for(unsigned int i=0;i<modelLines.size();i++)
	{
	  Soup X;
	  X.addToSoup(modelLines[i],file);

	  char no[6];
	  sprintf(no, "%d", i+1);
	  
	  X.name = name + string("_") + string(no);
	  soups.push_back(X);
	}
    }
  else
    {
      printf("1 model found in %s - creating 1 soup\n",file.c_str());
      Soup X;
      X.addToSoup(lines,file);
      X.name = name;
      soups.push_back(X);
    }
  //if models present make one soup for each model and push back in all
}



void Control::loadToSoup(int offset, string file, string name)
{
  /*** add contents of file to soups, beginning at soup offset ***/

  /*** read pdb file ***/
  vector<string> lines = readFile(file);
  vector<vector<string> > modelLines;
  /*** divide file into models ***/
  string::size_type found_pdb = file.find(".pdb",0);
  string::size_type found_mol2 = file.find(".mol2",0);
  
  if(found_pdb != string::npos)
    modelLines = divide_pdb(lines);
  else if(found_mol2 != string::npos)
    modelLines = divide_mol2(lines);
  else
    printf("Unknown file extension of file %s !", file.c_str());


  if(modelLines.size() != 0)
    {
      printf("%d models found in %s - adding to %d soups beginning at soup %d\n",
	     modelLines.size(),file.c_str(),modelLines.size(),offset);
      for(unsigned int i=0;i<modelLines.size();i++)
	{
	  if(offset <(int) soups.size()) //add to existing soup
	    {
	      soups.at(offset).addToSoup(modelLines[i],file);
	    }
	  else //make new soup
	    {
	      Soup X;
	      X.addToSoup(modelLines[i],file);

	      char no[6];
	      sprintf(no, "%d", i+1);
	      
	      X.name = name + string("_") + string(no);
	      soups.push_back(X);

	    }
	  offset++;

	}
    }
  else
    {
      printf("1 model found in %s - adding to soup %d\n",file.c_str(),offset);
      soups.at(offset).addToSoup(lines,file);
    }

}



vector<vector<string> > Control::divide_pdb(vector<string> lines)
{
  /*** divide file into models ***/
  
  vector<string> tempLines;
  vector<vector<string> > modelLines;
  bool go = false;
  for(unsigned int i=0;i<lines.size();i++)
    {
      string identifier(lines[i],0,6);
      if(go)
	tempLines.push_back(lines[i]);

      if(identifier == "MODEL ")
	go = true;

      if(identifier == "ENDMDL")
	{
	  go = false;
	  modelLines.push_back(tempLines);
	  tempLines.clear();
	}
    }
  /*** in case somebody forgot last ENDMDL tag ***/
  if(tempLines.size() != 0)
      modelLines.push_back(tempLines);

  return modelLines;
}

vector<vector<string> > Control::divide_mol2(vector<string> lines)
{
  /*** divide file into models ***/
  
  vector<string> tempLines;
  vector<vector<string> > modelLines;
  bool go = false;
  for(unsigned int i=0;i<lines.size();i++)
    {
      string identifier(lines[i],0,17);

      if(go)
	tempLines.push_back(lines[i]);


      if(identifier == "@<TRIPOS>MOLECULE" && go)
	{
	  modelLines.push_back(tempLines);
	  tempLines.clear();
	}
      if(identifier == "@<TRIPOS>MOLECULE" && !go)
	go = true;

    }
  /*** including the last molecule ***/
  if(tempLines.size() != 0)
      modelLines.push_back(tempLines);

  /*** if only one molecule present there is no reason to divide into more than one soup ***/
  if(modelLines.size() == 1)
    modelLines.clear(); 

  return modelLines;
}


void Control::stat()
{
  for(unsigned int i=0;i<soups.size();i++)
    soups.at(i).soupStat();
}


void Control::list()
{  
  for(unsigned int i=0;i<soups.size();i++)
    soups.at(i).listAtoms();
}

void Control::clear()
{
  printf("Clearing soup\n");
  soups.clear();
  //  for(int i=0;i<models.size();i++)
  //  models.at(i).clearSoup();
}

vector<string> Control::readFile(string filename)
{
  vector<string> lines;
  string line;

  ifstream file(filename.c_str());

  if(!file.is_open())
    {
      printf("File %s not found!\n",filename.c_str());
      exit(0);
    }
  
  while(!file.eof())
    {
      getline(file, line);
      lines.push_back(line);
    }

  return lines;
}
