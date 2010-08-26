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


#include "targets.h"


Targets::Targets(FILE * reporter, string dataset):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  T = 300;// K
  //R = 1.9872 / 1000;// kcal/(K mol)
  R = 8.314472/1000; // kJ/(K mol)
  log_conversion = 2.302585093;

  if (dataset == "wang")
    read_file_wang("/home/people/chresten/run/docksets/wang/all_original_dockset/complex_list_all");
  else if(dataset == "pdbbind")
    read_file_pdbbind("/home/people/chresten/run/docksets/pdbbind/original/INDEX.2004.refined.pkd");
  else if(dataset == "pdbbind2007")
    read_file_pdbbind("/rnase_local/chresten/run/docksets/pdbbind2007//INDEX.2007.refined.data");
  else if(dataset == "fictive")
    read_file_fictive("targets.txt");

}

Targets::~Targets(){}


float Targets::get_mlogKi(string key)
{
  int pos = find_key(key);
  if(pos == -1)
    {
      printf("Warning: Key '%s' not found in target keys!\n",key.c_str());
      return 0.0;   
    }
  else
    return mlogKis[pos];
}


float Targets::get_G(string key)
{
  int pos = find_key(key);
  if(pos == -1)
    {
      printf("Warning: Key '%s' not found in target keys!\n",key.c_str());
      return 0.0;   
    }
  else
    return Gs[pos];
}


float Targets::get_decoy_G(string key, float rmsd)
{

  //calculates an 'energy' using the energy of the experimentel biding mode
  //and a rmsd value of the decoy

  int pos = find_key(key);
  if(pos == -1)
    {
      printf("Warning: Key '%s' not found in target keys!\n",key.c_str());
      return 0.0;   
    }
  
  //  return Gs[pos] + 0.05*rmsd;
  if (rmsd < 0.00000001)
    return 1.0;
  else
    return 0.0;


}

float Targets::get_rmsd(string key, int decoy_no)
{
  //looks up the rmsd using a pdbid and a decoy no.

  string filename = "/home/people/chresten/run/docksets/wang/all_original_dockset/summary/"+key+"_result.all";
  float res = -1;

  ifstream file (filename.c_str());

  if(!file.is_open())
    {
      printf("Error: Could not open file '%s'!\n",filename.c_str());
      return 0.0;
    }

  string line;
  while(file)
    {
      getline(file, line);
      
      if(atoi(line.substr(0,3).c_str()) == decoy_no)
	res = atof(line.substr(4,10).c_str());

    }
  file.close();
  
  if(res == -1)
    printf("WARNING: Could not find decoy no '%d' of '%s'!\n",decoy_no, key.c_str());

  // printf("rmsd: %f\n",res);
  return res;

}




float Targets::get_resolution(string key)
{
  int pos = find_key(key);
  if(pos == -1)
    {
      printf("Warning: Key '%s' not found in target keys!\n",key.c_str());
      return 0.0;   
    }
  else
    return resolutions[pos];
}

float Targets::get_no_rotor(string key)
{
  int pos = find_key(key);
  if(pos == -1)
    {
      printf("Warning: Key '%s' not found in target keys!\n",key.c_str());
      return 0.0;   
    }
  else
    return no_rotors[pos];
}



int Targets::find_key(string key)
{
  int pos = -1;
  for(unsigned int i=0; i<keys.size(); i++)
    if(keys[i] == key)
      pos = i;

  return pos;
}




void Targets::read_file_wang(string filename)
{
  ifstream file (filename.c_str());

  if(!file.is_open())
    {
      printf("Error: Could not open file '%s'!\n",filename.c_str());
      return;
    }

  string line;
  while(file)
    {
      getline(file, line);
      
      if(line.substr(0,1) != "#" && line.size() > 4)
	{

	  string temp_key   = line.substr(0,4);
	  float temp_mlogKi = atof(line.substr(5,11).c_str());
	  float temp_resolution  = atof(line.substr(14,18).c_str());
	  int temp_no_rotor = atoi(line.substr(19,21).c_str());

	  float temp_G = -log_conversion*R*T*temp_mlogKi;

	  keys.push_back(temp_key);
	  mlogKis.push_back(temp_mlogKi);
	  resolutions.push_back(temp_resolution);
	  no_rotors.push_back(temp_no_rotor);
	  Gs.push_back(temp_G);

	}
    }//end while

  file.close();


  printf("key  -logKi  resolution  #rotor \n");
  for(unsigned int i=0; i<keys.size();i++)
    printf("%s %f %f %d\n",keys[i].c_str(),mlogKis[i],resolutions[i],no_rotors[i]);
  



  printf("Affinity data for use in database - wang (temp)\n");
 

  for(unsigned int i=0; i<keys.size();i++)
    printf("%-4s   -1 %8.3f\n",
	   keys[i].c_str(),
	   get_G(keys[i]));	
  

  


  /*  printf("Decoy Data for use in database (temp)\n");
 

  for(unsigned int i=0; i<keys.size();i++)
    for(unsigned int d=1; d<=101; d++)
      printf("%-4s %4d %8.3f\n",
	     keys[i].c_str(),
	     d,
	     get_decoy_G(keys[i], get_rmsd(keys[i], d)));	
  

  */
  /*      printf("%-4s %-4d %8f %8f\n",
	     keys[i].c_str(),
	     d,
	     get_decoy_G(keys[i], get_rmsd(keys[i], d)),
	     get_rmsd(keys[i], d));	
  
  */

}


void Targets::read_file_pdbbind(string filename)
{
  ifstream file (filename.c_str());

  if(!file.is_open())
    {
      printf("Error: Could not open file '%s'!\n",filename.c_str());
      return;
    }

  string line;
  while(file)
    {
      getline(file, line);
      
      if(line.substr(0,1) != "#" && line.size() > 4)
	{

	  string temp_key   = line.substr(0,4);
	  float temp_resolution = atof(line.substr(5,11).c_str());
	  float temp_mlogKi  = atof(line.substr(18,23).c_str());

	  float temp_G = -log_conversion*R*T*temp_mlogKi;
	  keys.push_back(temp_key);
	  mlogKis.push_back(temp_mlogKi);
	  resolutions.push_back(temp_resolution);
	  Gs.push_back(temp_G);

	}
    }//end while

  file.close();

  /*
  printf("Affinity data for use in database - pdbbind (temp)\n");
 

  for(unsigned int i=0; i<keys.size();i++)
    printf("%-4s   -1 %8.3f\n",
	   keys[i].c_str(),
	   get_G(keys[i]));	
  */
//   printf("key  -logKi  resolution  \n");
//   for(unsigned int i=0; i<keys.size();i++)
//     printf("%s %f %f\n",keys[i].c_str(),mlogKis[i],resolutions[i]);
  
//   printf("Data for use in database (temp)\n");
//   for(unsigned int i=0; i<keys.size();i++)
//     printf("%-10s %5.3f\n",keys[i].c_str(),Gs[i]);	
  
}


void Targets::read_file_fictive(string filename)
{
  ifstream file (filename.c_str());

  if(!file.is_open())
    {
      printf("Error: Could not open file '%s'!\n",filename.c_str());
      return;
    }

  string line;
  while(file)
    {
      getline(file, line);
      
      if(line.substr(0,1) != "#" && line.size() > 4)
	{

	  string temp_key   = line.substr(0,4);
	  float temp_G = atof(line.substr(5,16).c_str());

	  float temp_mlogKi = 0.0;
	  float temp_resolution  = 0.0;
	  int temp_no_rotor = 0;


	  keys.push_back(temp_key);
	  mlogKis.push_back(temp_mlogKi);
	  resolutions.push_back(temp_resolution);
	  no_rotors.push_back(temp_no_rotor);
	  Gs.push_back(temp_G);

	}
    }//end while

  file.close();


  printf("key  G\n");
  for(unsigned int i=0; i<keys.size();i++)
    printf("%s %f\n",keys[i].c_str(),Gs[i]);
  
}





void Targets::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: targets.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
