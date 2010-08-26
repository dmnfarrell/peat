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

/*
$Id: metal.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "metal.h"


Metal::Metal(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  
  set_parameters();
}


Metal::Metal(FILE * resultsFile, FILE * reporter, vector<Soup *> A, vector<Soup *> B):Method(reporter)
{

  /*** write out to out file **/
  writeDescription(reporter);
 
  set_parameters();

  if(A.size() != B.size())
    {
      cout << "WARNING: An equal amount of protein and ligand structures must be given to the Metal function!\n";
      return;
    }

  float result;
  vector<vector<float> > distances;

  Distances dist(reporter);


  for(unsigned int i=0; i<A.size(); i++)
    {
      dist.calculate(A[i], B[i], true);
      distances = dist.getResult();

      result = calculate(A[i],B[i],distances);

      string decoy_number;
      if(B[i]->name.find_last_of("_")== string::npos)
      	decoy_number = "-1";
      else
	decoy_number = B[i]->name.substr(B[i]->name.find_last_of("_")+1);

      fprintf(resultsFile,"%-4s %-4s %8f\n",
	      A[i]->name.substr(0,4).c_str(),
	      decoy_number.c_str(),
	      result);


    }  

}


Metal::~Metal(){}

float Metal::calculate(Soup * protein, Soup * ligand, vector<vector<float> > distances)
{

  vector<Atom*> protein_atoms = protein->getAtoms();
  vector<Atom*> ligand_atoms = ligand->getAtoms();
  vector<string>::iterator result;  
  float threshold = 10.0;
  float res = 0;

  for(unsigned int p=0; p<protein_atoms.size(); p++)
    {
      result = find(metals.begin(), metals.end(), protein_atoms[p]->element);
      if(result != metals.end()) 
	{
	  printf("Metal %s-%d has the interactions:\n",protein_atoms[p]->element.c_str(),p+1);
	  for(unsigned int i=0; i<distances[p].size(); i++)
	    {
	      if(distances[p][i]<threshold)
		printf(" %s-%d at the distance %f\n",ligand_atoms[i]->element.c_str(),i+1, distances[p][i]);
	    }
	  printf("--------------\n");

	}              


    }              
    
  return res;
}










void Metal::set_parameters()
{
  
  metals.push_back("F");
  metals.push_back("NA");
  metals.push_back("CA");

  metals.push_back("CL");
  metals.push_back("CA");
  metals.push_back("BR");
  metals.push_back("I");
  metals.push_back("P");


}







void Metal::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: metal.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
