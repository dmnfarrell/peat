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


#include "up_complementarity.h"


Up_Complementarity::Up_Complementarity(FILE * resultsFile, FILE * reporter, vector<Soup*> A, vector<Soup*> B):Method(reporter)
{
  set_parameters();

  double result =0;

  bool include = true;
  /*** write out to out file **/
  writeDescription(reporter);
  Distances calc(reporter);

  for(unsigned int i=0; i< A.size(); i++ )
    {
      calc.calculate(A[i], B[i], false);
      vector<vector<float> > d =calc.getResult();

      vector<Atom *> Aatoms = A[i]->getAtoms();
      vector<Atom *> Batoms = B[i]->getAtoms();

      for(unsigned int x=0; x<d.size(); x++)
	for(unsigned int y=0; y<d[0].size(); y++)
	  {
	    include = true;
	    
	    if(ignore_h == 1 && Aatoms[x]->element =="H")
	      include = false;
	    if(ignore_h == 1 && Batoms[y]->element =="H")
	      include = false;
	    
	    if(Aatoms[x]->element =="C")
	      include = false;
	    if(Batoms[y]->element =="C")
	      include = false;

	    //	    cout<<"A:"<<Aatoms[x]->element<<" B:"<<Batoms[y]->element<<" Dist:"<<sqrt(d[x][y]);
	    if(d[x][y] < tsq && include)
	      {
		result += complEnergy;
		//		cout<<" INCLUDED";
	      }
	    //	    cout<<endl;
	  }  

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

  //  fprintf(reporter,"The result is: %f kcal/mol\n", result);	
  // fprintf(resultsFile,"%-10s %5.3f\n",A[0]->identify().substr(0,4).c_str(),result);	
  
}




Up_Complementarity::~Up_Complementarity(){}


















void Up_Complementarity::set_parameters()
{

  string dummy;

  file.open("up_complementarity.cfg");

  if(file.is_open())
    {
      while(!file.eof())
	{
	  file >> dummy;
	  if(dummy == "THRESHOLD:")
	    file >>	threshold;
	  if(dummy == "COMPLENERGY:")
	    file >> complEnergy;
	  if(dummy == "IGNORE_H")
	    file >> ignore_h;
	}
    }
  else
    {
      threshold = 15;
      complEnergy = -0.00008; //-80 mCal/mol
      ignore_h = 0;
      printf("Using standard values:threshold is %f and complementary energy is %f, ignore_h is %d\n",threshold,complEnergy, ignore_h);
    }



  printf("Configurations for Complementarity:\n");
  printf("\n");
  printf("\n");

  printf("THRESHOLD:              %f\n", threshold);
  printf("COMPLENERGY:            %f\n", complEnergy);
  printf("IGNORE_H:               %d\n", ignore_h);
  printf("\n");

  printf("Version: $Id: up_complementarity.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");






  tsq = threshold*threshold;
}



void Up_Complementarity::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: up_complementarity.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
