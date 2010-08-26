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
$Id: nn_discrimination_preparation.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "nn_discrimination_preparation.h"


NN_DISCRIMINATION_PREPARATION::NN_DISCRIMINATION_PREPARATION(FILE * resultsFile,FILE * reporter)
{

  /*** write out to out file **/
  writeDescription(reporter,resultsFile);
 
  go();

}

NN_DISCRIMINATION_PREPARATION::~NN_DISCRIMINATION_PREPARATION(){}


void NN_DISCRIMINATION_PREPARATION::build_inputs()
{
}


void NN_DISCRIMINATION_PREPARATION::set_parameters()
{


  ifstream parameters;
  parameters.open("nn_discrimination_preparation.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'nn_discrimination_preparation.cfg' could not be found\n";
      exit(0);
    }

  printf("Setting parameters\n");

  string dummy;

  while (!parameters.eof())
    {

      parameters >> dummy;

      if(dummy == "DATASET:")
	parameters >> dataset;

    } //end read file

  //printing all configurations

  printf("Configurations for fitter NN:\n");
  printf("\n");
  printf("\n");

  printf("DATASET:              %s\n",dataset.c_str());

  printf("\n");
  printf("\n");
  printf("Version: $Id: nn_discrimination_preparation.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");

}


void NN_DISCRIMINATION_PREPARATION::writeDescription(FILE * reporter, FILE * resultsFile)
{

  rep=reporter;
  resultsF=resultsFile;
  /*** write out to out file **/
  version =string("$Id: nn_discrimination_preparation.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
