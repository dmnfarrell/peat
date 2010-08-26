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


#include "vdw.h"


Vdw::Vdw(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  


}


//constructor so that vdw can be started directly without using energy

Vdw::Vdw(FILE * resultsFile, FILE * reporter, vector<Soup *> A, vector<Soup *> B):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  /*** read parameters ***/
  PrepinReader prepA(reporter, A);
  PrepinReader prepB(reporter, B);



  Distances dist(reporter);
  double result = 0;
  vector<vector<float> > distances;  

  if(A.size() != B.size())
    {
      cout << "WARNING: An equal amount of protein and ligand structures must be given to the van der Waals function!\n";
      return;
    }
  
  for(unsigned int i=0; i<A.size(); i++)
    {
      dist.calculate(A[i], B[i], true);
      distances = dist.getResult();

      result = calculate(A[i],B[i],distances,0);

      string decoy_number;
      if(B[i]->name.find_last_of("_")== string::npos)
      	decoy_number = "-1";
      else
	decoy_number = B[i]->name.substr(B[i]->name.find_last_of("_")+1);

      fprintf(resultsFile,"%-4s %-4s %8f\n",
	      A[i]->name.substr(0,4).c_str(),
	      decoy_number.c_str(),
	      result);


      //      fprintf(resultsFile,"%-10s %5.3f\n",convert_soup_to_objects(A[i])[0]->identify().substr(0,4).c_str(),result);
    }  



}




Vdw::~Vdw(){}


double Vdw::calculate(Soup * A, Soup * B, vector<vector<float> > distances, double cutoff){
  vector<Atom *> atomsA = A->getAtoms();  
  vector<Atom *> atomsB = B->getAtoms();  

  result = 0;

  /*** No cutoff - set cutoff to a silly and high value ***/
  if(cutoff == 0)
    cutoff = 100000000;

  for(unsigned int ai=0; ai< atomsA.size(); ai++)
    for(unsigned int bi=0; bi< atomsB.size(); bi++) 
      if(distances[ai][bi] < cutoff)
	result += calculate_vdw(atomsA[ai], atomsB[bi], distances[ai][bi]);

   return result;

}


double Vdw::calculate_vdw(Atom * atom1, Atom * atom2, double dist)
{
  double R = atom1->vdw_dist + atom2->vdw_dist;
  double E = sqrt(atom1->vdw_energy * atom2->vdw_energy);

  double ratio = R/dist;
  double ratio2 = ratio*ratio;
  double ratio6 = ratio2*ratio2*ratio2;
  double ratio12 = ratio6*ratio6;

  double res = E*(ratio12 - 2*ratio6);
  printf("R: %f, E: %f, dist: %f, res: %f\n",R,E,dist,res);
  exit(0);

  return res;
}


void Vdw::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: vdw.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("Calculates van der Waals interactions");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
