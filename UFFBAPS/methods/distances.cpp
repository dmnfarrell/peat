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
$Id: distances.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "distances.h"

Distances::Distances(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
}

Distances::Distances():Method(){}

Distances::~Distances(){}

void Distances::calculate(SoupObject * A, bool doSqrt)
{
  //
  // Calculates intra-atomic distances for a soup object
  //
 

  calculate(A->getAtoms(), doSqrt);

}

void Distances::calculate(vector<Atom*> atoms, bool doSqrt)
{
  //
  // Calculates intra-atomic distances for a vector of atoms
  //

  int noAtoms = atoms.size();
  vector<float> temp;

  temp.assign(noAtoms, 0);
  result.assign(noAtoms, temp);

  

  /*** calculate sqrt'ed distances ***/
  for(int i=0; i< noAtoms; i++)
    for(int j=0; j< noAtoms; j++)
      {
	result[i][j] = 
	  (atoms[i]->x - atoms[j]->x)*(atoms[i]->x - atoms[j]->x)+
	  (atoms[i]->y - atoms[j]->y)*(atoms[i]->y - atoms[j]->y)+
	  (atoms[i]->z - atoms[j]->z)*(atoms[i]->z - atoms[j]->z);
      }

  /*** do sqrt ***/
  if(doSqrt)
    for(int i=0; i< noAtoms; i++)
      for(int j=0; j< noAtoms; j++)
	result[i][j] = sqrt(result[i][j]);



}


void Distances::calculate(SoupObject * A, SoupObject * B, bool doSqrt)
{
  //
  // Calculates inter-atomic distances for two soup objects
  //


  vector<Atom *> atomsA = A->getAtoms();
  vector<Atom *> atomsB = B->getAtoms();

  calculate(atomsA, atomsB, doSqrt);
}



void Distances::calculate(Soup * A, Soup * B, bool doSqrt)
{
  //
  // Calculates inter-atomic distances for two soups
  //


  vector<Atom *> atomsA = A->getAtoms();  
  vector<Atom *> atomsB = B->getAtoms();  

  calculate(atomsA, atomsB, doSqrt);
}




void Distances::calculate(vector<SoupObject *> A, vector<SoupObject *> B, bool doSqrt)
{
  //
  // Calculates inter-atomic distances for two vectors of soup objects
  //


  vector<float> temp;
  int noAtomsA = 0;
  int noAtomsB = 0;
  
  printf("Calculating distances this set of soup objects:\n");
  for(unsigned int i=0; i<A.size(); i++)
    {
      printf("   %s which has %d atoms\n",A[i]->identify().c_str(),A[i]->numberOfAtoms());
      noAtomsA += A[i]->getAtoms().size();
    }
  printf("and this set:\n");
  for(unsigned int i=0; i<B.size(); i++)
    {
      printf("   %s which has %d atoms\n",B[i]->identify().c_str(),B[i]->numberOfAtoms());
      noAtomsB += B[i]->getAtoms().size();
    }
  

  temp.assign(noAtomsB, 0);
  result.assign(noAtomsA, temp);

  int i=0;
  int j=0;

  /*** calculate sqrt'ed distances ***/
  for(unsigned int obA=0; obA<A.size(); obA++)
    {
      vector<Atom * > atomsA =  A[obA]->getAtoms();
      int number_atoms_A = atomsA.size();
      for(int noA=0; noA<number_atoms_A; noA++)
	{
	  for(unsigned int obB=0; obB<B.size(); obB++)
	    {
	      vector<Atom * > atomsB =  B[obB]->getAtoms();
	      int number_atoms_B = atomsB.size();
	      for(int noB=0; noB<number_atoms_B; noB++)
		{
		  result.at(i).at(j) = 
		    (atomsA[noA]->x - atomsB[noB]->x)*(atomsA[noA]->x - atomsB[noB]->x)+
		    (atomsA[noA]->y - atomsB[noB]->y)*(atomsA[noA]->y - atomsB[noB]->y)+
		    (atomsA[noA]->z - atomsB[noB]->z)*(atomsA[noA]->z - atomsB[noB]->z);

		  j++;
		}
	    }
	  i++;
	  j=0;
	}
    }

  /*** do sqrt ***/

  if(doSqrt)
    for(int i=0; i< noAtomsA; i++)
      for(int j=0; j< noAtomsB; j++)
	result[i][j] = sqrt(result[i][j]);
    
}



void Distances::calculate(vector<Atom *> atomsA, vector<Atom *> atomsB, bool doSqrt)
{
  //
  // Calculates inter-atomic distances for two vectors of atoms
  //




  vector<float> temp;
  int noAtomsA = atomsA.size();
  int noAtomsB = atomsB.size();


  //  cout <<"Calculating distances for "<< noAtomsA<<" X "<< noAtomsB<<endl;

  temp.assign(noAtomsB, 0);
  result.assign(noAtomsA, temp);

  /*** calculate sqrt'ed distances ***/
  for(int i=0; i< noAtomsA; i++)
    for(int j=0; j< noAtomsB; j++)
      {
	result[i][j] = 
	  (atomsA[i]->x - atomsB[j]->x)*(atomsA[i]->x - atomsB[j]->x)+
	  (atomsA[i]->y - atomsB[j]->y)*(atomsA[i]->y - atomsB[j]->y)+
	  (atomsA[i]->z - atomsB[j]->z)*(atomsA[i]->z - atomsB[j]->z);
      }

  /*** do sqrt ***/
  if(doSqrt)
    for(int i=0; i< noAtomsA; i++)
      for(int j=0; j< noAtomsB; j++)
	result[i][j] = sqrt(result[i][j]);

}

void Distances::calculate(Atom * atomA, vector<Atom *> atomsB, bool doSqrt)
{
  //
  // Calculates inter-atomic distances for an atom and a vector of atoms
  //

  vector<Atom*> temp;
  temp.push_back(atomA);

  calculate(temp,atomsB, doSqrt);

}


void Distances::print_result()
{
  int no_B = result.size();
  if(no_B > 0)
    {
      int no_A = result.at(0).size();

      for(int i=0; i< no_A; i++)
	{
	  for(int j=0; j< no_B; j++)
	    {
	      printf("%3.3f",result[i][j]); 
	    }
	  printf("\n");
	}
    }
}


float Distances::calculate(Atom * atomA, Atom * atomB, bool doSqrt)
{
  //
  // Calculates inter-atomic distance for two atoms
  //

  dx = (atomA->x - atomB->x);
  dy = (atomA->y - atomB->y);
  dz = (atomA->z - atomB->z);

  res = 
    dx*dx+
    dy*dy+
    dz*dz;
    
  /*** do sqrt ***/
  if(doSqrt)
    res = sqrt(res);

  return res;

}

vector<vector<float> > Distances::getResult(){return result;}
vector<vector<float> > * Distances::getResultPointer(){return &result;}

void Distances::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: distances.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("Calculates inter- or intra protein atom distances.\n");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
