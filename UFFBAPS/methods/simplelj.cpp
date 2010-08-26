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


#include "simplelj.h"


SimpleLJ::SimpleLJ(FILE * resultsFile, FILE * reporter,vector<SoupObject *> A, vector<SoupObject *> B):Method(reporter)
{

  /*** write out to out file **/
  writeDescription(reporter, A, B);
  /*
  float rmin[6][6] = {{0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
		      {0.00, 4.00, 3.75, 3.60, 4.00, 3.00},
		      {0.00, 3.75, 3.50, 3.35, 3.75, 2.75},
		      {0.00, 3.60, 3.35, 3.20, 3.60, 2.60},
		      {0.00, 4.00, 3.75, 3.60, 4.00, 3.00},
		      {0.00, 3.00, 2.75, 2.60, 3.00, 2.00}};


  float epsilon[6][6] = {{0.00, 0.00,  0.00,  0.00,  0.00,  0.00},
			 {0.00, 0.150, 0.155, 0.173, 0.173, 0.055},
			 {0.00, 0.155, 0.160, 0.179, 0.179, 0.057},
			 {0.00, 0.173, 0.179, 0.200, 0.200, 0.063},
			 {0.00, 0.173, 0.179, 0.200, 0.200, 0.063},
			 {0.00, 0.055, 0.057, 0.063, 0.063, 0.020}};

  */

 
  result =0;
  
  /*  
  printf("conversion[CL] = %d\n", conversion["Cl"]);
  printf("eps[CL] = %f\n", epsilon[conversion["Cl"]][conversion["Cl"]]);

  float eps, rm,distSq,rfrac2,rfrac6,rfrac12;

  for(int noA=0; noA<A.size();noA++)
    for(int i=0; i<A[noA]->numberOfAtoms(); i++)
      for(int noB=0; noB<B.size();noB++)
	for(int j=0; j<B[noB]->numberOfAtoms(); j++)
	  { 
	    try{
	      eps = epsilon[conversion[A[noA].getAtoms()[i]->element]][conversion[B[noB].getAtoms()[j].element]];
	      rm = rmin[conversion[A[noA]->getAtoms()[i]->element]][conversion[B[noB]->getAtoms()[j].element]];
	      
	      if(eps ==0 || rm == 0)
		throw "Unknown type";

	      distSq = 
		(A[noA]->getAtoms()[i]->x - B[noB]->getAtoms()[j].x)*(A[noA]->getAtoms()[i]->x - B[noB]->getAtoms()[j].x)+
		(A[noA]->getAtoms()[i]->y - B[noB]->getAtoms()[j].y)*(A[noA]->getAtoms()[i]->y - B[noB]->getAtoms()[j].y)+
		(A[noA]->getAtoms()[i]->z - B[noB]->getAtoms()[j].z)*(A[noA]->getAtoms()[i]->z - B[noB]->getAtoms()[j].z);
	      
	      rfrac2 = rm*rm/distSq;
	      rfrac6 = rfrac2*rfrac2*rfrac2;
	      rfrac12 = rfrac6*rfrac6;
	      
	      result += eps*(rfrac12 - rfrac6); ////BUG? 2*rfrac6 ?
	    }
	    catch(...)
	      {
		//		printf("found an exception!\n");
		printf("Unexpected atom type '%s' or '%s' \n",
		       A[noA]->getAtoms()[i]->element.c_str(),B[noB]->getAtoms()[j].element.c_str());
	      }
	}
  */
  fprintf(reporter,"The result is: %f kcal/mol\n", result);	
  fprintf(resultsFile,"%-10s %5.3f\n",A[0]->identify().substr(0,4).c_str(),result);	

}
SimpleLJ::~SimpleLJ(){}



void SimpleLJ::writeDescription(FILE * reporter,vector<SoupObject *> A, vector<SoupObject *> B)
{
  version =string("$Id: simplelj.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
    string("Doing simple element based Lennard-Jones poteinal.\n")+
    string("Source of parameters is\n")+
    string("http://www.scripps.edu/mb/olson/gmm/autodock/ad305/Using_AutoDock_305.a.html#29613");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
  
  fprintf(reporter,"Calculating interaction energy between { ");
  for(unsigned int i=0;i<A.size();i++)
    fprintf(reporter,"%s ",A[i]->identify().c_str());
  fprintf(reporter,"} and { ");
  for(unsigned int i=0;i<B.size();i++)
    fprintf(reporter,"%s ",B[i]->identify().c_str());
  fprintf(reporter,"}.\n");



  /*** doing a bit of setup ***/
  conversion[string("C")] = 1;
  conversion[string("N")] = 2;
  conversion[string("O")] = 3;
  conversion[string("S")] = 4;
  conversion[string("H")] = 5;

  

 /*
source:
http://www.scripps.edu/mb/olson/gmm/autodock/ad305/Using_AutoDock_305.a.html#29613

 C     N     O     S     H
C 4.00 3.75 3.60 4.00 3.00
N 3.75 3.50 3.35 3.75 2.75
O 3.60 3.35 3.20 3.60 2.60
S 4.00 3.75 3.60 4.00 3.00
H 3.00 2.75 2.60 3.00 2.00

 C      N      O      S      H
C 0.150 0.155 0.173 0.173 0.055
N 0.155 0.160 0.179 0.179 0.057
O 0.173 0.179 0.200 0.200 0.063
S 0.173 0.179 0.200 0.200 0.063
H 0.055 0.057 0.063 0.063 0.020


Atoms i-j  reqm,ij/ Å   eij/ kcal mol-1
C-C 4.00 0.150
C-N 3.75 0.155
C-O 3.60 0.173
C-S 4.00 0.173
C-H 3.00 0.055
N-C 3.75 0.155
N-N 3.50 0.160
N-O 3.35 0.179
N-S 3.75 0.179
N-H 2.75 0.057
O-C 3.60 0.173
O-N 3.35 0.179
O-O 3.20 0.200
O-S 3.60 0.200
O-H 2.60 0.063
S-C 4.00 0.173
S-N 3.75 0.179
S-O 3.60 0.200
S-S 4.00 0.200
S-H 3.00 0.063 
H-C 3.00 0.055
H-N 2.75 0.057
H-O 2.60 0.063
H-S 3.00 0.063
H-H 2.00 0.020
  */




}


double SimpleLJ::getResult()
{return result;}

