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


#include "peoe.h"


Peoe::Peoe(FILE * reporter,vector<SoupObject *> A):Method(reporter)
{

  /*** write out to out file **/
  version =string("$Id: peoe.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
    string("Doing Partial Equalization of Orbital Electronegativity\n")+
    string("Ref: Gasteiger and Marsili 1980, Tetrahedron V 36 pp 3219-3288");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
  
  fprintf(reporter,"Doing PEOE for ");
  for(unsigned int i=0;i<A.size();i++)
    fprintf(reporter,"%s ",A[i]->identify().c_str());
  fprintf(reporter,"\n");

  /*** setting up parameters ***/
  a[string("H")] =   7.17;
  b[string("H")] =   6.24;
  c[string("H")] =  -0.56;
  d[string("H")] =  12.85;

  a[string("C3")] =  7.98;
  b[string("C3")] =  9.18;
  c[string("C3")] =  1.88;
  d[string("C3")] = 19.04;

  a[string("C2")] =  8.79;
  b[string("C2")] =  9.32;
  c[string("C2")] =  1.51;
  d[string("C2")] = 19.62;

  a[string("C1")] = 10.39;
  b[string("C1")] =  9.45;
  c[string("C1")] =  0.73;
  d[string("C1")] = 20.57;

  a[string("N3")] = 11.54;
  b[string("N3")] = 10.82;
  c[string("N3")] =  1.36;
  d[string("N3")] = 23.72;

  a[string("N2")] = 12.87;
  b[string("N2")] = 11.15;
  c[string("N2")] =  0.85;
  d[string("N2")] = 24.87;

  a[string("N1")] = 15.68;
  b[string("N1")] = 11.70;
  c[string("N1")] = -0.27;
  d[string("N1")] = 27.11;

  a[string("O3")] = 14.18;
  b[string("O3")] = 12.92;
  c[string("O3")] =  1.39;
  d[string("O3")] = 28.49;

  a[string("O2")] = 17.07;
  b[string("O2")] = 13.79;
  c[string("O2")] =  0.47;
  d[string("O2")] = 31.33;

  a[string("F")] =  14.66;
  b[string("F")] =  13.85;
  c[string("F")] =   2.31;
  d[string("F")] =  30.82;

  a[string("Cl")] = 11.00;
  b[string("Cl")] =  9.69;
  c[string("Cl")] =  1.35;
  d[string("Cl")] = 22.04;

  a[string("Br")] = 10.08; 
  b[string("Br")] =  8.47;
  c[string("Br")] =  1.16;
  d[string("Br")] = 19.71;

  a[string("I")] =   9.90;
  b[string("I")] =   7.96;
  c[string("I")] =   0.96; 
  d[string("I")] =  18.82;

  a[string("S3")] = 10.14; 
  b[string("S3")] =  9.13;
  c[string("S3")] =  1.38;
  d[string("S3")] = 20.65;

 

  /*
Parameters from http://msdlocal.ebi.ac.uk/docs/Html/pages/appx_a.htm


; Ref. Tetrahedron, 36, 3219, 1980
;      Croat.Chem.Acta, 53, 601, 1980

#Title Gasteiger-Marsili charges

; Type        a        b        c        d
; =============================================
  H          7.17     6.24    -0.56    12.85
  C3         7.98     9.18     1.88    19.04
  C2         8.79     9.32     1.51    19.62
  C1        10.39     9.45     0.73    20.57
  N3        11.54    10.82     1.36    23.72
  N2        12.87    11.15     0.85    24.87
  N1        15.68    11.70    -0.27    27.11
  O3        14.18    12.92     1.39    28.49
  O2        17.07    13.79     0.47    31.33
  F         14.66    13.85     2.31    30.82
  Cl        11.00     9.69     1.35    22.04
  Br        10.08     8.47     1.16    19.71
  I          9.90     7.96     0.96    18.82
  S3        10.14     9.13     1.38    20.65
  */



  for(unsigned int m=0; m<A.size(); m++)
    {
      int noAtoms = A[m]->numberOfAtoms();

      /*** finding number of bonds ***/
      BondMatrix BM(reporter, A[m]);
      vector<vector<int> > bonds = BM.getResult();

      vector<int> noBonds;
      noBonds.assign(noAtoms,0);
      
      for(int i=0;i<noAtoms;i++)
	for(int j=0;j<noAtoms;j++)
	  if(bonds[i][j]>0)
	    noBonds[i] += 1;

      /*** determining hydrdization ***/
      vector<string> hybrid;
      hybrid.assign(noAtoms,"dummy");

      for(int i=0;i<noAtoms;i++)
	{
	  string elem = A[m]->getAtoms()[i]->element;
	  
	  if(elem == "H")
	    {
	      if(noBonds[i] == 0)
		printf("Waring: Hydrogen without bonds found! Atom no. %d\n",i+1);
	      else if(noBonds[i] == 1)
		hybrid[i] = "H";
	      else
		printf("Waring: Hydrogen with more than one bond found! Atom no. %d\n",i+1);
	    }
	  else if(elem == "C")
	    {
	      if(noBonds[i] == 0)
		printf("Waring: Carbon without bonds found! Atom no. %d\n",i+1);
	      else if(noBonds[i] == 1)
		printf("Waring: Carbon with one bond found! Atom no. %d\n",i+1);
	      else if(noBonds[i] == 2)
		hybrid[i] = "C1";
	      else if(noBonds[i] == 3)
		hybrid[i] = "C2";
	      else if(noBonds[i] == 4)
		hybrid[i] = "C3";
	      else
		printf("Waring: carbon with more than four bonds found! Atom no. %d\n",i+1);

	      
	    }
	  else if(elem == "O")
	    {
	      if(noBonds[i] == 0)
		printf("Waring: Oxygen without bonds found! Atom no. %d\n",i+1);
	      else if(noBonds[i] == 1)
		hybrid[i] = "O2";
	      else if(noBonds[i] == 2)
		hybrid[i] = "O3";
	      else
		printf("Waring: Oxygen with more than two bonds found! Atom no. %d\n",i+1);
	    }
	  else if(elem == "N")
	    {
	      if(noBonds[i] == 0)
		printf("Waring: Nitrogen without bonds found! Atom no. %d\n",i+1);
	      else if(noBonds[i] == 1)
		printf("Waring: Nitrogen with one bond found! Atom no. %d\n",i+1);
	      else if(noBonds[i] == 2)
		hybrid[i] = "N1";
	      else if(noBonds[i] == 3)
		hybrid[i] = "N2";
	      else if(noBonds[i] == 4)
		hybrid[i] = "N3";
	      else
		printf("Waring: Nitrogen with more than four bonds found! Atom no. %d\n",i+1);
      	    }
	  else if(elem == "S")
	    {
	      if(noBonds[i] == 0)
		printf("Waring: Sulphur without bonds found! Atom no. %d\n",i+1);
	      else if(noBonds[i] == 1)
		printf("Waring: Sulphur with one bond found! Atom no. %d\n",i+1);
	      else if(noBonds[i] == 2)
		hybrid[i] = "S3";
	      else
		printf("Waring: Sulphur with more than two bonds found! Atom no. %d\n",i+1);
	    }
	  else if(elem == "F")
	     {
	      if(noBonds[i] == 0)
		printf("Waring: Flourine without bonds found! Atom no. %d\n",i+1);
	      else if(noBonds[i] == 1)
		hybrid[i] = "F";
	      else
		printf("Waring: flourine with more than one bond found! Atom no. %d\n",i+1);
	    }

	  else if(elem == "Cl")
	     {
	      if(noBonds[i] == 0)
		printf("Waring: Chlorine without bonds found! Atom no. %d\n",i+1);
	      else if(noBonds[i] == 1)
		hybrid[i] = "Cl";
	      else
		printf("Waring: Chlorine with more than one bond found! Atom no. %d\n",i+1);
	    }

	  else if(elem == "Br")
	     {
	      if(noBonds[i] == 0)
		printf("Waring: Bromine without bonds found! Atom no. %d\n",i+1);
	      else if(noBonds[i] == 1)
		hybrid[i] = "Br";
	      else
		printf("Waring: Bromine with more than one bond found! Atom no. %d\n",i+1);
	     }

	  else if(elem == "I")
	     {
	      if(noBonds[i] == 0)
		printf("Waring: Iodine without bonds found! Atom no. %d\n",i+1);
	      else if(noBonds[i] == 1)
		hybrid[i] = "I";
	      else
		printf("Waring: Iodine with more than one bond found! Atom no. %d\n",i+1);
	     }

	}

      printf("The following atoms have 'unsual' bonding:\n");
      for(int i=0;i<noAtoms;i++)
	{
	  string elem = A[m]->getAtoms()[i]->element;
	  int bondds = 0;
	  for(int j=0;j<noAtoms;j++)
	    bondds += bonds[i][j];
	  
	  bool odd = false;
	  if(elem == "C" && bondds != 4)
	    odd =true;
	  if(elem == "N" && bondds != 3)
	    odd =true;
	  if(elem == "O" && bondds != 2)
	    odd =true;
	  if(elem == "S" && bondds != 2)
	    odd =true;
	  if(odd)
	    printf("Atom %d %s %s in res %s has %d bonds and hydrid %s\n",
		   i+1,
		   A[m]->getAtoms()[i]->element.c_str(),
		   A[m]->getAtoms()[i]->name.c_str(),
		   A[m]->getAtoms()[i]->residue.c_str(),
		   noBonds[i],
		   hybrid[i].c_str());
	  

	}
      
      /*      for(int i=0;i<noAtoms;i++)
	printf("Atom %d %s has %d bonds and hydrid %s\n",
	       i+1,
	       A[m]->atoms[i].element.c_str(),noBonds[i],
	       hybrid[i].c_str());

      */

      vector<double> electronegativities,electronegativitiesPlus1,charges,tempCharge;
      electronegativitiesPlus1.assign(noAtoms,0);
      electronegativities.assign(noAtoms,0);
      charges.assign(noAtoms,0);
      tempCharge.assign(noAtoms,0);

      /*** setting initial charges ***/
      for(int i=0;i<noAtoms;i++)
	{
	  //how to generelly set formal charges???
	}


      /*** setting electronegativities ***/
      for(int i=0;i<noAtoms;i++)
	{
	  electronegativitiesPlus1[i] = 
	    a[hybrid[i]] + 
	    b[hybrid[i]] + 
	    c[hybrid[i]];
	  electronegativities[i] = 
	    a[hybrid[i]] + 
	    b[hybrid[i]] * charges[i] + 
	    c[hybrid[i]] * charges[i] * charges[i];
	}


      for(int i=0;i<noAtoms;i++)
	printf("Atom %d %s has %d bonds and hydrid %s, electroneg %f, electroneg+1 %f, charges %f\n",
	       i+1,
	       A[m]->getAtoms()[i]->element.c_str(),noBonds[i],
	       hybrid[i].c_str(),
	       electronegativities[i],
	       electronegativitiesPlus1[i],
	       charges[i]);
      
      vector<vector<double> >  deltaCharge;

      /*** doing iterations ***/
      for(int iter=1;iter<=6;iter++)
       {
	 for(int i=0;i<noAtoms;i++)
	   for(int j=0;j<noAtoms;j++)
	     {
	       
	       if(bonds[i][j] == 1)
		 {
		   /*** neighbour more electroneg. => neg. charge trans. to neighbour ***/
		   if(electronegativities[i] < electronegativities[j]) 
		     {                         
		       tempCharge[i] += 
			 (electronegativities[j] - electronegativities[i])/
			 electronegativitiesPlus1[i] * pow(0.5,iter);
		     }
		   /*** neighbour less electroneg. => neg. charge trans. from neighbour ***/
		   else                                                
		     {
		       tempCharge[i] += 
			 (electronegativities[j] - electronegativities[i])/
			 electronegativitiesPlus1[j] * pow(0.5,iter);
		     
		     }
		 }

	     }

	 deltaCharge.push_back(tempCharge);
	 tempCharge.assign(noAtoms,0);

	 /*** update charges ***/
	 charges.assign(noAtoms,0);
	 for(int i=0;i<noAtoms;i++)
	   for(int j=0;j<iter;j++)
	     charges[i] += deltaCharge[j][i];
	 
	 /*** update electronegativities ***/
	 for(int i=0;i<noAtoms;i++)
	   electronegativities[i] = 
	     a[hybrid[i]] + 
	     b[hybrid[i]] * charges[i] + 
	     c[hybrid[i]] * charges[i] * charges[i];
	   

      for(int i=0;i<noAtoms;i++)
	printf("At iteration %d atom %d %s has %d bonds and hydrid %s, electroneg %f, electroneg+1 %f, charges %f\n",
	       iter,
	       i+1,
	       A[m]->getAtoms()[i]->element.c_str(),
	       noBonds[i],
	       hybrid[i].c_str(),
	       electronegativities[i],
	       electronegativitiesPlus1[i],
	       charges[i]);
       }


      double totalCharge =0;
      for(int i=0;i<noAtoms;i++)
	totalCharge+=charges[i];

      printf("Total charge is %f \n",totalCharge);

    } // loop over soupObjects vector


}
Peoe::~Peoe(){}



