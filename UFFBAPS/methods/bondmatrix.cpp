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
$Id: bondmatrix.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "bondmatrix.h"


BondMatrix::BondMatrix(FILE * reporter, SoupObject * A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  
  calculate(A);
}

BondMatrix::BondMatrix(FILE * reporter, vector<Atom*> A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  
  calculate(A);
}


BondMatrix::BondMatrix(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
}

BondMatrix::BondMatrix():Method(){}


BondMatrix::~BondMatrix(){}


void BondMatrix::calculate(SoupObject * A)
{
  calculate(A->getAtoms());
}

void BondMatrix::calculate(vector<Atom*> atoms) {
  /*** setting things up ***/
  result.clear();
  vector<int> temp;
  temp.assign(atoms.size(), 0);
  result.assign(atoms.size(), temp);
  
  double SSdist = 2.5;
  double Hdist = 1.5;
  double defaultDist = 2;
  double SSdistSquare = SSdist*SSdist;
  double HdistSquare = Hdist*Hdist;
  double defaultDistSquare=defaultDist*defaultDist;
  printf ("number of atoms: %d\n",atoms.size());
  //
  // Divide into boxes
  //
  Boxes BOX(atoms,static_cast<float>(SSdist)); // Set up the Box definition

  int no_atoms = atoms.size();
  vector<int> close_atoms;
  /*** finding bonds ***/
  for(int i=0; i<no_atoms; i++) {
      string elem1 =atoms[i]->element;
      close_atoms=BOX.get_close_atoms(atoms[i]); // Jens
      for (unsigned int jcount=0;jcount<close_atoms.size();jcount++) { //Jens
	int j=close_atoms[jcount]; // Jens
	/*** for every atom i ***/
	//for(int j=0; j<no_atoms; j++) { // Chresten
	/*** find bonded atoms to i ***/
	double dist=static_cast<double>(d.calculate(atoms[i],atoms[j],false)); // Jens
	//if(d.calculate(atoms[i],atoms[j],false) < SSdistSquare) //As SS dist is the greatest distance considered, it is used to get rid of all pairs not binding
	if (dist<SSdistSquare) {
	  if(i != j) {
	      string elem2 =atoms[j]->element;		  
	      /*** if di-sulphide ***/
	      if(elem1 =="S" && elem2 == "S") {
		//if(d.calculate(atoms[i],atoms[j],false) < SSdistSquare) // Chresten
		if (dist<SSdistSquare) { // Jens
		  result[i][j] = result[j][i] = 1;
		}
	      }
	      /*** if one hydrogen ***/
	      else if(elem1 =="H" || elem2 == "H") {
		if(!(elem1 =="H" && elem2 == "H"))
		  //if(d.calculate(atoms[i],atoms[j],false) < HdistSquare) // Chresten
		  if (dist<HdistSquare) { // Jens
		    result[i][j] = result[j][i] = 1;
		  }
	      }	
		/*** all other atom pairs ***/
	      else {
		//if(d.calculate(atoms[i],atoms[j],false) < defaultDistSquare)
		if (dist<defaultDistSquare) { // Jens
		  result[i][j] = result[j][i] = 1;
		}
	      } 
	  }
	} 
      }
  }
  printf ("This function is executed\n");
  /*** make sure hydrogens only have one bond ***/
  Boxes BOX2(atoms,Hdist);
  unsigned int minLocation =0;
  for(unsigned int i=0; i<atoms.size(); i++) {
    string elem1 =atoms[i]->element;
    double minDist = 0.0;
    if(elem1 == "H") {
	/*** find closest bond ***/ 
	//for(unsigned int j=0; j<atoms.size(); j++)
	close_atoms=BOX2.get_close_atoms(atoms[i]);
	for (unsigned int jcount=0;jcount<close_atoms.size();jcount++) {
	  unsigned int j=close_atoms[jcount];
	  if(result[i][j] == 1) {
	    if(minDist == 0.0) {
	      minDist = d.calculate(atoms[i],atoms[j],false);
	      minLocation = j;
	    }
	    else if(minDist > d.calculate(atoms[i],atoms[j],false)) {
	      minDist = d.calculate(atoms[i],atoms[j],false);
	      minLocation = j;
	    }
	  }
	  /*** and delete the rest ***/
	  for(unsigned int j=0; j<atoms.size(); j++)
	    if(j==minLocation) {
	      result[i][j] = result[j][i] = 1;
	    }
	    else {
	      result[i][j] = result[j][i] = 0;
	    }
	
      }
    }
  }    
  /*** set up double bonds ***/

  for(unsigned int i=0; i<atoms.size(); i++)
    {
      if(noBondsCompExp(i,atoms) < 0)
	{
	  int doubleBondAtom = -1;
	  for(unsigned int j=0; j<atoms.size(); j++)
	    if(result[i][j]==1 && noBondsCompExp(j,atoms) < 0)
	      doubleBondAtom = j;

	  if(doubleBondAtom != -1)
	    {
	      /*
	      printf("Found double bond bt %d and %d\n", i,doubleBondAtom );
	      atoms[i]->display();
	      atoms[doubleBondAtom]->display();
	      printf("Exps are %d and %d\n",noBondsCompExp(i,atoms),noBondsCompExp(doubleBondAtom,atoms));
	      */
	      result[i][doubleBondAtom] = result[doubleBondAtom][i] = 2;		
	    }
	}
    }




  /*  for(unsigned int i=0; i<atoms.size(); i++)
    {
      printf("Atom %d %s is bonded to:",i+1,atoms[i]->element.c_str());
      for(unsigned int j=0; j<atoms.size(); j++)
	{
	  if(result[i][j] ==1 && result[j][i] == 1)
	    printf(" %d %s ,",j+1,atoms[j]->element.c_str());
	  if(result[i][j] ==2 && result[j][i] == 2)
	    printf(" =%d %s ,",j+1,atoms[j]->element.c_str());


	}
      printf("\n");
    }

  */
}

void BondMatrix::assign_bonds(Soup * A)
{
  calculate_list(A->getAtoms());
  store_bonding_information_in_atoms(A->getAtoms());
  A->bonds_assigned = true;

  return;
}


void BondMatrix::calculate_list(vector<Atom*> atoms) {
  /*** setting things up ***/
  result.clear();
  //
  // Provide an empty list of bonds for each atom
  //
  _bondmatrix.clear();
  for (unsigned int j=0;j<atoms.size();j++) {
    vector<int> dummy;
    _bondmatrix.push_back(dummy);
  }

   
  double SSdist = 2.5;
  double Hdist = 1.5;
  double defaultDist = 2;
  double SSdistSquare = SSdist*SSdist;
  double HdistSquare = Hdist*Hdist;
  double defaultDistSquare=defaultDist*defaultDist;  

  int no_atoms = atoms.size();
  Boxes BOX(atoms,static_cast<float>(SSdist));
  vector<int> close_atoms;

  /*** finding bonds ***/
  for(int i=0; i<no_atoms; i++) {
    string elem1 =atoms[i]->element;
    close_atoms=BOX.get_close_atoms(atoms[i]);
    for (unsigned int jcount=0;jcount<close_atoms.size();jcount++) {
      int j=close_atoms[jcount];
      if (j<i) {
	  /*** for every atom i ***/
	  //for(int j=0; j<no_atoms; j++) {
	  /*** find bonded atoms to i ***/
	  double dist=static_cast<double>(d.calculate(atoms[i],atoms[j],false));
	if(dist < SSdistSquare) { //As SS dist is the greatest distance considered, it is used to get rid of all pairs not binding
	  if(i != j) {
	    string elem2 =atoms[j]->element;		  
	    /*** if di-sulphide ***/
	    if(elem1 =="S" && elem2 == "S") {
	      if(dist < SSdistSquare) {
		make_list_entry(i,j,1);
	      }
	    }
	    /*** if one hydrogen ***/
	    else if(elem1 =="H" || elem2 == "H") {
	      if(!(elem1 =="H" && elem2 == "H")) {
		if(dist < HdistSquare) {
		  make_list_entry(i,j,1);
		}
	      }
	    }	
	    /*** all other atom pairs ***/
	    else {
	      if(dist < defaultDistSquare) {
		make_list_entry(i,j,1);
	      }
	    } 
	  }
	}
      }
    }
  }

  /*** make sure hydrogens only have one bond ***/
  for(unsigned int i=0; i<atoms.size(); i++) {
    string elem1 =atoms[i]->element;
    if(elem1 == "H") {
      if (_bondmatrix[i].size()>1) {
	//printf ("H has more than one bond: %d\n",i);
	/*** find closest bond ***/ 
	//for(unsigned int j=0; j<atoms.size(); j++) {
	unsigned int minLocation =0;
	double minDist = 0.0;
	for (unsigned int jcount=0;jcount<_bondmatrix[i].size();jcount++) {
	  unsigned int j=_bondmatrix[i][jcount];
	  //if(is_in_list(i,j)) {
	  if(minDist == 0.0) {
	    minDist = d.calculate(atoms[i],atoms[j],false);
	    minLocation = j;
	  }
	  else if(minDist > d.calculate(atoms[i],atoms[j],false)) {
	    minDist = d.calculate(atoms[i],atoms[j],false);
	    minLocation = j;
	  }
	}
	only_bond(i,minLocation,1);
      }
    }
  }
  return;
}


void BondMatrix::make_list_entry(int i,int j,int type) {
  vector<int> temp;
  temp.assign(3, 0);
 
  if(!is_in_list(i,j))
    {
      temp[0] = i;
      temp[1] = j;
      temp[2] = type;
      
      result.push_back(temp);
      // Update the bond matrix
      _bondmatrix[i].push_back(j);
      //_bondmatrix[j].push_back(i);
    }
  if(!is_in_list(j,i))
    {
      temp[0] = j;
      temp[1] = i;
      temp[2] = type;
    
      result.push_back(temp);
      //_bondmatrix[i].push_back(j);
      _bondmatrix[j].push_back(i);
    }


}


void BondMatrix::only_bond(int i,int j,int type)
{

  // delete all bonds to i
  vector<vector<int> >::iterator tempIterator;
  
  for( tempIterator = result.begin(); tempIterator != result.end(); ) 
    {
      if( tempIterator->at(0)==i  || tempIterator->at(1)==i )
	tempIterator = result.erase(tempIterator);      
      else
	tempIterator++;
    }
  // Clear the bond matrix for atom i
  _bondmatrix[i].clear();
  // Remove i in bondmatrix for j
  vector<int> newj;
  for (unsigned int jj=0;jj<_bondmatrix[j].size();jj++) {
    if (_bondmatrix[j][jj]!=i) {
      newj.push_back(_bondmatrix[j][jj]);
    }
  }
  _bondmatrix[j].clear();
  _bondmatrix[j]=newj;
  //Remake the good bond
  make_list_entry(i,j,type);

  return;
}


int BondMatrix::no_bonds(int i) {
  int res = 0;
  for(unsigned int r=0; r<result.size(); r++)
    if(result[r][0]==i) {
      res += result[r][2];
    }
  return res;
}


bool BondMatrix::is_in_list(int i,int j) {
  bool res = false;
  ///printf ("Iterating over %d values\n",result.size());
  //for(unsigned int r=0; r<result.size(); r++) {
  //  if(result[r][0]==i && result[r][1]==j) {
  //    res = true;
  //    break;
  //  }
  //}
  //bool res2=false;
  
  //printf ("Iterating over %d records\n",_bondmatrix[i].size());
  for (unsigned int jbond=0;jbond<_bondmatrix[i].size();jbond++) {
   if (_bondmatrix[i][jbond]==j) {
      res=true;
      break;
    }
  }
  //if (res!=res2) {
  //  printf ("bools differ: %d %d, i: %d, j: %d\n",res,res2,i,j);
  //  for (unsigned int jbond=0;jbond<_bondmatrix[i].size();jbond++) {
  //    printf ("%d ",_bondmatrix[i][jbond]);
  //  }
  //  printf ("\n");
  // }
  return res;
}


void BondMatrix::change_type(int i,int j,int newtype)
{

  for(unsigned int r=0; r<result.size(); r++)
    {
      if(result[r][0]==i && result[r][1]==j)
	result[r][2]=newtype;
      if(result[r][0]==j && result[r][1]==i)
	result[r][2]=newtype;
      
    }
}


vector<vector<int> > BondMatrix::getResult()
{
  // for(unsigned int r=0; r<result.size(); r++)
  //    printf("Bond: %4d-%4d\n",result[r][0],result[r][1]);
  
  return result;
}

int BondMatrix::no_bonds_matrix(int i)
{

  int res = 0;
  for(unsigned int r=0; r<result.size(); r++)
    if(result[r][i]>0)
      {
	res += result[r][i];
	//printf(" -> %d and %d are bonded\n",i,r);
      }
  //printf("Found %d bonds for %d\n",res,i);

  return res;
}

int BondMatrix::noBondsCompExp(int a, vector<Atom*> atoms)
{

  string elem = atoms[a]->element;
  
  int bonds = no_bonds_matrix(a);
  int res = 0;

  if(elem == "C")
    res = bonds - 4;
 
  if(elem == "N")
    res = bonds - 3;

  if(elem == "O")
    res = bonds - 2;

  if(elem == "S")
    res = bonds - 2;
 
  return res;
}

bool BondMatrix::connection(vector<vector<int> > * bonds, int i, int j, int no_bonds, int cur_no_bonds, int last)
{
  /*
  printf("-->Connection started with i:%d,j:%d,no_bonds:%d,cur_no_bonds:%d, last:%d\n",
	 i,j,no_bonds,cur_no_bonds, last);
  */
  
  // finds out if atoms with indexes i and j are connected with no_bonds bonds or less
  bool res = false;

  // find partners and remove the previous atom
  vector<int> partners = find_bonding_atoms(bonds, i);
  remove(partners.begin(), partners.end(), last);

  /*  
  printf("    Partners are:");
  for(unsigned int p=0; p<partners.size();p++)
    printf("%d,",partners[p]);
  printf("\n");
  */

  if (count(partners.begin(),partners.end(),j)>0)
    res = true;
  else if(cur_no_bonds < no_bonds) 
    for(unsigned int p=0; p<partners.size();p++)
      if(connection(bonds,partners[p],j,no_bonds,cur_no_bonds+1,i))
	res=true;

  return res;
}


vector<int> BondMatrix::find_bonding_atoms(vector<vector<int> > * bonds, int a)
{
  // returns a list of indexes for atoms bonding to the atom with index a
  
  vector<int> res;
 
  // printf("Found these bonding atoms for atom %d:\n",a);
  for(unsigned int i=0; i<bonds->size(); i++)
    if(bonds->at(i).at(0)==a)
      {
	//printf("  %d\n",bonds->at(i).at(1));
	res.push_back(bonds->at(i)[1]);
      }
  
  return res;
}


void BondMatrix::store_bonding_information_in_atoms(vector<Atom*> atoms)
{
  //
  // Stores bonding atoms in atoms
  // Can only be called after calculate_list
  //
  

  vector<Atom*> bonding_atoms;

  vector<vector<int> >::iterator bond;

  for (unsigned int a=0; a<atoms.size(); a++)
    {
      bonding_atoms.clear();
      
      bond = result.begin();
      while(bond != result.end())
	{
	  if (bond->at(0) == (int) a)
	    bonding_atoms.push_back(atoms[bond->at(1)]);
      
	  bond++;
	}

      atoms[a]->bonded_atoms = bonding_atoms;
      //      cout<<atoms[a]->bonded_atoms.size()<<" are bonded to ";
      //      atoms[a]->display();
    }

  /*
  for (int a=0; a<atoms.size(); a++)
    {
      cout<<"Bonding atoms for ";
      atoms[a]->display();

      for(unsigned int b=0; b<atoms[a]->bonded_atoms.size(); b++)
	{
	  cout<<"    ";
	  atoms[a]->bonded_atoms[b]->display();

	}
    }
  */
  return;
}



void BondMatrix::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  rep = reporter;
  version =string("$Id: bondmatrix.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("Calculates bondmatrix");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
