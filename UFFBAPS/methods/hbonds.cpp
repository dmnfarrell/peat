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
$Id: hbonds.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "hbonds.h"


Hbonds::Hbonds():Method(){  set_parameters();} //Dummy constructor for use in ElementsNN2


Hbonds::Hbonds(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  set_parameters();

}

Hbonds::Hbonds(FILE * resultsFile,FILE * reporter, vector<Soup*> A, vector<Soup*> B):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  set_parameters();


  if(A.size() != B.size())
    {
      cout << "WARNING: An equal amount of protein and ligand structures must be given to the Hbonds function!\n";
      return;
    }

 
  float result = 0;


  for(unsigned int i=0; i<A.size(); i++)
    {
      //dist.calculate(A[i], B[i], true);
      // distances = dist.getResult();
      
      result = calculate(A[i],B[i]);
      
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

Hbonds::Hbonds(FILE * resultsFile,FILE * reporter, vector<Soup*> A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  set_parameters();

  float result = 0;

  for(unsigned int i=0; i<A.size(); i++)
    {
      result = calculate(A[i]);
      
      fprintf(resultsFile,"%-4s %8f\n",
 	      A[i]->name.substr(0,4).c_str(),
	      result);
    }  

}


Hbonds::~Hbonds(){}




float Hbonds::calculate(Soup * A, Soup * B)
{
  /*** bonds is used for bookkeeping of hydrogen bonds ***/
  bonds.clear();
  donor_count.clear();
  acceptor_count.clear();
  no_bounded_hydrogens.clear();
  float result = 0;

  atomsA = A->get_dry_atoms();  
  atomsB = B->get_dry_atoms();  

  BondMatrix bm;
  bm.calculate_list(atomsA);
  bondsA = bm.getResult();
  bm.calculate_list(atomsB);
  bondsB = bm.getResult();


  // make list of how many hydrogens are bonded to each heavy atom
  make_list_of_number_of_bonded_hydrogens();


  float factor,d;
  for(unsigned int ai=0; ai< atomsA.size(); ai++)
    for(unsigned int bi=0; bi< atomsB.size(); bi++) 
      {
	d = dist.calculate(atomsA[ai],atomsB[bi],false);
	if (is_hbond(atomsA[ai], atomsB[bi], d)||is_hbond(atomsB[bi], atomsA[ai], d)) { // This is wrong - I have to rewrite from scratch. Jens. 
	    d=sqrt(d);
	    factor = 1.0;
	    /*** factor - distance to water accessible surface ***/
	    if(include_surf_factor == 1)
	      factor = factor*(atomsA[ai]->desolvation_factor + atomsB[bi]->desolvation_factor)/2; //using average desolvation_factor factor of the two atoms
	    
	    /*** factor - distance between donor and acceptor ***/
	    if(include_dist_factor == 1)
	      factor = factor*distance_factor(d);
	    /*** factor - root-donor and root-acceptor distances ***/
	    if(include_angle_factor == 1)
	      factor = factor*angle_factors(ai,bi);
	    
	    if(atomsA[ai]->hydrogen_bond_could_be_replaced_by_water or 
	       atomsB[bi]->hydrogen_bond_could_be_replaced_by_water) {
	      factor = factor*replacement_by_water_factor;
	      //cout<<"Can be replaced "<<replacement_by_water_factor<<endl;
	    }
	    //else
	    //  cout<<"Can NOT be replaced"<<endl;

	    result = result + factor;
	    
	    /*
	      cout<<"Atom: ";
	      atomsA[ai]->display();
	      cout<<" --> ";
	      atomsB[bi]->display();
	      cout<<": "<<factor<<endl;
	    */
	  }
      }

  /**********************************************************************/
//   printf("TESTING SURFACE:\n");
//   float d = 0;
//   while(d<12)
//     {
//       printf("%6.2f  %6.2f\n",d,bf.surface_dist_factor(d,d+1));
//       d =d+0.1;
//     }

//   printf("TESTING DISTANCES:\n");
//   d = 0;
//   while(d<5)
//     {
//       printf("%6.2f  %6.2f\n",d,distance_factor(d));
//       d =d+0.1;
//     }
//   printf("TESTING ANGLES:\n");
//   d = 0;
//   while(d<370)
//     {
//       printf("%6.2f  %6.2f\n",d,angle_factor(d));
//       d =d+1;
//     }


  /*************************************************************************/

  //each H-bond contributes 1.4 kcal/mol  
  //ref: Guerois,Nielsen and Serrano: Predictig Changes in the stability of Proteins and Protein Complexes: 
  //                                  A study of more than 1000 mutations.
  result = -1.4 * 4.1868 * result; // kJ/mol
  return result;
}

float Hbonds::calculate(Soup * A)
{
  /*** bonds is used for bookkeeping of hydrogen bonds***/
  bonds.clear();
  donor_count.clear();
  acceptor_count.clear();
  no_bounded_hydrogens.clear();

  float result = 0;
  atomsA = A->get_dry_atoms();  
  atomsB.clear();
  BondMatrix bm;
  bm.calculate_list(atomsA);
  bondsA = bm.getResult();

  // make list of how many hydrogens are bonded to each heavy atom
  make_list_of_number_of_bonded_hydrogens();
  // Divide into boxes
  Boxes BOX(A,sqrt(max_dist));
  float d;
  unsigned int bi;
  vector<int> close_atoms;
  for(unsigned int ai=0; ai< atomsA.size(); ai++) {
    //for(unsigned int bi=0; bi< ai; bi++) {//Chresten org code
    if (atomsA[ai]->element!="H") {
      close_atoms=BOX.get_close_atoms(atomsA[ai]); //Jens
      for (unsigned int jcount=0;jcount<close_atoms.size();jcount++) { // Jens
	bi=close_atoms[jcount]; // Jens
	if (atomsA[bi]->element!="H") {
	  if (bi<ai) { //jens
	    d = dist.calculate(atomsA[ai],atomsA[bi],false);
	    result=result+hbond_factor(atomsA[ai],ai, atomsA[bi],bi, d) + hbond_factor(atomsA[bi],bi, atomsA[ai],ai, d);
	  } //Jens
	}
      }
    }
  }
  //printf ("Number of Hbonds: %f\n",last_result);
  //each H-bond contributes 1.4 kcal/mol  
  //ref: Guerois,Nielsen and Serrano: Predictig Changes in the stability of Proteins and Protein Complexes: 
  //                                  A study of more than 1000 mutations.

  /**********************************************************************/
//   printf("TESTING SURFACE:\n");
//   float d = 0;
//   while(d<12)
//     {
//       printf("%6.2f  %6.2f\n",d,bf.surface_dist_factor(d,d+1));
//       d =d+0.1;
//     }

//   printf("TESTING DISTANCES:\n");
//   d = 0;
//   while(d<5)
//     {
//       printf("%6.2f  %6.2f\n",d,distance_factor(d));
//       d =d+0.1;
//     }
//   printf("TESTING ANGLES:\n");
//   d = 0;
//   while(d<370)
//     {
//       printf("%6.2f  %6.2f\n",d,angle_factor(d));
//       d =d+1;
//     }


  /*************************************************************************/
  //  cout<<"-------------------"<<endl;
  //  cout<<"Total number of HBs: "<<bonds.size()<<endl;

  result = -1.4 * 4.1868 * result; // kJ/mol

  return result;
}

float Hbonds::hbond_factor(Atom * acceptor,int ai, Atom * donor,int bi, double dist) {
  //bool res = false;
  float factor=0.0;
  //check acceptor element
  if(count(possible_acceptors.begin(), possible_acceptors.end(),acceptor->element) > 0 && acceptor->name!="N") { // Jens excluded mainchain Ns as acceptors
    //check donor element
    if(count(possible_donors.begin(), possible_donors.end(),donor->element) > 0) {
      //check distance
      if(dist > min_dist && dist < max_dist) {
	//check number of hbonds 
	if(number_of_hbonds_ok(acceptor,donor)) {
	  // 
	  // The factor calculation moved here by Jens
	  // It does not make sense to have it anywhere else
	  //
	  dist=sqrt(dist);
	  factor=1;
	  if(include_surf_factor == 1) {
	    factor = factor*(atomsA[ai]->desolvation_factor + atomsA[bi]->desolvation_factor)/2; //using average desolvation_factor factor of the two atoms
	  }	    
	  /*** factor - distance between donor and acceptor ***/
	  if(include_dist_factor == 1) {
	    factor = factor*distance_factor(dist);
	  }
	  /*** factor - root-donor and root-acceptor distances ***/
	  if(include_angle_factor == 1) {
	    factor = factor*angle_factors_internal(ai,bi);
	  }
	  if (factor>0.01) {
	    //res = true;
	    register_bond(acceptor, donor); 
	  }
	}
      }
    }
  }
  return factor;
}

bool Hbonds::is_hbond(Atom * acceptor, Atom * donor, double dist) {
  bool res = false;
  //check acceptor element
  if(count(possible_acceptors.begin(), possible_acceptors.end(),acceptor->element) > 0 && acceptor->name!="N") { // Jens excluded mainchain Ns as acceptors
    //check donor element
    if(count(possible_donors.begin(), possible_donors.end(),donor->element) > 0) {    
      //check distance
      if(dist > min_dist && dist < max_dist) {
	//check number of hbonds 
	if(number_of_hbonds_ok(acceptor,donor)) {
	  res = true;
	  register_bond(acceptor, donor); 
	}
      }
    }
  }
  return res;
}

bool Hbonds::water_bridge_hbond(Atom * acceptor, Atom * donor, double dist) //used by the water bridge method
{
  bool res = false;
  
  //check acceptor element
  if(count(possible_acceptors.begin(), possible_acceptors.end(),acceptor->element) > 0)
    //check donor element
    if(count(possible_donors.begin(), possible_donors.end(),donor->element) > 0)
      //check distance
      if(dist > min_dist && dist < max_dist)
	res = true;
	
  return res;
}



bool Hbonds::is_potential_hbond(Atom * acceptor, Atom * donor)
{
  bool res = false;

  //check acceptor element
  if(count(possible_acceptors.begin(), possible_acceptors.end(),acceptor->element) > 0)
    //check donor element
    if(count(possible_donors.begin(), possible_donors.end(),donor->element) > 0)
      res = true;

  return res;
  
}


bool Hbonds::is_potential_acceptor(Atom * acceptor)
{
  bool res = false;

  //check acceptor element
  if(count(possible_acceptors.begin(), possible_acceptors.end(),acceptor->element) > 0)
    res = true;

  return res;
  
}

bool Hbonds::is_potential_donor(Atom * donor)
{
  bool res = false;

  //check donor element
  if(count(possible_donors.begin(), possible_donors.end(),donor->element) > 0)
    res = true;
  
  return res;
  
}


void Hbonds::register_bond(Atom * acceptor, Atom * donor)
{
  bonds[donor]=acceptor;

  if(donor_count.count(donor)==0)
    donor_count[donor]=1;
  else
    donor_count[donor]++;

  if(acceptor_count.count(acceptor)==0)
    acceptor_count[acceptor]=1;
  else
    acceptor_count[acceptor]++;

  /*** store hydrogen bond in atom objects ***/
  acceptor->hydrogen_bonds.push_back(donor);
  donor->hydrogen_bonds.push_back(acceptor);

  return;
}


bool Hbonds::number_of_hbonds_ok(Atom * acceptor, Atom * donor)
{
  /*** 
       Checks that the number of hydrogen bonds to an atom is reasonable

       Acceptor: number of hydrogen bonds must not exceed number of lone pairs
       Donor: number of hydrogen bonds must not exceed number of bound hydrogen atoms

  ***/
  bool res = true;

  if(donor_count.count(donor)>0)
    {
      if(donor_count[donor] == number_bounded_hydrogens(donor))
	{  
	  /*printf("++ Donor atom already has %d hydrogen bonds and has %d hydrogen atoms attached\n    ",
	    donor_count[donor],number_bounded_hydrogens(donor));
	  donor->display();*/
	  res = false;
	}
    }
  else if(number_bounded_hydrogens(donor)==0)
    res = false;

  if(acceptor_count.count(acceptor)>0)
    if(acceptor_count[acceptor] == no_lone_pairs[acceptor->element]) {
	/*	printf("++ Acceptor atom already has %d hydrogen bonds and has %d lone pairs\n    ",
	  acceptor_count[acceptor],no_lone_pairs[acceptor->element]);
	acceptor->display();*/
	res = false;
      }
  

  return res;

}




int Hbonds::number_bounded_hydrogens(Atom * a)
{
  return no_bounded_hydrogens[a];
}


void Hbonds::make_list_of_number_of_bonded_hydrogens() {
  for(unsigned int a=0; a<atomsA.size(); a++) {
      int res =0;
      for(unsigned int b=0; b<bondsA.size(); b++)
	if(bondsA[b][0]==(int)a && atomsA[bondsA[b][1]]->element=="H")
	  res++;
      no_bounded_hydrogens[atomsA[a]] = res;
      //      cout<<"A--> "<<res<<" hydrogens found for ";
      //      atomsA[a]->display();
    }

  for(unsigned int a=0; a<atomsB.size(); a++) {
      int res =0;
      for(unsigned int b=0; b<bondsB.size(); b++)
	if(bondsB[b][0]==(int)a && atomsB[bondsB[b][1]]->element=="H")
	  res++;
      no_bounded_hydrogens[atomsB[a]] = res;

      //      cout<<"B--> "<<res<<" hydrogens found for ";
      //      atomsB[a]->display();
    }
}


float Hbonds::distance_factor(float distance)
{
  float factor = 1.0;
  
  if(2.8<distance && distance<3.5)
    factor = (1/0.7)*(3.5-distance);
  else if(distance>= 3.5)
    factor = 0;
  
  return factor;
}

float Hbonds::angle_factors(int a,int b)
{
  float factor = 1.0;

  /*** positions of atoms ***/  
  vector<float> atomA,atomB;
  atomA.assign(3,0.0);
  atomB.assign(3,0.0);
  
  atomA[0] = atomsA[a]->x;
  atomA[1] = atomsA[a]->y;
  atomA[2] = atomsA[a]->z;

  atomB[0] = atomsB[b]->x;
  atomB[1] = atomsB[b]->y;
  atomB[2] = atomsB[b]->z;
  
  /*** position of roots ***/
  vector<float> rootA = find_root_A(a); 
  vector<float> rootB = find_root_B(b); 

  /*** angles ***/
  float angleA = find_angle(rootA,atomA,atomB);
  float angleB = find_angle(atomA,atomB,rootB);

  /*** factors ***/
  factor = factor*angle_factor(angleA);
  factor = factor*angle_factor(angleB);
  
  
  /*************************************/
//   printf("atomA: (%6f,%6f,%6f)\n",atomA[0],atomA[1],atomA[2]);
//   printf("atomB: (%6f,%6f,%6f)\n",atomB[0],atomB[1],atomB[2]);


//   printf("rootA: (%6f,%6f,%6f)\n",rootA[0],rootA[1],rootA[2]);
//   printf("rootB: (%6f,%6f,%6f)\n",rootB[0],rootB[1],rootB[2]);

//   printf("angleA: %6f\n",angleA);
//   printf("angleB: %6f\n",angleB);

//   printf("factorA: %6f\n",angle_factor(angleA));
//   printf("factorB: %6f\n",angle_factor(angleB));



  /***************************************/


  return factor;

}


float Hbonds::angle_factors_internal(int a,int b)
{
  float factor = 1.0;

  /*** positions of atoms ***/  
  vector<float> atomA,atomB;
  atomA.assign(3,0.0);
  atomB.assign(3,0.0);
  
  atomA[0] = atomsA[a]->x;
  atomA[1] = atomsA[a]->y;
  atomA[2] = atomsA[a]->z;

  atomB[0] = atomsA[b]->x;
  atomB[1] = atomsA[b]->y;
  atomB[2] = atomsA[b]->z;
  
  /*** position of roots ***/
  vector<float> rootA = find_root_A(a); 
  vector<float> rootB = find_root_A(b); 

  /*** angles ***/
  float angleA = find_angle(rootA,atomA,atomB);
  float angleB = find_angle(atomA,atomB,rootB);

  /*** factors ***/
  factor = factor*angle_factor(angleA);
  factor = factor*angle_factor(angleB);
  
  
  /*************************************/
//   printf("atomA: (%6f,%6f,%6f)\n",atomA[0],atomA[1],atomA[2]);
//   printf("atomB: (%6f,%6f,%6f)\n",atomB[0],atomB[1],atomB[2]);


//   printf("rootA: (%6f,%6f,%6f)\n",rootA[0],rootA[1],rootA[2]);
//   printf("rootB: (%6f,%6f,%6f)\n",rootB[0],rootB[1],rootB[2]);

//   printf("angleA: %6f\n",angleA);
//   printf("angleB: %6f\n",angleB);

//   printf("factorA: %6f\n",angle_factor(angleA));
//   printf("factorB: %6f\n",angle_factor(angleB));



  /***************************************/


  return factor;
}


float Hbonds::angle_factor(float angle)
{
  float factor = 1.0;

  if(60<angle && angle<120)
    factor = (1.0/60.0)*(angle-60.0);
  else if(angle <= 60)
    factor = 0;

  return factor;

}


float Hbonds::find_angle(vector<float> p,vector<float> q,vector<float> r)
{
  
  /*** find inter-point vectors ***/
  vector<float> vectorA,vectorB;
  vectorA.assign(3,0.0);
  vectorB.assign(3,0.0);
 
  for(unsigned int i=0; i<3; i++)
    {
      vectorA[i] = p[i] - q[i];
      vectorB[i] = r[i] - q[i];      
    }

  /*** find dot product ***/
  float dot = 0;
  for(unsigned int i=0; i<3; i++)
    dot += vectorA[i]*vectorB[i];

  /*** find lenghts ***/
  float lenghtA = find_lenght(vectorA);
  float lenghtB = find_lenght(vectorB);
  
  /*** find angle ***/
  float v = acos(dot/(lenghtA*lenghtB));

  /*** convert to degrees ***/
  v = v * piConv;

  return v;
}


float Hbonds::find_lenght(vector<float> a)
{
  float res = 0;
  for(unsigned int i=0; i<3; i++)
    res += a[i]*a[i];

  res = sqrt(res);
  
  return res;
}

vector<float> Hbonds::find_root_A(int c) {

  float count = 0;
  vector<float> res;
  res.assign(3,0.0);


  for(unsigned int i=0; i<bondsA.size(); i++)
    if(bondsA[i][1] == c)
      if(atomsA[bondsA[i][0]]->element != "H")
	{
	  count = count + 1;
	  res[0] += atomsA[bondsA[i][0]]->x;
	  res[1] += atomsA[bondsA[i][0]]->y;
	  res[2] += atomsA[bondsA[i][0]]->z;
	}



//   //find bonded atoms
//   for(unsigned int a=0; a<atomsA.size(); a++)
//     for(unsigned int i=0; i<bondsA.size(); i++)
//       {
// 	if(bondsA[i][0] == a && bondsA[i][1] == c)
// 	  if(atomsA[a]->element != "H")
// 	    {
// 	      count = count + 1;
// 	      res[0] += atomsA[a]->x;
// 	      res[1] += atomsA[a]->y;
// 	      res[2] += atomsA[a]->z;
// // 	    printf("Including in rootA: %s-%d at (%6f,%6f,%6f)\n",
// // 		  atomsA[a]->element.c_str(),
// // 		  i+1,
// // 		  atomsA[a]->x,
// // 		  atomsA[a]->y,
// // 		  atomsA[a]->z);
// 	  }
//     }

  for(unsigned int i=0;i<res.size();i++)
    res[i]=res[i]/count;

  return res;
}

vector<float> Hbonds::find_root_B(int c)
{

  float count = 0;
  vector<float> res;
  res.assign(3,0.0);

  for(unsigned int i=0; i<bondsB.size(); i++)
    if(bondsB[i][1] == c)
      if(atomsB[bondsB[i][0]]->element != "H")
	{
	  count = count + 1;
	  res[0] += atomsB[bondsB[i][0]]->x;
	  res[1] += atomsB[bondsB[i][0]]->y;
	  res[2] += atomsB[bondsB[i][0]]->z;
	}
  
//   //find bonded atoms
//   for(unsigned int b=0; b<atomsB.size(); b++)
//     for(unsigned int i=0; i<bondsB.size(); i++)
//     {
//       if(bondsB[i][0] == b && bondsB[i][1] == c)
// 	if(atomsB[b]->element != "H")
// 	  {
// 	    count = count + 1;
// 	    res[0] += atomsB[b]->x;
// 	    res[1] += atomsB[b]->y;
// 	    res[2] += atomsB[b]->z;
// 	    // 	    printf("Including in rootB: %s-%d at (%6f,%6f,%6f)\n",
// 	    // 		   atomsB[b]->element.c_str(),
// 	    // 		   i+1,
// 	    // 		   atomsB[b]->x,
// 	    // 		  atomsB[b]->y,
// // 		  atomsB[b]->z);

// 	  }
//     }

  for(unsigned int i=0;i<res.size();i++)
    res[i]=res[i]/count;

  return res;
}



void Hbonds::set_parameters()
{
  no_lone_pairs["N"] = 1;
  no_lone_pairs["O"] = 2;
  no_lone_pairs["F"] = 3;

  no_lone_pairs["P"] = 1;
  no_lone_pairs["S"] = 2;
  no_lone_pairs["Cl"] = 3;


  piConv = 90/acos(0.0);


  ifstream parameters;
  parameters.open("hbonds.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'hbonds.cfg' could not be found\n";
      exit(0);
    }


  string dummy;

  while (!parameters.eof())
    {

      parameters >> dummy;

      if(dummy == "POSSIBLE_ACCEPTORS:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      possible_acceptors.push_back(dummy);
	      parameters >> dummy;
	    }
	}

      if(dummy == "POSSIBLE_DONORS:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      possible_donors.push_back(dummy);
	      parameters >> dummy;
	    }
	}


      if(dummy == "MIN_DIST:")
	parameters >> min_dist;
      if(dummy == "MAX_DIST:")
	parameters >> max_dist;
   
      if(dummy == "INCLUDE_SURF_FACTOR:")
	parameters >> include_surf_factor;
      if(dummy == "INCLUDE_ANGLE_FACTOR:")
	parameters >> include_angle_factor;
      if(dummy == "INCLUDE_DIST_FACTOR:")
	parameters >> include_dist_factor;
      
      if(dummy == "REPLACEMENT_BY_WATER_FACTOR:")
	parameters >>replacement_by_water_factor;


    } //end read file

  //printing all configurations
  /*
  printf("Configurations for Hbonds:\n");
  printf("\n");
  printf("\n");

  printf("POSSIBLE_ACCEPTORS:     ");
  for(unsigned int i=0; i<possible_acceptors.size(); i++)
    printf("%s ",possible_acceptors[i].c_str());
  printf("end\n");

  printf("POSSIBLE_DONORS:        ");
  for(unsigned int i=0; i<possible_donors.size(); i++)
    printf("%s ",possible_donors[i].c_str());
  printf("end\n");
  printf("\n");

  printf("MIN_DIST:               %f\n", min_dist);
  printf("MAX_DIST:               %f\n", max_dist);
  printf("\n");

  printf("INCLUDE_SURF_FACTOR:     %d\n", include_surf_factor);
  printf("\n");
  printf("INCLUDE_ANGLE_FACTOR:    %d\n", include_angle_factor);
  printf("\n");
  printf("INCLUDE_DIST_FACTOR:     %d\n", include_dist_factor);
  printf("\n");


  printf("Version: $Id: hbonds.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");
  */


  min_dist = min_dist*min_dist;
  max_dist = max_dist*max_dist;
  
}


void Hbonds::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: hbonds.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
  
  rep = reporter;
}
