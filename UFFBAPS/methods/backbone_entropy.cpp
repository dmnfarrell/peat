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
$Id: backboneEntropy.cpp 18 2005-11-15 09:56:08Z chresten $
*/

#include "backbone_entropy.h"

BackboneEntropy::BackboneEntropy():Method()
{
  read_in_maps();
  precalculate_energies();
  
}

BackboneEntropy::BackboneEntropy(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  read_in_maps();
  precalculate_energies();
}

BackboneEntropy::BackboneEntropy(FILE * reporter, vector<Soup*> A, bool train)
{
  /*** write out to out file **/
  writeDescription(reporter);
  
  if(train)
    do_statistics(A);
  else
    {
      read_in_maps();
      precalculate_energies();
      calculate(A);
    }
}


BackboneEntropy::~BackboneEntropy(){}


void BackboneEntropy::calculate(vector<Soup*> A)
{


  for(vector<Soup*>::iterator soup = A.begin(); soup != A.end(); soup++)
    calculate(*soup);

  return;
}

float BackboneEntropy::calculate(Soup * A)
{

  vector<ProteinChain * > protein_chains = A->getProteinChains();
  float res = 0.0;

  for(vector<ProteinChain*>::iterator protein_chain = protein_chains.begin(); protein_chain != protein_chains.end(); protein_chain++)
    res += calculate(*protein_chain);
  
  //cout<<"Soup "<<A->name<<" has backbone entropy "<<res<<endl;
  return res;
  
}

float BackboneEntropy::calculate(ProteinChain * protein_chain)
{

  vector<Residue*> residues = protein_chain->getResidues();
  float res = 0.0;

  if(residues.size()>2)
    for(vector<Residue*>::iterator residue = residues.begin()+1; residue!=residues.end()-1; residue++)
      if(set_phi_and_psi(*(residue-1),*residue,*(residue+1)))
	{
	  res += propensity_maps[(*residue)->residueType].energies
	    [make_index((*residue)->phi)]
	    [make_index((*residue)->psi)];
	  
	}
  return res;
}




void BackboneEntropy::read_in_maps()
{

  string array_residues[] = 
    {"ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};

  list<string> residues (array_residues, array_residues + 20);  
  
  // for each residue
  for(list<string>::iterator res= residues.begin(); res != residues.end(); res++)
    {
      // construct filename
      string file = "parameters/"+*res+"_propensity_map.m";
      //cout<<file<<endl;

      // make propensity map for residue type
      propensity_maps[*res] = new_propensity_map();
      
      
      // read in each line
      ifstream ind(file.c_str()); 
      char line[200];
      vector<string> words;
      while(!ind.eof())
	{
	  ind.getline(line,200);
	  // and break it into words
	  words.clear();
	  char temp[100];
	  istringstream stream(line);
	  while(stream>>temp)
	    words.push_back(removeWhite(string(temp)));

	  // if we get two words (phi and psi angles)
	  if(words.size() == 2)
	    {
	      // cout<<words[0]<<" "<<words[1]<<endl;
	      // insert the phi and psi angles
	      insert_angles(&propensity_maps[*res],atof(words[0].c_str()), atof(words[1].c_str()));
	    }
	  
	}
      
    }
  


  /// checkiing
  /*
  for(map<string,propensity_map>::iterator propensity_map = propensity_maps.begin();
      propensity_map != propensity_maps.end(); propensity_map++)
    {
      print_map(&propensity_map->second, propensity_map->first+"_propensity_map.m");

      cout<<propensity_map->second.phis.size()<<" angles found for "<<propensity_map->first<<endl;
    }
  print_all_maps();
  */



}


void BackboneEntropy::precalculate_energies()
{
  float upper_limit = 20; // kJ/mol
  float R = 8.31451e-3; // kJ/(K mol)
  float T = 300; // K

  float offset_energy = -10; //kJ/mol  -> Offset added to all energy values to avoid too big penalties

  // for each propensity map
  for(map<string,propensity_map>::iterator propensity_map = propensity_maps.begin();
      propensity_map != propensity_maps.end(); propensity_map++)
    {
      // get number of keys in the map
      int n_total = propensity_map->second.phis.size();
      //cout<<propensity_map->first<<" has "<<n_total<<" angles "<<endl;

      // make an energy map
      for(unsigned int i=0; i<propensity_map->second.angles.size();i++)
	for(unsigned int j=0; j<propensity_map->second.angles[i].size();j++)
	  {
	    if (propensity_map->second.angles[i][j] == 0)
	      propensity_map->second.energies[i][j] = upper_limit;
	    else
	      propensity_map->second.energies[i][j] = -R*T*log(((float) propensity_map->second.angles[i][j])/((float) n_total));
	    
	    if (propensity_map->second.energies[i][j]>upper_limit)
	      propensity_map->second.energies[i][j] = upper_limit;

	    propensity_map->second.energies[i][j] += offset_energy; // Add offset!


	    //cout<<"   angles["<<i<<"]["<<j<<"]: "<<propensity_map->second.angles[i][j]
	    //<<" energies["<<i<<"]["<<j<<"]: "<<propensity_map->second.energies[i][j]<<endl;
	  }
      
    }


}



void BackboneEntropy::do_statistics(vector<Soup*> A)
{
  
  for(vector<Soup*>::iterator soup = A.begin(); soup != A.end(); soup ++)
    {
      vector<ProteinChain*> protein_chains = (*soup)->getProteinChains();
      for(vector<ProteinChain*>::iterator protein_chain = protein_chains.begin(); protein_chain != protein_chains.end(); protein_chain++)
	{
	  vector<Residue*> residues = (*protein_chain)->getResidues();
	  if(residues.size()>2)
	    for(vector<Residue*>::iterator residue = residues.begin()+1; residue != residues.end()-1; residue++)
	      {
		// calculate phi and psi angles if posible 
		if(set_phi_and_psi(*(residue-1),*residue,*(residue+1)))
		  {
		    printf("%3s  %5.2f  %5.2f\n",
			   (*residue)->residueType.c_str(),
			   (*residue)->phi,
			   (*residue)->psi);
		    
		    // make new propensity map, if one does not exist for the current reisdue type
		    if (propensity_maps.count((*residue)->residueType)==0)
		      propensity_maps[(*residue)->residueType] = new_propensity_map();
		    
		    // insert the phi and psi angles
		    insert_angles(&propensity_maps[(*residue)->residueType],(*residue)->phi, (*residue)->psi);

		  }
		
	    }
	}
    }

  // write out stats
  for(map<string,propensity_map>::iterator propensity_map = propensity_maps.begin();
      propensity_map != propensity_maps.end(); propensity_map++)
    {
      print_map(&propensity_map->second, propensity_map->first+"_propensity_map.m");

      cout<<propensity_map->second.phis.size()<<" angles found for "<<propensity_map->first<<endl;
    }
  print_all_maps();

  

  return;
}


void BackboneEntropy::print_all_maps()
{
  ofstream ud("all_propensity_map.m",fstream::app); 

  //  for(int i=0; i<map->angles.size(); i++)
  //    for(int j=0; j<map->angles[i].size(); j++)
  //      ud<<setw(5)<<i*18+9<<" "<<setw(5)<<j*18+9<<" "<<setw(10)<<map->angles[i][j] <<endl;

  for(map<string,propensity_map>::iterator propensity_map = propensity_maps.begin();
      propensity_map != propensity_maps.end(); propensity_map++)
    {
      vector<float>::iterator phi = (propensity_map->second).phis.begin();
      vector<float>::iterator psi = (propensity_map->second).psis.begin();
      
      while(phi!=propensity_map->second.phis.end() && psi!=propensity_map->second.psis.end())
	{
	  ud<<setw(6)<<*phi<<" "<<setw(6)<<*psi<<endl;
	  phi++;
	  psi++;
	}
    }
  ud.close();
  return;

}

void BackboneEntropy::print_map(propensity_map * map, string filename)
{

  ofstream ud(filename.c_str(),fstream::app); 

  //  for(int i=0; i<map->angles.size(); i++)
  //    for(int j=0; j<map->angles[i].size(); j++)
  //      ud<<setw(5)<<i*18+9<<" "<<setw(5)<<j*18+9<<" "<<setw(10)<<map->matrix[i][j] <<endl;
  
  vector<float>::iterator phi = map->phis.begin();
  vector<float>::iterator psi = map->psis.begin();

  while(phi!=map->phis.end() && psi!=map->psis.end())
    {
      ud<<setw(6)<<*phi<<" "<<setw(6)<<*psi<<endl;
      phi++;
      psi++;
    }

  ud.close();
  return;
}

void BackboneEntropy::insert_angles(propensity_map * map, float phi, float psi)
{
  map->phis.push_back(phi);
  map->psis.push_back(psi);

  int i = make_index(phi), j = make_index(psi);
  //cout<<"indexes: "<<i<<" "<<j<<endl;

  map->angles[i][j] ++;
  
  return;
}

int BackboneEntropy::make_index(float angle)
{

  int res = (int) floor(angle/18)+10;
  
  if(res>19) //in case angle = 180
    res=19;
  if(res<0) //in case angle = -180
    res=0;

  return res;
}


propensity_map BackboneEntropy::new_propensity_map()
{
  propensity_map res;

  vector<int> atemp(20,0);
  vector<vector<int> > atemp2(20,atemp);
  res.angles = atemp2;

  vector<float> etemp(20,0);
  vector<vector<float> > etemp2(20,etemp);
  res.energies = etemp2;

  return res;
}



bool BackboneEntropy::set_phi_and_psi(Residue * prev_residue, Residue * residue, Residue * next_residue)
{
  // set default values in case something goes wrong

  residue->phi = 1000;
  residue->psi = 1000;

  //Find atoms
  Atom 
    * Cp = new Atom,
    * N  = new Atom,
    * CA = new Atom, 
    * C  = new Atom,
    * Nn = new Atom;

  prev_residue->find_atom_named("C",Cp);
  residue->find_atom_named("N",N);
  residue->find_atom_named("CA",CA);
  residue->find_atom_named("C",C);
  next_residue->find_atom_named("N",Nn);

  //check that all atoms were found
  if(Cp->name=="Dummy" || N->name=="Dummy")
    return false;
  if(CA->name=="Dummy" || C->name=="Dummy")
    return false;
  if(Nn->name=="Dummy")
    return false;

  // make vectors
  //
  //      H
  //      |
  //      Ca
  //     /  \
  // Cp-N    C-Nn
  //    |    ||
  //    H    O

  la::vector<float> CpN(3);
  CpN[0] = (N)->x - (Cp)->x;
  CpN[1] = (N)->y - (Cp)->y;
  CpN[2] = (N)->z - (Cp)->z;
 
  la::vector<float> NCa(3);
  NCa[0] = (CA)->x - (N)->x;
  NCa[1] = (CA)->y - (N)->y;
  NCa[2] = (CA)->z - (N)->z;
 
  la::vector<float> CaC(3);
  CaC[0] = (C)->x - (CA)->x;
  CaC[1] = (C)->y - (CA)->y;
  CaC[2] = (C)->z - (CA)->z;
 
  la::vector<float> CNn(3);
  CNn[0] = (Nn)->x - (C)->x;
  CNn[1] = (Nn)->y - (C)->y;
  CNn[2] = (Nn)->z - (C)->z;

  // make normal vectors
  la::vector<float> n_phi_1 = cross(CpN,NCa);
  la::vector<float> n_phi_2 = cross(NCa,CaC);

  la::vector<float> n_psi_1 = cross(NCa,CaC);
  la::vector<float> n_psi_2 = cross(CaC,CNn);

  // find dihedral angles
  float phi_prod = dot(n_phi_1,n_phi_2)/(norm(n_phi_1)*norm(n_phi_2));
  float phi_sign = dot(cross(n_phi_1,n_phi_2),NCa);
  phi_sign = phi_sign/fabs(phi_sign);

  residue->phi = phi_sign*acos(phi_prod)*90/acos(0);
  
  float psi_prod = dot(n_psi_1,n_psi_2)/(norm(n_psi_1)*norm(n_psi_2));
  float psi_sign = dot(cross(n_psi_1,n_psi_2),CaC);
  psi_sign = psi_sign/fabs(psi_sign);
  
  residue->psi = psi_sign*acos(psi_prod)*90/acos(0);

  return true;
}


la::vector<float> BackboneEntropy::cross(la::vector<float> a, la::vector<float> b)
{
  la::vector<float> res(3);
  res[0] = a[1]*b[2] - a[2]*b[1];// (a2b3 - a3b2);
  res[1] = a[2]*b[0] - a[0]*b[2];// (a3b1 - a1b3);
  res[2] = a[0]*b[1] - a[1]*b[0];// (a1b2 - a2b1)

  return res;
}

float BackboneEntropy::dot(la::vector<float> a, la::vector<float> b)
{
  float res = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return res;
}

float BackboneEntropy::norm(la::vector<float> a)
{
  float res = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  return res;
}


void BackboneEntropy::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: backboneEntropy.cpp 18 2005-11-15 09:56:08Z chresten $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
