/*
 #
 # FFF - Flexible Force Field
 # Copyright (C) 2010 Jens Erik Nielsen
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
#include "Model.h"

using namespace std;

vector<double> model_class::Mutate(const std::string chainname,const std::string pdbresnumber, const std::string newresidue,int mode,double max_clash) {
    for (unsigned int chain=0;chain<_P.chains.size();chain++) {
      if (strip(chainname)==strip(_P.chains[chain].name)) {
            for (unsigned int residue=0;residue<_P.chains[chain].residues.size();residue++) {
                //printf ("prsnumber: %s \n",strip(_P.chains[chain].residues[residue].pdbnum).c_str());
                if (strip(pdbresnumber)==strip(_P.chains[chain].residues[residue].pdbnum)) {
                    return Mutate_2(chain,residue,newresidue,mode,max_clash);
                }
            }
        }
    }
    printf ("Could not find residue to mutate: '%s' '%s'\n",chainname.c_str(),pdbresnumber.c_str());
    vector<double> dummy;
    return dummy;
}

vector<int> model_class::get_chain_and_residue(const std::string resid) {
  // Given a Protool residue ID, get the corresponding FFF residuenumber and chainnumber
  vector<int> numbers;
  string chainname(resid,0,1);
  string pdbresnumber(resid,2,6);
  if (chainname==":") {
    chainname="";
    string pdbresnumber(resid,1,5);
  }
  printf ("Extracted chain name: %s, residue number: %s",chainname.c_str(),pdbresnumber.c_str());
  //
  for (unsigned int chain=0;chain<_P.chains.size();chain++) {
    if (strip(chainname)==strip(_P.chains[chain].name)) {
      numbers.push_back(chain);
      for (unsigned int residue=0;residue<_P.chains[chain].residues.size();residue++) {
	//printf ("prsnumber: %s \n",strip(_P.chains[chain].residues[residue].pdbnum).c_str());
	if (strip(pdbresnumber)==strip(_P.chains[chain].residues[residue].pdbnum)) {
	  numbers.push_back(residue);
	}
      }
    }
  }
  return numbers;
}

//
// -------------
//


vector<double> model_class::Mutate_2(int chainnumber,int resnumber,const std::string newresidue,int mode,double max_clash) {
  //
  // Mutate a single residue
  // mode=0: Just place the new side chain in the protein
  // mode=1:  do + maintain original chi angles
  // mode=2: Explore all possible combinations of side chain chi-angles
  // mode=3: Loop through all the conformations in the rotamer library and select 
  // the one with the lowest energy
  //
    vector<double> energy;
  switch (mode) {
  case 0: {
    _mutate(chainnumber,resnumber,newresidue);
    break;
  }
  case 1: {
    _save_chiangles(chainnumber,resnumber);
    _mutate(chainnumber,resnumber,newresidue);
    _restore_chiangles(chainnumber,resnumber);
    break;
  }
      case 2: {
    //Brute force exhaustive search
        double chistep=30.0;
      _save_chiangles(chainnumber,resnumber);
      _mutate(chainnumber,resnumber,newresidue);
      _restore_chiangles(chainnumber,resnumber);
          double chinums=static_cast<double>(getnumchi(chainnumber,resnumber));
          chistep=fmax(0.1,360.0/pow(1000.0,1.0/chinums));
      energy=_optimise1res_exhaustive(chainnumber,resnumber,chistep,max_clash);
          break;
  }
      case 3: {
          // Use rotamer library but debump only
          _save_chiangles(chainnumber,resnumber);
          _mutate(chainnumber,resnumber,newresidue);
          _restore_chiangles(chainnumber,resnumber);
          energy=_optimise1res(chainnumber,resnumber,max_clash);
          break;
      }
  }
  return energy;
}

//
// -----------------------------------------------------------
//

void model_class::_mutate(int chainnumber,int resnumber,string newresidue) {
    //
    // Mutate residue <resnumber> to <newresidue> - Takes the coordinates from the general
    // amino acid definition file
    //
    // This is the routine that actually changes the atom names and the coordinates
    // Don't call this function, but call Mutate instead
    //
  
    //
    // Make backup of original residue
    //
    string oldresiduename=_P.chains[chainnumber].residues[resnumber].name;
    vector<atom_class> oldatoms;
    for (unsigned int oldatom=0;oldatom<_P.chains[chainnumber].residues[resnumber].atoms.size();oldatom++) {
        oldatoms.push_back(_P.chains[chainnumber].residues[resnumber].atoms[oldatom]);
        // Put in delete atoms
	//delete_atoms.push_back(
      
    }
  _org_atoms.push_back(oldatoms);
  //
  // Deal with the simple cases
  // Asp <-> Asn, Gln <-> Glu and so on ..
  //
  if (simple_mutate(chainnumber,resnumber,oldresiduename,newresidue)) {
    return;
  }
  //
  // Prepare the new coordinates
  //
  int respointer=-1;
  //printf ("Number of AA definitions: %d \n",static_cast<int>(AAsets[0].aadefs.size()));
  for (int tempres=0;tempres<static_cast<int>(AAsets[0].aadefs.size());tempres++) {
    if (newresidue==AAsets[0].aadefs[tempres].name) {
      respointer=tempres;
      break;
    }
  }
  if (respointer==-1) {
    printf ("New residue %s not found.\n",newresidue.c_str());
    throw FFFError();
  }
  //
  // Mutate_special takes care of waters, ions etc. The rest
  // is handled here in this routine
  //
  vector<atom_class> newcoords;
  //    
  if (!_P.chains[chainnumber].residues[resnumber].is_aa) {
    newcoords=mutate_special(chainnumber,resnumber,newresidue);
  } 
  else {
    //
    // Find the N, CA & C coordinates in the protein and in the definition and store them in a vector
    //
    //
    // Coordinates from the protein
    //
    string old_pdbresnum;
    vector<atom_class> refcoords;
    atom_class N("N",
		 _P.chains[chainnumber].residues[resnumber].N->x,
		 _P.chains[chainnumber].residues[resnumber].N->y,
		 _P.chains[chainnumber].residues[resnumber].N->z);
    refcoords.push_back(N);
    old_pdbresnum=_P.chains[chainnumber].residues[resnumber].N->pdbresnum;
    //
    atom_class CA("CA",
		  _P.chains[chainnumber].residues[resnumber].CA->x,
		  _P.chains[chainnumber].residues[resnumber].CA->y,
		  _P.chains[chainnumber].residues[resnumber].CA->z);
    refcoords.push_back(CA);
    //
    atom_class C("C",
		 _P.chains[chainnumber].residues[resnumber].C->x,
		 _P.chains[chainnumber].residues[resnumber].C->y,
		 _P.chains[chainnumber].residues[resnumber].C->z);
    refcoords.push_back(C);
    //
    // Coordinates from the defintion
    //
    vector<atom_class> fitcoords;
    fitcoords.push_back(AAsets[0].aadefs[respointer].N());
    fitcoords.push_back(AAsets[0].aadefs[respointer].CA());
    fitcoords.push_back(AAsets[0].aadefs[respointer].C());
    //
    // Superpose the N, CA and C
    //
    quatfit_class Q;
    Q.fit(3,refcoords,fitcoords);
    //
    // Now transform the coordinates of the entire amino acid definition
    //
    fitcoords.resize(0);
    int numcoords=static_cast<int>(AAsets[0].aadefs[respointer].atoms.size());
    for (int atom=0;atom<numcoords;atom++) {
      fitcoords.push_back(AAsets[0].aadefs[respointer].atoms[atom]);
    }
    Q.transform(numcoords,fitcoords);
    for (int atom=0;atom<numcoords;atom++) {
      AAsets[0].aadefs[respointer].atoms[atom].x=fitcoords[atom].x;
      AAsets[0].aadefs[respointer].atoms[atom].y=fitcoords[atom].y;
      AAsets[0].aadefs[respointer].atoms[atom].z=fitcoords[atom].z;
    }
    //
    // Preserve the coordinates for the backbone atoms (identical names required)
    //
    for (int atom=0;atom<numcoords;atom++) {
      atom_class NEWATOM=AAsets[0].aadefs[respointer].atoms[atom];
      bool found=false;
      //
      // See which atoms to copy in
      //
      for (unsigned int oldatom=0;oldatom<_P.chains[chainnumber].residues[resnumber].atoms.size();oldatom++) {
	string oldname(_P.chains[chainnumber].residues[resnumber].atoms[oldatom].name);
	if (strip(NEWATOM.name)==strip(oldname)) {
	  //
	  // Identical names. Find out which atom to choose
	  //
	  found=true;
	  if (NEWATOM.is_backbone()) {	
	    atom_class newatom(oldname,
			       _P.chains[chainnumber].residues[resnumber].atoms[oldatom].x,
			       _P.chains[chainnumber].residues[resnumber].atoms[oldatom].y,
			       _P.chains[chainnumber].residues[resnumber].atoms[oldatom].z);
	    newcoords.push_back(newatom);
	  }
	  else {
	    //
	    // No, it is not a backbone atom, so we want the new coordinates
	    //
	    atom_class newatom(NEWATOM.name,NEWATOM.x,NEWATOM.y,NEWATOM.z);
	    if (_use_hydrogens) {
	      if (not NEWATOM.is_hydrogen()) {
		newcoords.push_back(newatom);
	      }
	    } else {
	      newcoords.push_back(newatom);
	    }
	  }
	}
      }   
      //
      // If we did not find an atom with the same name, then add the new atom
      //
      if (!found) {
	if (!NEWATOM.is_hydrogen() || NEWATOM.name=="HN" || _use_hydrogens==true) {
	  atom_class newatom(NEWATOM.name,NEWATOM.x,NEWATOM.y,NEWATOM.z);
	  newcoords.push_back(newatom);
	}
      }
    }
  }
  //
  // Fill in the new coordinates for this residue
  //
  _P.chains[chainnumber].residues[resnumber].atoms.resize(0); //Remove all old atoms
  for (int atom=0;atom<static_cast<int>(newcoords.size());atom++) {
    _P.chains[chainnumber].residues[resnumber].atoms.push_back(newcoords[atom]);
  } 
  //
  // Change the resiue name
  //
  _P.chains[chainnumber].residues[resnumber].name=newresidue;
  //
  // Calculate the new chiangles
  //
  _all_chi_angles[chainnumber][resnumber]=_calc_residue_chiangles(chainnumber,resnumber);
  //
  // Call the residue update function
  //
  _P.chains[chainnumber].residues[resnumber].update();
    update_bonds();
   return;
};

//
// ----------
//

vector<atom_class> model_class::mutate_special(int chainnumber,int resnumber,string newresidue) {
  // 
  // Mutate groups that are not amino acids
  //
  unsigned int respointer=_get_AA_template(newresidue);
  //
  // Find the atoms that are needed to superimpose coordinates
  //
  vector<atom_class> add_atoms;
  string resname=_P.chains[chainnumber].residues[resnumber].name;
  if (resname=="WAT") {
    int atnum=_P.chains[chainnumber].residues[resnumber].get_atom_number("O");
    atom_class refcoord("O",
			_P.chains[chainnumber].residues[resnumber].atoms[atnum].x,
			_P.chains[chainnumber].residues[resnumber].atoms[atnum].y,
			_P.chains[chainnumber].residues[resnumber].atoms[atnum].z);
    //
    // Get the definition coords
    //
    atom_class fitcoord=AAsets[0].aadefs[respointer].get_atom("O");
    // 
    // Calculate offset
    //
    vector3d offset=refcoord-fitcoord;
    //
    // Add all the atoms in the definition
    //
    for (unsigned atom=0;atom<AAsets[0].aadefs[respointer].atoms.size();atom++) {
      atom_class NEWATOM=AAsets[0].aadefs[respointer].atoms[atom];
      NEWATOM.x=NEWATOM.x+offset.x;
      NEWATOM.y=NEWATOM.y+offset.y;
      NEWATOM.z=NEWATOM.z+offset.z;
      add_atoms.push_back(NEWATOM);
    }
  }
  return add_atoms;
}
    


//
// -----------
//

void model_class::undo_mutation() {
  // 
  // Undo the last mutation
  //
  unsigned int last_record=_org_atoms.size()-1;
  vector<atom_class> restore_atoms=_org_atoms[last_record];
  //
  // Get the residue(s)
  //
  for (unsigned int atom=0;atom<restore_atoms.size();atom++) {
    int resnumber=restore_atoms[atom].inresidue;
    int chainnumber=restore_atoms[atom].inchain;
    //
    // Get rid of all atoms in those residues 
    //
    _P.chains[chainnumber].residues[resnumber].atoms.resize(0); //Remove all old atoms
  }
  // 
  // Reinstate old atoms
  //
  for (unsigned int atom=0;atom<restore_atoms.size();atom++) {
    int resnumber=restore_atoms[atom].inresidue;
    int chainnumber=restore_atoms[atom].inchain;
    _P.chains[chainnumber].residues[resnumber].atoms.push_back(restore_atoms[atom]);
    // change the residue name
    _P.chains[chainnumber].residues[resnumber].name=restore_atoms[atom].inresiduename;
  }
  //printf ("----end---\n");
  //
  // Update all the residues
  //
  return;
}

//
// ----------------
//

bool model_class::simple_mutate(int chain,int resnum,string oldres,string newres) {
  //
  // Perform simple mutations quickly
  //
  bool mutate_done=false;
  bool found=false;
  if (oldres=="ASP" && newres=="ASN") {
    _P.chains[chain].residues[resnum].name="ASN";

    for (unsigned int atom=0;atom<_P.chains[chain].residues[resnum].atoms.size();atom++) {
      if (_P.chains[chain].residues[resnum].atoms[atom].name=="OD2") {
	_P.chains[chain].residues[resnum].atoms[atom].name="ND2";
	found=true;
      }
    }
    mutate_done=true;
    _P.chains[chain].residues[resnum].update();
  }
  // 
  // Could we do it?
  //
  if (not found && mutate_done) {
    printf ("Could not perform simple mutation: %d %d %s %s\n",chain,resnum,oldres.c_str(),newres.c_str());
    throw FFFError();
  }
  return mutate_done;
}

//
// ------
//

vector<double> model_class::_optimise1res_exhaustive(int chainnumber, int resnumber, double chistep,double max_clash) {
    _P.update_BOXLJ(); // We must update the LJ boxes before calculating energies
    // We must update these boxes every time we move, remove or add an atom
    //
    if (_P.chains[chainnumber].residues[resnumber].name=="PRO") {
      printf ("Cannot do brute force optimisation for PRO\n");
      return _ENERGY.get_external_energy(chainnumber,resnumber);
    }
    printf ("Brute force optimisation using chistep of %5.1f degrees\n",chistep);
    vector<double> energy=_ENERGY.get_energy(chainnumber,resnumber);
    printf ("Starting to optimise with clash= %5.2f energy=%5.2f\n",energy[0],energy[1]);
    int optchis=getnumchi(chainnumber,resnumber);
    iteration_count=0;
    if (optchis>0) {
        energy=opt_chi(chainnumber,resnumber,1,chistep,max_clash);
    }
    energy=_ENERGY.get_external_energy(chainnumber,resnumber);
    printf ("Final energy after optimisation: clash %5.2f energy %5.2f\n-----------------------\n",energy[0],energy[1]);
  return energy;
}

//
// -----
//
    
vector<double> model_class::opt_chi(int chainnumber, int resnumber, int chinum, double chistep,double max_clash) {
    // Optimise the current chinum, and call self if we're not the last chinum
    //printf ("Optimising chinum %3d in chain: %3d residue: %5d\n",chinum,chainnumber,resnumber)
    vector<double> new_energy;
    double startchi=getchi(chainnumber,resnumber,chinum);
    double nowchi=startchi;
    int steps=static_cast<int>(360.0/chistep)+1;
    //
    // Array for storing best result until now
    //
    vector<double> best_energy=_ENERGY.get_energy(chainnumber,resnumber);
    vector<double> best_chis;
    for (unsigned int j=static_cast<unsigned int>(chinum);j<=static_cast<unsigned int>(getnumchi(chainnumber,resnumber));j++) {
        best_chis.push_back(getchi(chainnumber,resnumber,j));
    }
    // 
    // Start looping
    //
    for (int step=1;step<steps;step++) {
        nowchi=nowchi+chistep;
        iteration_count++;
        if (nowchi>180.0) nowchi=nowchi-360.0;
        setchi(chainnumber,resnumber,chinum,nowchi);
        if (chinum<getnumchi(chainnumber,resnumber)) {
            // Optimise the next chi angle for this value
            new_energy=opt_chi(chainnumber,resnumber,chinum+1,chistep,max_clash);
        } else {
            new_energy=_ENERGY.get_energy(chainnumber,resnumber);
        }
        // Continue in the loop
        if (new_energy[1]<best_energy[1] && new_energy[0]<=max_clash) {
            best_energy=new_energy;
            best_chis.resize(0);
            for (unsigned int j=chinum;j<=static_cast<unsigned int>(getnumchi(chainnumber,resnumber));j++) {
                best_chis.push_back(getchi(chainnumber,resnumber,j));
            }
        }
    }
    //printf ("Best energy for chinum %2d at E= %7.2f\n",chinum,best_energy);
    int lchi=chinum;
    for (unsigned int j=0;j<best_chis.size();j++) {
        setchi(chainnumber,resnumber,lchi,best_chis[j]);
        lchi++;
    }
    // Check for consistency
    //double energy=_ENERGY.get_energy(chainnumber,resnumber);
    //if (fabs(energy-best_energy)>0.01) {
    //    printf ("Energy is differenct from best_energy: %6.2f %6.2f %6.2f\n",energy,best_energy,fabs(energy-best_energy));
    //}
    return best_energy;
}
        
        

vector<double> model_class::_optimise1res(int chainnumber, int resnumber,double max_clash) {
    // 
    // Optimise one residue using the rotamer library
    //
    _P.update_BOXLJ(); // We must update the LJ boxes before calculating energies
    printf ("Optimizing using rotamer library\n");
    vector<double> best_energy=_ENERGY.get_energy(chainnumber,resnumber);
    vector<double> start_chis;
    for (unsigned int j=0;j<static_cast<unsigned int>(getnumchi(chainnumber,resnumber));j++) {
      start_chis.push_back(getchi(chainnumber,resnumber,j+1));
    }
    //
    printf ("Starting clash: %7.3f energy: %7.3f\n",best_energy[0],best_energy[1]);
    vector<double> new_energy;
    unsigned int best_rotamer=30000;
    // Get the residue name
    string resname=_P.chains[chainnumber].residues[resnumber].name;
    if (resname=="ALA" || resname=="GLY") {
        printf ("Residue with no freedom: %s. Skipping optimisation\n",resname.c_str());
        return _ENERGY.get_external_energy(chainnumber,resnumber);
    }
    // Find the phi/psi angle
    double phi, psi;
    phi=getphi(chainnumber,resnumber);
    psi=getpsi(chainnumber,resnumber);
    printf ("Phi: %5.2f, psi: %5.2f\n",getphi(chainnumber,resnumber),getpsi(chainnumber,resnumber));
    vector<vector<float> > use_rotamers=_ROT.get_rotamers(resname,getphi(chainnumber,resnumber),getpsi(chainnumber,resnumber));

    bool brute=false;
    
    if (phi>180.0 || phi<-180.0 || psi>180.0 || psi<-180.0) {
        printf ("Phi or Psi undefined, defaulting to brute force optimisation");
        brute=true;
    }
    if (use_rotamers.size()==0) {
        printf ("No rotamers found for chain: %d residue: %d with phi: %5.1f, psi: %5.1f. Defaulting to brute force optimisation",chainnumber,resnumber,
                phi,psi);
        brute=true;
    }
    if (brute) {
        double chinums=static_cast<double>(getnumchi(chainnumber,resnumber));
        double chistep=fmax(1.0,360.0/pow(1000.0,1.0/chinums));
        return _optimise1res_exhaustive(chainnumber,resnumber,chistep,max_clash);
    }
    //
    // All ok, do the real rotamer optimisation
    //
    for (unsigned int rotamer=0;rotamer<use_rotamers.size();rotamer++) {
        unsigned int array_index=7;
        //printf ("new rotamer number: %d, number of chis: %d\n",rotamer,getnumchi(chainnumber,resnumber));
        for (int chinum=1;chinum<getnumchi(chainnumber,resnumber)+1;chinum++) {
            setchi(chainnumber,resnumber,chinum,use_rotamers[rotamer][array_index]);
            array_index++;
            //printf ("C: %d, R: %d, Chinum: %d, angle: %7.2f\n",chainnumber,resnumber,chinum,use_rotamers[rotamer][array_index]);
        }
        new_energy=_ENERGY.get_energy(chainnumber,resnumber);
        //printf ("Rotamer: %d, Clash: %7.2f, Energy: %7.3f\n",rotamer,new_energy[0],new_energy[1]);
        if (new_energy[1]<best_energy[1] && new_energy[0]<=max_clash) {
            best_energy=new_energy;
            best_rotamer=rotamer;
        }
    }
    //
    // Get the best rotamer
    //
    if (best_rotamer>29999 || best_energy[0]>1.0) {
      // Restore starting state or best state and do brute force optimisation
      brute=true;
      if (best_rotamer>29999) {
          //printf ("Rotamer optimisation only made things worse. Switching to brute force\n");
          for (unsigned int chinum=0;chinum<start_chis.size();chinum++) {
              setchi(chainnumber,resnumber,chinum+1,start_chis[chinum]);
          }
	}

    }
    //
    // Rotamer opt worked
    //
    if (best_rotamer<29999) {
      unsigned int array_index=7;
      for (int chinum=1;chinum<getnumchi(chainnumber,resnumber)+1;chinum++) {
        setchi(chainnumber,resnumber,chinum,use_rotamers[best_rotamer][array_index]);
        array_index++;
        //printf ("C: %d, R: %d, Chinum: %d, angle: %7.2f\n",chainnumber,resnumber,chinum,use_rotamers[rotamer][array_index]-180.0);
      }
    }
    //
    if (brute) {
      double chinums=static_cast<double>(getnumchi(chainnumber,resnumber));
      double chistep=fmax(1.0,360.0/pow(1000.0,1.0/chinums));
      return _optimise1res_exhaustive(chainnumber,resnumber,chistep,max_clash);
    }
    vector<double> ext_energy=_ENERGY.get_external_energy(chainnumber,resnumber);
    printf ("Best rotamer was # %d, with clash: %7.3f, energy: %7.3f. External clash: %7.3f, energy: %7.3f\n",best_rotamer,best_energy[0],best_energy[1],ext_energy[0],ext_energy[1]);
    return ext_energy;
}
        
    
//
// ----------------------------------------------------------
//

void model_class::_calc_all_angles() {
  //
  // Calculate all chi angles and store them in a vector .. 
  //
  cout << "Calculating all chiangles...\n";
  _all_chi_angles.resize(0);
  for (unsigned int chain=0;chain<_P.chains.size();chain++) {
    vector< vector <double> > _all_chain_angles;
    for (unsigned int residue=0;residue<_P.chains[chain].residues.size();residue++) {
      _all_chain_angles.push_back(_calc_residue_chiangles(chain,residue));
      // Calculate Phi and Psi
        calc_phi(chain,residue);
        calc_psi(chain,residue);
    }
    _all_chi_angles.push_back(_all_chain_angles);
  }
  cout << "done.\n\n";
  return;
}

void model_class::calc_phi(int chain,int residue) {
    vector<int> prev_res;
    double phi=-0.0;
    try {
        prev_res=_get_prev_residue(chain,residue);
           phi=Dihedral(*(_get_atom(chain,prev_res[1],"C")),
                        *(_get_atom(chain,residue,"N")),
                         *(_get_atom(chain,residue,"CA")),
                        *(_get_atom(chain,residue,"C")));
    } catch (FFFError) {
        phi=phi;
    }    
    _P.chains[chain].residues[residue].phi=phi;
    return;
}

void model_class::calc_psi(int chain,int residue) {
    vector<int> next_res;
    double psi=0.0;
    try {
        next_res=_get_next_residue(chain,residue);
        //
        // Get the psi angle for the residue
        //
        // N, CA, C, N+1 
        // 
        psi=Dihedral(*(_get_atom(chain,residue,"N")),
                     *(_get_atom(chain,residue,"CA")),
                     *(_get_atom(chain,residue,"C")),
                     *(_get_atom(chain,next_res[1],"N")));
    } catch (FFFError) {
        psi=-9999.9;
    }
    _P.chains[chain].residues[residue].psi=psi; 
    return;
}

//
// ------------------------------------------------
//
vector<double> model_class::_calc_residue_chiangles(int chainnumber,int residuenumber) {
  vector<double> chis;
  for (int chinum=1;chinum<=getnumchi(chainnumber,residuenumber);chinum++) {
    chis.push_back(calc_chi(chainnumber,residuenumber,chinum));
  }
  return chis;
}

//
// ------------------------------------------------
// 

double model_class::getchi(int chainnumber,int residuenumber, int chinum) {
  //
  // Return the chiangle from the internal array
  //
  // First chiangle is number 1 when calling the routine
  //
  if (chinum>getnumchi(chainnumber,residuenumber) || chinum==0) {
    cout << "This residue: c/r " << chainnumber << "/" << residuenumber << " has only " 
	 << getnumchi(chainnumber,residuenumber) << " Chi angles\n";
    cout << "You requested: " << chinum << endl;
    exit(0);
  }
  return _all_chi_angles[chainnumber][residuenumber][chinum-1];
}

//
// -----------------------------------------------------------
//

atom_class* model_class::_get_atom(int chainnumber,int residuenumber,string searchname) throw(FFFError) {
  //
  // Given a residue number return an instance of atom class
  // that holds the atom with atomname
  //
  //
  // Clean the search name
  //
  searchname=strip(searchname);
  //
  // Go find the atom
  //
  for (unsigned int atom=0;atom<_P.chains[chainnumber].residues[residuenumber].atoms.size();atom++) {
    //
    // See if we have the right name
    //
    if (strip(_P.chains[chainnumber].residues[residuenumber].atoms[atom].name)==searchname) {
      return &_P.chains[chainnumber].residues[residuenumber].atoms[atom];
    }
  }
  //printf ("Atom not found: '%s' in chain: %d residue: %d %s \n",searchname.c_str(),
	//  chainnumber,residuenumber,_P.chains[chainnumber].residues[residuenumber].name.c_str());
  throw FFFError();
}

vector<int> model_class::_get_prev_residue(int chainnumber, int residuenumber) throw(FFFError) {
    //
    // The the previous residue
    //
    vector<int> prevres;
    if (residuenumber==0) {
            throw FFFError();
    }
    residuenumber=residuenumber-1;
    prevres.push_back(chainnumber);
    prevres.push_back(residuenumber);
    return prevres;
}

//
// ------------------------------------------------
//

vector<int> model_class::_get_next_residue(int chainnumber, int residuenumber) throw(FFFError) {
    //
    // Get the next residue
    //
    vector<int> nextres;
    if (residuenumber+1>static_cast<int>(_P.chains[chainnumber].residues.size()-1)) {
        throw FFFError();
    }
    nextres.push_back(chainnumber);
    nextres.push_back(residuenumber+1);
    return nextres;
}

//
// ------------------------------------------------
//


double model_class::getpsi(int chainnumber, int residuenumber) {
    //
    // Get the psi angle
    //
    return _P.chains[chainnumber].residues[residuenumber].psi;
}

//
// -----
//

double model_class::getphi(int chainnumber,int residuenumber) {
    //
    // Get the phi angle
    //
    return _P.chains[chainnumber].residues[residuenumber].phi;
}



//
// -----------------------------------------------------------
//

double model_class::calc_chi(int chainnumber,int residuenumber,int chinum) {
  //
  // Get the chi angle[in degrees] <chinum> of residue <residue>
  //
  // When you call this routines, then the chinumber starts with 1. 
  // One (1) is subtracted from the chinum in this routine.
  //
  // Check that the residue has that chinumber
  //
  if (chinum>getnumchi(chainnumber,residuenumber) || chinum==0) {
    cout << "This residue: c/r " << chainnumber << "/" << residuenumber << " has only " 
	 << getnumchi(chainnumber,residuenumber) << " Chi angles\n";
    cout << "You requested: " << chinum << endl;
    exit(0);
  }
  //
  // Find the amino acid definition
  //
  int defres=-1;
  defres=_get_AA_template(_P.chains[chainnumber].residues[residuenumber].name);
  //
  // OK, find the atoms and get the angle
  //
  vector<atom_class*> tors_atoms;
  string defname;
  for (int n=0;n<4;n++) {
    // Get the name of the atom from the definition of the torsion angles
    defname=AAsets[0].aadefs[defres].dihedral_atoms[(chinum-1)*4+n].name;
    try {
      tors_atoms.push_back(_get_atom(chainnumber,residuenumber,defname));
    }
    catch (FFFError) {
      //
      // Atom not found, so we return a large torsion angle
      //
      return 1000.0;
    }
  }
  //
  // Get the angle
  //
  return Dihedral(*tors_atoms[0],*tors_atoms[1],*tors_atoms[2],*tors_atoms[3]);
}

//
// ------------------------------------------------
//

void model_class::setchi(int chain,int residue,int chinum,double chiang) throw(exception) {
  //
  // Set the chi angle <chinum> of residue <residue> to angle[in degrees] <chiang>
  //
  // When you call this routines, then the chinumber starts with 1. 
  // One (1) is subtracted from the chinum in this routine.
  //
  double oldchi;
  oldchi=getchi(chain,residue,chinum);
  if (oldchi>180.0) {
    throw exception();
  }
  double difchi=chiang-oldchi;
  //if (DEBUG>0) {
  //printf("In setchi. Old chiangle is: %f, desired chiangle is: %f\n",oldchi,chiang);
  //printf ("\nRotating Chi angle # %d for chain %d residue #%d by %f degrees.\n",chinum,chain,residue,difchi);
    //}
  //
  // Find the amino acid definition
  //
  int defres=-1;
  defres=_get_AA_template(_P.chains[chain].residues[residue].name);
  //
  // Find the atoms and construct the vector
  //
  vector<atom_class*> tors_atoms;
  string defname;
  for (int n=0;n<4;n++) {
    // Get the name of the atom from the definition of the torsion angles
    defname=AAsets[0].aadefs[defres].dihedral_atoms[(chinum-1)*4+n].name;
    try {
      tors_atoms.push_back(_get_atom(chain,residue,defname));
    } 
    catch (FFFError) {
      printf ("Cannot set chiangle# %d for chain: %d residue: %d %s\n",
	      chain,chinum,residue,
	      _P.chains[chain].residues[residue].name.c_str());
    }
  }
  vector3d V;
  int basis=1;
  int point=2;
  V=*(tors_atoms[point])-*(tors_atoms[basis]);
  quatfit_class Q;
  Q.rotate_around_vector(V,difchi);
  //
  // Change the coordinates in the sidechain
  //
    //printf ("Chinum: %d\n",chinum);
  vector<atom_class*> moveatoms=_get_move_atoms(chain,residue,(*(tors_atoms[2])).name);
  vector<atom_class> transform_coords;
  //
  //
  int numpoints=0;
  for (unsigned int n=0;n<moveatoms.size();n++) {
    numpoints++;
    atom_class X("DUM",
    		 (*(moveatoms[n])).x-(*(tors_atoms[basis])).x,
    		 (*(moveatoms[n])).y-(*(tors_atoms[basis])).y,
    		 (*(moveatoms[n])).z-(*(tors_atoms[basis])).z);
    transform_coords.push_back(X);
  }
  //
  // Transform the coordinates
  //
  Q.transform(numpoints,transform_coords);
  //
  // Put the new coordinates in the array
  //
  int counter=0;
  for (int n=0;n<static_cast<int>(moveatoms.size());n++) {
    (*(moveatoms[n])).x=transform_coords[counter].x+(*(tors_atoms[basis])).x;
    (*(moveatoms[n])).y=transform_coords[counter].y+(*(tors_atoms[basis])).y;
    (*(moveatoms[n])).z=transform_coords[counter].z+(*(tors_atoms[basis])).z;
    counter++;
  }
  //
  // Put the new chiangle in the internal array
  //
  _all_chi_angles[chain][residue][chinum-1]=calc_chi(chain,residue,chinum);
  return;
}

//
// ------------------------------------------------
//

vector<atom_class *> model_class::_get_move_atoms(int chainnumber,int residuenumber,string atomname) {
  //
  // Given a residue number and an atom name, this function returns
  // all atoms that need to have their coordinates changed if the
  // atom given is the first one to be changed
  //
  // Get the amino acid definition for this residue
  //
  int aadef=_get_AA_template(_P.chains[chainnumber].residues[residuenumber].name);
  vector<atom_class*> moveatoms;
  //
  // Get the distance from the CA for the first atom to be moved
  //
  int first_distance=AAsets[0].aadefs[aadef].bonds_from_CA[strip(atomname)];
  //
  // Loop over all atoms and see if they are the same distance
  // or farther away from the CA
  //
    //printf ("Moving these atoms\n");
  for (unsigned int atom_local=0;atom_local<_P.chains[chainnumber].residues[residuenumber].atoms.size();atom_local++) {
    if (AAsets[0].aadefs[aadef].bonds_from_CA
	[strip(_P.chains[chainnumber].residues[residuenumber].atoms[atom_local].name)]
	>first_distance) {
        //_P.chains[chainnumber].residues[residuenumber].atoms[atom_local].print_detail();
      moveatoms.push_back(&_P.chains[chainnumber].residues[residuenumber].atoms[atom_local]);
    }
  }
    //printf ("----end----\n");
  return moveatoms;
}

//
// ------------------------------------------------
//

int model_class::getnumchi(int chainnumber,int residuenumber) {
  //
  // Get the number of chi angles for residue <residue>
  //
  int aadef=-1;
  if (residuenumber<0 || residuenumber>static_cast<int>(_P.chains[chainnumber].residues.size())) {
    printf ("Error: Residue number outside range: %d\n",residuenumber);
    exit(0);
  }
  string S(_P.chains[chainnumber].residues[residuenumber].name);
  aadef=_get_AA_template(S);
  if (aadef==-1) {
    return -1;
  }
  return AAsets[0].aadefs[aadef].numchi;
}

//
// ----------------------------------------------------------------------------
//
// Functions for handling the amino acid definitions and the rotamer library
//

void model_class::_read_aa_defs(const std::string aadef_dir) {
  //
  // This function reads the Amino Acid Definition file and the hydrogen definition file
  //
  // The aa_def file holds info on the standard conformations of all amino acids whereas
  // HYDROGENS.DAT holds information on all the conformations of hydrogens
  //
  //if (DEBUG>=3) {
  cout << "Reading AA.DAT\n";
  //}
  string aadef_file=aadef_dir+"/AA.DAT";
  ifstream file(aadef_file.c_str());
  if (!file) {
    cout << "ERROR: File AA.DAT was not found!\n";
    exit(0);
  }
  string line;
  vector<string> lines;
  while(getline(file,line)) {
    lines.push_back(line);
  }
  file.close();
  cout << "AA.DAT has been read\n";
  //
  // Read OTHER.DAT
  //
  aadef_file=aadef_dir+"/OTHER.DAT";
  ifstream file2(aadef_file.c_str());
  if (!file2) {
    printf ("ERROR: File %s was not found!\n",aadef_file.c_str());
    exit(0);
  }
  //string line;
  //vector<string> lines;
  while(getline(file2,line)) {
    lines.push_back(line);
  }
  file.close();
  cout << "OTHER.DAT has been read\n";
  //
  // Parse the aa def file
  //
  parse_aadefs X(lines);
  printf ("returned from parsing\n");
  AAsets.push_back(X);
  // Copy the information on amino acids to the array in _P
  for (unsigned int count=0;count<X.amino_acids.size();count++) {
    _P.amino_acids.push_back(X.amino_acids[count]);
  }
  // ----------------------------
  //
  // Read HYDROGENS.DAT
  //
  printf ("Reading HYDROGENS.DAT\n");
  string Hdef_file=aadef_dir+"/HYDROGENS.DAT";
  ifstream Hfile(Hdef_file.c_str());
  if (!Hfile) {
    cout << "ERROR: File HYDROGENS.DAT was not found!\n";
    exit(0);
  }
  string Hline;
  vector<string> Hlines;
  while(getline(Hfile,Hline)) {
    Hlines.push_back(Hline);
  }
  Hfile.close();
  printf ("Done reading HYDROGENS.DAT\n");
  //
  // Parse the conformations
  //
  _Hydrogen_defs = new hydrogens_defclass(Hlines);
  printf ("Parsed HYDROGENS.DAT\n");
  return;
 }

//
// ----------------------------------------------------
//

int model_class::_get_AA_template(string residuename) {
  //
  // Return the number in AAsets[0] that corresponds to the residue name
  //
  int respointer=-1;
  for (int tempres=0;tempres<static_cast<int>(AAsets[0].aadefs.size());tempres++) {
    if (residuename==AAsets[0].aadefs[tempres].name) {
      respointer=tempres;
      break;
    }
  }
  return respointer;
}

void model_class::update_bonds(atom_class& atom) {
    // Update the bonds for this atom
    //printf ("Getting bonds for: %s in residue %s %d %d\n",atom.name.c_str(),atom.inresiduename.c_str(),atom.inchain,atom.inresidue);
    int aadef=_get_AA_template(atom.inresiduename);
    // Fix for OXT
    string usename=atom.name;
    if (atom.name=="OXT" || atom.name=="O2" || atom.name=="O''") {
        usename="O";
    }
    // Get the atoms bonded to this atom as per the topology
if (aadef<0) {
return;
}
    vector<string> bonded_atoms=AAsets[0].aadefs[aadef].get_bonded_atoms(usename);
    for (unsigned int i=0;i<bonded_atoms.size();i++) {
        // Assign bonds 
        try {
            atom.covalent_bonds.push_back(_get_atom(atom.inchain,atom.inresidue,bonded_atoms[i]));
        } catch (FFFError) {
            //printf ("Could not find: %s in chain: %d, residue: %d\n",bonded_atoms[i].c_str(),atom.inchain,atom.inresidue);
        }
    }
    // Look for connections between amino acids
    if (atom.name=="N") {
        try {
            vector<int> prev_res=_get_prev_residue(atom.inchain,atom.inresidue);
            string Cname("C");
            atom.covalent_bonds.push_back(_get_atom(prev_res[0],prev_res[1],Cname));
            //printf ("Added N-C bond between residues %d %d - %d %d\n",prev_res[0],prev_res[1],atom.inchain,atom.inresidue);
        } catch (FFFError) {
            // No, this is the terminal
        }
    }
    if (atom.name=="C") {
        try {
            vector<int> next_res=_get_next_residue(atom.inchain,atom.inresidue);
            string Nname("N");
            atom.covalent_bonds.push_back(_get_atom(next_res[0],next_res[1],Nname));
            //printf ("Added C-N bond between residues %d %d - %d %d\n",next_res[0],next_res[1],atom.inchain,atom.inresidue);
        } catch (FFFError) {
            // No, this is the C-terminal
        }
        // Also check if there's an OXT or O2 in the residue
        try {
            atom.covalent_bonds.push_back(_get_atom(atom.inchain,atom.inresidue,"OXT"));
            atom.covalent_bonds.push_back(_get_atom(atom.inchain,atom.inresidue,"O2"));
            atom.covalent_bonds.push_back(_get_atom(atom.inchain,atom.inresidue,"O''"));
        } catch (FFFError) {
            // No
        }
    }
        
}

//
// ----
//

void model_class::update_bonds() {
    // 
    // Update all the bonding information
    //
    _P.update_all_atoms();
    //
    // first assign bonds based on topology
    //
    for (unsigned int atom=0;atom<_P.all_atoms.size();atom++) {
        (*(_P.all_atoms[atom])).covalent_bonds.resize(0);
        update_bonds(*(_P.all_atoms[atom]));
    }
    //
    // The code below should be updated to look for bonds in ligands, covalently attached groups etc
    // for now we just skip that.
    return;
    printf ("We never get here\n");
    exit(0);
    Boxes* BOXB=new Boxes(_P.all_atoms,2.0); // 2.0 A box
    for (unsigned int atom=0;atom<_P.all_atoms.size();atom++) {
        _P.all_atoms[atom]->covalent_bonds.resize(0);
        vector<atom_class *> close_atoms=(*BOXB).get_close_atoms(_P.all_atoms[atom]);
        for (unsigned int atom2=0;atom2<close_atoms.size();atom2++) {
            double distance=Dist(*_P.all_atoms[atom],*close_atoms[atom2]);
            if (distance<2.0) {
                _P.all_atoms[atom]->covalent_bonds.push_back(close_atoms[atom2]);
                //close_atoms[atom2]->covalent_bonds.push_back(all_atoms[atom]);
            }
            //printf ("%d %d: %5.2f\n",atom,atom2,distance);
        }
    }
    return;
}

//
// ----
//

void model_class::_save_chiangles(int chain,int residue) {
  //
  // Just push all the chiangles into _saved_angles
  //
  int num_saved=0;
  for (int chinum=1;chinum<=getnumchi(chain,residue);chinum++) {
    _saved_angles.push_back(getchi(chain,residue,chinum));
    num_saved=num_saved+1;
  }
  _saved_angles_count.push_back(num_saved);
  return;
}

void model_class::_restore_chiangles(int chain,int residue) {
  //
  // Get the last angles of _saved_angles
  //
  // This works now only for mutating to the same residue - i.e. only
  // for the repair function.
  // I still have to write code for e.g. Val -> Thr etc.
  //
  int records=static_cast<int>(_saved_angles_count.size());
  if (records<1) {
    printf ("Error in internal book-keeping. Number of saved chiangles is too small: %d\n",static_cast<int>(records));
    exit(0);
  }
  int num_saved=_saved_angles_count[records-1];
  int length=static_cast<int>(_saved_angles.size());
  int max_restore=getnumchi(chain,residue);
  //printf ("length of _saved_angles: %d\n",length);
  int chinum=1;
  for (int val=length-num_saved;val<length&&chinum<=max_restore;val++) {
    if (val<0) {
      printf ("something is seriously wrong with the number of saved chiangles\n");
      exit(0);
    }
    //printf ("val: %d length: %d angle: %f\n",val,length,_saved_angles[val]);
    if (_saved_angles[val]<=180.0) {
      setchi(chain,residue,chinum,_saved_angles[val]);
    }
    chinum++;
  }
  //printf ("resizing _saved_angles to :%d\n",length-num_saved);
  _saved_angles.resize(length-num_saved);
  _saved_angles_count.resize(records-1);
  return;
}


//
// ----
//

vector<double> model_class::get_energy(const std::string chain,const std::string residue) {
  // Get the energy for a single residue
  update_bonds();
  _P.update_BOXLJ(); // We must update the LJ boxes before calculating energies
  vector<int> resnum=_P.find_residue(chain,residue);
    vector<double> energy;
    energy.push_back(-9999999.9);
    energy.push_back(-9999999.9);
  if (resnum.size()==2) {
    energy=_ENERGY.get_external_energy(resnum[0],resnum[1]);
  }
  return energy;
}

vector<double> model_class::get_soup_energy() {
  // Get the energy for a single residue
  update_bonds();
  _P.update_BOXLJ(); // We must update the LJ boxes before calculating energies
  return _ENERGY.get_external_energy();
}

double model_class::get_accessibility(int chainnumber, int residuenumber) {
  // Get the accessibility for a single residue
  update_bonds();
  _P.update_BOX10A();
  double access=_ENERGY.get_accessibility(chainnumber,residuenumber);
  printf ("Access: %5.3f\n",access);
  return access;
}
