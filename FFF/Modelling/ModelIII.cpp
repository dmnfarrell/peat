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

void model_class::repair_all() {
  //
  // First strip all hydrogens
  //
  delete_all_hydrogens();
  //
  // Look for missing heavy atoms
  //
  printf ("==========================================\n");
  printf ("Checking that all residues are complete\n");
  find_bad_residues();
  printf ("done\n\n");
}

//
// ---------------------------------------
//

void model_class::delete_all_hydrogens() {
  //
  // Loop over all residue
  //
  for (unsigned int chain=0;chain<_P.chains.size();chain++) {
    for (unsigned int residue=0;residue<_P.chains[chain].residues.size();residue++) {
      delete_hydrogens(_P.chains[chain].residues[residue]);
    }
  }
  return;
}

//
// ---------------------------------------
//

void model_class::delete_hydrogens(residue_class& residue) {
  //
  // Loop over all atoms and see if they are hydrogens
  //
  for (unsigned int atom=0;atom<residue.atoms.size();atom++) {
    if (residue.atoms[atom].is_hydrogen()) {
      //printf ("Deleting: %s residue: %d %s\n",residue.atoms[atom].name.c_str(),residue.number,residue.name.c_str());
      residue.atoms[atom].name="DELETE";
    }
  }
  residue.update();
}

//
// ---------------
//

void model_class::build_hydrogens() {
    //
    // Build all hydrogens and set the flags for using hydrogens from now on
    //
    _use_hydrogens=true;
    for (unsigned int chain=0;chain<_P.chains.size();chain++) {
        for (unsigned int residue=0;residue<_P.chains[chain].residues.size();residue++) {
	  Mutate_2(chain,residue,_P.chains[chain].residues[residue].name,1,10.0);
        }
    }
    return;
}

void model_class::assign_all_charges_and_radii() {
    //
    // Assign charges and radii according to the current protonation state
    //
    for (unsigned int chain=0;chain<_P.chains.size();chain++) {
        for (unsigned int residue=0;residue<_P.chains[chain].residues.size();residue++) {
            _crgrad->assign_all(_P.chains[chain].residues[residue]);
        }
    }
    return;
}
    

//
// ---------------------------------------
//

void model_class::find_bad_residues() {
  //
  // Loop through all residues and see if they have all the
  // heavy atoms that they should have
  //
  vector<residue_class*> residue_miss;
  vector<vector <string> > atom_miss;
  vector<string> temp;
  for (unsigned int chain=0;chain<_P.chains.size();chain++) {
    for (unsigned int residue=0;residue<_P.chains[chain].residues.size();residue++) {
      temp=find_missing_atoms(_P.chains[chain].residues[residue]);
      if (temp.size()>0) {
	atom_miss.push_back(temp);
	residue_miss.push_back(&_P.chains[chain].residues[residue]);
      }
    }
  }
  //
  // List the atoms that should be repaired
  //
  printf ("\nI found missing atoms in %d residues:\n\n",static_cast<int>(residue_miss.size()));
  for (unsigned int x=0;x<residue_miss.size();x++) {
    printf ("%4d %3s: %3d atoms: ",static_cast<int>((*(residue_miss[x])).number),(*(residue_miss[x])).name.c_str(),
	    static_cast<int>(atom_miss[x].size()));
    for (unsigned int at=0;at<atom_miss[x].size();at++) {
      printf ("%-4s",strip(atom_miss[x][at]).c_str());
    }
    printf ("\n");
  }
  printf ("----------------------------------------------\n");
  printf ("\nNow starting to repair those residues\n");
  for (unsigned int x=0;x<residue_miss.size();x++) {
    printf ("Repairing: %4d %3s.....",(*(residue_miss[x])).number,(*(residue_miss[x])).name.c_str());
    cout << flush;
    //
    // Simply mutate the residue to itself while preserving the chiangles
    //
      // This does not work in a lot of cases. Especially when it is just a single residue missing in a long chain
      //
      //
    Mutate_2( (*(residue_miss[x])).inchain,(*(residue_miss[x])).number,(*(residue_miss[x])).name,1,10.0);
    printf ("done\n");
  }
  return;
}

//
// ---------------------------------------
//

vector<string> model_class::find_missing_atoms(residue_class& residue) {
  //
  // Check that this residue has all the heavy atoms it should have
  //
  //
  // Find the template
  //
  int respointer=-1;
  vector<string> missing_atoms;
  for (unsigned int tempres=0;tempres<AAsets[0].aadefs.size();tempres++) {
    if (residue.name==AAsets[0].aadefs[tempres].name) {
      respointer=tempres;
      break;
    }
  }
  if (respointer==-1) {
    //
    // Residue not found, this means that we have a ligand?
    // Put this residue in the unknown residues vector and return
    //
    unknown_residues.push_back(residue);
    printf ("Residue without topology %s\n",residue.name.c_str());
    return missing_atoms;
  }
  //
  // Are all the atoms from the definition here?
  //
  for (unsigned int defatom=0;defatom<AAsets[0].aadefs[respointer].atoms.size();defatom++) {
    if (!AAsets[0].aadefs[respointer].atoms[defatom].is_hydrogen()) {
      bool found=false;
      for (unsigned int atom=0;atom<residue.atoms.size();atom++) { 
	if (strip(AAsets[0].aadefs[respointer].atoms[defatom].name)==
	    strip(residue.atoms[atom].name)) {
	  found=true;
	  break;
	}
      }
      if (!found) {
	//printf ("In residue: %d %s - Missing ",residue.number,residue.name.c_str());
	//AAsets[0].aadefs[respointer].atoms[defatom].print();
	missing_atoms.push_back(AAsets[0].aadefs[respointer].atoms[defatom].name);
      } 
    }
  }
  //
  // Are all the atoms in the residue found in the definition
  //
  for (unsigned int atom=0;atom<residue.atoms.size();atom++) {
    bool found=false;
    for (unsigned int defatom=0;defatom<AAsets[0].aadefs[respointer].atoms.size();defatom++) {
      if (strip(AAsets[0].aadefs[respointer].atoms[defatom].name)==
	  strip(residue.atoms[atom].name)) {
	found=true;
	break;
      }
    }
    if (!found) {
      //
      // Catch exception: C-terminal oxygen
      //
      if (residue.atoms[atom].name!="OXT") {
	printf ("Atom in residue is not in definition: ");
	residue.atoms[atom].print();
	//
	// This means that we have to recalculate charges for the entire residue
	//
	unknown_residues.push_back(residue);
      }
    }
  }
  return missing_atoms;
}
