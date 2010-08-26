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
#ifndef SOUPMANIP_H
#define SOUPMANIP_H
#include <stdlib.h>
#include <string>
#include <vector>
#include "soup.h"
#include "proteinchain.h"
#include "water.h"
#include "molecule.h"
#include "atom.h"

using namespace std;
                     
vector<SoupObject *> convert_soup_to_objects(Soup * soup);
vector<SoupObject *> convert_soup_to_objects(Soup * soup, bool incl_protein, bool incl_molecules, bool incl_water);
vector<SoupObject *> convert_soups_to_objects(vector<Soup *> soups);
vector<Atom *> convert_objects_to_atoms(vector<SoupObject *> objects);

vector<Residue> cut_residues_from_protein_chains_using_residue_number(vector<ProteinChain*> pc, int number);
vector<Atom*> get_atoms_from_protein_chains(vector<ProteinChain*> pcs);

vector<float> setup_boxes (vector<Atom*> atoms, float box_size);



#endif
