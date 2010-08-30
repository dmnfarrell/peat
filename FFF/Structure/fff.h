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
#ifndef FFF_H
#define FFF_H
#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <unistd.h>

//
// Local include files
//
#include "chain_class.h"
#include "residue_class.h"
#include "atom_class.h"
#include "TitGroup_class.h"
#include "ostools.h"
#include "Boxes.h"


using namespace std;

class FFF {
public:
    FFF() {
        //
        // Intialize environment
        //
        _pH=7.0;
    };
  void read_pdb(const std::string pdbfilename);
  void write_pdb(const std::string pdbfilename);
  void writepdb(const std::string pdbfilename);
    void parse_lines(const std::vector<std::string> lines);
    void remove_atoms_with_tag(const std::string);
    void remove_atom(atom_class& atom);
  void soup_stat(); // Prints a summary of the soup
  void check_soup(); // Check the integrety of the soup
    void update_all_atoms(); //Update the all_atoms array
        void update_BOXLJ(); //Update the boxes for LJ and Steric Clash
    //
    // 
    //
    void set_pH(double pH) {
        _pH=pH;
        return;
    }
    //
    // Utility functions
    //
    bool are_bonded(const atom_class& atom1, const atom_class& atom2);
    bool are_1_3_bonded(const atom_class& atom1,const atom_class atom2);
    bool are_1_4_bonded(const atom_class& atom1,const atom_class atom2);
    //bool is_acceptor(const atom_class& atom);
    //bool is_donor(const atom_class& atom);
    bool could_be_Hbonded(const atom_class& atom1, const atom_class& atom2);
    vector<int> find_residue(const std::string chainname,const std::string pdbresnumber);
  //
  // variables
  //
  vector<chain_class> chains;
  string pdbname;
  int nresidues;
  int natoms;
  int atomcounter, rescounter;
  string last_resnum;
  vector<atom_class *> all_atoms;
    vector <TitGroup> _titratable_groups;
    //Boxes* BOXB; 
    Boxes* BOXLJ;
    //
    // Environment
    //
    double _pH;
};
#endif
