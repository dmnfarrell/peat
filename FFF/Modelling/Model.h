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
#ifndef MODEL_H
#define MODEL_H
#include <vector>
#include <string>
#include <iostream>
#include <exception>
//
// Local header files
//
#include "quatfit.h"
#include "parse_aadef.h"
#include "assign.h"
#include "fff.h"
#include "residue_class.h"
#include "Rotamers.h"
#include "energy.h"
#include "parse_hydrogens.h"

//
// Class for throwing errors
//
class FFFError {};

// Forward declaration

//
// Class definition
//
class model_class {
 public:
  //
  // Public functions
  //
  // Constructor
  //
    model_class(FFF &P,Rotamer_class &ROT,const std::string aadef_dir): _P(P),_ROT(ROT),_ENERGY(P) {
    //
    // Read AAdefinition file
    //
    printf ("Instantiating model_class\n\n");
    _read_aa_defs(aadef_dir);
    //
    // Read the charge and radius files, names are derived from aadef_file
    //
    printf ("reading radius and charge files\n");
    _crgrad = new crgrad_class(aadef_dir);
    //
    // Calculate all chiangles
    //
    _calc_all_angles();
    //
    // Update all bonds
    //
    update_bonds();
    //
    // set standard flags
    //
    _use_hydrogens=false;
    printf ("Number of residue types in rotamer library: %d\n",static_cast<int>(_ROT._names.size()));
    cout << "Model_class instantiated\n" << endl;
  };
  //
  // Destructor
  //
  ~model_class() {
    cout << "Destructing model_class\n";
  }
  //
  // Normal functions
  //
  vector<int> get_chain_and_residue(const std::string resid) throw(exception);
  vector<double> Mutate(const std::string resid, const std::string newresidue,int mode,double max_clash);
  vector<double> Mutate_2(int chainnumber,int resnumber,const std::string newresidue,int mode,double max_clash);
  vector<atom_class> mutate_special(int chainnumber,int resnumber,const std::string newresidue);
  void undo_mutation();
  bool simple_mutate(int chainnumber,int resnumber,string oldres,string newres);
  void setchi(int chainnumber,int residuenumber,int chinum,double chiang) throw(exception);
  double getchi(int chainnumber,int residuenumber,int chinum);
  int getnumchi(int chainnumber,int residuenumber);
  //int get_atom(int chainnumber,int residuenumber,string atomname) ;
  atom_class* _get_atom(int chainnumber,int residuenumber,string atomname) throw(FFFError);
  //
  // Global modelling functions
  //
  double resolve_clashes(); // Resolves all clashes. Does not touch clash-free residues
  //
  // Backbone functions
  //
  //void setphi(int chainnumber,int residuenumber,double angle);
  vector<int> _get_prev_residue(int chainnumber,int residuenumber) throw(FFFError);
  vector<int> _get_next_residue(int chainnumber,int residuenumber) throw(FFFError);
  double getphi(int chainnumber,int residuenumber);
  double getpsi(int chainnumber,int residuenumber);
  void calc_phi(int chain,int residue);
  void calc_psi(int chain,int residue);
  double calc_chi(int chainnumber,int residuenumber,int chinum);
  //
  // Bonding information
  //
  void update_bonds(); //Update the boxes for bonding
  void update_bonds(atom_class& atom);
  
  //double getomega(int chainnumber,int residuenumber);
  void repair_all();
  void find_bad_residues();
  vector<string> find_missing_atoms(residue_class& residue);
  vector<residue_class> unknown_residues;
  //
  // Hydrogens
  //
  void build_hydrogens();
  void delete_all_hydrogens();
  void delete_hydrogens(residue_class& residue);
  // 
  // Energy utility functions
  //
  vector<double> get_energy(const std::string chain,const std::string residue);
  vector<double> get_soup_energy();
  double get_accessibility(const std::string resid);
  //
  // Parameters
  //
  void assign_all_charges_and_radii();
  //
  // Public variables
  //
  FFF& _P;
  crgrad_class* _crgrad;
 private:
    //
    // Private functions
    //
    void _mutate(int chainnumber,int residuenumber,string newresidue);   // Mutate a single residue
    vector<double> _optimise1res_exhaustive(int chainnumber,int resnumber,double chistep,double max_clash); // Optimise 1 residue using brute force search
    vector<double> opt_chi(int chainnumber, int resnumber, int chinum, double chistep,double max_clash); // Optimise this chiangle and any subsequent ones
    vector<double> _optimise1res(int chainnumber,int resnumber,double max_clash); // Optimise 1 residue using the rotamer library 
    void _save_chiangles(int chainnumber,int residuenumber);
    void _restore_chiangles(int chainnumber,int residuenumber);
    vector<double> _saved_angles;
    vector<int> _saved_angles_count;

    void _calc_all_angles(); // Calculates all angles and stores them locally
 
  
    vector<atom_class *> _get_move_atoms(int chainnumber,int residuenumber,string atomname);
    vector<int> _get_first_last(int chainnumber,int residuenumber);
    // Calculates all angles for a single residue
    vector<double> _calc_residue_chiangles(int chainnumber,int residuenumber);   
    // 
    //

    //
    // For reading the rotamer library
    //
    void _read_aa_defs(const std::string);
    int _get_AA_template(string residuename);
    //
    // Variables for parameters
    //
    vector<parse_aadefs> AAsets;
    vector<vector <vector <double> > > _all_chi_angles; //Holds all chiangles for a protein
    hydrogens_defclass* _Hydrogen_defs; // Holds information on how to change protonation states by removing/adding hydrogens
    //
    // Variables for mutate operations
    //
    vector<string> delete_atoms;
    vector<string> add_atoms;
 
  //  
  // Private variables
  //
  Rotamer_class& _ROT;
  vector<vector<atom_class> > _org_atoms;
  energy_class _ENERGY;
  int iteration_count;

  // Flag controlling hydrogens
public:
  bool _use_hydrogens;
};



#endif
