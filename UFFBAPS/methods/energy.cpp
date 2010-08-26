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
$Id: energy.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "energy.h"


Energy::Energy(FILE * reporter, vector<Soup *> proteins_in, vector<Soup *> ligands_in):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);


  /*** do A and B have equal lengths? ***/
  if(proteins.size() != ligands.size())
    {
      cout << "WARNING: An equal amount of protein and ligand structures must be given to the energy function!\n";
      return;
    }

  are_coefficients_set = false;

  /*** copy structures and read in setup ***/
  proteins = proteins_in;ligands = ligands_in;
  read_parameters();

  /*** set parameters ***/
  //  Generalparameters prepA(rep,proteins);
  //  Generalparameters prepB(rep,ligands);
  simple_parameters.prepare(proteins);
  simple_parameters.prepare(ligands);


  /*** set targets ***/
  generate_target_vector();

  test();

  if(mode == "train")
    train();


  /*** write out matrices for octave ***/
  Matrix_IO M(rep);
  M.open("binding.m");
  M.write_matrix("A", A);
  M.write_vector("b",b);
  M.update_lists(names,A,b);
  M.write_string_vector("Names",names);
  M.close();



  if(are_coefficients_set)
    stream_results();
}

Energy::Energy():Method(){}

Energy::~Energy(){}

void Energy::initiate()
{
  FILE * reporter = fopen("out.txt","w");

  // write out to out file
  writeDescription(reporter);
  are_coefficients_set = false;

  // read parameters from file
  read_parameters();

  return;
}

void Energy::generate_target_vector()
{
  Targets T(rep,dataset);

  for(vector<Soup*>::iterator s = proteins.begin(); s != proteins.end(); s++)
    {
      b.resize(b.size()+1,true);
      b(b.size()-1) = T.get_G((*s)->name.substr(0,4));
    } 
}


void Energy::test()
{

  for(unsigned int s=0; s<proteins.size(); s++)
    {
      cout<<s+1<<" / "<<proteins.size()<<" ";
      
      calculate_binding(proteins[s], ligands[s]);

      // force field components 
      cout<<proteins[s]->name<<":     "<<flush;
      cout<<" de:"<<de_component <<flush;
      cout<<" hb:"<<hb_component <<flush;
      cout<<" vdw:"<<vdw_component<<flush;
      cout<<" es:"<<es_component <<flush;
      cout<<" pe:"<<pe_component <<flush;
      cout<<" le:"<<le_component<<flush;
      cout<<" wb:"<<wb_component <<flush;
      cout<<endl;
	  
 

      //cout<<"Target: "<<b(count)<<" kJ/mol"<<endl;
      // vdW     ES     HB     HP     PE     LE     WE 
      // Insert in matrices
      A(s, 0) = vdw_component;
      A(s, 1) = es_component;
      A(s, 2) = hb_component;
      A(s, 3) = de_component;
      A(s, 4) = bb_component;
      A(s, 5) = pe_component;
      A(s, 6) = le_component;
      A(s, 7) = wb_component; 
      
      names.push_back(proteins[s]->name);
    }
}



void Energy::train()
{
  try{
    cout<<"A:\n"<<A<<endl;
    cout<<"b:\n"<<b<<endl;


    /*** linear least squares ***/
    la::matrix<double> AT = la::trans(A);
    la::matrix<double> ATA = la::prod(AT,A);
    la::vector<double> ATb = la::prod(AT,b);
    la::permutation_matrix<int> pm(ATA.size1());
    la::lu_factorize(ATA,pm);
    la::vector<double> c(ATb);
    lu_substitute(ATA, pm, c); 
    la::vector<double> est = la::prod(A,c);

    /*** set coeffficients  ***/
    set_coefficients(c);
    
    /*  cout<<"A:\n"<<print_matrix(A)<<endl;
	cout<<"b:\n"<<print_matrix(b)<<endl;
	cout<<"AT:\n"<<print_matrix(AT)<<endl;
	cout<<"ATA:\n"<<print_matrix(ATA)<<endl;
	cout<<"ATb:\n"<<print_matrix(ATb)<<endl;
	cout<<"ATA - after:\n"<<print_matrix(ATA)<<endl;
	cout<<"c:\n"<<print_matrix(c)<<endl;
	cout<<"est:\n"<<print_matrix(est)<<endl;
    */
    cout<<"TRAINING DONE"<<endl;
  }
  catch(...)
    {
      cout<<"WARNING: Training fail, will stil write matices to octave file"<<endl;
      
    }


  cout<<"Done writting to octave file"<<endl;

}

void Energy::set_coefficients(la::vector<double> c)
{
  coefficients = c;
  are_coefficients_set = true;
  return;
}

double Energy::calc_total_binding_energy(int i)
{
  la::vector<double> res = la::prod(coefficients, la::trans(la::project(A,la::range(i,i+1),la::range(0,no_components))));
  return res(0);
}

void Energy::stream_results()
{
  /*** write output header to stream ***/
  printf("All values below are in kJ/mol\n");
  printf("          |                  COMPONENTS                              |       DG total       \n");
  printf("----------+----------------------------------------------------------+-----------------------\n");
  printf("Structure |     vdW     ES     HB     DE     BB     PE     LE     WB |   Target   Prediction \n");  
  printf("----------+----------------------------------------------------------+-----------------------\n");
  
  for(unsigned int i=0; i<A.size1(); i++)
    {
      printf("%-10s|  %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f |  %6.3f       %6.3f \n",
	     names[i].c_str(), 
	     A(i,0),A(i,1),A(i,2),A(i,3),A(i,4),A(i,5),A(i,6),A(i,7),
	     b(i),calc_total_binding_energy(i));
    }
}


double Energy::calculate_binding(Soup * protein, Soup * ligand)
{

  simple_parameters.prepare(protein);
  simple_parameters.prepare(ligand);

  
  BM.assign_bonds(ligand);

  // force field components 

  de_apo = de.calculate(protein);
  de_holo = de.calculate(protein, ligand); //Desolvation_Factor must be calculated first, as desolvation_factor is needed in hb, vdw, es, wb
  de_component = de_holo - de_apo;

  wb.find_water_contacts(protein); //mark which acceptors and donors can be replaced by water
  //cout<<"  de:"<<de_component<<flush;
  
  hb_component  = hb.calculate(protein, ligand); //Hydrogen bonds must be calculated before es and vdw, as info on hydrogen bonds is needed in es and vdw 
  //cout<<"  hb:"<<hb_component<<flush;  

  vdw_component = vdw_line.calculate(protein, ligand, 10);
  //cout<<"  vw:"<<vdw_component<<flush;  

  es_component  = es.calculate(protein, ligand, 0);
  //cout<<"  es:"<<es_component<<flush;

  bb_component = 0.0;


  pe_component = 0.0;//pe.calculate(proteins[s], ligands[s]);
  //cout<<"  pe:"<<pe_component<<flush;  

  le_component = le.calculate_using_desolvation(protein, ligand);
  //  le_component = le.calculate(ligand);
  //cout<<"  le:"<<le_component<<flush;

  wb_component = wb.calculate(protein, ligand);
  //cout<<"  wb:"<<wb_component<<flush;
  //cout<<" wb:"<<r_wb <<flush;
  //cout<<endl;
  // Insert in matrices 
  components(0) = vdw_component;
  components(1) = es_component;
  components(2) = hb_component;
  components(3) = de_component;
  components(4) = bb_component;
  components(5) = pe_component;
  components(6) = le_component;
  components(7) = wb_component;


  if(are_coefficients_set)
    {
      calculate_result();
      return total_binding;
    }

  return 0.0;
}

void Energy::calculate_result()
{
  //  cout<<"size of components   "<<components.size()<<endl;
  //  cout<<"size of coefficients "<<coefficients.size()<<endl;
  total_binding = la::inner_prod(coefficients,components);
  return;
}

vector<double> Energy::get_components()
{
  vector<double> res;
  la::vector<double>::iterator it = components.begin();

  while(it != components.end())
    {
      res.push_back(*it);
      it++;
    }
  
  return res;

}

int Energy::get_no_components(){return no_components;}


void Energy::read_parameters()
{
  ifstream parameters;
  parameters.open("energy.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'energy.cfg' could not be found\n";
      exit(0);
    }
  
  printf("Setting parameters\n");

  string dummy;
  
  la::vector<double> coef;

  while (!parameters.eof())
    {
      parameters >> dummy;

      if(dummy == "MODE:")
	parameters >> mode;

      if(dummy == "DATASET:")
	parameters >> dataset;

      if(dummy == "COEFFICIENTS:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      coef.resize(coef.size()+1,true);
	      coef(coef.size()-1) = (double) atof(dummy.c_str());

	      parameters >> dummy;
	    }
	}

    } //end read file

  // set coefficients if given in parameter files
  if(coef.size()>0)
    set_coefficients(coef);

  //printing all configurations
  /*
  printf("# configuration file for Energy\n");
  printf("\n");
  printf("\n");

  printf("MODE:                 %s\n",mode.c_str());
  printf("\n");

  printf("DATASET:              %s\n",dataset.c_str());
  printf("\n");

  printf("COEFFICIENTS:         ");
  for(unsigned int i=0; i<coef.size(); i++)
    printf("%f ",coef[i]);
  printf("end\n");

  printf("\n");
  printf("\n");
  printf("Version: $Id: energy.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");
  */


  no_components = 8;
  A.resize(proteins.size(),no_components);
  components.resize(no_components);
  


}



void Energy::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: energy.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("Calculates energy terms");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());

  rep = reporter;
}
