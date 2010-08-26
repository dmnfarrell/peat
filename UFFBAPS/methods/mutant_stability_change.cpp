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
$Id: mutant_stability_change.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "mutant_stability_change.h"

MSC::MSC(FILE * reporter, FILE * results, vector<Soup *> wts_in, vector<Soup *> mutants_in):Method(reporter)
{

  /**********************************************************************
                      Mutant stability change

Input: Proberly protonated wild type and mutant structures


  **********************************************************************/

  wts = wts_in;
  mts = mutants_in;

  // write out to out file
  writeDescription(reporter);

  coefficients_are_set = 0;
  // does wts and mutants have equal lengths?
  if(wts.size() != mts.size())
    {
      cout << "ERROR: An equal amount of wild type and mutant structures must be given!\n";
      return;
    }

  // read parameters from file
  read_parameters();

  // set parameters
  if(include_ff_components)
    {
      //correct_names.go(wts);
      //correct_names.go(mts);
      
      //prepin.go(wts);
      //prepin.go(mts);
    
      simple_parameters.prepare(wts);
      simple_parameters.prepare(mts);
    }

  test();

  if(mode == "train")
    train();



  /*** write out matrices for octave ***/
  Matrix_IO M(rep);
  M.open("stability.m");
  M.write_matrix("Astability", A);
  M.write_vector("bstability",b);
  M.write_string_vector("Names",names);
  M.close();
  cout<<"Done writting to octave file"<<endl;





  if(coefficients_are_set)
    stream_results();
}

MSC::MSC():Method(){} // Constructor for scripting



MSC::~MSC(){}


void MSC::initiate()
{
  FILE * reporter = fopen("out.txt","w");
  //resultsFile = fopen(res.c_str(),"w");


  // write out to out file
  writeDescription(reporter);
  coefficients_are_set = 0;

  // read parameters from file
  read_parameters();


}

void MSC::test()
{

  vector<Soup*>::iterator 
    wt = wts.begin(), 
    mt = mts.begin();
  
  vector<vector<float> > wt_distances, mut_distances;
  cout<<setprecision(2)<<setw(8)<<fixed;
  while(wt != wts.end())
    {
      cout<<"******* Now doing no "<<count+1<<" of "<<wts.size()<<": "<<(*mt)->name<<endl;

      if(include_ff_components)
	{

	  A.resize(A.size1()+1,A.size2(),true);
	  /*** Insert mutant values in matrices ***/
	  calculate_stability(*mt);
	  for(int i=0; i<components.size();i++)
	    A(count,i) = components(i);

	  /*** Insert wild type values in matrices ***/
	  calculate_stability(*wt);
	  for(int i=0; i<components.size();i++)
	      A(count,i) = A(count,i) - components(i);


	  wt_comps.resize(wt_comps.size1()+1,wt_comps.size2(),true);
	  mut_comps.resize(mut_comps.size1()+1,mut_comps.size2(),true);
	  
	  wt_comps(count, 0) = 0.0;
	  wt_comps(count, 1) = 0.0;
	  wt_comps(count, 2) = 0.0;
	  wt_comps(count, 3) = 0.0;
	  wt_comps(count, 4) = 0.0;
	  
	  mut_comps(count, 0) = 0.0;
	  mut_comps(count, 1) = 0.0;
	  mut_comps(count, 2) = 0.0;
	  mut_comps(count, 3) = 0.0;
	  mut_comps(count, 4) = 0.0;
	  
	}
      if(include_aa_components)
	{
	  /*** insert unfolded terms ***/
	  for(int i=0;i<20;i++)
	    A(count,no_components-20+i) = 0.0;

	  A(count,no_components-21+aa_map[(*mt)->name[5]]) = -1.0;
	  A(count,no_components-21+aa_map[(*mt)->name[(*mt)->name.size()-1]]) =  1.0;
	}

      names.push_back((*mt)->name);
      
      count++;
      
      /*** iterate iterators ***/
      wt++;
      mt++;
    }
}

void MSC::train()
{

  /*** linear least squares ***/
  try{
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
    
    cout<<"A:\n"<<print_matrix(A)<<endl;
    cout<<"b:\n"<<print_matrix(b)<<endl;
    cout<<"AT:\n"<<print_matrix(AT)<<endl;
    cout<<"ATA:\n"<<print_matrix(ATA)<<endl;
    cout<<"ATb:\n"<<print_matrix(ATb)<<endl;
    cout<<"ATA - after:\n"<<print_matrix(ATA)<<endl;
    cout<<"c:\n"<<print_matrix(c)<<endl;
    cout<<"est:\n"<<print_matrix(est)<<endl;
    
    cout<<"TRAINING DONE"<<endl;
    
  }
  catch(...)
    {
      cout<<"WARNING: Training failed! Writting matrices to octave file\n";
    }
  
  
  return;

}

void MSC::set_coefficients(la::vector<double> in)
{
  coefficients_are_set = 1;
  coefficients = in;
  return;
}

string MSC::print_matrix(la::vector<double> v)
{
  ostringstream res;
  for(unsigned int i=0; i<v.size(); i++)
    {
      res<<" ";
      res<<v(i);
    }
  res<<"\n";
  return res.str();

}

string MSC::print_matrix(la::matrix<double> m)
{
  ostringstream res;
  for(unsigned int i=0; i<m.size1(); i++)
    {
      for(unsigned int j=0; j<m.size2(); j++)
	{
	  res<<" ";
	  res<<m(i,j);
	}
      res<<"\n";
    }
  return res.str();

}

double MSC::calc_total_DDG(int i)
{
  float res = 0;
  if (coefficients.size()>0)
    {
      la::vector<double> resv = la::prod(coefficients, la::trans(la::project(A,la::range(i,i+1),la::range(0,no_components))));
      res = resv(0);
    }

  return res;
}

void MSC::stream_results()
{

  /*** write output header to stream ***/
  printf("All values below are in kJ/mol\n");
  printf("          |                       WILDTYPE                   |                        MUTANT                    |            DG components            |       DG total       \n");
  printf("----------+--------------------------------------------------+--------------------------------------------------+-------------------------------------+----------------------\n");
  printf("Structure |        vdW        ES        HB       DE       BB |        vdW        ES        HB       DE       BB |     vdW     ES     HB     DE     BB |  Target   Prediction \n");
  printf("----------+--------------------------------------------------+--------------------------------------------------+-------------------------------------+----------------------\n");

  fprintf(output,"name,wt_vdW,wt_ES,wt_HB,wt_DE,wt_BB,mut_vdW,mut_ES,mut_HB,mut_DE,mut_BB,D_vdW,D_ES,D_HB,D_DE,D_BB,Target,Prediction\n");


  for(unsigned int i=0; i<A.size1(); i++)
    {
      
      printf("%-10s|  %6.2e %6.2e %6.2e %6.2e %6.2e |  %6.2e %6.2e %6.2e %6.2e %6.2e |  %6.2f %6.2f %6.2f %6.2f %6.2f |  %6.3f       %6.3f \n",
	     names[i].c_str(), 
	     wt_comps(i,0),wt_comps(i,1),wt_comps(i,2),wt_comps(i,3),wt_comps(i,4),
	     mut_comps(i,0),mut_comps(i,1),mut_comps(i,2),mut_comps(i,3),mut_comps(i,4),
	     A(i,0)*coefficients(0),A(i,1)*coefficients(1),A(i,2)*coefficients(2),A(i,3)*coefficients(3),A(i,4)*coefficients(4),
	     0.0,calc_total_DDG(i));
      
      fprintf(output,"%-10s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %%8.4f \n",
	      names[i].c_str(), 
	      wt_comps(i,0),wt_comps(i,1),wt_comps(i,2),wt_comps(i,3),wt_comps(i,4),
	      mut_comps(i,0),mut_comps(i,1),mut_comps(i,2),mut_comps(i,3),mut_comps(i,4),
	      A(i,0),A(i,1),A(i,2),A(i,3),A(i,4),
	      0.0,calc_total_DDG(i));
      
   }


  fclose(output);

}


double MSC::calculate_stability(Soup * soup)
{
  if(include_ff_components)
    simple_parameters.prepare(soup);
  
  cout<<soup->name<<flush;
  
  de_component  = de.calculate(soup);//Desolvation_Factor must be calculated first, as desolvation_factor is needed in hb, vdw, es, wb
    cout<<" de:"<<de_component <<flush;

  hb_component  = hb.calculate(soup); //Hydrogen bonds must be calculated before es and vdw, as info on hydrogen bonds is needed in es and vdw 
    cout<<" hb:"<<hb_component <<flush;

  vdw_component = vdw_line.calculate(soup, 10);
    cout<<" vdw:"<<vdw_component<<flush;

  es_component  = es.calculate(soup, 0);
    cout<<" es:"<<es_component <<flush;

  bb_component = bb.calculate(soup);
    cout<<" bb:"<<bb_component <<flush;

  pe_component  = 0.0;//pe.calculate(soup);
    cout<<" pe:"<<pe_component <<flush;	  

  le_component  = 0.0;//le.calculate(soup);
    cout<<" le:"<<le_component <<flush;
  
  wb_component = wb.calculate(soup);
    cout<<endl;
  
  cout << "calc is done" <<flush;
  // Insert in matrices 
  components(0) = vdw_component;
  components(1) = es_component;
  components(2) = hb_component;
  components(3) = de_component;
  components(4) = bb_component;
  components(5) = pe_component;
  components(6) = le_component;
  components(7) = wb_component;

  /*
  A.resize(A.size1()+1,A.size2(),true);
  for(int i=0; i<components.size();i++)
    A(A.size1()-1,i) = components(i);

  names.push_back(soup->name);
      
 
  count++;
  */

  if(coefficients_are_set)
    {
      calculate_result();
      return total_stability;
    }

  return 0.0;
}



void MSC::calculate_result()
{
  total_stability = la::inner_prod(coefficients,components);
  return;
}

vector<double> MSC::get_components()
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

int MSC::get_no_components(){return no_components;}



void MSC::read_parameters()
{
  ifstream parameters;
  parameters.open("msc.cfg");

  if(! parameters.is_open())
    {
      cout << "Error:File 'msc.cfg' could not be found\n";
      exit(0);
    }
  

  string dummy;
  la::vector<double> coef;


  while (!parameters.eof())
    {

      parameters >> dummy;

      if(dummy == "MODE:")
	parameters >> mode;

      if(dummy == "INCLUDE_FF:")
	parameters >> include_ff_components;

      if(dummy == "INCLUDE_AA:")
	parameters >> include_aa_components;

      if(dummy == "TARGETS:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      b.resize(b.size()+1,true);
	      b(b.size()-1) = (double) atof(dummy.c_str());
	  
	      parameters >> dummy;
	    }
	}

      if(dummy == "COEFFICIENTS:")
	{
	  parameters >> dummy;
	  while(dummy != "end")
	    {
	      coef.resize(coef.size()+1,true);
	      coef(coef.size()-1) = (double) atof(dummy.c_str());
	      //cout<<"COEF: "<<dummy.c_str()<<endl;
	      parameters >> dummy;
	    }
	}


    } //end read file

  // set coefficients if given in parameter files
  if(coef.size()>0)
    set_coefficients(coef);

  //printing all configurations
  /*
  printf("# configuration file for Mutant Stability Change\n");
  printf("\n");
  printf("\n");

  printf("MODE:                 %s\n",mode.c_str());
  printf("\n");

  printf("INCLUDE_FF:                 %d\n",include_ff_components);
  printf("\n");

  printf("INCLUDE_AA:                 %d\n",include_aa_components);
  printf("\n");
  
  printf("COEFFICIENTS:               ");
  for(int i=0; i<coef.size(); i++)
    printf("%f ",coef[i]);
  printf("end\n");
  printf("\n");
  
  printf("TARGETS:              ");
  for(int i=0; i<b.size(); i++)
    printf("%f ",b(i));
  printf("end\n");

  printf("\n");
  printf("\n");
  printf("Version: $Id: mutant_stability_change.cpp 6326 2010-08-26 16:08:21Z nielsen $\n");
  printf("\n");
  */




  aa_map['A'] = 1;
  aa_map['R'] = 2;
  aa_map['N'] = 3;
  aa_map['D'] = 4;
  aa_map['C'] = 5;
  aa_map['E'] = 6;
  aa_map['Q'] = 7;
  aa_map['G'] = 8;
  aa_map['H'] = 9;
  aa_map['I'] = 10;
  aa_map['L'] = 11;
  aa_map['K'] = 12;
  aa_map['M'] = 13;
  aa_map['F'] = 14;
  aa_map['P'] = 15;
  aa_map['S'] = 16;
  aa_map['T'] = 17;
  aa_map['W'] = 18;
  aa_map['Y'] = 19;
  aa_map['V'] = 20;

  // find number of components
  no_components = 0;
  if (include_ff_components == 1)
    no_components += 8;
  if (include_aa_components == 1)
    no_components += 20;

  // resizing matrices
  A.resize(0,no_components);
  wt_comps.resize(0,no_components);
  mut_comps.resize(0,no_components);

  components.resize(no_components);
 
  count=0;
  

}

void MSC::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: mutant_stability_change.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("Calculates MSC terms");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());

  rep = reporter;
  output = fopen( "output.txt", "w" );
}
