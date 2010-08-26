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
$Id: matrix_io.cpp,v 1.1 2005/11/15 09:56:07 chresten Exp $
*/

#include "matrix_io.h"


Matrix_IO::Matrix_IO(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);

  //lista = fopen("A.txt","a");
  //listb = fopen("b.txt","a");
}

Matrix_IO::Matrix_IO():Method(){}

Matrix_IO::~Matrix_IO(){}

int Matrix_IO::open(string filename)
{
  out = fopen(filename.c_str(),"w");
  int res = 0;
  if (out != NULL)
    res=1;

  return res;
}

void Matrix_IO::write_string_vector(string name, vector<string> v)
{

  fprintf(out,"%s = [\n",name.c_str());
  
  for(unsigned int i=0; i<v.size(); i++)
    {
      fprintf(out,"\"%s\"\n",v[i].c_str());
    }
  fprintf(out,"];\n\n\n");



}

void Matrix_IO::write_matrix(string name, la::matrix<double> m)
{
  fprintf(out,"%s = [\n",name.c_str());
  
  for(unsigned int i=0; i<m.size1(); i++)
    {
      for(unsigned int j=0; j<m.size2(); j++)
	fprintf(out,"%8.2f ",m(i,j));
      fprintf(out,"\n");
    }
  fprintf(out,"];\n\n\n");

}

void Matrix_IO::write_vector(string name, la::vector<double> v)
{
  fprintf(out,"%s = [\n",name.c_str());
  
  for(unsigned int i=0; i<v.size(); i++)
    fprintf(out,"%8.2f\n",v(i));
  fprintf(out,"];\n\n\n");

}



void Matrix_IO::write_matrix(string name, vector<vector<double> > m)
{
  fprintf(out,"%s = [\n",name.c_str());
  
  for(unsigned int i=0; i<m.size(); i++)
    {
      for(unsigned int j=0; j<m[i].size(); j++)
	fprintf(out,"%8.2f ",m[i][j]);
      fprintf(out,"\n");
    }
  fprintf(out,"];\n\n\n");

}

void Matrix_IO::write_vector(string name, vector<double> v)
{
  fprintf(out,"%s = [\n",name.c_str());
  
  for(unsigned int i=0; i<v.size(); i++)
    fprintf(out,"%8.2f\n",v[i]);
  fprintf(out,"];\n\n\n");

}




void Matrix_IO::update_lists(vector<string> names, la::matrix<double> to_a, la::vector<double> to_b)
{
  /*
  // print to A 
  for(unsigned int i=0; i<to_a.size1(); i++)
    {
      fprintf(lista,"%14s ",names[i].c_str());
      for(unsigned int j=0; j<to_a.size2(); j++)
	fprintf(lista,"%14f ",to_a(i,j));
      fprintf(lista,"\n");
    }
  
  // print to b 
  for(unsigned int j=0; j<to_b.size(); j++)
    {
      fprintf(listb,"%14s ",names[j].c_str());  
      fprintf(listb,"%14f ",to_b(j));
    }
  fprintf(listb,"\n");

  */
}

void Matrix_IO::close()
{
  fclose(out);
  //  fclose(lista);
  //  fclose(listb);
}

void Matrix_IO::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: matrix_io.cpp,v 1.1 2005/11/15 09:56:07 chresten Exp $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
