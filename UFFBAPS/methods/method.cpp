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

#include "method.h"


Method::Method()
{
  //to avoid to have to pass the reporter along all the time...
}

Method::Method(FILE * reporter)
{
  version =string("$Id: method.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
    string(" ");

  rep=reporter;
  fprintf(rep,"-------------------------------------------------------------------------------\n");	
  fprintf(rep,"Method:");	
}
Method::~Method()
{
  //fprintf(rep,"-------------------------------------------------------------------------------\n");	
}


void Method::setReporter(FILE * reporter)
{
  rep=reporter;
}

string Method::getDescription(){return description;}
string Method::getVersion(){return version;}
