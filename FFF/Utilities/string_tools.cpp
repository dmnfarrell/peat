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
#include "string_tools.h"
#include <string>
#include <iostream>
//
//
using namespace std;
//
// --------------------------------------
//
vector<string> split(string s) {
  //
  // Split the string s into sub-strings.
  // Separator is ' '
  //
  vector<string> sub_strings;
  for (int i = 0; i< static_cast<int>(s.size()); i++) {
    if (strncmp(&s[i]," ",1)) {
      string sub_string(" ");
      while (strncmp(&s[i]," ",1) && strncmp(&s[i],"\n",1)  && i<static_cast<int>(s.size())) {
	sub_string=sub_string+s[i];
	i++;
      }
      sub_strings.push_back(strip(sub_string));
      //printf ("substring: %s\n",strip(sub_string).c_str());
    }
  }
  return sub_strings;
}

//
// ---------------
//

vector<string> split(string str,string sep) {
  //
  // Split the string s into sub-strings.
  // Separator is sep
  //
  string::size_type lastPos = 0;
  // Find first instance of sep
  string::size_type pos = str.find_first_of(sep, lastPos);
  vector<string> sub_strings;
  while (string::npos!=pos || string::npos!=lastPos) {
    sub_strings.push_back(str.substr(lastPos,pos-lastPos));
    //printf ("word: %s\n",str.substr(lastPos,pos-lastPos).c_str());
    lastPos=str.find_first_not_of(sep,pos);
    pos=str.find_first_of(sep,lastPos);
  }
  //printf ("done with split\n");
  return sub_strings;
}

//
// ----
//

string join(vector<string> words,string joiner) {
    // Join the words
    string final;
    for (unsigned int jj=0;jj<words.size();jj++) {
        final=final+words[jj];
        if (jj<words.size()-1) {
            final=final+joiner;
        }
    }
    return final;
}
    

//
// -------------------------------------------
//
string strip(string s) {
  if (s.length()==0) {
    return s;
  }
  int b = s.find_first_not_of(" \t");
  int e = s.find_last_not_of(" \t");
  if (b==-1) {
    return "";
  }
  string sn;
  sn=s.substr(b,e-b+1);
  return sn;
}
//
// --------------------------------------------
//
vector<int> getints(string s) {
  //
  // Return all integers from <s> in a vector of ints
  //
  vector<int> numbers;
  char buffer[120];
  for (int i = 0; i< static_cast<int>(s.size()); i++) {
    if (strncmp(&s[i]," ",1)) {
      string number("   ");
      while (strncmp(&s[i]," ",1) && i<static_cast<int>(s.size())) {
	number=number+s[i];
	i++;
      }
      strcpy(buffer,"                 ");
      number.copy(buffer,number.length());
      numbers.push_back(atoi(buffer));   
    }
  }
  return numbers;
}
//
// ---------------------------------
//
char * s_to_char(string s) {
  //
  // Convert a string to a char *
  //
  s=strip(s);
  char* name = new char[s.size()];
  s.copy(name,s.size());
  name[s.size()]=0;
  return name;
}
