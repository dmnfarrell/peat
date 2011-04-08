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

#include "parse_hydrogens.h"

// 
// Parse_hydrogens_class functions
//
hydrogens_defclass::hydrogens_defclass(vector<string> lines) {
  //
  // Split into sections and pass to titgroup_conformations
  //
  //printf ("what the flying fuck!\n");
  unsigned int count=0;
  while (count<lines.size()) {
    if (strip(lines[count]).compare("//END OF FILE")==0) { 
      break;
    }
    // Skip comments
    if (lines[count].substr(0,2).compare("//")==0 || lines[count].substr(0,1).compare("#")==0)  {
      count++;
      continue;
    }
    //
    if (lines[count].substr(0,1).compare("*")==0) {
      count++;
      //printf ("thisline: %s\n",lines[count].c_str());
      vector<string> words=split(lines[count]);
      string titgroupname=words[0];
      string acidbase=words[1];
      count++;
      vector<string> titgrouplines;
      //
      // Catch comments and end of file
      //
      if (strip(lines[count]).compare("END OF FILE")==0) { 
          break;
      }
      while (lines[count].substr(0,1).compare("*")!=0 && count<lines.size()) {
          printf ("LL: %s\n",lines[count].c_str());
          if (!lines[count].substr(0,2).compare("//") || !lines[count].substr(0,1).compare("#"))  {
              count++;
              if (count>=lines.size()) {
                  break;
              }
              continue;
          }
          if (strip(lines[count]).compare("END OF FILE")==0) { 
              break;
          }
	//printf ("Adding titgroupline: %s \n",lines[count].c_str());
	titgrouplines.push_back(lines[count]);
	count++;
      }
      //printf ("SEnding section down\n");
      Hbuild_titgroups.push_back(titgroup_conformations(titgroupname,acidbase,titgrouplines));
      continue;
    }
    //
    // Empty line. Increment
    //
    count++;
  }
  return;
}
    
// ----------------------------------------------
//
// Titgroup_conformations functions
//
titgroup_conformations::titgroup_conformations(string groupname,string acidbase,vector<string> lines) {
    //
    // Parse the lines for a single titgroup definition
    //
    name=groupname;
    string line;
    unsigned int count=0;
    while (count<lines.size()) {
      //printf ("HERE: %d %s\n",count,lines[count].c_str());
      // Skip comments
      if (!lines[count].substr(0,2).compare("//") || !lines[count].substr(0,1).compare("#"))  {
	count++;
	continue;
      }
      // Find conformations
      if (lines[count].substr(0,1).compare(">")==0) {
	string confname=strip(lines[count]).substr(1); // Returns rest of the string
	vector<string> conflines;
	count++;
	while (lines[count].substr(0,1).compare(">")!=0 && lines[count].substr(0,1).compare("*")!=0 && count<lines.size()) {
	  conflines.push_back(strip(lines[count]));
	  //printf ("LL: %s\n",lines[count].c_str());
	  count++;
	  if (count>=lines.size()) {
	    break;
	  }
	}
	conformations.push_back(hydrogen_conformation_class(confname,conflines)); // Instantiate this conformation
	continue;
      }
      count++;
    }
    return;
}

// -------------------------------------------------
//
// Hydrogen_conformation_class functions
//
hydrogen_conformation_class::hydrogen_conformation_class(string name,vector<string> lines) {
    // Parse the lines in a single conformation
    //printf ("Parsing this conformation: %s\n",name.c_str());
    //printf ("size of lines: %d\n",static_cast<int>(lines.size()));
    for (unsigned int count=0;count<lines.size();count++) {
      //printf ("SUB %d %s\n", count,lines[count].c_str());
    }
    //printf ("returning\n");
    return;
}
        
        

