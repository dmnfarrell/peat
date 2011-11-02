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
#ifndef STRING_TOOLS_H
#define STRING_TOOLS_H
#include <string>
#include <vector>
#include <string>
#include <string.h>
#include <cstdlib>
//
using namespace std;
//
// Functions
//
string strip(string s); // Remove leading and trailing whitespace
string zstrip(string s); //Remove leading zeros and trailing whitespace
//
// Get all integers from a string
//
vector<int> getints(string s);
vector<float> getfloats(string s);
//
// Split string into substrings.
// Separator is whitespace
//
vector<string> split(string s);
vector<string> split(string s,string sep);
string join(vector<string> words,string joiner); // Join strings together
//
// Convert a string to a char *
//
char * s_to_char(string s);
#endif
