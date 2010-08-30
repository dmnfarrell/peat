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
#include "ostools.h"
using namespace std;

vector<string> splitpath(string path) {
    //Split a path to get the root and last part
    vector<string> result;
    string separator("/");
    vector<string> words=split(path,separator);
    string filename=words[words.size()-1];
    words.resize(words.size()-1);
    string dirname=join(words,"/");
    result.push_back(dirname);
    result.push_back(filename);
    return result;
};

bool isdirectory(string dirname) {
    // Check if dirname exists
    char pwd[200];
    bool result=true;
	getcwd(pwd, sizeof(pwd));
    if(chdir(dirname.c_str())!=0) {
        result=false;
    }
    chdir(pwd);
    return result;
};
