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
#ifndef TITGROUP_H
#define TITGROUP_H
#include <vector>
#include <string>
#include <iostream>
#include <exception>

class TitGroup {
public:
    string name;
    double intpKa;
    int acidbase; // acidbase is -1 for acids and +1 for bases
    vector<string> start_states;
    vector<string> end_states;
};

class chargestate_class {
    vector<vector<double> > charges;
    vector<double> pHs;
    vector<TitGroup> titgroups;
};

#endif
