#!/usr/bin/env python
#
# pKa - various programs and scripts for pKa value analysis, calculation and redesign
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

def parse_arguments(args,extern_defaults):
    #
    # Given a list of commandline arguments (from sys.argv) and
    # the defaults for the application, this routine updates
    # defaults with the new parameters
    # defaults is a dictionary of the form: defaults={<parm name>:[<default value>,<type>]}
    # where <type> can be 'number', 'T/F' or any other value. Only 'T/F' and 'number' types are
    # changed
    #
    defaults=extern_defaults.copy()
    import string
    args=string.join(args[1:])
    args=string.split(args,'-')
    for arg in args:
        split=string.split(string.strip(arg))
        if split==[]:
            continue
        parm_name=split[0]
        if not defaults.has_key(parm_name):
            raise 'Unknown parameter: ',parm_name
        #
        # Deal with T/F
        #
        if len(split)==1:
            if defaults[parm_name][1]=='T/F':
                if defaults[parm_name][0]:
                    defaults[parm_name][0]=None
                else:
                    defaults[parm_name][0]=1
        #
        # Deal with all the other cases
        #
        elif len(split)==2:
            if defaults[parm_name][1]=='number':
                defaults[parm_name][0]=string.atof(split[1])
            else:
                defaults[parm_name][0]=split[1]
        else:
            raise 'Incorrect usage'
    return defaults
