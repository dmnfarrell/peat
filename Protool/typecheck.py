#!/usr/bin/env python
#
# # Protool - Python class for manipulating protein structures
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

import types, string, errors

def isnumber(s):
    # Returns 1 if s is a string that holds a number,
    # and None if not
    if not type(s) is types.StringType:
        raise errors.ProtoolError, "Non-string passed to isanumber %s" % s
    try:
        a=string.atof(s)
    except ValueError:
        return None
    return 1

def isint(s):
    # Returns 1 if s is a string that holds an integer, and
    # None if not
    if not type(s) is types.StringType:
        raise errors.ProtoolError, "Non-string passed to isanumber %s" % s
    try:
        a=string.atoi(s)
    except ValueError:
        return None
    return 1        

def isstring(s):
    # Return 1 if s is a string, and None if not
    if not type(s) is types.StringType:
        return None
    return 1

def containsletter(s):
    #Returns 1 if the string s contains a letter, None if not
    if not isstring(s):
        raise errors.TypeCheckError, 'Non-string passed to containsletter'
    ok=None
    for letter in s:
        if letter in string.letters:
            ok=1
    return ok
